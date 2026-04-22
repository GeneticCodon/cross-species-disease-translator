import streamlit as st
import pandas as pd
import requests
import time
from typing import List, Dict, Any, Tuple, Optional
import json
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import numpy as np

st.set_page_config(page_title="Cross-Species Disease Model Translator", layout="wide")

ENSEMBL_REST = "https://rest.ensembl.org"
OPENTARGETS_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

SUPPORTED_SPECIES = {
    "Mouse (Mus musculus)": "mus_musculus",
    "Rat (Rattus norvegicus)": "rattus_norvegicus",
    "Zebrafish (Danio rerio)": "danio_rerio",
    "Fruit fly (Drosophila melanogaster)": "drosophila_melanogaster",
    "C. elegans (Caenorhabditis elegans)": "caenorhabditis_elegans",
    "Human (Homo sapiens)": "homo_sapiens",
}

ID_PREFIX_TO_SPECIES = {
    "ENSMUSG": "mus_musculus",
    "ENSRNOG": "rattus_norvegicus",
    "ENSDARG": "danio_rerio",
    "FBgn": "drosophila_melanogaster",
    "WBGene": "caenorhabditis_elegans",
    "ENSG": "homo_sapiens",
}

# -------------- HTTP helpers --------------
def http_get_json(url: str, headers: Dict[str, Any] = None, params: Dict[str, Any] = None, retries: int = 3, backoff: float = 1.0, debug: bool = False):
    base_headers = {"Accept": "application/json", "User-Agent": "GeneticCodon-MVP/5.1"}
    if headers:
        base_headers.update(headers)
    if params is None:
        params = {}
    params.setdefault("content-type", "application/json")
    last_err = None
    
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=base_headers, params=params, timeout=25)
            if r.status_code == 200:
                if debug: st.write(f"✅ GET {r.url}")
                return r.json()
            elif r.status_code == 429:  # Rate limited
                wait_time = backoff * (2 ** attempt) + 0.5  # Exponential backoff
                last_err = f"Rate limited (429), retrying in {wait_time:.1f}s"
                if debug: st.warning(f"Attempt {attempt+1}: {last_err}")
                time.sleep(wait_time)
            elif r.status_code in (500, 502, 503, 504):  # Server errors
                wait_time = backoff * (1.5 ** attempt)
                last_err = f"Server error ({r.status_code}), retrying in {wait_time:.1f}s"
                if debug: st.warning(f"Attempt {attempt+1}: {last_err}")
                time.sleep(wait_time)
            elif r.status_code == 404:
                last_err = f"Not found (404): {url.split('/')[-1]}"
                break  # Don't retry 404s
            else:
                last_err = f"HTTP {r.status_code}: {r.text[:200]}"
                break
        except requests.Timeout:
            last_err = f"Timeout after 25s"
            time.sleep(backoff * (1.2 ** attempt))
        except requests.RequestException as e:
            last_err = f"Connection error: {str(e)[:100]}"
            time.sleep(backoff * (1.2 ** attempt))
    
    if debug and last_err: 
        st.error(f"❌ GET {url.split('/')[-2:]}: {last_err}")
    return None

def http_post_json(url: str, json_payload: Dict[str, Any], headers: Dict[str, Any] = None, retries: int = 3, backoff: float = 0.8, debug: bool = False):
    base_headers = {"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "GeneticCodon-MVP/5.0"}
    if headers:
        base_headers.update(headers)
    last_err = None
    for i in range(retries):
        try:
            r = requests.post(url, json=json_payload, headers=base_headers, timeout=25)
            if r.status_code == 200:
                if debug: st.write(f"✅ POST {url}")
                return r.json()
            elif r.status_code in (429, 500, 502, 503):
                last_err = f"HTTP {r.status_code}: {r.text[:200]}"
                time.sleep(backoff * (i + 1))
            else:
                last_err = f"HTTP {r.status_code}: {r.text[:200]}"
                break
        except requests.RequestException as e:
            last_err = str(e)
            time.sleep(backoff * (i + 1))
    if debug: 
        st.warning(f"POST {url} failed → {last_err}")
        st.code(json.dumps(json_payload, indent=2))
    return None

# -------------- Ensembl helpers --------------
@st.cache_data(show_spinner=False, ttl=24*3600)
def ensembl_symbol_to_id(species_code: str, symbol: str, debug: bool = False) -> Optional[Tuple[str, str]]:
    url = f"{ENSEMBL_REST}/xrefs/symbol/{species_code}/{symbol}"
    data = http_get_json(url, debug=debug)
    if not data:
        return None
    for entry in data:
        if str(entry.get("type", "")).lower() == "gene":
            return entry.get("id"), entry.get("display_id") or symbol
    entry = data[0] if data else None
    if entry:
        return entry.get("id"), entry.get("display_id") or symbol
    return None

@st.cache_data(show_spinner=False, ttl=24*3600)
def ensembl_archive_latest_id(ensembl_id: str, debug: bool = False) -> Optional[str]:
    url = f"{ENSEMBL_REST}/archive/id/{ensembl_id}"
    data = http_get_json(url, debug=debug)
    if not data:
        return None
    mapped = data.get("mapped")
    if mapped:
        return mapped
    mappings = data.get("mappings")
    if isinstance(mappings, list) and mappings:
        m = mappings[0].get("mapped") or mappings[0].get("target")
        if m:
            return m
    return data.get("id") or ensembl_id

@st.cache_data(show_spinner=False, ttl=24*3600)
def ensembl_id_to_symbol(ensembl_id: str, debug: bool = False) -> Optional[str]:
    url_lookup = f"{ENSEMBL_REST}/lookup/id/{ensembl_id}"
    data = http_get_json(url_lookup, debug=debug)
    if data and data.get("display_name"):
        return data["display_name"]
    url = f"{ENSEMBL_REST}/xrefs/id/{ensembl_id}"
    data = http_get_json(url, debug=debug)
    if not data:
        return None
    for e in data:
        if (e.get("dbname") or "").upper() == "HGNC" and e.get("display_id"):
            return e["display_id"]
    for e in data:
        if str(e.get("type","")).lower() == "gene" and e.get("display_id"):
            return e["display_id"]
    return next((e.get("display_id") for e in data if e.get("display_id")), None)

def _pick_best(homs, debug: bool = False):
    if not homs:
        return None
    candidates = []
    for h in homs:
        sp = (h.get("target", {}) or {}).get("species") or h.get("species") or ""
        if sp:
            sp_norm = sp.replace(" ", "_").lower()
            if sp_norm not in ("homo_sapiens", "human"):
                continue
        candidates.append(h)
    if not candidates:
        candidates = homs
    def score(h):
        t = (h.get("type") or h.get("homology_type") or "").lower()
        # Fixed: use correct patterns with underscores
        if "one_to_one" in t or "ortholog_one2one" in t:
            return 3
        if "one_to_many" in t or "many_to_one" in t:
            return 2
        if "ortholog" in t:
            return 1
        return 0
    best = sorted(candidates, key=score, reverse=True)[0]
    if debug:
        st.write("🔎 Selected homology:", f"type={best.get('type')}, score={score(best)}")
        st.write("🔎 Full object:", json.dumps(best, indent=2)[:800])
    return best

@st.cache_data(show_spinner=False, ttl=24*3600)
def ensembl_homology_to_human(source_gene: str, source_species_code: str = "mus_musculus", debug: bool = False) -> Optional[Dict[str, Any]]:
    symbol_try = source_gene
    if any(source_gene.startswith(p) for p in ID_PREFIX_TO_SPECIES):
        latest = ensembl_archive_latest_id(source_gene, debug=debug) or source_gene
        sym = ensembl_id_to_symbol(latest, debug=debug)
        if sym:
            symbol_try = sym
        else:
            symbol_try = None
            source_gene = latest
    if symbol_try:
        url_sym = f"{ENSEMBL_REST}/homology/symbol/{source_species_code}/{symbol_try}"
        data = http_get_json(url_sym, params={"target_species": "homo_sapiens", "format": "condensed"}, debug=debug)
        if data and data.get("data"):
            best = _pick_best(data["data"][0].get("homologies"), debug=debug)
            if best:
                return best
    if any((source_gene or "").startswith(p) for p in ID_PREFIX_TO_SPECIES):
        url_id = f"{ENSEMBL_REST}/homology/id/{source_gene}"
        data = http_get_json(url_id, params={"target_species": "homo_sapiens", "format": "condensed"}, debug=debug)
        if data and data.get("data"):
            best = _pick_best(data["data"][0].get("homologies"), debug=debug)
            if best:
                return best
    return None

# -------------- Open Targets (v5) --------------
OT_QUERY = """
query TargetByEnsembl($ensemblId: String!){
  target(ensemblId: $ensemblId){
    id
    approvedSymbol
    tractability{ modality label value }
    associatedDiseases{
      count
      rows{
        score
        disease{ id name }
      }
    }
  }
}
"""

@st.cache_data(show_spinner=False, ttl=24*3600)
def opentargets_for_ensembl(human_ensembl_id: str, debug: bool = False) -> Dict[str, Any]:
    payload = {"query": OT_QUERY, "variables": {"ensemblId": human_ensembl_id}}
    data = http_post_json(OPENTARGETS_GQL, payload, debug=debug)
    out = {"diseases": [], "tractability": "None"}
    try:
        t = (data or {}).get("data", {}).get("target")
        if not t:
            return out
        rows = (t.get("associatedDiseases") or {}).get("rows", []) or []
        diseases = [{"name": r.get("disease", {}).get("name"), "score": r.get("score")} for r in rows if r]
        # derive coarse druggability
        tract_list = t.get("tractability") or []
        high_labels = {"Approved Drug", "Advanced Clinical", "Phase 1 Clinical"}
        any_true = any(bool(x.get("value")) for x in tract_list)
        high = any(x.get("label") in high_labels and x.get("value") for x in tract_list)
        if high:
            tract = "High"
        elif any_true:
            tract = "Medium"
        else:
            tract = "None"
        out = {"diseases": diseases, "tractability": tract}
    except Exception:
        pass
    return out

# -------------- Core pipeline --------------
def norm_symbol(x: str) -> str:
    return (x or "").strip().replace('"', "").replace("'", "")

def validate_gene_input(gene: str) -> Tuple[bool, str]:
    """Validate a single gene symbol/ID. Returns (is_valid, error_msg)"""
    if not gene or len(gene.strip()) == 0:
        return False, "Empty gene symbol"
    if len(gene) > 50:
        return False, "Gene symbol too long (>50 chars)"
    # Allow alphanumeric, underscores, hyphens, dots
    import re
    if not re.match(r'^[a-zA-Z0-9._-]+$', gene):
        return False, "Invalid characters (only letters, numbers, ._- allowed)"
    return True, ""

def parse_input_genes(text: str) -> Tuple[List[str], List[str]]:
    """Parse input genes and return (valid_genes, errors)"""
    raw = [norm_symbol(x) for x in text.replace(",", "\n").splitlines()]
    seen = set()
    valid_genes = []
    errors = []
    
    for g in raw:
        if not g:
            continue
        if g in seen:
            continue
        seen.add(g)
        
        is_valid, error_msg = validate_gene_input(g)
        if is_valid:
            valid_genes.append(g)
        else:
            errors.append(f"{g}: {error_msg}")
    
    return valid_genes, errors

def load_genes_from_file(uploaded) -> Tuple[List[str], List[str]]:
    """Load genes from file and return (valid_genes, errors)"""
    MAX_GENES_FROM_FILE = 1000
    
    try:
        name = uploaded.name.lower()
        if name.endswith(".csv"):
            df = pd.read_csv(uploaded)
        elif name.endswith(".tsv") or name.endswith(".txt"):
            df = pd.read_csv(uploaded, sep="\t")
        elif name.endswith(".xlsx") or name.endswith(".xls"):
            df = pd.read_excel(uploaded)
        else:
            return [], ["Unsupported file format. Use CSV, TSV/TXT, or Excel."]
        
        # Find gene column
        gene_col = None
        for col in df.columns:
            if col.lower() in ("gene", "gene_symbol", "symbol", "genes", "ensembl", "gene_id"):
                gene_col = col
                break
        
        if gene_col is None:
            gene_col = df.columns[0]  # Use first column as fallback
            
        raw_genes = [norm_symbol(str(x)) for x in df[gene_col].dropna().tolist()]
        
        if len(raw_genes) > MAX_GENES_FROM_FILE:
            return [], [f"Too many genes in file ({len(raw_genes)}). Maximum allowed: {MAX_GENES_FROM_FILE}"]
        
        # Validate each gene
        valid_genes = []
        errors = []
        seen = set()
        
        for g in raw_genes:
            if not g or g in seen:
                continue
            seen.add(g)
            
            is_valid, error_msg = validate_gene_input(g)
            if is_valid:
                valid_genes.append(g)
            else:
                errors.append(f"{g}: {error_msg}")
                
        return valid_genes, errors[:20]  # Limit error display
        
    except Exception as e:
        return [], [f"Failed to read file: {str(e)}"]


def run_pipeline(input_genes: List[str], source_species_code: str, debug: bool = False) -> pd.DataFrame:
    MAX_BATCH_SIZE = 200
    
    if len(input_genes) > MAX_BATCH_SIZE:
        st.warning(f"⚠️ Too many genes ({len(input_genes)}). Processing first {MAX_BATCH_SIZE} only.")
        input_genes = input_genes[:MAX_BATCH_SIZE]
    
    rows = []
    unresolved = []
    api_failures = []
    
    # Create progress bar
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, g in enumerate(input_genes):
        # Update progress
        progress = (i + 1) / len(input_genes)
        progress_bar.progress(progress)
        status_text.text(f"Processing gene {i+1}/{len(input_genes)}: {g}")
        
        try:
            if any(g.startswith(pfx) for pfx in ID_PREFIX_TO_SPECIES):
                ens_id = g
                sym = ensembl_id_to_symbol(g, debug=debug) or g
            else:
                sym = g
                r = ensembl_symbol_to_id(source_species_code, g, debug=debug)
                ens_id = r[0] if r else None
                if not r:
                    unresolved.append(g)

            h = ensembl_homology_to_human(sym if sym else (ens_id or g), source_species_code, debug=debug)
            if not h:
                rows.append({
                    "source_symbol": sym,
                    "source_gene_id": ens_id,
                    "human_symbol": None,
                    "human_gene_id": None,
                    "orthology_type": "No homolog found",
                    "diseases": [],
                    "disease_count": 0,
                    "druggability": "None",
                })
                continue

            human_id = h.get("target", {}).get("id") or h.get("id")
            human_symbol = ensembl_id_to_symbol(human_id, debug=debug) if human_id else None
            orth_type = h.get("type") or h.get("homology_type") or "Unknown"

            # Try to get disease data with error handling
            try:
                ot = opentargets_for_ensembl(human_id, debug=debug) if human_id else {"diseases": [], "tractability": "None"}
                diseases = ot.get("diseases", []) or []
                druggability = ot.get("tractability", "None")
            except Exception as e:
                if debug: st.warning(f"Open Targets failed for {human_id}: {e}")
                diseases = []
                druggability = "API Error"
                api_failures.append(f"{g} (Open Targets)")

            rows.append({
                "source_symbol": sym,
                "source_gene_id": ens_id,
                "human_symbol": human_symbol,
                "human_gene_id": human_id,
                "orthology_type": orth_type,
                "diseases": diseases,
                "disease_count": len(diseases),
                "druggability": druggability,
            })
            
        except Exception as e:
            if debug: st.error(f"Pipeline failed for {g}: {e}")
            api_failures.append(f"{g} (Pipeline error)")
            # Add empty row for failed gene
            rows.append({
                "source_symbol": g,
                "source_gene_id": None,
                "human_symbol": None, 
                "human_gene_id": None,
                "orthology_type": "Processing failed",
                "diseases": [],
                "disease_count": 0,
                "druggability": "Error",
            })
        
        # Rate limiting
        time.sleep(0.08)
    
    # Clear progress indicators
    progress_bar.empty()
    status_text.empty()
    
    df = pd.DataFrame(rows)
    
    # Show summary of issues
    if unresolved:
        st.info(f"🔍 **Unresolved symbols** ({len(unresolved)}): {', '.join(unresolved[:10])}" + (" ..." if len(unresolved)>10 else ""))
    if api_failures:
        st.warning(f"⚠️ **API failures** ({len(api_failures)}): {', '.join(api_failures[:5])}" + (" ..." if len(api_failures)>5 else ""))
    
    return df

def diseases_to_str(dlist: List[Dict[str, Any]]) -> str:
    if not dlist:
        return ""
    tops = sorted([d for d in dlist if d.get("name")], key=lambda x: -(x.get("score") or 0))[:5]
    return "; ".join([f"{d['name']}" + (f" ({d['score']:.2f})" if d.get("score") is not None else "") for d in tops])

# -------------- Visualization functions --------------
def create_orthology_pie_chart(df: pd.DataFrame) -> go.Figure:
    """Create pie chart showing distribution of orthology types"""
    orthology_counts = df['orthology_type'].fillna('No ortholog').value_counts()
    
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD']
    
    fig = go.Figure(data=[go.Pie(
        labels=orthology_counts.index,
        values=orthology_counts.values,
        hole=0.3,
        marker_colors=colors[:len(orthology_counts)]
    )])
    
    fig.update_traces(textposition='inside', textinfo='percent+label')
    fig.update_layout(
        title="Distribution of Orthology Types",
        height=400,
        showlegend=True,
        legend=dict(orientation="v", yanchor="middle", y=0.5)
    )
    return fig

def create_druggability_bar_chart(df: pd.DataFrame) -> go.Figure:
    """Create bar chart showing druggability distribution"""
    drug_counts = df['druggability'].fillna('None').value_counts()
    
    # Define colors for druggability levels
    color_map = {'High': '#2ECC71', 'Medium': '#F39C12', 'None': '#E74C3C', 'API Error': '#95A5A6', 'Error': '#95A5A6'}
    colors = [color_map.get(label, '#BDC3C7') for label in drug_counts.index]
    
    fig = go.Figure(data=[go.Bar(
        x=drug_counts.index,
        y=drug_counts.values,
        marker_color=colors,
        text=drug_counts.values,
        textposition='auto'
    )])
    
    fig.update_layout(
        title="Druggability Assessment Distribution",
        xaxis_title="Druggability Level",
        yaxis_title="Number of Genes",
        height=400
    )
    return fig

def create_disease_histogram(df: pd.DataFrame) -> go.Figure:
    """Create histogram of disease counts per gene"""
    disease_counts = df['disease_count'].fillna(0)
    
    fig = go.Figure(data=[go.Histogram(
        x=disease_counts,
        nbinsx=min(20, disease_counts.max() + 1),
        marker_color='#3498DB',
        opacity=0.7
    )])
    
    fig.update_layout(
        title="Distribution of Disease Associations per Gene",
        xaxis_title="Number of Associated Diseases",
        yaxis_title="Number of Genes",
        height=400
    )
    return fig

def create_top_diseases_chart(df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    """Create bar chart of most common diseases across all genes"""
    all_diseases = []
    for diseases_list in df['diseases']:
        if diseases_list:
            for disease in diseases_list:
                if disease.get('name'):
                    all_diseases.append(disease['name'])
    
    if not all_diseases:
        # Return empty chart if no diseases found
        fig = go.Figure()
        fig.add_annotation(
            text="No disease associations found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16, color="gray")
        )
        fig.update_layout(title="Top Associated Diseases", height=400)
        return fig
    
    disease_counts = Counter(all_diseases)
    top_diseases = disease_counts.most_common(top_n)
    
    diseases, counts = zip(*top_diseases) if top_diseases else ([], [])
    
    fig = go.Figure(data=[go.Bar(
        y=list(diseases)[::-1],  # Reverse for horizontal bar chart
        x=list(counts)[::-1],
        orientation='h',
        marker_color='#E74C3C',
        text=list(counts)[::-1],
        textposition='auto'
    )])
    
    fig.update_layout(
        title=f"Top {len(diseases)} Most Common Disease Associations",
        xaxis_title="Number of Genes Associated",
        yaxis_title="Disease",
        height=max(400, len(diseases) * 25),
        margin=dict(l=200)  # More space for disease names
    )
    return fig

def create_species_detection_chart(genes: List[str]) -> go.Figure:
    """Create chart showing detected species from gene IDs"""
    species_counts = {}
    unidentified = 0
    
    for gene in genes:
        detected = False
        for prefix, species in ID_PREFIX_TO_SPECIES.items():
            if gene.startswith(prefix):
                species_name = species.replace('_', ' ').title()
                species_counts[species_name] = species_counts.get(species_name, 0) + 1
                detected = True
                break
        if not detected:
            unidentified += 1
    
    if unidentified > 0:
        species_counts['Gene Symbols (Species TBD)'] = unidentified
    
    if not species_counts:
        fig = go.Figure()
        fig.add_annotation(
            text="No genes to analyze",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        fig.update_layout(title="Input Gene Species Detection", height=300)
        return fig
    
    fig = go.Figure(data=[go.Bar(
        x=list(species_counts.keys()),
        y=list(species_counts.values()),
        marker_color=['#3498DB' if 'TBD' not in k else '#95A5A6' for k in species_counts.keys()],
        text=list(species_counts.values()),
        textposition='auto'
    )])
    
    fig.update_layout(
        title="Detected Species from Input Genes",
        xaxis_title="Species",
        yaxis_title="Number of Genes",
        height=350,
        xaxis_tickangle=-45
    )
    return fig

def create_success_metrics_chart(df: pd.DataFrame, total_input: int) -> go.Figure:
    """Create funnel chart showing mapping success rates"""
    total_genes = total_input
    mapped_to_human = int(df['human_symbol'].notna().sum())
    with_diseases = int((df['disease_count'] > 0).sum())
    high_druggability = int((df['druggability'] == 'High').sum())
    
    stages = ['Input Genes', 'Mapped to Human', 'With Disease Links', 'High Druggability']
    values = [total_genes, mapped_to_human, with_diseases, high_druggability]
    
    fig = go.Figure(go.Funnel(
        y=stages,
        x=values,
        textinfo="value+percent initial",
        marker_color=['#3498DB', '#2ECC71', '#F39C12', '#E74C3C']
    ))
    
    fig.update_layout(
        title="Gene Mapping Success Funnel",
        height=400
    )
    return fig

def create_disease_score_scatter(df: pd.DataFrame) -> go.Figure:
    """Create scatter plot of disease count vs average disease scores"""
    scatter_data = []
    
    for _, row in df.iterrows():
        if row['diseases'] and len(row['diseases']) > 0:
            disease_count = len(row['diseases'])
            avg_score = np.mean([d.get('score', 0) or 0 for d in row['diseases']])
            scatter_data.append({
                'gene': row['source_symbol'] or row['source_gene_id'],
                'disease_count': disease_count,
                'avg_score': avg_score,
                'druggability': row['druggability']
            })
    
    if not scatter_data:
        fig = go.Figure()
        fig.add_annotation(
            text="No disease score data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        fig.update_layout(title="Disease Count vs Average Score", height=400)
        return fig
    
    scatter_df = pd.DataFrame(scatter_data)
    
    fig = px.scatter(
        scatter_df,
        x='disease_count',
        y='avg_score',
        color='druggability',
        hover_data=['gene'],
        title="Disease Associations: Count vs Average Score",
        labels={
            'disease_count': 'Number of Disease Associations',
            'avg_score': 'Average Disease Association Score'
        },
        color_discrete_map={'High': '#2ECC71', 'Medium': '#F39C12', 'None': '#E74C3C'}
    )
    
    fig.update_layout(height=400)
    return fig

# -------------- UI --------------
st.title("🧬 Cross-Species Disease Model Translator")

with st.expander("ℹ️ What is this tool? (click to expand)"):
    st.markdown("""
    Map **model organism genes** to **human orthologs** and attach disease + druggability hints (Open Targets).
    **Research-only MVP** — not for clinical use.
    """)

with st.sidebar:
    st.header("Input")
    species_label = st.selectbox("Source species", list(SUPPORTED_SPECIES.keys()), index=0)
    species_code = SUPPORTED_SPECIES[species_label]
    st.caption("👆 Choose the organism of your input genes (the tool maps to **human**).")
    demo = st.toggle("Use demo mouse genes", value=False)
    uploaded = st.file_uploader("Upload gene list (CSV / TSV / Excel)", type=["csv", "tsv", "txt", "xlsx", "xls"])
    pasted = st.text_area("Or paste gene symbols / Ensembl IDs (one per line or comma-separated)", height=140, placeholder="Akt1\nBraf\nGfap\nMapt\nSnca\nApp")
    debug_mode = st.toggle("Debug API calls (dev only)", value=True)
    run = st.button("Run Mapping →", type="primary", use_container_width=True)

genes: List[str] = []
all_errors: List[str] = []

if demo:
    genes = ["Akt1", "Braf", "Gfap", "Mapt", "Snca", "App", "Tp53", "Pten", "Egfr", "Kras"]
elif uploaded is not None:
    genes, file_errors = load_genes_from_file(uploaded)
    all_errors.extend(file_errors)
else:
    genes, parse_errors = parse_input_genes(pasted)
    all_errors.extend(parse_errors)

# Show validation errors
if all_errors:
    with st.expander(f"⚠️ Input validation issues ({len(all_errors)})", expanded=len(genes)==0):
        for error in all_errors[:20]:  # Limit display
            st.text(f"• {error}")
        if len(all_errors) > 20:
            st.text(f"... and {len(all_errors)-20} more")

# Show gene count
if genes:
    st.info(f"📋 Ready to process **{len(genes)}** valid genes")

if run:
    if not genes:
        st.warning("⚠️ No valid genes to process. Please upload or paste a gene list first.")
        if all_errors:
            st.error("Fix the input validation issues shown above.")
    else:
        with st.spinner("🔄 Mapping orthologs and fetching disease relevance…"):
            # Auto-detect species from Ensembl IDs
            if species_code == "homo_sapiens":
                for g in genes[:5]:  # Check first few genes only
                    for pfx, sp in ID_PREFIX_TO_SPECIES.items():
                        if g.startswith(pfx) and sp != "homo_sapiens":
                            st.info(f"🔍 Auto-detected {pfx} IDs → switching to **{sp.replace('_',' ').title()}**")
                            species_code = sp
                            break
                    if species_code != "homo_sapiens":
                        break
            
            df = run_pipeline(genes, species_code, debug=debug_mode)

        if df.empty:
            st.warning("No results to show.")
        else:
            st.subheader("📊 Results Summary")
            
            # Key metrics
            c1, c2, c3, c4 = st.columns(4)
            with c1: st.metric("Input genes", len(genes))
            with c2: st.metric("Mapped to human", int(df["human_symbol"].notna().sum()))
            with c3: st.metric("With ≥1 disease link", int((df["disease_count"]>0).sum()))
            with c4: st.metric("High druggability", int((df["druggability"]=="High").sum()))
            
            # Create visualization tabs
            tab1, tab2, tab3, tab4 = st.tabs(["📈 Overview Charts", "🎯 Disease Analysis", "🧬 Gene Details", "📋 Data Table"])
            
            with tab1:
                st.markdown("### Overview Visualizations")
                
                # First row of charts
                col1, col2 = st.columns(2)
                with col1:
                    if len(genes) > 0:
                        fig_species = create_species_detection_chart(genes)
                        st.plotly_chart(fig_species, use_container_width=True)
                    
                    fig_success = create_success_metrics_chart(df, len(genes))
                    st.plotly_chart(fig_success, use_container_width=True)
                    
                with col2:
                    fig_orthology = create_orthology_pie_chart(df)
                    st.plotly_chart(fig_orthology, use_container_width=True)
                    
                    fig_druggability = create_druggability_bar_chart(df)
                    st.plotly_chart(fig_druggability, use_container_width=True)
            
            with tab2:
                st.markdown("### Disease Association Analysis")
                
                col1, col2 = st.columns(2)
                with col1:
                    fig_disease_hist = create_disease_histogram(df)
                    st.plotly_chart(fig_disease_hist, use_container_width=True)
                    
                with col2:
                    fig_scatter = create_disease_score_scatter(df)
                    st.plotly_chart(fig_scatter, use_container_width=True)
                
                # Full-width chart for top diseases
                fig_top_diseases = create_top_diseases_chart(df, top_n=20)
                st.plotly_chart(fig_top_diseases, use_container_width=True)
                
                # Disease insights
                total_diseases = sum(len(d) for d in df['diseases'] if d)
                unique_diseases = len(set(d['name'] for diseases in df['diseases'] if diseases for d in diseases if d.get('name')))
                
                if total_diseases > 0:
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total Disease Links", total_diseases)
                    with col2:
                        st.metric("Unique Diseases", unique_diseases)
                    with col3:
                        avg_diseases = total_diseases / len([d for d in df['diseases'] if d])
                        st.metric("Avg Diseases/Gene", f"{avg_diseases:.1f}")
            
            with tab3:
                st.markdown("### Individual Gene Details")
                
                # Filter options
                col1, col2, col3 = st.columns(3)
                with col1:
                    orthology_filter = st.selectbox(
                        "Filter by Orthology Type:",
                        ["All"] + list(df['orthology_type'].dropna().unique())
                    )
                with col2:
                    druggability_filter = st.selectbox(
                        "Filter by Druggability:",
                        ["All"] + list(df['druggability'].dropna().unique())
                    )
                with col3:
                    min_diseases = st.number_input(
                        "Min Disease Count:",
                        min_value=0,
                        max_value=int(df['disease_count'].max()) if not df['disease_count'].empty else 0,
                        value=0
                    )
                
                # Apply filters
                filtered_df = df.copy()
                if orthology_filter != "All":
                    filtered_df = filtered_df[filtered_df['orthology_type'] == orthology_filter]
                if druggability_filter != "All":
                    filtered_df = filtered_df[filtered_df['druggability'] == druggability_filter]
                if min_diseases > 0:
                    filtered_df = filtered_df[filtered_df['disease_count'] >= min_diseases]
                
                st.write(f"Showing {len(filtered_df)} genes (filtered from {len(df)})")
                
                # Gene details with expanded disease info
                if not filtered_df.empty:
                    for _, row in filtered_df.head(20).iterrows():  # Limit to first 20 for performance
                        with st.expander(f"🧬 {row['source_symbol'] or row['source_gene_id']} → {row['human_symbol'] or 'No human ortholog'}"):
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write(f"**Source ID:** {row['source_gene_id'] or 'N/A'}")
                                st.write(f"**Human ID:** {row['human_gene_id'] or 'N/A'}")
                                st.write(f"**Orthology:** {row['orthology_type'] or 'N/A'}")
                                st.write(f"**Druggability:** {row['druggability']}")
                            with col2:
                                if row['diseases']:
                                    st.write(f"**Disease Associations ({len(row['diseases'])}):**")
                                    for i, disease in enumerate(row['diseases'][:5]):
                                        score_str = f" (score: {disease['score']:.3f})" if disease.get('score') else ""
                                        st.write(f"  {i+1}. {disease['name']}{score_str}")
                                    if len(row['diseases']) > 5:
                                        st.write(f"  ... and {len(row['diseases']) - 5} more")
                                else:
                                    st.write("**No disease associations found**")
            
            with tab4:
                st.markdown("### Complete Data Table")
                
                df_view = df.copy()
                df_view["Diseases (top)"] = df_view["diseases"].apply(diseases_to_str)
                df_view = df_view.drop(columns=["diseases"])
                
                # Make the table more readable
                display_columns = [
                    "source_symbol", "source_gene_id", "human_symbol", "human_gene_id",
                    "orthology_type", "disease_count", "Diseases (top)", "druggability"
                ]
                
                st.dataframe(
                    df_view[display_columns],
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "source_symbol": "Source Gene",
                        "source_gene_id": "Source ID",
                        "human_symbol": "Human Gene",
                        "human_gene_id": "Human ID",
                        "orthology_type": "Orthology Type",
                        "disease_count": st.column_config.NumberColumn("Disease Count", format="%d"),
                        "Diseases (top)": "Top Diseases",
                        "druggability": "Druggability"
                    }
                )
                
                # Download options
                col1, col2 = st.columns(2)
                with col1:
                    csv_bytes = df.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        "⬇️ Download Complete Results (CSV)",
                        csv_bytes,
                        file_name=f"cross_species_results_{len(genes)}_genes.csv",
                        mime="text/csv"
                    )
                with col2:
                    # Create a summary CSV with just key metrics
                    summary_df = df_view[display_columns].copy()
                    summary_csv = summary_df.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        "⬇️ Download Summary Table (CSV)",
                        summary_csv,
                        file_name=f"cross_species_summary_{len(genes)}_genes.csv",
                        mime="text/csv"
                    )
                
                with st.expander("🔧 Debug: Raw Data (JSON)"):
                    st.code(json.dumps(df.head(5).to_dict(orient="records"), indent=2, default=str))
else:
    st.info("Upload or paste a gene list, choose the source species, then click **Run Mapping →**. Tip: enable **demo** in the sidebar to try it instantly.")