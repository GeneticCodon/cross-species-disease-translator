# Cross-Species Disease Model Translator — Full Documentation

**Version:** 5.0  
**Maintained by:** Genetic Codon ([geneticcodon.com](https://geneticcodon.com))  
**License:** MIT

---

## Table of Contents

1. [Background and Motivation](#1-background-and-motivation)
2. [Biological Concepts Explained](#2-biological-concepts-explained)
3. [How the App Works — Step by Step](#3-how-the-app-works--step-by-step)
4. [Input Guide](#4-input-guide)
5. [Understanding the Results](#5-understanding-the-results)
6. [Understanding the Visualizations](#6-understanding-the-visualizations)
7. [External APIs and Data Sources](#7-external-apis-and-data-sources)
8. [Technical Architecture](#8-technical-architecture)
9. [Known Limitations](#9-known-limitations)
10. [Worked Example](#10-worked-example)
11. [FAQ](#11-faq)

---

## 1. Background and Motivation

Most of what we know about human disease at the molecular level was first discovered in a model organism — a mouse, a zebrafish, a fruit fly, or a nematode worm. These organisms are used because:

- Experiments that would be unethical in humans are feasible in animals
- Their genomes can be precisely edited (knockouts, knockins, overexpression)
- Their short lifespans allow observation of disease progression over weeks rather than decades
- Large numbers of individuals can be studied at low cost

The critical challenge is **translational relevance**: when you find that knocking out gene X in a mouse causes neurodegeneration, you need to know whether gene X has a human equivalent, whether that human gene is already known to be involved in human neurological disease, and whether any drugs or drug candidates already target it.

Manually answering those three questions requires searching Ensembl, browsing Open Targets, reading literature, and manually compiling results — often taking hours per gene and days for a list of 50–100 genes.

This app automates that entire process. You provide a gene list from your model organism experiment; the app returns human orthologs with disease associations and druggability ratings in minutes.

---

## 2. Biological Concepts Explained

### Orthologs

Orthologs are genes in different species that descended from the **same gene in a common ancestor**. They typically perform the same biological function. For example, the mouse gene *Tp53* and the human gene *TP53* are orthologs — both encode the p53 tumour suppressor protein.

This is different from **paralogs**, which are genes related by duplication within the same genome. Paralogs often have diverged functions and are less reliable for cross-species translation.

### Orthology Types

| Type | Meaning | Translational reliability |
|------|---------|--------------------------|
| one-to-one | One gene in species A maps to exactly one gene in species B | Highest — function is most likely conserved |
| one-to-many | One gene in species A maps to multiple genes in species B | Moderate — one of the human copies may be relevant |
| many-to-one | Multiple genes in species A map to one human gene | Moderate — the human gene may consolidate functions |
| many-to-many | Multiple genes on both sides | Lowest confidence — interpret cautiously |

The app scores one-to-one orthologs highest and selects the best available match when multiple homologs exist.

### Evolutionary Conservation

Genes that are conserved across very distant species (e.g. a fruit fly gene that has a clear human ortholog) tend to be involved in **fundamental cellular processes** — cell division, DNA repair, protein folding, synaptic signalling. These are often the most disease-relevant genes because they are also the most mutation-intolerant.

The fact that the same gene has been maintained for 600 million years of evolution (from worm to human) is itself strong evidence that the gene performs a critical, non-redundant function.

### Disease Association Scores

Open Targets aggregates evidence from multiple independent sources and calculates an association score between 0 and 1 for every gene-disease pair:

| Score range | Interpretation |
|-------------|---------------|
| 0.00 – 0.09 | No or very weak evidence |
| 0.10 – 0.24 | Limited evidence, single source |
| 0.25 – 0.49 | Moderate evidence, some consistency |
| 0.50 – 0.74 | Good evidence from multiple sources |
| 0.75 – 1.00 | Strong evidence, well-validated association |

Evidence sources include: GWAS, rare variant studies, somatic mutations, gene expression changes in disease tissue, animal model phenotypes, drug clinical trial data, and curated literature.

### Druggability

Druggability refers to whether a protein target can be modulated by a drug molecule. Open Targets classifies targets into tractability tiers:

- **High**: An approved drug or advanced clinical candidate (Phase 2/3) already targets this protein. The protein is known to be accessible to small molecules or biologics.
- **Medium**: Some tractability evidence exists (early clinical trials, structural data suggesting a binding site, tool compounds in research).
- **None**: No tractability evidence currently available. The target may still be scientifically important but is harder to pursue therapeutically.

---

## 3. How the App Works — Step by Step

```
User input (gene symbols or Ensembl IDs)
            |
            v
    Input validation
    - Remove duplicates
    - Check for invalid characters
    - Enforce 200-gene batch limit
            |
            v
    For each gene:
            |
    [If gene symbol]                    [If Ensembl ID]
    Ensembl xrefs/symbol API           Ensembl archive/id API
    -> resolve to Ensembl ID           -> get latest stable ID
            |                                   |
            +-----------------------------------+
                            |
                            v
              Ensembl homology API
              -> find human ortholog
              -> score and select best match
              (one-to-one > one-to-many > many-to-many)
                            |
                            v
              Resolve human Ensembl ID -> human gene symbol
                            |
                            v
              Open Targets GraphQL query
              -> disease associations (name + score)
              -> tractability assessments
                            |
                            v
              Aggregate results into DataFrame
                            |
                            v
              Generate visualizations
              Display results in tabbed interface
              Offer CSV download
```

### API Retry Logic

All API calls use exponential backoff:
- On HTTP 429 (rate limit): waits 2× longer on each retry
- On HTTP 500/502/503/504 (server error): retries up to 3 times
- On HTTP 404 (not found): fails immediately, no retry
- Timeout: 25 seconds per request

Results are cached for 24 hours using Streamlit's `@st.cache_data` decorator, so repeated queries for the same genes are instant.

---

## 4. Input Guide

### Accepted Gene Formats

| Format | Example | Notes |
|--------|---------|-------|
| Gene symbol | `Akt1`, `Tp53`, `Braf` | Case-insensitive for most species |
| Mouse Ensembl ID | `ENSMUSG00000001156` | Auto-detected by prefix |
| Human Ensembl ID | `ENSG00000139618` | Auto-detected |
| Rat Ensembl ID | `ENSRNOG00000004540` | Auto-detected |
| Zebrafish Ensembl ID | `ENSDARG00000037912` | Auto-detected |
| Fly gene ID | `FBgn0003721` | Auto-detected |
| Worm gene ID | `WBGene00000915` | Auto-detected |

### Accepted File Formats

| Format | Required column name (any of) |
|--------|-------------------------------|
| CSV | `gene`, `gene_symbol`, `symbol`, `genes`, `ensembl`, `gene_id` |
| TSV / TXT | Same as above |
| Excel (.xlsx, .xls) | Same as above |

If no recognised column name is found, the first column is used automatically.

### Batch Size Limits

- **Pasted input or file upload**: up to 200 genes per run
- **File upload**: up to 1000 genes in the file (first 200 are processed)

These limits exist to respect Ensembl and Open Targets API rate limits and to keep run times reasonable.

---

## 5. Understanding the Results

### Results Table Columns

| Column | Description |
|--------|-------------|
| `source_symbol` | The gene symbol as you entered it |
| `source_gene_id` | The resolved Ensembl ID in the source species |
| `human_symbol` | The human ortholog gene symbol (e.g. BRAF, TP53) |
| `human_gene_id` | The human Ensembl ID (ENSG...) |
| `orthology_type` | Type of homology relationship (see section 2) |
| `disease_count` | Number of human diseases associated with this gene in Open Targets |
| `Diseases (top)` | Top 5 diseases by association score, with scores in brackets |
| `druggability` | High / Medium / None (see section 2) |

### What "No homolog found" Means

If the orthology column shows "No homolog found", it means either:
1. The gene is species-specific (no human equivalent exists — common in some immune or olfactory genes)
2. The gene symbol could not be resolved to an Ensembl ID in the source species (check spelling)
3. The Ensembl API returned no homology data at the time of query

### What Empty Disease Fields Mean

An empty disease field does not mean the gene is unimportant — it may mean:
- The human ortholog has not been studied enough to appear in GWAS or clinical datasets
- The gene is relevant to biology but not yet linked to a named disease in Open Targets
- The Open Targets API call failed for that gene (check the API failures warning if shown)

---

## 6. Understanding the Visualizations

### Overview Tab

**Species Detection Chart**  
Shows how many of your input genes were recognised as Ensembl IDs from specific species vs. plain gene symbols. Useful for confirming the app interpreted your input correctly.

**Gene Mapping Success Funnel**  
A step-by-step funnel showing how many genes survived each stage: input → mapped to human → have disease links → high druggability. A steep drop at any stage tells you something specific — e.g. a large drop at "mapped to human" suggests many species-specific genes or symbol resolution failures.

**Orthology Type Distribution**  
Pie chart of the orthology types across your gene set. A high proportion of one-to-one orthologs is a good sign for translational confidence.

**Druggability Assessment**  
Bar chart of the High/Medium/None breakdown. If you are doing drug target prioritisation, the ratio of High to None is a key metric.

### Disease Analysis Tab

**Disease Count Histogram**  
Distribution of how many diseases each gene is associated with. Genes with very high disease counts (right tail) are often well-studied "hub" genes like TP53 or EGFR — they may not be the most specific targets. Genes with moderate counts (3–15) are often more actionable.

**Disease Count vs. Average Score Scatter Plot**  
Each point is one gene. X-axis is the number of disease associations; Y-axis is the average confidence score. Points in the top-right are the most broadly and confidently disease-associated genes. Colour represents druggability. This is the most information-dense chart in the app — use it for prioritisation decisions.

**Top Diseases Chart**  
Which diseases appear most frequently across your entire gene set. If you are working on a specific disease (e.g. Parkinson's), you should see it dominate this chart for relevant gene sets. Unexpected diseases appearing here can reveal off-target biology worth investigating.

### Gene Details Tab

Interactive per-gene explorer with three filters:
- **Orthology type**: restrict to high-confidence (one-to-one) orthologs only
- **Druggability**: show only High or Medium druggability targets
- **Minimum disease count**: filter out genes with few associations

Each gene card expands to show all top disease associations with their individual scores, the source and human gene IDs, and orthology type.

---

## 7. External APIs and Data Sources

### Ensembl REST API

- **Base URL**: `https://rest.ensembl.org`
- **Endpoints used**:
  - `/xrefs/symbol/{species}/{symbol}` — resolve gene symbol to Ensembl ID
  - `/lookup/id/{id}` — get gene metadata including display name
  - `/xrefs/id/{id}` — cross-references for an Ensembl ID
  - `/archive/id/{id}` — resolve retired IDs to current stable IDs
  - `/homology/symbol/{species}/{symbol}` — find orthologs by symbol
  - `/homology/id/{id}` — find orthologs by Ensembl ID
- **Rate limits**: ~15 requests/second; the app enforces 0.08s delay between calls
- **No API key required**

### Open Targets Platform API

- **Base URL**: `https://api.platform.opentargets.org/api/v4/graphql`
- **Query**: GraphQL query fetching disease associations and tractability assessments per human Ensembl ID
- **Data includes**: disease name, EFO ID, association score, tractability modality, label, and value
- **No API key required**
- **Data is updated** quarterly by the Open Targets consortium

### Data Freshness

Both APIs are queried live and results are cached locally for 24 hours per session. The underlying databases (Ensembl gene annotations, Open Targets associations) are updated on their own release schedules — Ensembl releases roughly every 2–3 months, Open Targets every quarter.

---

## 8. Technical Architecture

```
app.py
  |
  |-- HTTP helpers
  |     |-- http_get_json()     GET with retry + exponential backoff
  |     |-- http_post_json()    POST with retry + exponential backoff
  |
  |-- Ensembl helpers (all @st.cache_data, TTL 24h)
  |     |-- ensembl_symbol_to_id()        symbol -> Ensembl ID
  |     |-- ensembl_id_to_symbol()        Ensembl ID -> symbol
  |     |-- ensembl_archive_latest_id()   retired ID -> current ID
  |     |-- ensembl_homology_to_human()   find best human ortholog
  |     |-- _pick_best()                  score and select best homolog
  |
  |-- Open Targets helper (cached)
  |     |-- opentargets_for_ensembl()     diseases + druggability
  |
  |-- Input processing
  |     |-- validate_gene_input()         per-gene validation
  |     |-- norm_symbol()                 strip quotes/whitespace
  |     |-- parse_input_genes()           text -> validated list
  |     |-- load_genes_from_file()        file upload -> validated list
  |
  |-- Core pipeline
  |     |-- run_pipeline()                orchestrates all API calls, builds DataFrame
  |
  |-- Visualization functions
  |     |-- create_orthology_pie_chart()
  |     |-- create_druggability_bar_chart()
  |     |-- create_disease_histogram()
  |     |-- create_top_diseases_chart()
  |     |-- create_species_detection_chart()
  |     |-- create_success_metrics_chart()
  |     |-- create_disease_score_scatter()
  |
  |-- Streamlit UI
        |-- Sidebar: species, file upload, text input, debug toggle
        |-- Main: metrics, 4 tabs (Overview / Disease / Gene Details / Table)
```

### Caching Strategy

Expensive API calls are wrapped in `@st.cache_data(ttl=24*3600)`. This means:
- The first run for a given gene is slow (live API calls)
- Subsequent runs within 24 hours are instant (served from memory cache)
- The cache is per-session (not shared across users on the same server)

---

## 9. Known Limitations

| Limitation | Impact | Workaround |
|-----------|--------|------------|
| Only maps to human | Cannot compare mouse vs. zebrafish directly | Run two separate queries, compare manually |
| 200-gene batch limit | Large screens require multiple runs | Split your list and combine the CSV downloads |
| 24h cache TTL | Very recent Ensembl updates not reflected immediately | Clear Streamlit cache or restart app |
| one-to-one scoring relies on Ensembl annotation | Occasionally misses functional orthologs in families with complex duplication history | Cross-check important hits manually in Ensembl |
| Open Targets data is quarterly | Very recent publications not included | Check literature directly for cutting-edge findings |
| No statistical correction | Disease scores are not corrected for gene length or study bias | Treat scores as evidence strength, not p-values |
| Research use only | Not validated for clinical or diagnostic decisions | — |

---

## 10. Worked Example

**Scenario**: You have run an RNA-seq experiment in a mouse model of Parkinson's disease and identified 6 significantly dysregulated genes. You want to know which of these are relevant to human Parkinson's and whether any are druggable.

**Input** (paste into the app with Mouse selected as source species):
```
Snca
Mapt
Lrrk2
Pink1
Prkn
Gba
```

**Expected output summary**:

| Source | Human ortholog | Orthology | Top disease | Score | Druggability |
|--------|---------------|-----------|-------------|-------|-------------|
| Snca | SNCA | one-to-one | Parkinson's disease | 0.98 | High |
| Mapt | MAPT | one-to-one | Alzheimer's disease | 0.95 | Medium |
| Lrrk2 | LRRK2 | one-to-one | Parkinson's disease | 0.96 | High |
| Pink1 | PINK1 | one-to-one | Parkinson's disease | 0.91 | Medium |
| Prkn | PRKN | one-to-one | Parkinson's disease | 0.93 | Medium |
| Gba | GBA | one-to-one | Gaucher disease | 0.99 | High |

**Interpretation**:
- All 6 genes have clear one-to-one human orthologs — high translational confidence
- 5 of 6 are directly linked to Parkinson's disease in humans with very high scores
- SNCA, LRRK2, and GBA are High druggability — clinical programmes already exist
- MAPT's top association is Alzheimer's but is also a known Parkinson's risk gene — worth checking its full disease list in the Gene Details tab
- This result strongly validates your mouse model and prioritises SNCA and LRRK2 as the most actionable targets

---

## 11. FAQ

**Q: My gene symbol is not being found. What should I do?**  
A: Check the exact capitalisation used in the source species database (mouse genes are typically title case: `Akt1`; human genes are all caps: `AKT1`). Try the Ensembl ID instead of the symbol. Also confirm you have selected the correct source species in the sidebar.

**Q: Why does a gene show "No homolog found" when I know a human ortholog exists?**  
A: This can happen if the gene is very recently annotated and Ensembl's homology pipeline has not yet processed it, or if the gene belongs to a rapidly evolving family where automated ortholog assignment is uncertain. Check Ensembl directly for that gene.

**Q: The disease list seems short or empty for a gene I know is disease-relevant.**  
A: Open Targets requires substantial curated evidence to include a gene-disease association. Newly published findings may not yet be integrated. Also check that the human ortholog was correctly resolved — a wrong human gene ID would return wrong or empty disease data.

**Q: Can I use this for human genes directly?**  
A: Yes. Select "Human (Homo sapiens)" as the source species and enter human gene symbols or ENSG IDs. The app will query Open Targets directly for those genes without needing ortholog mapping.

**Q: How do I interpret a gene with 200+ disease associations?**  
A: Genes like TP53, EGFR, and KRAS accumulate large numbers of disease associations because they are central to many pathways and have been extensively studied in cancer. High disease counts do not make them better targets — they reflect research bias as much as biology. Focus on the top association scores and the specific disease you are interested in.

**Q: Is this app suitable for clinical use?**  
A: No. This is a research tool for hypothesis generation and target prioritisation. It is not validated for clinical decision-making, diagnostic use, or patient management. All findings should be independently validated through appropriate experimental and regulatory processes.

---

*For questions, feedback, or collaboration enquiries, visit [geneticcodon.com](https://geneticcodon.com).*
