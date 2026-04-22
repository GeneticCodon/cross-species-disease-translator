# 🧬 Cross-Species Disease Model Translator

> Map model organism genes to human orthologs and instantly surface disease associations and druggability insights — powered by Ensembl REST API & Open Targets.

![Python](https://img.shields.io/badge/Python-3.10+-blue?logo=python)
![Streamlit](https://img.shields.io/badge/Streamlit-1.37+-red?logo=streamlit)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Research%20MVP-orange)

---

## 🌍 Who Is This For?

| User | Use Case |
|------|----------|
| **Wet lab biologists** | Translate mouse/rat/zebrafish gene hits from experiments to human relevance |
| **Drug target scientists** | Quickly assess if a model organism finding is druggable in humans |
| **Bioinformaticians** | Batch-process gene lists from RNA-seq, GWAS, or proteomics studies |
| **Academic researchers** | Cross-reference model organism data with human disease databases |
| **Biotech/pharma analysts** | Prioritise candidates based on orthology confidence and druggability |
| **Students & educators** | Explore how model organisms relate to human biology interactively |

---

## ✨ Features

### 🔬 Core Biology
- **6 supported species**: Mouse, Rat, Zebrafish, Fruit fly, *C. elegans*, Human
- **Ensembl REST API** for ortholog mapping (one-to-one preferred, fallback to one-to-many)
- **Open Targets GraphQL** for disease associations & tractability scores
- Accepts gene **symbols** (e.g. `Akt1`) or **Ensembl IDs** (e.g. `ENSMUSG00000001`)
- Auto-detects species from Ensembl ID prefix

### 📊 Interactive Visualizations (Plotly)
- **Funnel chart** — gene mapping success pipeline
- **Pie chart** — orthology type distribution
- **Bar charts** — druggability & top diseases
- **Histogram** — disease count distribution per gene
- **Scatter plot** — disease count vs. average association score
- **Species detection chart** — shows breakdown of input gene origins

### 🧠 Gene Explorer
- Filter by orthology type, druggability, or minimum disease count
- Expandable cards per gene with full disease + score details

### 📁 Input / Output
- Upload CSV, TSV, TXT, Excel files
- Paste gene lists (comma or newline separated)
- Download full results or summary as CSV

---

## 🚀 Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/YOUR_USERNAME/cross-species-translator.git
cd cross-species-translator

# 2. Create a virtual environment
python -m venv .venv
source .venv/bin/activate       # Linux/macOS
# .venv\Scripts\activate        # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the app
streamlit run app.py
```

Then open **http://localhost:8501** in your browser.

---

## 📦 Dependencies

| Package | Purpose |
|---------|---------|
| `streamlit` | Web app framework |
| `pandas` | Data manipulation |
| `requests` | API calls (Ensembl, Open Targets) |
| `plotly` | Interactive charts |
| `matplotlib` / `seaborn` | Supporting visualization |
| `numpy` | Numerical computation |
| `openpyxl` | Excel file support |

---

## 🌐 External APIs Used

| API | Purpose | Docs |
|-----|---------|------|
| [Ensembl REST API](https://rest.ensembl.org) | Gene symbol/ID lookup, homology mapping | [Docs](https://rest.ensembl.org/documentation) |
| [Open Targets Platform](https://platform.opentargets.org) | Disease associations, tractability | [Docs](https://platform-docs.opentargets.org) |

> **No API key required** — both APIs are free and open access.

---

## 📂 Project Structure

```
cross_species_translator_mvp_v5/
├── app.py              # Main Streamlit application
├── requirements.txt    # Python dependencies
├── sample_data/        # Example gene list files
└── README.md           # This file
```

---

## 🖥️ Demo

Enable **"Use demo mouse genes"** in the sidebar to instantly try the app with 10 well-known mouse genes:

```
Akt1, Braf, Gfap, Mapt, Snca, App, Tp53, Pten, Egfr, Kras
```

---

## ⚠️ Limitations & Disclaimer

- **Research use only** — not for clinical diagnosis or treatment decisions
- Results depend on Ensembl & Open Targets data freshness (cached 24h)
- Maximum **200 genes** per run (API rate limits)
- Complex gene families (paralogs) may return unexpected orthologs

---

## 🗺️ Roadmap

- [ ] Gene-disease network graph (NetworkX / Cytoscape)
- [ ] PDF export of full analysis report
- [ ] Support for additional species (pig, marmoset, dog)
- [ ] OMIM / ClinVar integration
- [ ] Pathway enrichment analysis (KEGG / Reactome)
- [ ] Docker image for easy deployment
- [ ] Streamlit Community Cloud one-click deploy

---

## 👩‍🔬 Built By

**Genetic Codon** — bioinformatics tools for researchers.

---

## 📄 License

MIT — free to use, modify and share. Please cite if used in publications.
```