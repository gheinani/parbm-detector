# ParylationPredictor

**ParylationPredictor** is a Python package for predicting PARylation (poly-ADP-ribosylation) sites in proteins and detecting PAR-binding domains. It accepts UniProt IDs, raw amino-acid sequences, or FASTA files as input and produces publication-quality figures, formatted data exports, and an interactive HTML report — all in a single command.

> Developed by the [Hashemi Gheinani Lab](https://www.dbmr.unibe.ch/research/research_programs/hashemi_gheinani_lab/index_eng.html), Department of BioMedical Research, University of Bern.

---

## What does it do?

| Capability | Details |
|---|---|
| **PAR-binding domain detection** | Queries the InterPro API for 11 domain families: WWE, Macro domain, MacroD-type, BRCT, PBZ zinc finger, RRM, CCCH zinc finger, PARP catalytic & regulatory |
| **PARylation site prediction** | Scores every residue using SxxE, SxxD, TxxE, ExxE, and acidic-cluster motifs with position-weighted confidence |
| **Disorder prediction** | Per-residue IUPred2A scores to identify modification-accessible flexible regions |
| **Cross-database enrichment** | Known PTMs from PHOSPHO.ELM · pathways from Reactome & KEGG · literature from EuropePMC |
| **Publication figures** | 8-panel composite PNG/PDF/SVG + 8 individual panel files |
| **Interactive HTML report** | 9-tab Plotly + Bootstrap 5 self-contained report |
| **Structured data export** | Formatted Excel workbook (4 sheets) · JSON · TSV |

---

## Requirements

- **Python** ≥ 3.8
- **pip** ≥ 21

Python dependencies (installed automatically):

| Package | Version | Purpose |
|---|---|---|
| `requests` | ≥ 2.28 | UniProt, InterPro, IUPred2A, PHOSPHO.ELM, EuropePMC, Reactome, KEGG API calls |
| `matplotlib` | ≥ 3.5 | 8-panel publication figures |
| `numpy` | ≥ 1.21 | Sliding-window profiles and scoring |
| `openpyxl` | ≥ 3.0 | Formatted Excel workbook export |
| `plotly` | ≥ 5.0 | Interactive HTML report |

An active internet connection is required for UniProt, InterPro, and enrichment API calls.

---

## Installation

### 1 — Clone the repository

```bash
git clone https://github.com/gheinani/parbm-detector.git
cd parbm-detector
```

### 2 — Install in editable mode

```bash
pip install -e .
```

This installs the package and all dependencies automatically.

### 3 — Verify

```python
import parbm_detector_pkg as pp
print(pp.__version__)   # → 3.1.0
```

---

## Quick start

### Analyze a protein by UniProt ID

```python
from parbm_detector_pkg import PARBMDetector

detector = PARBMDetector()

# Analyze RNF146 — a PAR-binding E3 ubiquitin ligase
result = detector.analyze("Q9NTX7")

# Print formatted text report
detector.print_analysis(result)
```

**Example output:**
```
────────────────────────────────────────────────────────
           ParylationPredictor — ANALYSIS REPORT
────────────────────────────────────────────────────────
Protein                          E3 ubiquitin-protein ligase RNF146
Gene                             RNF146
UniProt ID                       Q9NTX7
Organism                         Homo sapiens (Human)
Amino acids                      359

PAR-binding domains              WWE domain, WWE domain subgroup
High-conf sites (≥0.70)          37 sites
Top site                         E286 via SxxE motif, score 0.95
────────────────────────────────────────────────────────
```

---

### Full pipeline — enrich + export everything

```python
from parbm_detector_pkg import PARBMDetector

detector = PARBMDetector()

# Step 1: Analyze
result = detector.analyze("Q9NTX7")

# Step 2: Enrich with disorder, known PTMs, pathways, literature
result = detector.enrich_result(result)

# Step 3: Export all outputs to a named folder
folder = detector.export_to_folder(result, base_dir=".")
# → Creates ./Q9NTX7_RNF146/ with all files listed below
```

---

### Other input types

```python
# Raw amino-acid sequence
result = detector.analyze("MSEQAAKG...LNRK")

# FASTA string
result = detector.analyze(">MyProtein\nMSEQAAKG...LNRK")

# FASTA file — returns list of results
results = detector.analyze_fasta_file("proteins.fasta")
for r in results:
    detector.export_to_folder(r, base_dir="output/")
```

---

## Output files

All files are written to a folder named `{UniProtID}_{gene}/`:

```
Q9NTX7_RNF146/
├── Q9NTX7_RNF146_analysis.png      ← 8-panel composite figure (200 dpi)
├── Q9NTX7_RNF146_analysis.pdf      ← Vector PDF for journal submission
├── Q9NTX7_RNF146_analysis.svg      ← Scalable SVG
├── Q9NTX7_RNF146_report.xlsx       ← Formatted Excel (4 sheets)
├── Q9NTX7_RNF146_report.html       ← Interactive 9-tab HTML report
├── Q9NTX7_RNF146_data.json         ← Structured JSON (schema v3.0)
├── Q9NTX7_RNF146_sites.tsv         ← Site table with disorder & PTM columns
└── panels/                         ← 8 panels × PNG + PDF + SVG = 24 files
    ├── *_A_domain_architecture.png
    ├── *_B_confidence_landscape.png
    ├── *_C_hydrophobicity.png
    ├── *_D_charge_profile.png
    ├── *_E_motif_breakdown.png
    ├── *_F_residue_composition.png
    ├── *_G_disorder_profile.png
    └── *_H_residue_type_breakdown.png
```

---

## Figure panels

| Panel | Content |
|---|---|
| **A** | Domain architecture — backbone with PAR-binding domain blocks and PARylation site ticks |
| **B** | Confidence landscape — sliding-window PARylation score with threshold bands |
| **C** | Hydrophobicity profile — Kyte–Doolittle, sliding window |
| **D** | Charge profile — local net charge, highlighting acidic clusters |
| **E** | Motif breakdown — site counts per motif type (SxxE, SxxD, TxxE, ExxE, acidic) |
| **F** | Residue composition — amino-acid pie chart + summary table |
| **G** | Intrinsic disorder — IUPred2A per-residue scores |
| **H** | Residue-type breakdown — all vs. high-confidence sites by amino acid (S/T/E/D/Y/K/R) |

---

## API reference

```python
detector = PARBMDetector()

detector.analyze(input_data)                      # → result dict
detector.enrich_result(result, skip=[])           # → enriched result dict
detector.print_analysis(result)                   # formatted text to stdout
detector.visualize_protein(result, output_file)   # → PNG + PDF + SVG
detector.export_to_folder(result, base_dir=".")   # → folder path (all files)
detector.export_excel(result, output_file)        # → .xlsx path
detector.export_html(result, output_file)         # → .html path
detector.export_json_structured(result, file)     # → .json path
detector.export_tsv(result, output_file)          # → .tsv path
```

**`enrich_result` skip options:**

```python
# Skip specific modules (useful for offline use or speed)
result = detector.enrich_result(result, skip=["kegg", "literature"])
# Available: "disorder", "ptm", "pathways", "literature"
```

---

## Running the tests

```bash
pytest tests/ -v
```

---

## Contact

**Dr. Ali Hashemi Gheinani**  
Group Leader  
Department of BioMedical Research, University of Bern  
✉️ [ali.hashemi@unibe.ch](mailto:ali.hashemi@unibe.ch)  
🔬 [Hashemi Gheinani Lab](https://www.dbmr.unibe.ch/research/research_programs/hashemi_gheinani_lab/index_eng.html)

For bug reports and feature requests, please use the [GitHub Issues](https://github.com/gheinani/parbm-detector/issues) tracker.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
