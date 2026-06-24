# ParylationPredictor

**ParylationPredictor** predicts PARylation (poly-ADP-ribosylation) sites in proteins and detects PAR-binding domains. You give it a UniProt ID — it gives you publication-quality figures, a formatted Excel file, an interactive HTML report, and structured data exports, all in one command.

> Developed by the [Hashemi Gheinani Lab](https://www.dbmr.unibe.ch/research/research_programs/hashemi_gheinani_lab/index_eng.html), Department of BioMedical Research, University of Bern.

---

## The quickest way to try it — no installation needed

Go to **[gheinani.github.io/parbm-detector/tool.html](https://gheinani.github.io/parbm-detector/tool.html)**, type your UniProt ID (for example `Q9NTX7`), and click **Analyze**. All 8 analysis panels appear in the browser and you can download a ZIP with the full results. No Python, no installation, no command line.

**Use the Python package when you need:**
- Batch analysis of many proteins at once
- Full Excel workbook output with formatted sheets
- The complete 8-panel composite figure as PNG/PDF/SVG for a paper

---

## Installing and running the Python package

### What you need before starting

- A Mac or Windows PC
- Python 3.8 or newer already installed (see below if you are not sure)
- An internet connection

You do **not** need to know how to program. You only need to type a few commands exactly as written.

---

### Step 0 — Check that Python is installed

**On Mac:** open the **Terminal** application.
- Press `Command + Space`, type `Terminal`, press Enter.

**On Windows:** open **Command Prompt**.
- Press the Windows key, type `cmd`, press Enter.

In the window that opens, type the following and press Enter:

```
python --version
```

You should see something like `Python 3.11.2`. If you see a version number starting with 3, you are ready. If you see an error or a version starting with 2, visit [python.org/downloads](https://www.python.org/downloads/) and install the latest version before continuing.

---

### Step 1 — Install the package

In the Terminal / Command Prompt window, type the following and press Enter:

```
pip install parylation-predictor
```

This downloads and installs the package and all its dependencies automatically (requests, matplotlib, numpy, openpyxl, plotly). It takes about 1–2 minutes. When the prompt reappears with no error message, installation is complete.

> **If you see "pip: command not found"**, try `pip3 install parylation-predictor` instead.

You do not need to download anything from GitHub manually. The package is published on PyPI at [pypi.org/project/parylation-predictor](https://pypi.org/project/parylation-predictor/).

---

### Step 2 — Run your first analysis

You will write a short script — a plain text file with the instructions for the tool. Do not worry, you only need to change one line.

1. Open **Notepad** (Windows) or **TextEdit** (Mac) — any plain text editor
   - On Mac: if TextEdit opens in rich-text mode, go to Format > Make Plain Text
2. Copy and paste the following text exactly:

```python
from parbm_detector_pkg import PARBMDetector

detector = PARBMDetector()
result = detector.analyze("Q9NTX7")        # <-- replace Q9NTX7 with your UniProt ID
result = detector.enrich_result(result)
folder = detector.export_to_folder(result, base_dir=".")
print("Done! Your results are in:", folder)
```

3. Replace `Q9NTX7` with your UniProt accession ID (for example `P18887`)
4. Save the file as `run_analysis.py` on your Desktop or anywhere you like
   - On Windows: in the Save dialog, change "Save as type" to "All Files" and type the filename as `run_analysis.py`
   - On Mac: save as `run_analysis.py` and make sure TextEdit does not add `.txt`

Now go back to the Terminal / Command Prompt and type:

```
python run_analysis.py
```

Press Enter. The tool will connect to UniProt, fetch your protein, run the analysis, and print a progress log. When it finishes, it prints:

```
Done! Your results are in: ./Q9NTX7_RNF146
```

> **If you see "python: command not found"**, try `python3 run_analysis.py` instead.

---

### Step 3 — Find your results

Open the folder where you saved `run_analysis.py` in Finder (Mac) or File Explorer (Windows). You will see a new folder named after your protein, for example `Q9NTX7_RNF146`. Open it — it contains:

```
Q9NTX7_RNF146/
├── Q9NTX7_RNF146_analysis.png      <- 8-panel figure (open with any image viewer)
├── Q9NTX7_RNF146_analysis.pdf      <- same figure as PDF, ready for a journal
├── Q9NTX7_RNF146_analysis.svg      <- scalable vector version
├── Q9NTX7_RNF146_report.xlsx       <- Excel workbook with 4 data sheets
├── Q9NTX7_RNF146_report.html       <- interactive report (open in any browser)
├── Q9NTX7_RNF146_data.json         <- raw data in JSON format
├── Q9NTX7_RNF146_sites.tsv         <- table of predicted sites (open in Excel)
└── panels/                          <- each of the 8 panels as separate files
    ├── *_A_domain_architecture.png
    ├── *_B_confidence_landscape.png
    ├── *_C_hydrophobicity.png
    ├── *_D_charge_profile.png
    ├── *_E_motif_breakdown.png
    ├── *_F_residue_composition.png
    ├── *_G_disorder_profile.png
    └── *_H_residue_type_breakdown.png
```

Double-click `Q9NTX7_RNF146_report.html` to open the interactive report in your browser. Double-click the `.xlsx` file to open it in Excel. The `.png` and `.pdf` files open in any image or PDF viewer.

---

### Analyzing multiple proteins

To analyze several proteins at once, replace the analysis lines in `run_analysis.py` with a list:

```python
from parbm_detector_pkg import PARBMDetector

detector = PARBMDetector()

proteins = ["Q9NTX7", "P18887", "Q9Y6K9"]   # add as many UniProt IDs as you like

for uid in proteins:
    result = detector.analyze(uid)
    result = detector.enrich_result(result)
    detector.export_to_folder(result, base_dir=".")
    print("Finished:", uid)
```

Each protein gets its own output folder.

---

## What does the tool analyse?

| Capability | Details |
|---|---|
| PAR-binding domain detection | Queries InterPro for 11 domain families: WWE, Macro, MacroD-type, BRCT, PBZ, RRM, CCCH, PARP catalytic and regulatory |
| PARylation site prediction | Scores every residue for SxxE, SxxD, TxxE, ExxE, and acidic-cluster motifs |
| Disorder prediction | Per-residue IUPred2A scores to find modification-accessible flexible regions |
| Cross-database enrichment | Known PTMs from PHOSPHO.ELM, pathways from Reactome and KEGG, literature from EuropePMC |
| Publication figures | 8-panel composite PNG/PDF/SVG plus 8 individual panel files |
| Interactive report | Self-contained HTML report viewable in any browser |
| Structured exports | Formatted Excel (4 sheets), JSON, and TSV |

---

## Figure panels explained

| Panel | What it shows | How to interpret it |
|---|---|---|
| A | Domain architecture | Protein backbone with PAR-binding domain blocks. Sites inside domain boxes are highest-priority candidates |
| B | Confidence landscape | Sliding-window score across the sequence. Broad peaks indicate modification hotspots |
| C | Hydrophobicity | Kyte-Doolittle profile. Sites in hydrophilic dips are most solvent-accessible |
| D | Charge profile | Local net charge. Deep acidic (negative) regions are primary modification zones |
| E | Motif breakdown | Counts per motif type (SxxE, SxxD, TxxE, ExxE, acidic cluster) |
| F | Residue composition | Amino-acid pie chart with key metrics summary |
| G | Disorder profile | IUPred2A scores. Sites in disordered regions (score > 0.5) are most accessible |
| H | Residue-type breakdown | All predicted vs. high-confidence sites by amino acid (S, T, E, D, Y, K, R) |

---

## Troubleshooting

**"pip is not recognized" or "python is not recognized"**
Python may not be added to your system PATH. On Windows, reinstall Python from [python.org](https://www.python.org/downloads/) and make sure to tick the box **"Add Python to PATH"** during installation.

**"No module named parbm_detector_pkg"**
You are running `python` from the wrong folder. Make sure your Terminal is inside the `parbm-detector-main` folder (Step 2) before running the install and analysis commands.

**The script runs but I get a network error**
The tool needs internet access to contact UniProt, InterPro, and other databases. Check your connection and try again. Corporate or university firewalls occasionally block API calls — try from a different network if the problem persists.

**The analysis takes a long time**
Each protein analysis makes several API requests. 30–60 seconds per protein is normal. Enrichment (`enrich_result`) adds another 30–60 seconds because it queries EuropePMC, Reactome, and KEGG.

---

## Requirements

- Python 3.8 or newer
- pip (comes with Python)

Python packages installed automatically by `pip install .`:

| Package | Purpose |
|---|---|
| requests | API calls to UniProt, InterPro, IUPred2A, PHOSPHO.ELM, EuropePMC, Reactome, KEGG |
| matplotlib | 8-panel publication figures |
| numpy | Sliding-window scoring |
| openpyxl | Formatted Excel workbook |
| plotly | Interactive HTML report |

---

## Contact

**Dr. Ali Hashemi Gheinani**
Group Leader
Department of BioMedical Research, University of Bern
Email: [ali.hashemi@unibe.ch](mailto:ali.hashemi@unibe.ch)
Lab: [Hashemi Gheinani Lab](https://www.dbmr.unibe.ch/research/research_programs/hashemi_gheinani_lab/index_eng.html)

For bug reports and feature requests, please use the [GitHub Issues](https://github.com/gheinani/parbm-detector/issues) tracker.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
