#!/usr/bin/env python3
"""
PARBM Detector v3.0 – Example usage
Run: python example_usage.py
"""

from parbm_detector_pkg import PARBMDetector

def main():
    detector = PARBMDetector()

    # ── 1. UniProt lookup (PARP1) ────────────────────────────────────────────
    print("=" * 70)
    print("1. Fetching PARP1 (P09874) from UniProt + InterPro …")
    result = detector.analyze("P09874")
    detector.print_analysis(result)
    detector.visualize_protein(result, output_file="parbm_PARP1.png")

    # ── 2. Raw sequence ──────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("2. Raw sequence analysis …")
    seq = (
        "MSVSVKKLQARVKKENGQSTSDLTLKESLEELKEKLKEVNEELKKELRELEKDLNQKLNELE"
        "KKLKEQNESLDSEDKSTSSDSLTSDSSVSSDSVNTDSDPESTSSDSLTSDSSVSSDSVNTDS"
    )
    result2 = detector.analyze(seq)
    detector.print_analysis(result2)
    detector.visualize_protein(result2, output_file="parbm_custom_seq.png")

    # ── 3. Batch + comparison figure ─────────────────────────────────────────
    print("\n" + "=" * 70)
    print("3. Batch analysis …")
    results = detector.batch_analyze(["P09874", "Q13547"])
    for r in results:
        detector.print_analysis(r)
    detector.visualize_comparison(results, output_file="parbm_comparison.png")
    detector.export_results(results, "parbm_results.json")
    detector.generate_csv_report(results, "parbm_sites.csv")

    print("\n✓ Done.  PNG files, JSON and CSV written to current directory.")

if __name__ == "__main__":
    main()
