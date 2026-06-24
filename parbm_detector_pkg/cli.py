"""
Command-line interface for PARBM Detector.

Usage:
  parbm-detect Q13547
  parbm-detect --visualize Q13547
  parbm-detect --batch Q13547 P38398 --csv out.csv
"""

import argparse
import sys
from .detector import PARBMDetector


def main():
    parser = argparse.ArgumentParser(
        prog="parbm-detect",
        description="PAR-binding Motif Detector – find PAR-binding domains and PARylation sites",
    )
    parser.add_argument("inputs", nargs="+",
                        help="UniProt IDs, sequences, or FASTA file paths")
    parser.add_argument("--visualize", "-v", action="store_true",
                        help="Generate visualisation PNG for each result")
    parser.add_argument("--compare", action="store_true",
                        help="Generate a side-by-side comparison figure (batch mode)")
    parser.add_argument("--json", metavar="FILE", default=None,
                        help="Export full results to JSON")
    parser.add_argument("--csv", metavar="FILE", default=None,
                        help="Export site table to CSV")
    parser.add_argument("--out-dir", metavar="DIR", default=".",
                        help="Output directory for PNG files (default: current dir)")
    args = parser.parse_args()

    detector = PARBMDetector()
    results  = detector.batch_analyze(args.inputs) if len(args.inputs) > 1 else []

    if len(args.inputs) == 1:
        result = detector.analyze(args.inputs[0])
        results = result if isinstance(result, list) else [result]

    for r in results:
        detector.print_analysis(r)
        if args.visualize and "error" not in r:
            import os
            safe = r["protein_name"][:30].replace(" ", "_")
            out  = os.path.join(args.out_dir, f"parbm_{safe}.png")
            detector.visualize_protein(r, output_file=out)

    if args.compare and len(results) > 1:
        import os
        detector.visualize_comparison(results,
            output_file=os.path.join(args.out_dir, "parbm_comparison.png"))

    if args.json:
        detector.export_results(results, args.json)
    if args.csv:
        detector.generate_csv_report(results, args.csv)


if __name__ == "__main__":
    main()
