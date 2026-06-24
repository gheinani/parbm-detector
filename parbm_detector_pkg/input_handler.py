"""
Input handling and detection for multiple file formats
"""

import re
from typing import List, Tuple


class InputHandler:
    """Intelligently handles different input types."""

    @staticmethod
    def detect_input_type(input_str: str) -> str:
        """Detect input type: uniprot_id, sequence, fasta, or fasta_file"""
        input_str = input_str.strip()

        if (
            input_str.startswith("/")
            or input_str.startswith(".")
            or input_str.endswith(".fasta")
            or input_str.endswith(".fa")
        ):
            return "fasta_file"

        if input_str.startswith(">"):
            return "fasta"

        if re.match(r"^[A-Z0-9]{6,10}$", input_str):
            return "uniprot_id"

        seq = input_str.upper().replace("-", "").replace("*", "").replace(" ", "")
        if len(seq) > 20:
            valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
            valid_chars = sum(1 for c in seq if c in valid_aa)
            if valid_chars / len(seq) > 0.95:
                return "sequence"

        if len(input_str) > 30:
            return "sequence"
        else:
            return "uniprot_id"

    @staticmethod
    def parse_fasta(fasta_text: str) -> List[Tuple[str, str]]:
        """Parse FASTA format text. Returns list of (header, sequence) tuples."""
        entries = []
        current_header = None
        current_seq = []

        for line in fasta_text.split("\n"):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header and current_seq:
                    entries.append((current_header, "".join(current_seq)))
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

        if current_header and current_seq:
            entries.append((current_header, "".join(current_seq)))

        return entries

    @staticmethod
    def read_fasta_file(filepath: str) -> List[Tuple[str, str]]:
        """Read FASTA file from disk."""
        try:
            with open(filepath, "r") as f:
                content = f.read()
            return InputHandler.parse_fasta(content)
        except FileNotFoundError:
            raise FileNotFoundError(f"FASTA file not found: {filepath}")
        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")
