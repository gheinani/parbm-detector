"""
Tests for PARBM Detector
"""

import pytest
from parbm_detector_pkg import PARBMDetector, InputHandler, PARylationSitePredictor


TEST_SEQ = "MKVLRVDVDMALDALDSLDESXXEMTXXE" * 3   # contains motifs


class TestInputHandler:
    def test_uniprot_id(self):
        assert InputHandler.detect_input_type("Q13547") == "uniprot_id"

    def test_sequence(self):
        assert InputHandler.detect_input_type("MKVLRVDVDMALDALDALD" * 2) == "sequence"

    def test_fasta_text(self):
        assert InputHandler.detect_input_type(">header\nMKVL") == "fasta"

    def test_fasta_file(self):
        assert InputHandler.detect_input_type("proteins.fasta") == "fasta_file"

    def test_parse_fasta(self):
        fasta = ">prot1\nMKVL\nVDVD\n>prot2\nACDE"
        entries = InputHandler.parse_fasta(fasta)
        assert len(entries) == 2
        assert entries[0] == ("prot1", "MKVLVDVD")
        assert entries[1] == ("prot2", "ACDE")


class TestPredictor:
    def test_predict_returns_list(self):
        sites = PARylationSitePredictor.predict_sites(TEST_SEQ)
        assert isinstance(sites, list)

    def test_scores_in_range(self):
        sites = PARylationSitePredictor.predict_sites(TEST_SEQ)
        for s in sites:
            assert 0.0 <= s.score <= 1.0

    def test_sorted_descending(self):
        sites = PARylationSitePredictor.predict_sites(TEST_SEQ)
        scores = [s.score for s in sites]
        assert scores == sorted(scores, reverse=True)

    def test_local_charge(self):
        c = PARylationSitePredictor.calculate_local_charge("KDEKDE", 3, 2)
        assert isinstance(c, float)

    def test_local_hydrophobicity(self):
        h = PARylationSitePredictor.calculate_local_hydrophobicity("ALILAV", 2, 2)
        assert isinstance(h, float)


class TestDetector:
    def test_init(self):
        assert PARBMDetector() is not None

    def test_analyze_sequence(self):
        d = PARBMDetector()
        r = d.analyze("MKVLRVDVDMALDALDSLDE" * 3)
        assert "error" not in r
        assert r["input_type"] == "sequence"
        assert r["sequence_length"] > 0

    def test_analyze_fasta_text(self):
        d = PARBMDetector()
        fasta = ">prot_A\nMKVLSEDESXXEMTXXE\n>prot_B\nACDEFGHIKLMNP"
        rs = d.analyze(fasta)
        assert isinstance(rs, list)
        assert len(rs) == 2

    def test_batch_analyze(self):
        d  = PARBMDetector()
        rs = d.batch_analyze(["MKVLRVDVDMALDALDS" * 2, "ACDEFGHIKLMNPQRST" * 2])
        assert len(rs) == 2

    def test_export_csv(self, tmp_path):
        d  = PARBMDetector()
        r  = d.analyze("MKVLRVDVDMALDALDSLDESXXE" * 3)
        out = str(tmp_path / "out.csv")
        d.generate_csv_report([r], out)
        import os
        assert os.path.exists(out)

    def test_visualize_no_crash(self, tmp_path):
        d   = PARBMDetector()
        r   = d.analyze("MKVLRVDVDMALDALDSLDESXXE" * 5)
        out = str(tmp_path / "fig.png")
        path = d.visualize_protein(r, output_file=out)
        import os
        assert path is not None
        assert os.path.exists(out)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
