"""
ParylationPredictor v3.1 — PARylation Site Prediction Tool

Identifies PAR-binding domains and predicts PARylation sites from
UniProt IDs, raw sequences, or FASTA files.
"""

from .detector import PARBMDetector, PARylationSitePredictor, DomainHit, PARylationSite
from .input_handler import InputHandler
from .enrichment import EnrichmentEngine, DisorderPredictor

__version__ = "3.1.0"
__all__ = [
    "PARBMDetector",
    "PARylationSitePredictor",
    "InputHandler",
    "DomainHit",
    "PARylationSite",
    "EnrichmentEngine",
    "DisorderPredictor",
]
