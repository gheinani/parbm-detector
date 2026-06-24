"""
External data enrichment for PARBM Detector.

Modules:
  DisorderPredictor  – IUPred2A REST API (fallback: AA-propensity)
  PTMDatabaseClient  – PHOSPHO.ELM experimentally validated sites
  LiteratureSearcher – EuropePMC full-text search
  PathwayClient      – Reactome + KEGG pathways
  EnrichmentEngine   – orchestrates all modules, caches per UniProt ID

All methods fail gracefully – empty / None when API is unreachable.
"""

import re
import time
import json
import requests
from typing import Dict, List, Optional


# ── Amino-acid disorder propensity (fallback when IUPred2A is down) ───────────

_DISORDER_SCORE = {
    "A": 0.35, "R": 0.70, "N": 0.55, "D": 0.65, "C": 0.20,
    "Q": 0.60, "E": 0.75, "G": 0.45, "H": 0.40, "I": 0.15,
    "L": 0.15, "K": 0.75, "M": 0.25, "F": 0.15, "P": 0.65,
    "S": 0.60, "T": 0.50, "W": 0.15, "Y": 0.25, "V": 0.15,
}

def _smooth(values: List[float], window: int = 13) -> List[float]:
    n = len(values)
    out = []
    for i in range(n):
        lo, hi = max(0, i - window // 2), min(n, i + window // 2 + 1)
        out.append(sum(values[lo:hi]) / (hi - lo))
    return out


# ── Disorder predictor ────────────────────────────────────────────────────────

class DisorderPredictor:
    """Predict per-residue disorder using IUPred2A; AA-propensity as fallback."""

    IUPRED_BASE = "https://iupred2a.elte.hu"

    def __init__(self, session: requests.Session = None):
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": "PARBMDetector/3.0"})

    def predict(self, sequence: str, uniprot_id: str = None) -> Dict:
        scores = self._iupred_api(sequence, uniprot_id)
        method = "IUPred2A"
        if not scores or len(scores) != len(sequence):
            scores = _smooth([_DISORDER_SCORE.get(aa, 0.45) for aa in sequence])
            method = "AA-propensity (fallback)"

        mean = sum(scores) / len(scores) if scores else 0.0
        frac = sum(1 for s in scores if s > 0.5) / len(scores) if scores else 0.0
        return {
            "scores": scores,
            "mean_disorder": round(mean, 3),
            "disordered_fraction": round(frac, 3),
            "method": method,
        }

    def _iupred_api(self, sequence: str, uniprot_id: str = None) -> List[float]:
        # Try by UniProt ID first (faster server cache)
        if uniprot_id:
            try:
                r = self.session.get(
                    f"{self.IUPRED_BASE}/iupred2a/long/{uniprot_id}",
                    timeout=12,
                )
                if r.status_code == 200:
                    scores = self._parse(r.text)
                    if scores:
                        return scores
            except Exception:
                pass

        # Fall back to sequence POST
        try:
            r = self.session.post(
                f"{self.IUPRED_BASE}/iupred2a/long",
                data={"sequence": sequence},
                timeout=20,
            )
            if r.status_code == 200:
                return self._parse(r.text)
        except Exception:
            pass
        return []

    @staticmethod
    def _parse(text: str) -> List[float]:
        scores = []
        for line in text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    scores.append(float(parts[2]))
                except ValueError:
                    continue
        return scores


# ── Known PTM sites ───────────────────────────────────────────────────────────

class PTMDatabaseClient:
    """Fetch experimentally validated PTM sites from PHOSPHO.ELM."""

    BASE = "http://phospho.elm.eu.org/api"

    def __init__(self, session: requests.Session = None):
        self.session = session or requests.Session()
        self.session.headers.update({
            "User-Agent": "PARBMDetector/3.0",
            "Accept": "application/json",
        })

    def get_known_sites(self, uniprot_id: str) -> List[Dict]:
        try:
            r = self.session.get(
                f"{self.BASE}/byprotein/{uniprot_id}",
                timeout=12,
            )
            if r.status_code == 200:
                raw = r.json()
                sites = []
                for entry in raw if isinstance(raw, list) else raw.get("objects", []):
                    sites.append({
                        "position":     entry.get("position"),
                        "residue":      entry.get("code", ""),
                        "modification": entry.get("modification", ""),
                        "kinase":       entry.get("kinase", ""),
                        "source":       "PHOSPHO.ELM",
                        "pmid":         entry.get("pmid", ""),
                    })
                return sites
        except Exception:
            pass
        return []


# ── Literature search ─────────────────────────────────────────────────────────

class LiteratureSearcher:
    """Search EuropePMC for PARylation-relevant publications."""

    BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    def __init__(self, session: requests.Session = None):
        self.session = session or requests.Session()

    def search(self, gene_name: str, protein_name: str,
               max_results: int = 6) -> List[Dict]:
        try:
            terms = []
            if gene_name and gene_name not in ("N/A", "Unknown"):
                terms.append(f'"{gene_name}"')
            if protein_name and protein_name != "Unknown":
                terms.append(f'"{protein_name}"')
            q_protein = " OR ".join(terms) if terms else "PARylation"
            query = (
                f"({q_protein}) AND "
                f'(PARylation OR "ADP-ribosylation" OR "PAR-binding" OR PARP)'
            )
            params = {
                "query":      query,
                "format":     "json",
                "pageSize":   max_results,
                "sort":       "CITED desc",
                "resultType": "core",
            }
            r = self.session.get(self.BASE, params=params, timeout=12)
            r.raise_for_status()
            results = r.json().get("resultList", {}).get("result", [])
            return [
                {
                    "title":    pub.get("title", ""),
                    "authors":  pub.get("authorString", ""),
                    "journal":  pub.get("journalTitle", ""),
                    "year":     pub.get("pubYear", ""),
                    "pmid":     pub.get("pmid", ""),
                    "doi":      pub.get("doi", ""),
                    "cited_by": pub.get("citedByCount", 0),
                }
                for pub in results
            ]
        except Exception:
            return []


# ── Pathway annotation ────────────────────────────────────────────────────────

class PathwayClient:
    """Fetch biological pathway associations from Reactome and KEGG."""

    REACTOME = "https://reactome.org/ContentService"
    KEGG     = "https://rest.kegg.jp"

    def __init__(self, session: requests.Session = None):
        self.session = session or requests.Session()

    def get_reactome(self, uniprot_id: str) -> List[Dict]:
        try:
            url = (f"{self.REACTOME}/data/mapping/UniProt"
                   f"/{uniprot_id}/pathways")
            r = self.session.get(url, params={"species": "9606"}, timeout=12)
            if r.status_code == 200:
                return [
                    {"id": p.get("stId", ""),
                     "name": p.get("displayName", ""),
                     "source": "Reactome"}
                    for p in r.json()[:8]
                ]
        except Exception:
            pass
        return []

    def get_kegg(self, gene_name: str) -> List[Dict]:
        if not gene_name or gene_name in ("N/A", "Unknown"):
            return []
        try:
            r = self.session.get(
                f"{self.KEGG}/find/genes/hsa:{gene_name}", timeout=10
            )
            if r.status_code != 200 or not r.text.strip():
                return []
            gene_id = r.text.strip().split("\t")[0]

            r2 = self.session.get(
                f"{self.KEGG}/link/pathway/{gene_id}", timeout=10
            )
            if r2.status_code != 200:
                return []

            pathways = []
            for line in r2.text.strip().splitlines():
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                pw_id = parts[1].strip()
                r3 = self.session.get(f"{self.KEGG}/get/{pw_id}", timeout=8)
                pw_name = pw_id
                if r3.status_code == 200:
                    for ln in r3.text.splitlines()[:5]:
                        if ln.startswith("NAME"):
                            pw_name = ln.split(None, 1)[1].strip()
                            break
                pathways.append({"id": pw_id, "name": pw_name, "source": "KEGG"})
                if len(pathways) >= 5:
                    break
            return pathways
        except Exception:
            return []


# ── Enrichment engine ─────────────────────────────────────────────────────────

class EnrichmentEngine:
    """
    Orchestrates all enrichment modules.
    Results cached per UniProt ID within the session.
    """

    def __init__(self):
        self.session  = requests.Session()
        self.session.headers.update({"User-Agent": "PARBMDetector/3.0"})
        self.disorder   = DisorderPredictor(self.session)
        self.ptm        = PTMDatabaseClient(self.session)
        self.literature = LiteratureSearcher(self.session)
        self.pathways   = PathwayClient(self.session)
        self._cache: Dict[str, Dict] = {}

    def enrich(self, result: Dict, skip: List[str] = None) -> Dict:
        """
        Add an 'enrichment' block to a detector result dict.

        Parameters
        ----------
        result : dict  – output from PARBMDetector.analyze()
        skip   : list  – modules to skip: 'disorder', 'ptm', 'literature', 'pathways'

        Returns a shallow copy of result with 'enrichment' key added.
        """
        skip = set(skip or [])
        uid   = result.get("uniprot_id")
        seq   = result.get("sequence", "")
        gene  = result.get("gene_name", "")
        pname = result.get("protein_name", "")

        # serve from cache if same protein already enriched
        cache_key = uid or pname[:30]
        if cache_key in self._cache:
            enriched = dict(result)
            enriched["enrichment"] = self._cache[cache_key]
            return enriched

        enrichment: Dict = {}

        # ── Disorder ──────────────────────────────────────────────────────────
        if "disorder" not in skip and seq:
            print("    · Disorder prediction (IUPred2A) …", end=" ", flush=True)
            enrichment["disorder"] = self.disorder.predict(seq, uid)
            print(f"✓  ({enrichment['disorder']['method']})")

        # ── Known PTM sites ───────────────────────────────────────────────────
        if "ptm" not in skip and uid:
            print("    · Known PTM sites (PHOSPHO.ELM) …", end=" ", flush=True)
            enrichment["known_ptm_sites"] = self.ptm.get_known_sites(uid)
            print(f"✓  ({len(enrichment['known_ptm_sites'])} sites)")

        # ── Literature ────────────────────────────────────────────────────────
        if "literature" not in skip:
            print("    · Literature (EuropePMC) …", end=" ", flush=True)
            enrichment["literature"] = self.literature.search(gene, pname)
            print(f"✓  ({len(enrichment['literature'])} refs)")

        # ── Pathways ──────────────────────────────────────────────────────────
        if "pathways" not in skip:
            print("    · Pathways (Reactome + KEGG) …", end=" ", flush=True)
            pws: Dict[str, List] = {}
            if uid:
                pws["reactome"] = self.pathways.get_reactome(uid)
            if gene and gene not in ("N/A", "Unknown"):
                pws["kegg"] = self.pathways.get_kegg(gene)
            enrichment["pathways"] = pws
            total_pw = sum(len(v) for v in pws.values())
            print(f"✓  ({total_pw} pathways)")

        # ── Motif context ─────────────────────────────────────────────────────
        # Compute ±4 aa context around each high-confidence site
        if seq and result.get("high_confidence_sites"):
            contexts = []
            for site in result["high_confidence_sites"]:
                pos = site.position - 1
                lo  = max(0, pos - 4)
                hi  = min(len(seq), pos + 5)
                contexts.append({
                    "position": site.position,
                    "residue":  site.residue,
                    "context":  seq[lo:hi],
                    "score":    site.score,
                })
            enrichment["site_contexts"] = contexts

        self._cache[cache_key] = enrichment
        enriched = dict(result)
        enriched["enrichment"] = enrichment
        return enriched

    # ── Residue-type breakdown (no API needed) ────────────────────────────────

    @staticmethod
    def residue_type_stats(result: Dict) -> Dict[str, Dict]:
        """Count PARylation sites per residue type and confidence tier."""
        tiers = {"high": [], "medium": [], "low": []}
        by_aa: Dict[str, List] = {}

        for site in result.get("paryation_sites", []):
            aa = site.residue
            by_aa.setdefault(aa, []).append(site)
            if site.score >= 0.7:
                tiers["high"].append(site)
            elif site.score >= 0.5:
                tiers["medium"].append(site)
            else:
                tiers["low"].append(site)

        total = sum(len(v) for v in by_aa.values())
        residue_stats = {}
        for aa in ["S", "T", "E", "D", "Y", "K", "R"]:
            sites = by_aa.get(aa, [])
            residue_stats[aa] = {
                "count":   len(sites),
                "percent": round(100 * len(sites) / total, 1) if total else 0,
                "mean_score": (round(sum(s.score for s in sites) / len(sites), 3)
                               if sites else 0),
            }

        return {
            "by_residue": residue_stats,
            "tiers": {k: len(v) for k, v in tiers.items()},
            "total": total,
        }

    # ── Prediction quality metrics (against known PTM sites) ─────────────────

    @staticmethod
    def quality_metrics(result: Dict, window: int = 3) -> Dict:
        """
        Compare predicted sites to PHOSPHO.ELM known sites.
        Counts a predicted site as correct if it falls within ±window of a known site.
        Returns sensitivity, precision, F1, and counts.
        """
        enrichment   = result.get("enrichment", {})
        known_sites  = enrichment.get("known_ptm_sites", [])
        pred_sites   = result.get("paryation_sites", [])

        if not known_sites or not pred_sites:
            return {"available": False,
                    "note": "Requires known PTM sites from PHOSPHO.ELM"}

        known_pos = {s["position"] for s in known_sites if s.get("position")}
        pred_pos  = [s.position for s in pred_sites]

        # TP: predicted positions within ±window of a known site
        tp = sum(
            1 for p in pred_pos
            if any(abs(p - k) <= window for k in known_pos)
        )
        fp = len(pred_pos) - tp
        fn = len(known_pos) - sum(
            1 for k in known_pos
            if any(abs(p - k) <= window for p in pred_pos)
        )

        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        precision   = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1          = (2 * sensitivity * precision / (sensitivity + precision)
                       if (sensitivity + precision) > 0 else 0)

        return {
            "available":   True,
            "tp":          tp,
            "fp":          fp,
            "fn":          fn,
            "known_sites": len(known_pos),
            "pred_sites":  len(pred_pos),
            "sensitivity": round(sensitivity * 100, 1),
            "precision":   round(precision   * 100, 1),
            "f1_score":    round(f1, 3),
            "window_aa":   window,
            "note": f"Site counted correct if within ±{window} aa of a known site",
        }
