#!/usr/bin/env python3
"""
PPI Cell-Type Expression Validator
====================================
Validates whether each protein-protein interaction (PPI) pair in BFS Beam Search pathway data
can be realized in a specific cell type.

Validation items:
  1. Expression status of each protein (HPA single-cell expression data)
  2. HPA Tau Score (cell-type specificity)
  3. Intracellular organelle co-localization (subcellular location compatibility)
  4. STRING PPI confidence score + evidence channels
  5. Overall pathway score (validation grade)

APIs used:
  - MyGene.info : Gene Symbol → Ensembl ID mapping
  - Human Protein Atlas REST API : expression / Tau / subcellular
  - STRING REST API v12.x : combined PPI confidence + channel evidence

Output: interactive HTML report
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import json
import re
import os
import time
import urllib.request
import urllib.parse
from datetime import datetime
from collections import defaultdict

# ─── Constants ────────────────────────────────────────────────────────────────────
APP_VERSION = "2.1.0"
STRING_API = "https://version-12-0.string-db.org/api/json"

# HPA API — correct URL and column structure
# rnascm = RNA single cell type specific nTPM (core column)
# rnasctm = RNA single cell type Tau score
HPA_BASE = "https://www.proteinatlas.org"
HPA_JSON_ENSG = "https://www.proteinatlas.org/{ensg}.json"
HPA_JSON_SYM = "https://www.proteinatlas.org/{symbol}.json"
HPA_SEARCH = ("https://www.proteinatlas.org/api/search_download.php"
              "?search={symbol}&format=json"
              "&columns=g,eg,rnascm,rnasctm,subcell_location"
              "&compress=no")

# MyGene.info (Gene Symbol → Ensembl ID)
MYGENE_API = "https://mygene.info/v3/query?q=symbol:{sym}&species=human&fields=ensembl.gene,name"

STRING_SCORE_THRESHOLD = 400  # medium confidence

# ─── Built-in ENSG mapping (major signaling proteins — fallback for network errors) ─────────
BUILTIN_ENSG = {
    "EGFR": "ENSG00000146648", "ERBB2": "ENSG00000141736",
    "ERBB3": "ENSG00000065361", "ERBB4": "ENSG00000178568",
    "SRC": "ENSG00000197122", "CAV1": "ENSG00000105974",
    "PIK3CA": "ENSG00000121879", "PIK3R1": "ENSG00000145675",
    "PIK3R2": "ENSG00000105647", "PIK3R3": "ENSG00000117461",
    "PIK3CB": "ENSG00000051382", "PIK3CD": "ENSG00000171608",
    "PIK3CG": "ENSG00000105851",
    "AKT1": "ENSG00000142208", "AKT2": "ENSG00000105221",
    "PTEN": "ENSG00000171862", "MTOR": "ENSG00000198793",
    "KRAS": "ENSG00000133703", "NRAS": "ENSG00000213281",
    "HRAS": "ENSG00000174775", "RAF1": "ENSG00000132155",
    "BRAF": "ENSG00000157764", "MAP2K1": "ENSG00000169032",
    "MAP2K2": "ENSG00000126934", "MAPK1": "ENSG00000100030",
    "MAPK3": "ENSG00000102882",
    "SHC1": "ENSG00000197634", "GRB2": "ENSG00000177885",
    "SOS1": "ENSG00000115085", "SOS2": "ENSG00000100485",
    "GAB1": "ENSG00000109458", "GAB2": "ENSG00000033327",
    "PTPN11": "ENSG00000179295", "CBL": "ENSG00000110395",
    "CRK": "ENSG00000167193", "CRKL": "ENSG00000099942",
    "IRS1": "ENSG00000169047", "IRS2": "ENSG00000185950",
    "NRG1": "ENSG00000157168", "ERBIN": "ENSG00000187210",
    "VAV1": "ENSG00000141968", "SYK": "ENSG00000165025",
    "LCP2": "ENSG00000043462", "ABL1": "ENSG00000097007",
    "BCAR1": "ENSG00000050820", "PDGFRB": "ENSG00000113721",
    "PLCG1": "ENSG00000124181", "RAPGEF1": "ENSG00000107263",
    "CAVIN1": "ENSG00000177469",
    # HeLa pathway
    "PXN": "ENSG00000089159", "ILK": "ENSG00000166033",
    "TLN1": "ENSG00000137076", "TLN2": "ENSG00000171914",
    "LIMS1": "ENSG00000169756", "ZYX": "ENSG00000159111",
    "PARVB": "ENSG00000168243", "PARVA": "ENSG00000137309",
    "HSP90AA1": "ENSG00000080824", "HSP90AB1": "ENSG00000096384",
    "STUB1": "ENSG00000103266", "STIP1": "ENSG00000135914",
    "BAG3": "ENSG00000151929", "HSPA4": "ENSG00000170606",
    "HSPA8": "ENSG00000109971", "ESR1": "ENSG00000091831",
    "CDC42": "ENSG00000070831", "IQGAP1": "ENSG00000140612",
    "IQGAP2": "ENSG00000109189", "RAP1A": "ENSG00000116473",
    # Thrombin/PAR1 pathway
    "F2R": "ENSG00000181104", "GNAQ": "ENSG00000156052",
    "GNA12": "ENSG00000196587", "GNA13": "ENSG00000064687",
    "PLCB3": "ENSG00000149782", "GRK2": "ENSG00000108953",
    "ARRB1": "ENSG00000137486", "ARRB2": "ENSG00000141480",
    "ADRB2": "ENSG00000169252", "RHOA": "ENSG00000067560",
    "ROCK2": "ENSG00000134318", "RTKN": "ENSG00000117399",
    "ARHGEF11": "ENSG00000107960", "WASL": "ENSG00000106299",
    "MAPK8": "ENSG00000107643", "MAP2K4": "ENSG00000065559",
    "HRH1": "ENSG00000196639", "ARHGDIA": "ENSG00000141522",
    "ANLN": "ENSG00000011426", "ARHGAP1": "ENSG00000112183",
    "RTKN": "ENSG00000117399", "CDC42": "ENSG00000070831",
    "SLC9A3R1": "ENSG00000109062", "RDX": "ENSG00000088247",
    "EZR": "ENSG00000092820", "CD44": "ENSG00000026508",
    "FN1": "ENSG00000115414", "ITGA5": "ENSG00000161638",
    "ITGB1": "ENSG00000150093",
}

# Major cell types (HPA reference)
CELL_TYPES = [
    "cardiomyocytes", "skeletal muscle cells", "smooth muscle cells",
    "hepatocytes", "cholangiocytes", "neurons", "astrocytes", "oligodendrocytes",
    "microglia", "endothelial cells", "fibroblasts", "adipocytes",
    "T cells", "B cells", "NK cells", "monocytes", "macrophages",
    "dendritic cells", "neutrophils", "erythroid cells",
    "keratinocytes", "melanocytes",
    "proximal tubular cells", "podocytes",
    "alveolar cells", "club cells",
    "intestinal epithelial cells", "goblet cells",
    "pancreatic acinar cells", "beta cells",
    "breast glandular cells", "cervical cells",
    "spermatocytes", "oocytes",
    "HeLa (cancer cell line)", "MDA-MB-468 (cancer cell line)",
]

# Subcellular location compatibility matrix (same group = interaction possible)
COMPARTMENT_GROUPS = {
    "cytosol": ["cytosol", "cytoplasm", "cytoskeleton"],
    "nucleus": ["nucleus", "nucleoplasm", "nucleoli", "nuclear membrane", "nuclear bodies", "nuclear speckles"],
    "membrane": ["plasma membrane", "cell membrane", "vesicles", "endosomes", "lysosomes"],
    "er": ["endoplasmic reticulum", "golgi apparatus"],
    "mitochondria": ["mitochondria", "mitochondrial matrix"],
    "secreted": ["extracellular", "secreted"],
}


def get_compartment_group(loc: str) -> str:
    loc_lower = loc.lower()
    for grp, terms in COMPARTMENT_GROUPS.items():
        if any(t in loc_lower for t in terms):
            return grp
    return "other"


# ─── API utilities ─────────────────────────────────────────────────────────────
def fetch_json(url: str, timeout: int = 15) -> dict | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PPIValidator/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8"))
    except Exception as e:
        return None


def gene_symbol_to_ensembl(symbol: str, cache: dict) -> str | None:
    """Convert Gene Symbol to Ensembl ID (3-step fallback)."""
    if symbol in cache:
        return cache[symbol]

    # Step 1: built-in mapping (instant, no network required)
    builtin = BUILTIN_ENSG.get(symbol.upper())
    if builtin:
        cache[symbol] = builtin
        return builtin

    # Step 2: MyGene.info API
    try:
        url = MYGENE_API.format(sym=urllib.parse.quote(symbol))
        data = fetch_json(url, timeout=10)
        if data and data.get("hits"):
            hit = data["hits"][0]
            ens = hit.get("ensembl", {})
            ensg = None
            if isinstance(ens, dict):
                ensg = ens.get("gene")
            elif isinstance(ens, list) and ens:
                ensg = ens[0].get("gene")
            if ensg:
                cache[symbol] = ensg
                return ensg
    except Exception:
        pass

    # Step 3: query HPA directly by gene symbol without ENSG ID (fallback marker)
    cache[symbol] = f"SYM:{symbol}"  # special marker: use symbol directly
    return f"SYM:{symbol}"


def ntpm_to_level(ntpm) -> str:
    """Convert HPA nTPM value to an expression level."""
    try:
        v = float(ntpm)
        if v >= 100: return "High"
        if v >= 10:  return "Medium"
        if v >= 1:   return "Low"
        return "Not detected"
    except Exception:
        return "Not detected"


def _parse_hpa_json(data, symbol: str = "") -> dict:
    """Parse HPA JSON response — based on the actual structure

    Actual top-level HPA JSON fields:
      "RNA single cell type specific nTPM" : {cell_type: nTPM_float, ...}
      "RNA single cell type specificity"   : "Cell type enhanced"
      "RNA single cell type Tau score"     : 0.62
      "Subcellular location"               : {...}
    """
    if isinstance(data, list):
        data = data[0] if data else {}
    if not isinstance(data, dict) or not data:
        return {}

    result = {
        "name": data.get("Gene name", symbol),
        "ensg": data.get("Ensembl", ""),
        "tau_score": None,
        "cell_specificity": data.get("RNA single cell type specificity",
                                     data.get("RNA cell type specificity", "Unknown")),
        "subcellular_locations": [],
        "single_cell_expression": {},  # {cell_type_lower: "High"/"Medium"/"Low"/"Not detected"}
        "single_cell_ntpm": {},  # {cell_type_lower: nTPM_float}  (original numeric value)
        "raw": data,
    }

    # ── Tau score ────────────────────────────────────────────────────
    tau_raw = (data.get("RNA single cell type Tau score") or
               data.get("RNA single cell type tau score") or
               data.get("Tau score"))
    if tau_raw is not None:
        try:
            result["tau_score"] = float(tau_raw)
        except:
            pass

    # ── Single-cell nTPM → level conversion ──────────────────────────────────
    # Actual HPA key: "RNA single cell type specific nTPM"
    sc_ntpm = (data.get("RNA single cell type specific nTPM") or
               data.get("RNA single cell type nTPM") or
               data.get("Single cell type") or {})
    if isinstance(sc_ntpm, dict):
        for ct, val in sc_ntpm.items():
            ct_lower = ct.lower()
            if isinstance(val, dict):
                # Legacy JSON: mixed format {"Level": "High", "nTPM": 234.5}
                if "nTPM" in val:
                    ntpm = val["nTPM"]
                    result["single_cell_expression"][ct_lower] = ntpm_to_level(ntpm)
                    result["single_cell_ntpm"][ct_lower] = float(ntpm)
                elif "Level" in val:
                    result["single_cell_expression"][ct_lower] = val["Level"]
                    result["single_cell_ntpm"][ct_lower] = None
            else:
                # New JSON: direct nTPM numeric value
                level = ntpm_to_level(val)
                result["single_cell_expression"][ct_lower] = level
                try:
                    result["single_cell_ntpm"][ct_lower] = float(val)
                except:
                    pass

    # ── Subcellular locations ───────────────────────────────────────────────
    sub = data.get("Subcellular location", {})
    if isinstance(sub, dict):
        locs = sub.get("Approved", []) or sub.get("Predicted", []) or []
        result["subcellular_locations"] = locs if isinstance(locs, list) else [locs]
    elif isinstance(sub, list):
        result["subcellular_locations"] = sub
    elif isinstance(sub, str) and sub:
        result["subcellular_locations"] = [s.strip() for s in sub.split(";") if s.strip()]

    return result


def fetch_hpa_data(ensg_or_marker: str, cache: dict, symbol: str = "") -> dict:
    """HPA JSON API — 3-step URL fallback strategy
    ensg_or_marker: ENSG ID or "SYM:GENENAME" format
    """
    cache_key = ensg_or_marker
    if cache_key in cache:
        return cache[cache_key]

    # Handle SYM: marker — try direct symbol query without ENSG
    if ensg_or_marker.startswith("SYM:"):
        sym = ensg_or_marker[4:]
    else:
        sym = symbol or ensg_or_marker

    empty_result = {
        "ensg": ensg_or_marker if not ensg_or_marker.startswith("SYM:") else "",
        "name": sym, "tau_score": None,
        "cell_specificity": "Unknown",
        "subcellular_locations": [],
        "single_cell_expression": {},
        "raw": None,
    }

    urls_to_try = []
    if not ensg_or_marker.startswith("SYM:"):
        urls_to_try.append(("ENSG", HPA_JSON_ENSG.format(ensg=ensg_or_marker)))
    if sym:
        urls_to_try.append(("SYM", HPA_JSON_SYM.format(symbol=sym)))
        urls_to_try.append(("SRCH", HPA_SEARCH.format(symbol=urllib.parse.quote(sym))))

    parsed = {}
    for tag, url in urls_to_try:
        data = fetch_json(url, timeout=12)
        if data:
            parsed = _parse_hpa_json(data, sym)
            if parsed.get("name") or parsed.get("single_cell_expression"):
                break

    if not parsed:
        cache[cache_key] = empty_result
        return empty_result

    # Merge
    final = {**empty_result, **parsed}
    if not final["ensg"] and not ensg_or_marker.startswith("SYM:"):
        final["ensg"] = ensg_or_marker
    cache[cache_key] = final
    return final


def fetch_string_scores(gene_a: str, gene_b: str, species: int = 9606) -> dict | None:
    """
    Fetch STRING evidence scores for a specific protein pair.

    Returns a dictionary with the combined score and channel-specific scores
    on the standard STRING JSON scale (0.0–1.0):
      score   = combined score
      nscore  = gene neighborhood
      fscore  = gene fusion
      pscore  = phylogenetic profile
      ascore  = coexpression
      escore  = experimental
      dscore  = database
      tscore  = text mining
    """
    url = (f"{STRING_API}/network?"
           f"identifiers={urllib.parse.quote(gene_a)}%0D{urllib.parse.quote(gene_b)}"
           f"&species={species}&required_score=0&network_type=functional"
           f"&caller_identity=PPIValidator")
    data = fetch_json(url)
    if data and isinstance(data, list):
        best_payload = None
        best_score = None
        for interaction in data:
            a_name = (interaction.get("preferredName_A", "") or "").upper()
            b_name = (interaction.get("preferredName_B", "") or "").upper()
            a_id = (interaction.get("stringId_A", "") or "").upper()
            b_id = (interaction.get("stringId_B", "") or "").upper()
            combined = interaction.get("score", None)
            if combined is None:
                continue

            matched = False
            if {a_name, b_name} == {gene_a.upper(), gene_b.upper()}:
                matched = True
            elif gene_a.upper() in a_id and gene_b.upper() in b_id:
                matched = True
            elif gene_b.upper() in a_id and gene_a.upper() in b_id:
                matched = True

            if not matched:
                continue

            try:
                combined = float(combined)
            except Exception:
                continue

            payload = {
                "score": combined,
                "nscore": float(interaction.get("nscore", 0.0) or 0.0),
                "fscore": float(interaction.get("fscore", 0.0) or 0.0),
                "pscore": float(interaction.get("pscore", 0.0) or 0.0),
                "ascore": float(interaction.get("ascore", 0.0) or 0.0),
                "escore": float(interaction.get("escore", 0.0) or 0.0),
                "dscore": float(interaction.get("dscore", 0.0) or 0.0),
                "tscore": float(interaction.get("tscore", 0.0) or 0.0),
            }

            if best_score is None or combined > best_score:
                best_score = combined
                best_payload = payload

        return best_payload
    return None


# ─── Path file parser ──────────────────────────────────────────────────────────
def parse_pathway_file(filepath: str):
    """
    Parse pathway files in multiple formats:
    - "Path #N: A -> B -> C"
    - "Unique Path #N (original #M): A -> B -> C"
    """
    paths = []
    seen = set()
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            # Original pathway formats
            m = re.search(r'(?:Path|Unique Path)\s*#\d+.*?:\s*(.+)', line)
            if m:
                seq = m.group(1).strip()
                nodes = [n.strip() for n in seq.split("->")]
                nodes = [n for n in nodes if n]
                # Exclude single-node paths
                if len(nodes) >= 2:
                    key = tuple(nodes)
                    if key not in seen:
                        seen.add(key)
                        paths.append(nodes)
    return paths


# ─── Validation engine ────────────────────────────────────────────────────────────────
class ValidationEngine:
    def __init__(self, log_fn):
        self.log = log_fn
        self.ensembl_cache = {}
        self.hpa_cache = {}
        self.string_cache = {}

    def get_expression_level(self, hpa_data: dict, cell_type: str) -> tuple:
        """Get expression level — based on HPA nTPM (3-step fallback)
        Returns: (level, source_description)
        """
        sc = hpa_data.get("single_cell_expression", {})
        ntpm = hpa_data.get("single_cell_ntpm", {})
        spec = hpa_data.get("cell_specificity", "Unknown")
        ct_lower = cell_type.lower()

        def ntpm_str(key):
            v = ntpm.get(key)
            return f" (nTPM={v:.1f})" if v is not None else ""

        # Step 1: exact key match
        if ct_lower in sc:
            return sc[ct_lower], f"HPA nTPM exact match: '{ct_lower}'{ntpm_str(ct_lower)}"

        # Step 2: partial word-set matching
        ct_words = set(ct_lower.replace("(", "").replace(")", "").split())
        best_key, best_ov = None, 0
        for k in sc:
            kw = set(k.replace("(", "").replace(")", "").split())
            ov = len(ct_words & kw)
            if ov > best_ov and ov >= 1:
                best_ov, best_key = ov, k
        if best_key:
            return sc[best_key], f"HPA nTPM partial match: '{best_key}'{ntpm_str(best_key)}"

        # Step 3: specificity fallback (when nTPM is unavailable)
        spec_map = {
            "Cell type enriched": "High",
            "Tissue enriched": "High",
            "Cell type enhanced": "Medium",
            "Tissue enhanced": "Medium",
            "Group enriched": "Medium",
            "Low cell type specificity": "Low",
            "Low tissue specificity": "Low",
            "Not detected": "Not detected",
        }
        fb_level = "Not detected"
        for sk, lv in spec_map.items():
            if sk.lower() in spec.lower():
                fb_level = lv
                break

        avail_sample = ", ".join(list(sc.keys())[:3]) if sc else "none"
        if sc:
            return fb_level, f"⚠️ Cell type mismatch (e.g. {avail_sample}) → specificity fallback({spec})"
        else:
            return fb_level, f"⚠️ Single cell nTPM No data → specificity fallback({spec})"

    def check_colocalization(self, hpa_a: dict, hpa_b: dict) -> tuple[bool, str]:
        """Check subcellular co-localization between two proteins."""
        locs_a = set(get_compartment_group(l) for l in hpa_a.get("subcellular_locations", []))
        locs_b = set(get_compartment_group(l) for l in hpa_b.get("subcellular_locations", []))
        if not locs_a or not locs_b:
            return None, "Unknown (no localization data)"
        overlap = locs_a & locs_b
        if overlap:
            return True, ", ".join(overlap)
        return False, f"A:{','.join(locs_a)} vs B:{','.join(locs_b)}"

    def get_string_score(self, gene_a: str, gene_b: str) -> dict | None:
        key = tuple(sorted([gene_a, gene_b]))
        if key in self.string_cache:
            return self.string_cache[key]
        score_info = fetch_string_scores(gene_a, gene_b)
        self.string_cache[key] = score_info
        return score_info

    def score_ppi(self, expr_a, expr_b, coloc: bool | None,
                  tau_a: float | None, tau_b: float | None,
                  string_score: dict | float | None, cell_type: str) -> tuple[str, int, list]:
        """
        Calculate the overall validation score
        Returns: (grade, score, reason list)
        """
        score = 0
        reasons = []

        # 1. Expression score (both sides must be expressed)
        # If a tuple is provided (level, source), extract only level
        if isinstance(expr_a, tuple): expr_a = expr_a[0]
        if isinstance(expr_b, tuple): expr_b = expr_b[0]
        level_map = {"High": 3, "Medium": 2, "Low": 1, "Not detected": 0}
        lv_a = level_map.get(expr_a, 0)
        lv_b = level_map.get(expr_b, 0)

        if lv_a == 0 or lv_b == 0:
            reasons.append(f"⚠️ Not expressed: A={expr_a}, B={expr_b}")
        elif lv_a >= 2 and lv_b >= 2:
            score += 3
            reasons.append(f"✅ Both highly expressed: A={expr_a}, B={expr_b}")
        else:
            score += 1
            reasons.append(f"🔶 Low expression: A={expr_a}, B={expr_b}")

        # 2. Co-localization score
        if coloc is True:
            score += 2
            reasons.append("✅ Subcellular co-localization confirmed")
        elif coloc is False:
            reasons.append("❌ Subcellular location mismatch")
        else:
            reasons.append("❓ No localization data")

        # 3. STRING confidence
        if string_score is not None:
            if isinstance(string_score, dict):
                ss = float(string_score.get("score", 0.0) or 0.0)
                esc = float(string_score.get("escore", 0.0) or 0.0)
                dsc = float(string_score.get("dscore", 0.0) or 0.0)
                tsc = float(string_score.get("tscore", 0.0) or 0.0)
                asc = float(string_score.get("ascore", 0.0) or 0.0)
            else:
                ss = float(string_score)
                esc = dsc = tsc = asc = 0.0

            if ss >= 0.9:
                score += 3
                reasons.append(f"✅ High STRING confidence: {ss:.3f}")
            elif ss >= 0.7:
                score += 2
                reasons.append(f"🔶 Medium STRING confidence: {ss:.3f}")
            elif ss >= 0.4:
                score += 1
                reasons.append(f"🔶 Low STRING confidence: {ss:.3f}")
            else:
                reasons.append(f"⚠️ Very low STRING confidence: {ss:.3f}")

            channel_pairs = [("database", dsc), ("experimental", esc), ("textmining", tsc), ("coexpression", asc)]
            best_label, best_value = max(channel_pairs, key=lambda x: x[1])
            if best_value > 0:
                reasons.append(f"ℹ️ Top STRING evidence: {best_label}={best_value:.3f}")
        else:
            reasons.append("❓ STRING No data")

        # 4. Tau score (cell-type specificity)
        if tau_a is not None and tau_b is not None:
            avg_tau = (tau_a + tau_b) / 2
            if avg_tau >= 0.7:
                score += 2
                reasons.append(f"✅ High cell-type specificity Tau={avg_tau:.2f}")
            elif avg_tau >= 0.4:
                score += 1
                reasons.append(f"🔶 Medium cell-type specificity Tau={avg_tau:.2f}")
            else:
                reasons.append(f"🔶 Low cell-type specificity Tau={avg_tau:.2f}")

        # Grade assignment (maximum 10 points)
        if lv_a == 0 or lv_b == 0:
            grade = "NOT_EXPRESSED"
        elif score >= 7:
            grade = "HIGH"
        elif score >= 4:
            grade = "MEDIUM"
        elif score >= 2:
            grade = "LOW"
        else:
            grade = "UNLIKELY"

        return grade, score, reasons

    def validate_paths(self, paths: list, cell_type: str,
                       progress_fn=None, stop_flag=None) -> list:
        """Validate the full list of pathways."""
        results = []
        total_edges = sum(len(p) - 1 for p in paths)
        done = 0

        for path_idx, path in enumerate(paths):
            if stop_flag and stop_flag():
                break

            path_result = {
                "path_idx": path_idx + 1,
                "path": path,
                "edges": [],
                "proteins": {},
                "path_grade": "UNKNOWN",
                "path_score": 0,
            }

            # Collect HPA data for each protein
            for gene in path:
                if gene not in path_result["proteins"]:
                    self.log(f"🔍 [{path_idx + 1}/{len(paths)}] HPA query: {gene}")
                    ensg = gene_symbol_to_ensembl(gene, self.ensembl_cache)
                    if ensg:
                        hpa = fetch_hpa_data(ensg, self.hpa_cache, symbol=gene)
                        self.log(
                            f"  [HPA] {'✅' if hpa.get('raw') else '⚠️'} {gene}: {hpa.get('cell_specificity', '?')}")
                    else:
                        hpa = {"ensg": None, "name": gene,
                               "tau_score": None, "cell_specificity": "Unknown",
                               "subcellular_locations": [], "single_cell_expression": {}}
                    path_result["proteins"][gene] = hpa
                    time.sleep(0.3)  # API rate limit

            # Edge-by-edge validation
            edge_grades = []
            for i in range(len(path) - 1):
                gene_a, gene_b = path[i], path[i + 1]
                hpa_a = path_result["proteins"][gene_a]
                hpa_b = path_result["proteins"][gene_b]

                expr_a, expr_src_a = self.get_expression_level(hpa_a, cell_type)
                expr_b, expr_src_b = self.get_expression_level(hpa_b, cell_type)
                coloc, coloc_detail = self.check_colocalization(hpa_a, hpa_b)

                self.log(f"  → STRING: {gene_a} — {gene_b}")
                string_score = self.get_string_score(gene_a, gene_b)
                if isinstance(string_score, dict):
                    self.log(
                        "    [STRING] "
                        f"combined={float(string_score.get('score', 0.0) or 0.0):.3f} | "
                        f"exp={float(string_score.get('escore', 0.0) or 0.0):.3f}, "
                        f"db={float(string_score.get('dscore', 0.0) or 0.0):.3f}, "
                        f"text={float(string_score.get('tscore', 0.0) or 0.0):.3f}, "
                        f"coexp={float(string_score.get('ascore', 0.0) or 0.0):.3f}"
                    )
                time.sleep(0.2)

                grade, score, reasons = self.score_ppi(
                    expr_a, expr_b, coloc,
                    hpa_a.get("tau_score"), hpa_b.get("tau_score"),
                    string_score, cell_type
                )

                edge_result = {
                    "gene_a": gene_a, "gene_b": gene_b,
                    "expr_a": expr_a, "expr_b": expr_b,
                    "expr_src_a": expr_src_a, "expr_src_b": expr_src_b,
                    "coloc": coloc, "coloc_detail": coloc_detail,
                    "string_score": (string_score.get("score") if isinstance(string_score, dict) else string_score),
                    "string_channels": (string_score if isinstance(string_score, dict) else None),
                    "grade": grade, "score": score, "reasons": reasons,
                    "locs_a": hpa_a.get("subcellular_locations", []),
                    "locs_b": hpa_b.get("subcellular_locations", []),
                    "tau_a": hpa_a.get("tau_score"),
                    "tau_b": hpa_b.get("tau_score"),
                    "sc_keys_a": list(hpa_a.get("single_cell_expression", {}).keys())[:5],
                    "sc_keys_b": list(hpa_b.get("single_cell_expression", {}).keys())[:5],
                }
                path_result["edges"].append(edge_result)
                edge_grades.append(grade)

                done += 1
                if progress_fn:
                    progress_fn(done, total_edges)

            # Overall pathway grade (based on the worst edge)
            grade_order = ["HIGH", "MEDIUM", "LOW", "UNLIKELY", "NOT_EXPRESSED", "UNKNOWN"]
            worst = min(edge_grades, key=lambda g: grade_order.index(g) if g in grade_order else 99,
                        default="UNKNOWN")
            path_result["path_grade"] = worst
            path_result["path_score"] = sum(e["score"] for e in path_result["edges"])
            results.append(path_result)

        return results


# ─── HTML generator ──────────────────────────────────────────────────────────────
def grade_color(grade: str) -> str:
    return {
        "HIGH": "#22c55e",
        "MEDIUM": "#f59e0b",
        "LOW": "#f97316",
        "UNLIKELY": "#ef4444",
        "NOT_EXPRESSED": "#991b1b",
        "UNKNOWN": "#94a3b8",
    }.get(grade, "#94a3b8")


def grade_bg(grade: str) -> str:
    return {
        "HIGH": "#dcfce7",
        "MEDIUM": "#fef3c7",
        "LOW": "#ffedd5",
        "UNLIKELY": "#fee2e2",
        "NOT_EXPRESSED": "#fecaca",
        "UNKNOWN": "#f1f5f9",
    }.get(grade, "#f1f5f9")


def grade_label_kr(grade: str) -> str:
    return {
        "HIGH": "✅ High",
        "MEDIUM": "🔶 Medium",
        "LOW": "⚠️ Low",
        "UNLIKELY": "❌ Unlikely",
        "NOT_EXPRESSED": "🚫 Not expressed",
        "UNKNOWN": "❓ Unknown",
    }.get(grade, grade)


def expr_badge(level) -> str:
    if isinstance(level, tuple): level = level[0]  # tuple guard
    colors = {"High": "#166534", "Medium": "#92400e", "Low": "#1e40af", "Not detected": "#6b7280"}
    bgs = {"High": "#dcfce7", "Medium": "#fef3c7", "Low": "#dbeafe", "Not detected": "#f3f4f6"}
    c = colors.get(level, "#6b7280")
    b = bgs.get(level, "#f3f4f6")
    return f'<span style="background:{b};color:{c};padding:2px 8px;border-radius:12px;font-size:12px;font-weight:600">{level}</span>'


def generate_html(results: list, cell_type: str, filepath: str,
                  source_file: str = "") -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Summary statistics
    grade_counts = defaultdict(int)
    for r in results:
        grade_counts[r["path_grade"]] += 1
    total = len(results)
    high = grade_counts["HIGH"]
    medium = grade_counts["MEDIUM"]
    low = grade_counts["LOW"] + grade_counts["UNLIKELY"] + grade_counts["NOT_EXPRESSED"]

    # Collect all proteins
    all_proteins = {}
    for r in results:
        for g, hpa in r["proteins"].items():
            if g not in all_proteins:
                all_proteins[g] = hpa

    # ── Build pathway card HTML ─────────────────────────────────────────
    path_cards = []
    for r in results:
        grade = r["path_grade"]
        path_str = " → ".join(r["path"])
        edge_rows = []
        for e in r["edges"]:
            ss = f"{float(e['string_score']):.3f}" if e['string_score'] is not None else "N/A"
            channels = e.get("string_channels") or {}
            string_sub = (
                f"Exp {float(channels.get('escore', 0.0) or 0.0):.3f} | "
                f"DB {float(channels.get('dscore', 0.0) or 0.0):.3f} | "
                f"Text {float(channels.get('tscore', 0.0) or 0.0):.3f} | "
                f"Coexp {float(channels.get('ascore', 0.0) or 0.0):.3f}"
            ) if channels else "No channel data"
            locs_a = ", ".join(e["locs_a"][:2]) or "Unknown"
            locs_b = ", ".join(e["locs_b"][:2]) or "Unknown"
            tau_a = f"{e['tau_a']:.2f}" if e['tau_a'] else "N/A"
            tau_b = f"{e['tau_b']:.2f}" if e['tau_b'] else "N/A"
            coloc_icon = "✅" if e["coloc"] else ("❌" if e["coloc"] is False else "❓")
            reasons_html = "<br>".join(e["reasons"])
            src_a = e.get("expr_src_a", "")
            src_b = e.get("expr_src_b", "")
            sc_a = ", ".join(e.get("sc_keys_a", [])[:3]) or "none"
            sc_b = ", ".join(e.get("sc_keys_b", [])[:3]) or "none"
            edge_rows.append(f"""
            <tr style="border-bottom:1px solid #e2e8f0">
              <td style="padding:8px;font-weight:600">{e['gene_a']}</td>
              <td style="padding:8px;text-align:center">→</td>
              <td style="padding:8px;font-weight:600">{e['gene_b']}</td>
              <td style="padding:8px;text-align:center">{expr_badge(e['expr_a'])}<br><small style="color:#94a3b8;font-size:10px;display:block;max-width:130px;overflow:hidden">{src_a[:40]}</small></td>
              <td style="padding:8px;text-align:center">{expr_badge(e['expr_b'])}<br><small style="color:#94a3b8;font-size:10px;display:block;max-width:130px;overflow:hidden">{src_b[:40]}</small></td>
              <td style="padding:8px;text-align:center;font-size:12px">{locs_a}<br><span style="color:#64748b">{tau_a}</span></td>
              <td style="padding:8px;text-align:center;font-size:12px">{locs_b}<br><span style="color:#64748b">{tau_b}</span></td>
              <td style="padding:8px;text-align:center">{coloc_icon} {e.get('coloc_detail', '')[:30]}</td>
              <td style="padding:8px;text-align:center;font-weight:600">{ss}<br><span style="color:#64748b;font-size:10px">{string_sub}</span></td>
              <td style="padding:8px;background:{grade_bg(e['grade'])};color:{grade_color(e['grade'])};font-weight:700;text-align:center">{grade_label_kr(e['grade'])}</td>
              <td style="padding:8px;font-size:11px;color:#475569">{reasons_html}</td>
            </tr>""")

        path_cards.append(f"""
        <div class="path-card" style="margin:16px 0;border:2px solid {grade_color(grade)};border-radius:12px;overflow:hidden">
          <div style="background:{grade_color(grade)};padding:10px 16px;display:flex;justify-content:space-between;align-items:center">
            <span style="color:white;font-weight:700;font-size:14px">Path #{r['path_idx']}: {path_str}</span>
            <span style="background:white;color:{grade_color(grade)};padding:3px 12px;border-radius:16px;font-weight:800;font-size:13px">{grade_label_kr(grade)} (Score: {r['path_score']})</span>
          </div>
          <div style="overflow-x:auto">
            <table style="width:100%;border-collapse:collapse;font-size:13px">
              <thead>
                <tr style="background:#f8fafc;font-size:12px;color:#64748b">
                  <th style="padding:8px;text-align:left">Protein A</th>
                  <th style="padding:8px"></th>
                  <th style="padding:8px;text-align:left">Protein B</th>
                  <th style="padding:8px;text-align:center">Expression (A)</th>
                  <th style="padding:8px;text-align:center">Expression (B)</th>
                  <th style="padding:8px;text-align:center">Location (A) / Tau</th>
                  <th style="padding:8px;text-align:center">Location (B) / Tau</th>
                  <th style="padding:8px;text-align:center">Co-localization</th>
                  <th style="padding:8px;text-align:center">STRING combined / channel scores</th>
                  <th style="padding:8px;text-align:center">Grade</th>
                  <th style="padding:8px;text-align:left">Evidence</th>
                </tr>
              </thead>
              <tbody>{"".join(edge_rows)}</tbody>
            </table>
          </div>
        </div>""")

    # ── Protein summary table ───────────────────────────────────────────
    prot_rows = []
    for gene, hpa in sorted(all_proteins.items()):
        ensg = hpa.get("ensg", "N/A") or "N/A"
        tau = f"{hpa['tau_score']:.2f}" if hpa.get("tau_score") else "N/A"
        spec = hpa.get("cell_specificity", "Unknown")
        locs = ", ".join(hpa.get("subcellular_locations", [])[:3]) or "N/A"
        hpa_link = f'<a href="https://www.proteinatlas.org/{ensg}" target="_blank" style="color:#3b82f6">{ensg}</a>' if ensg != "N/A" else "N/A"
        prot_rows.append(f"""
        <tr style="border-bottom:1px solid #e2e8f0">
          <td style="padding:8px;font-weight:700">{gene}</td>
          <td style="padding:8px">{hpa_link}</td>
          <td style="padding:8px">{tau}</td>
          <td style="padding:8px;font-size:12px">{spec}</td>
          <td style="padding:8px;font-size:12px">{locs}</td>
        </tr>""")

    # ── Bar chart data
    # ─ Expression matching diagnostics table
    seen_genes_debug = set()
    debug_rows_list = []
    spec_map = {"Cell type enriched": "High", "Tissue enriched": "High",
                "Cell type enhanced": "Medium", "Tissue enhanced": "Medium",
                "Group enriched": "Medium", "Low cell type specificity": "Low",
                "Low tissue specificity": "Low", "Not detected": "Not detected"}
    for r in results:
        for gene, hpa in r["proteins"].items():
            if gene in seen_genes_debug: continue
            seen_genes_debug.add(gene)
            sc = hpa.get("single_cell_expression", {})
            spec = hpa.get("cell_specificity", "Unknown")
            ct_lower = cell_type.lower()
            ct_words = set(ct_lower.replace("(", "").replace(")", "").split())
            level_found = "Not detected";
            src_label = "No data"
            if ct_lower in sc:
                level_found = sc[ct_lower];
                src_label = f"✅ Exact: '{cell_type}'"
            else:
                best_key = None;
                best_ov = 0
                for k in sc:
                    kw = set(k.lower().replace("(", "").replace(")", "").split())
                    ov = len(ct_words & kw)
                    if ov > best_ov and ov >= 1: best_ov = ov; best_key = k
                if best_key:
                    level_found = sc[best_key];
                    src_label = f"🔶 Partial match: '{best_key}'"
                else:
                    for sk, lv in spec_map.items():
                        if sk.lower() in spec.lower(): level_found = lv; break
                    src_label = f"⚠️ Cell-type mismatch → fallback({spec})"
            avail = ", ".join(list(sc.keys())[:5]) or "none"
            bg = {"High": "#dcfce7", "Medium": "#fef3c7", "Low": "#dbeafe", "Not detected": "#fce4ec"}.get(level_found,
                                                                                                           "#f1f5f9")
            debug_rows_list.append(
                f'<tr style="border-bottom:1px solid #e2e8f0">'
                f'<td style="padding:8px;font-weight:700">{gene}</td>'
                f'<td style="padding:8px"><span style="background:{bg};padding:2px 8px;border-radius:10px;font-size:12px;font-weight:600">{level_found}</span></td>'
                f'<td style="padding:8px;font-size:12px;color:#475569">{src_label}</td>'
                f'<td style="padding:8px;font-size:11px;color:#64748b">{avail}</td></tr>')
    debug_rows = "".join(
        debug_rows_list) if debug_rows_list else "<tr><td colspan=4 style='padding:12px;color:#94a3b8'>No diagnostic data</td></tr>"

    grade_labels = ["HIGH", "MEDIUM", "LOW", "UNLIKELY", "NOT_EXPRESSED", "UNKNOWN"]
    grade_data = [grade_counts[g] for g in grade_labels]
    grade_colors_js = [grade_color(g) for g in grade_labels]
    grade_labels_kr = [grade_label_kr(g) for g in grade_labels]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PPI Cell-Type Feasibility Validation Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4"></script>
<style>
  * {{ box-sizing:border-box; margin:0; padding:0; }}
  body {{ font-family:Arial, sans-serif;background:#f0f4f8;color:#1e293b }}
  .container {{ max-width:1400px;margin:0 auto;padding:24px }}
  .header {{ background:linear-gradient(135deg,#0d2137 0%,#1e40af 100%);color:white;padding:32px;border-radius:16px;margin-bottom:24px }}
  .header h1 {{ font-size:24px;font-weight:800;margin-bottom:8px }}
  .header .sub {{ font-size:14px;opacity:0.8 }}
  .stat-grid {{ display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:16px;margin-bottom:24px }}
  .stat-card {{ background:white;border-radius:12px;padding:20px;text-align:center;box-shadow:0 1px 4px rgba(0,0,0,.08) }}
  .stat-card .num {{ font-size:36px;font-weight:800 }}
  .stat-card .lbl {{ font-size:13px;color:#64748b;margin-top:4px }}
  .section {{ background:white;border-radius:12px;padding:24px;margin-bottom:24px;box-shadow:0 1px 4px rgba(0,0,0,.08) }}
  .section h2 {{ font-size:18px;font-weight:700;margin-bottom:16px;padding-bottom:8px;border-bottom:2px solid #e2e8f0 }}
  .path-card {{ transition:box-shadow .2s }}
  .path-card:hover {{ box-shadow:0 4px 16px rgba(0,0,0,.12) }}
  .filter-bar {{ display:flex;gap:12px;margin-bottom:16px;flex-wrap:wrap }}
  .filter-btn {{ padding:6px 16px;border:2px solid #e2e8f0;border-radius:20px;cursor:pointer;font-size:13px;background:white;transition:all .2s }}
  .filter-btn.active {{ border-color:#3b82f6;background:#dbeafe;color:#1e40af;font-weight:600 }}
  .legend {{ display:flex;gap:16px;flex-wrap:wrap;margin-bottom:12px }}
  .legend-item {{ display:flex;align-items:center;gap:6px;font-size:13px }}
  .legend-dot {{ width:12px;height:12px;border-radius:50% }}
  table {{ width:100%;border-collapse:collapse }}
  th {{ padding:10px 8px;text-align:left;font-size:12px;color:#64748b;font-weight:600 }}
  tr:hover td {{ background:#f8fafc }}
  @media print {{ .filter-bar {{ display:none }} }}
</style>
</head>
<body>
<div class="container">

<!-- HEADER -->
<div class="header">
  <h1>🔬 PPI Cell-Type Feasibility Validation Report</h1>
  <div class="sub">
    Target cell type: <strong>{cell_type}</strong> &nbsp;|&nbsp;
    Number of analyzed pathways: <strong>{total}</strong> &nbsp;|&nbsp;
    Analysis time: {now}
    {f'<br>Source file: {os.path.basename(source_file)}' if source_file else ''}
  </div>
</div>

<!-- Summary statistics -->
<div class="stat-grid">
  <div class="stat-card"><div class="num" style="color:#22c55e">{high}</div><div class="lbl">✅ High feasibility (HIGH)</div></div>
  <div class="stat-card"><div class="num" style="color:#f59e0b">{medium}</div><div class="lbl">🔶 Medium feasibility (MEDIUM)</div></div>
  <div class="stat-card"><div class="num" style="color:#ef4444">{low}</div><div class="lbl">❌ Low / not expressed</div></div>
  <div class="stat-card"><div class="num">{total}</div><div class="lbl">Total validated pathways</div></div>
  <div class="stat-card"><div class="num">{len(all_proteins)}</div><div class="lbl">Number of unique proteins</div></div>
</div>

<!-- Charts -->
<div class="section">
  <h2>📊 Validation grade distribution</h2>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:24px">
    <canvas id="gradeChart" height="200"></canvas>
    <canvas id="pieChart" height="200"></canvas>
  </div>
</div>

<!-- Validation legend -->
<div class="section">
  <h2>🎯 Validation methodology</h2>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px;font-size:13px">
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>① Expression validation</strong><br>Based on HPA single-cell expression data<br>High/Medium/Low/Not detected</div>
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>② Tau score</strong><br>Cell-type specificity metric (0~1)<br>≥0.7 = high specificity, ≥0.4 = medium specificity</div>
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>③ Subcellular co-localization</strong><br>Comparison of HPA subcellular locations<br>Same organelle group = co-localization possible</div>
    <div style="background:#f8fafc;padding:12px;border-radius:8px"><strong>④ STRING confidence</strong><br>STRING v11.5 combined PPI score<br>≥0.9 = high confidence, ≥0.7 = medium confidence</div>
  </div>
  <div class="legend" style="margin-top:16px">
    {"".join(f'<div class="legend-item"><div class="legend-dot" style="background:{grade_color(g)}"></div><span>{grade_label_kr(g)}</span></div>' for g in grade_labels)}
  </div>
</div>

<!-- Path filters -->
<div class="section">
  <h2>🧬 PPI validation results by pathway</h2>
  <div class="filter-bar">
    <button class="filter-btn active" onclick="filterPaths('ALL')">Show all</button>
    {"".join(f'<button class="filter-btn" onclick="filterPaths(&quot;{g}&quot;)" style="border-color:{grade_color(g)}">{grade_label_kr(g)}</button>' for g in grade_labels if grade_counts[g] > 0)}
  </div>
  <div id="path-container">
    {"".join(f'<div class="path-card grade-{r["path_grade"]}">{card}</div>' for r, card in zip(results, path_cards))}
  </div>
</div>

<!-- Protein summary -->
<div class="section">
  <h2>🔍 Expression matching diagnostics (selected cell type:  {cell_type})</h2>
  <p style="font-size:13px;color:#64748b;margin-bottom:12px">
    ✅ Exact match: exact match to an HPA Single Cell Type key &nbsp;|&nbsp;
    🔶 Partial match: word overlap &nbsp;|&nbsp;
    ⚠️ Fallback: no matching cell type, using HPA specificity field
  </p>
  <div style="overflow-x:auto">
    <table><thead><tr style="background:#f8fafc">
      <th style="padding:8px">Gene</th>
      <th style="padding:8px">Expression level</th>
      <th style="padding:8px">Matching evidence</th>
      <th style="padding:8px">Available HPA cell types (top 5)</th>
    </tr></thead>
    <tbody>{debug_rows}</tbody></table>
  </div>
</div>

<div class="section">
  <h2>🔬 Protein expression profile summary</h2>
  <div style="overflow-x:auto">
    <table>
      <thead>
        <tr style="background:#f8fafc">
          <th>Gene</th>
          <th>Ensembl ID</th>
          <th>Tau Score</th>
          <th>Cell specificity</th>
          <th>Subcellular location</th>
        </tr>
      </thead>
      <tbody>{"".join(prot_rows)}</tbody>
    </table>
  </div>
</div>

<!-- Footer -->
<div style="text-align:center;color:#94a3b8;font-size:12px;padding:16px 0">
  PPI Cell-Type Validator v{APP_VERSION} &nbsp;|&nbsp; AI4Emotion &nbsp;|&nbsp;
  HPA API · STRING v12 API · MyGene.info based on &nbsp;|&nbsp; Generated: {now}
</div>
</div>

<script>
// ─── Charts
const grades  = {json.dumps(grade_labels_kr)};
const counts  = {json.dumps(grade_data)};
const colors  = {json.dumps(grade_colors_js)};

new Chart(document.getElementById('gradeChart'), {{
  type:'bar',
  data:{{ labels:grades, datasets:[{{ label:'Path count', data:counts,
    backgroundColor:colors, borderRadius:8 }}]}},
  options:{{ responsive:true,plugins:{{legend:{{display:false}},
    title:{{display:true,text:'Path count by grade'}}}} }}
}});

new Chart(document.getElementById('pieChart'), {{
  type:'doughnut',
  data:{{ labels:grades, datasets:[{{ data:counts, backgroundColor:colors }}]}},
  options:{{ responsive:true,plugins:{{legend:{{position:'right'}},
    title:{{display:true,text:'Proportion by grade'}}}} }}
}});

// ─── Filters
function filterPaths(grade) {{
  document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
  document.querySelectorAll('#path-container .path-card').forEach(card => {{
    card.style.display = (grade === 'ALL' || card.classList.contains('grade-' + grade)) ? 'block' : 'none';
  }});
}}
</script>
</body>
</html>"""

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(html)
    return filepath


# ─── Tkinter GUI ──────────────────────────────────────────────────────────────
class PPIValidatorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"PPI Cell-Type Validator v{APP_VERSION}")
        self.geometry("860x800")
        self.resizable(True, True)
        self.configure(bg="#f0f4f8")

        self._stop = False
        self._running = False
        self._paths = []
        self._source_file = ""

        self._build_ui()

    def _build_ui(self):
        # ── Top header
        hdr = tk.Frame(self, bg="#0d2137", pady=14)
        hdr.pack(fill="x")
        tk.Label(hdr, text="🔬  PPI Cell-Type Validator",
                 font=("Arial", 16, "bold"), bg="#0d2137", fg="white").pack()
        tk.Label(hdr, text="Validate cell-type feasibility of PPIs in BFS pathway data",
                 font=("Arial", 10), bg="#0d2137", fg="#90caf9").pack()

        body = tk.Frame(self, bg="#f0f4f8", padx=20, pady=16)
        body.pack(fill="both", expand=True)

        # ── File selection
        f1 = ttk.LabelFrame(body, text="① Select pathway file", padding=12)
        f1.pack(fill="x", pady=6)
        self.file_var = tk.StringVar(value="Please select a file...")
        tk.Label(f1, textvariable=self.file_var, width=60, anchor="w",
                 relief="sunken", bg="white", font=("Arial", 9)).pack(side="left", padx=4)
        ttk.Button(f1, text="Open", command=self._choose_file).pack(side="left", padx=4)
        self.path_count_lbl = tk.Label(f1, text="", fg="#16a34a", font=("Arial", 9))
        self.path_count_lbl.pack(side="left", padx=8)

        # ── Cell type selection
        f2 = ttk.LabelFrame(body, text="② Select target cell type", padding=12)
        f2.pack(fill="x", pady=6)
        tk.Label(f2, text="Cell type:", font=("Arial", 9)).pack(side="left")
        self.cell_var = tk.StringVar(value=CELL_TYPES[0])
        cb = ttk.Combobox(f2, textvariable=self.cell_var, values=CELL_TYPES,
                          state="readonly", width=40)
        cb.pack(side="left", padx=8)
        tk.Label(f2, text="Or enter manually:", font=("Arial", 9)).pack(side="left")
        self.cell_custom = ttk.Entry(f2, width=20)
        self.cell_custom.pack(side="left", padx=4)

        # ── Parameters
        f3 = ttk.LabelFrame(body, text="③ Analysis parameters", padding=12)
        f3.pack(fill="x", pady=6)
        tk.Label(f3, text="Maximum path count:", font=("Arial", 9)).grid(row=0, column=0, sticky="w")
        self.max_paths = ttk.Spinbox(f3, from_=1, to=999, width=6)
        self.max_paths.set(50)
        self.max_paths.grid(row=0, column=1, padx=8, sticky="w")
        tk.Label(f3, text="Minimum STRING score:", font=("Arial", 9)).grid(row=0, column=2, padx=16, sticky="w")
        self.min_score = ttk.Spinbox(f3, from_=0.0, to=1.0, increment=0.1, width=6)
        self.min_score.set(0.4)
        self.min_score.grid(row=0, column=3, padx=8, sticky="w")

        # ── Output file
        f4 = ttk.LabelFrame(body, text="④ HTML output file", padding=12)
        f4.pack(fill="x", pady=6)
        self.out_var = tk.StringVar(value="ppi_validation_result.html")
        tk.Label(f4, textvariable=self.out_var, width=55, anchor="w",
                 relief="sunken", bg="white", font=("Arial", 9)).pack(side="left", padx=4)
        ttk.Button(f4, text="Choose location", command=self._choose_output).pack(side="left", padx=4)

        # ── Progress bar
        self.progress = ttk.Progressbar(body, length=800, mode="determinate")
        self.progress.pack(fill="x", pady=8)
        self.status_lbl = tk.Label(body, text="Ready", font=("Arial", 9),
                                   fg="#475569", bg="#f0f4f8")
        self.status_lbl.pack()

        # ── Log
        f5 = ttk.LabelFrame(body, text="Execution log", padding=6)
        f5.pack(fill="both", expand=True, pady=6)
        self.log_box = scrolledtext.ScrolledText(f5, height=10, font=("Courier", 9),
                                                 state="disabled", bg="#0d1b2a", fg="#90ee90")
        self.log_box.pack(fill="both", expand=True)

        # ── Buttons
        btn_frame = tk.Frame(body, bg="#f0f4f8")
        btn_frame.pack(pady=8)
        self.run_btn = ttk.Button(btn_frame, text="▶  Start validation",
                                  command=self._start, style="Accent.TButton")
        self.run_btn.pack(side="left", padx=8, ipadx=12, ipady=4)
        ttk.Button(btn_frame, text="⏹  Stop",
                   command=self._stop_run).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="🌐  Open result",
                   command=self._open_result).pack(side="left", padx=4)

    def _log(self, msg: str):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", msg + "\n")
        self.log_box.see("end")
        self.log_box.configure(state="disabled")
        self.update_idletasks()

    def _choose_file(self):
        path = filedialog.askopenfilename(
            title="Select pathway file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if path:
            self._source_file = path
            self.file_var.set(path)
            self._paths = parse_pathway_file(path)
            self.path_count_lbl.config(text=f"→ {len(self._paths)}pathways loaded")
            self._log(f"✅ File loaded: {os.path.basename(path)}  ({len(self._paths)} Path)")

    def _choose_output(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".html",
            filetypes=[("HTML", "*.html"), ("All files", "*.*")],
            initialfile="ppi_validation_result.html"
        )
        if path:
            self.out_var.set(path)

    def _start(self):
        if self._running:
            return
        if not self._paths:
            messagebox.showwarning("Warning", "Please select a pathway file first.")
            return
        cell_type = self.cell_custom.get().strip() or self.cell_var.get()
        if not cell_type:
            messagebox.showwarning("Warning", "Please select or enter a cell type.")
            return
        self._running = True
        self._stop = False
        self.run_btn.config(state="disabled")
        self.progress["value"] = 0

        t = threading.Thread(target=self._run_validation,
                             args=(cell_type,), daemon=True)
        t.start()

    def _stop_run(self):
        self._stop = True
        self._log("⏹  Stop requested...")

    def _open_result(self):
        out = self.out_var.get()
        if os.path.exists(out):
            import webbrowser
            webbrowser.open(f"file://{os.path.abspath(out)}")
        else:
            messagebox.showinfo("Notice", "Run the validation first to create the result file.")

    def _update_progress(self, done: int, total: int):
        pct = (done / total * 100) if total else 0
        self.progress["value"] = pct
        self.status_lbl.config(text=f"Progress: {done}/{total} edges being validated...")
        self.update_idletasks()

    def _run_validation(self, cell_type: str):
        try:
            max_p = int(self.max_paths.get())
            paths = self._paths[:max_p]
            self._log(f"\n{'=' * 50}")
            self._log(f"🚀 Validation started | Cell type: {cell_type}")
            self._log(f"   Path count: {len(paths)} | API: HPA + STRING + MyGene.info")
            self._log(f"{'=' * 50}\n")

            engine = ValidationEngine(self._log)
            results = engine.validate_paths(
                paths, cell_type,
                progress_fn=self._update_progress,
                stop_flag=lambda: self._stop
            )

            if results:
                out_path = self.out_var.get()
                generate_html(results, cell_type, out_path, self._source_file)
                self._log(f"\n✅ HTML report saved: {out_path}")
                self._log(f"   Validated pathways: {len(results)}")

                grade_counts = defaultdict(int)
                for r in results:
                    grade_counts[r["path_grade"]] += 1
                for g, c in sorted(grade_counts.items()):
                    self._log(f"   {grade_label_kr(g)}: {c}")

                if messagebox.askyesno("Completed", f"✅ Validation completed!\n\nResult file: {out_path}\n\nOpen it in the browser?"):
                    import webbrowser
                    webbrowser.open(f"file://{os.path.abspath(out_path)}")
            else:
                self._log("⚠️ No results were generated.")

        except Exception as e:
            self._log(f"❌ Error: {e}")
            import traceback
            self._log(traceback.format_exc())
        finally:
            self._running = False
            self.run_btn.config(state="normal")
            self.status_lbl.config(text="Completed")
            self.progress["value"] = 100


# ─── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = PPIValidatorApp()
    app.mainloop()