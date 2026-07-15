#!/usr/bin/env python3
r"""
Build the normalized MALAT1 A549 decision matrix v0.3 from the raw v0.2 TSV.

Run from repository root:
    .\.venv\Scripts\python.exe scripts\build_normalized_decision_matrix_v0_3.py
"""
from pathlib import Path
import csv
from collections import OrderedDict

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "inputs" / "malat1_a549_raw_decision_matrix_v0_2.tsv"
OUTPUT = ROOT / "generated" / "normalized" / "malat1_a549_normalized_decision_matrix_v0_3.tsv"

PROFILES = OrderedDict([
    ("balanced", (0.30, 0.30, 0.15, 0.25)),
    ("efficacy_focused", (0.45, 0.20, 0.10, 0.25)),
    ("safety_focused", (0.20, 0.35, 0.15, 0.30)),
])

def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def b(v):
    return str(v).strip().lower() in {"true", "1", "yes"}

def clamp(v):
    return max(0.0, min(1.0, float(v)))

def position_score(d):
    d = float(d)
    if -50 <= d <= 300:
        return 1.0
    if d < -50:
        return clamp((d + 500) / 450)
    return clamp((500 - d) / 200)

def ranks_desc(values):
    order = sorted(enumerate(values), key=lambda x: (-x[1], x[0]))
    ranks = [0] * len(values)
    last = None
    rank = 0
    for pos, (idx, value) in enumerate(order, 1):
        if last is None or value != last:
            rank = pos
            last = value
        ranks[idx] = rank
    return ranks

def main():
    rows = read_tsv(INPUT)
    if len(rows) != 86:
        raise RuntimeError("Expected 86 candidates.")
    out = []
    for raw in rows:
        c1 = position_score(raw["midpoint_distance_to_tss_bp"])
        c2 = int(b(raw["overlaps_atac_peak"]))
        c3 = int(not b(raw["overlaps_talam1_annotation"]))
        c4 = int(not b(raw["overlaps_nearby_lncRNA_promoter_window"]))
        c5 = int(b(raw["local_sequence_qc_pass"]))
        c6 = clamp(raw["jost_crispri_specificity"])
        eff = (c1 + c2) / 2
        local = (c3 + c4) / 2
        seq = float(c5)
        off = c6
        row = OrderedDict(raw)
        row.update({
            "C1_position_soft_support": round(c1, 6),
            "position_core_minus50_plus300_flag": int(
                b(raw["midpoint_in_minus50_plus300"])
            ),
            "C2_A549_ATAC_support": c2,
            "C3_TALAM1_safety": c3,
            "C4_nearby_lncRNA_promoter_safety": c4,
            "C5_sequence_feasibility": c5,
            "C6_Jost_CRISPRi_specificity": round(c6, 6),
            "efficacy_domain_score": round(eff, 6),
            "local_safety_domain_score": round(local, 6),
            "sequence_domain_score": round(seq, 6),
            "whole_genome_safety_domain_score": round(off, 6),
            "one_mismatch_safeguard_pass": int(
                int(raw["one_mismatch_hits"]) == 0
            ),
            "normalization_status":
                "normalized_hierarchical_profiles_exploratory_not_SMAA",
        })
        for name, (we, wl, ws, wo) in PROFILES.items():
            row[name + "_score"] = round(
                we * eff + wl * local + ws * seq + wo * off, 8
            )
        out.append(row)
    for name in PROFILES:
        rr = ranks_desc([float(x[name + "_score"]) for x in out])
        for row, rank in zip(out, rr):
            row[name + "_rank"] = rank
    write_tsv(OUTPUT, out)
    print("Wrote {} candidates -> {}".format(len(out), OUTPUT))

if __name__ == "__main__":
    raise SystemExit(main())
