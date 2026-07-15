#!/usr/bin/env python3
"""Compute SMAA-2-type confidence-factor analogues at profile-specific central weights."""
from pathlib import Path
import csv
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NORM = ROOT / "generated" / "normalized" / "malat1_a549_normalized_decision_matrix_v0_3.tsv"
CENTRAL = ROOT / "generated" / "profile_smaa" / "malat1_a549_smaa_central_weights_v0_1.tsv"
OUT = ROOT / "generated" / "profile_smaa" / "malat1_a549_smaa_confidence_factors_v0_1.tsv"
N_ITER = 20000
BASE_SEED = 20261215

WEIGHT_COLS = [
    "central_weight_C1_position",
    "central_weight_C2_ATAC",
    "central_weight_C3_TALAM1_safety",
    "central_weight_C4_nearby_promoter_safety",
    "central_weight_C5_sequence_feasibility",
    "central_weight_C6_Jost_CRISPRi_specificity",
]

def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

def main():
    rows = read_tsv(NORM)
    central = read_tsv(CENTRAL)
    ids = np.array([r["candidate_id"] for r in rows], dtype=object)
    id_to_i = {cid:i for i,cid in enumerate(ids)}
    d = np.array([float(r["midpoint_distance_to_tss_bp"]) for r in rows])
    atac = np.array([float(r["C2_A549_ATAC_support"]) for r in rows])
    talam_overlap = 1.0 - np.array([float(r["C3_TALAM1_safety"]) for r in rows])
    prom_overlap = 1.0 - np.array([float(r["C4_nearby_lncRNA_promoter_safety"]) for r in rows])
    seq = np.array([float(r["C5_sequence_feasibility"]) for r in rows])
    jost = np.array([float(r["C6_Jost_CRISPRi_specificity"]) for r in rows])
    profiles = sorted(set(r["profile"] for r in central))
    output = []
    for pidx, profile in enumerate(profiles):
        eligible = [r for r in central if r["profile"] == profile and r[WEIGHT_COLS[0]] != ""]
        rng = np.random.default_rng(BASE_SEED + pidx)
        rank1 = {r["candidate_id"]:0 for r in eligible}
        top5 = {r["candidate_id"]:0 for r in eligible}
        done = 0
        while done < N_ITER:
            batch = min(1000, N_ITER-done)
            left_core = rng.triangular(-75.0, -50.0, -25.0, size=batch)
            right_core = rng.triangular(250.0, 300.0, 350.0, size=batch)
            shift = np.clip(rng.normal(0.0, 2.0, size=batch), -6.0, 6.0)
            adj = d[None,:] - shift[:,None]
            c1 = np.ones((batch,len(ids)))
            left = adj < left_core[:,None]; right = adj > right_core[:,None]
            lv = (adj+500.0)/(left_core[:,None]+500.0)
            rv = (500.0-adj)/(500.0-right_core[:,None])
            c1[left]=lv[left]; c1[right]=rv[right]; c1=np.clip(c1,0,1)
            rel = rng.beta(19.0,1.0,size=batch)
            c2 = np.where(atac[None,:]>.5, rel[:,None], 1-rel[:,None])
            tp = rng.triangular(.10,.30,.60,size=batch)
            pp = rng.triangular(.40,.70,.95,size=batch)
            c3 = 1-tp[:,None]*talam_overlap[None,:]
            c4 = 1-pp[:,None]*prom_overlap[None,:]
            c5 = np.broadcast_to(seq[None,:], (batch,len(ids)))
            a=np.maximum(jost*40,.01); b=np.maximum((1-jost)*40,.01)
            c6=rng.beta(a[None,:],b[None,:],size=(batch,len(ids)))
            criteria=(c1,c2,c3,c4,c5,c6)
            for record in eligible:
                cid=record["candidate_id"]; idx=id_to_i[cid]
                w=np.array([float(record[c]) for c in WEIGHT_COLS])
                if not np.isclose(w.sum(),1.0,atol=1e-6):
                    w=w/w.sum()
                scores=sum(w[k]*criteria[k] for k in range(6))
                target=scores[:,idx]
                better=(scores > target[:,None]).sum(axis=1)
                rank1[cid] += int(np.sum(better==0))
                top5[cid] += int(np.sum(better<5))
            done += batch
        for record in eligible:
            cid=record["candidate_id"]
            row={
                "profile":profile,
                "candidate_id":cid,
                "original_p_rank_1":record["p_rank_1"],
                "central_weight_confidence_factor_rank1":round(rank1[cid]/N_ITER,8),
                "central_weight_top5_probability":round(top5[cid]/N_ITER,8),
                "iterations":N_ITER,
                "interpretation":"probability candidate ranks first/top5 under criterion uncertainty with its profile-specific central criterion-weight vector fixed",
            }
            for c in WEIGHT_COLS: row[c]=record[c]
            output.append(row)
    write_tsv(OUT, output)
    print(f"Wrote {len(output)} confidence-factor rows -> {OUT}")

if __name__ == "__main__":
    main()
