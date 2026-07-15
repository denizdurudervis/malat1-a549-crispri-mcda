#!/usr/bin/env python3
r"""
MALAT1 A549 CRISPRi — SMAA-like Monte Carlo rank acceptability analysis v0.1

Place in scripts/ and run from repository root:

    .\.venv\Scripts\python.exe scripts\run_smaa_monte_carlo_v0_1.py

Input:
    data/decision_matrix/malat1_a549_normalized_decision_matrix_v0_3.tsv

Outputs:
    data/decision_matrix/smaa/
"""
from pathlib import Path
from collections import OrderedDict
import csv
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "generated" / "normalized"
    / "malat1_a549_normalized_decision_matrix_v0_3.tsv"
)
OUT = ROOT / "generated" / "profile_smaa"
OUT.mkdir(parents=True, exist_ok=True)

N_ITER = 30000
WEIGHT_CONCENTRATION = 60.0
JOST_CONCENTRATION = 40.0
BASE_SEED = 20260715

PROFILES = OrderedDict([
    ("balanced", np.array([0.30, 0.30, 0.15, 0.25], dtype=float)),
    ("efficacy_focused", np.array([0.45, 0.20, 0.10, 0.25], dtype=float)),
    ("safety_focused", np.array([0.20, 0.35, 0.15, 0.30], dtype=float)),
])

CRITERIA = [
    "C1_position",
    "C2_ATAC",
    "C3_TALAM1_safety",
    "C4_nearby_promoter_safety",
    "C5_sequence_feasibility",
    "C6_Jost_CRISPRi_specificity",
]


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def median_rank(counts):
    return int(np.searchsorted(np.cumsum(counts), counts.sum() / 2.0) + 1)


def simulate(anchor, seed, ids, d, atac, talam_overlap, prom_overlap,
             seq, jost, one_mm):
    rng = np.random.default_rng(seed)
    n = len(ids)
    rank_counts = np.zeros((n, n), dtype=np.int32)
    strict_counts = np.zeros((n, n), dtype=np.int32)
    score_sum = np.zeros(n)
    score_sq_sum = np.zeros(n)
    winner_count = np.zeros(n, dtype=np.int32)
    central_sum = np.zeros((n, 6))
    alpha = anchor * WEIGHT_CONCENTRATION

    done = 0
    while done < N_ITER:
        batch = min(2000, N_ITER - done)
        wd = rng.dirichlet(alpha, size=batch)
        es = rng.beta(20.0, 20.0, size=batch)
        ls = rng.beta(20.0, 20.0, size=batch)

        left_core = rng.triangular(-75.0, -50.0, -25.0, size=batch)
        right_core = rng.triangular(250.0, 300.0, 350.0, size=batch)
        tss_shift = np.clip(rng.normal(0.0, 2.0, size=batch), -6.0, 6.0)
        adj = d[None, :] - tss_shift[:, None]

        c1 = np.ones((batch, n))
        left = adj < left_core[:, None]
        right = adj > right_core[:, None]
        lv = (adj + 500.0) / (left_core[:, None] + 500.0)
        rv = (500.0 - adj) / (500.0 - right_core[:, None])
        c1[left] = lv[left]
        c1[right] = rv[right]
        c1 = np.clip(c1, 0.0, 1.0)

        reliability = rng.beta(19.0, 1.0, size=batch)
        c2 = np.where(
            atac[None, :] > 0.5,
            reliability[:, None],
            1.0 - reliability[:, None],
        )

        talam_penalty = rng.triangular(0.10, 0.30, 0.60, size=batch)
        prom_penalty = rng.triangular(0.40, 0.70, 0.95, size=batch)
        c3 = 1.0 - talam_penalty[:, None] * talam_overlap[None, :]
        c4 = 1.0 - prom_penalty[:, None] * prom_overlap[None, :]
        c5 = np.broadcast_to(seq[None, :], (batch, n))

        aj = np.maximum(jost * JOST_CONCENTRATION, 0.01)
        bj = np.maximum((1.0 - jost) * JOST_CONCENTRATION, 0.01)
        c6 = rng.beta(aj[None, :], bj[None, :], size=(batch, n))

        eff = es[:, None] * c1 + (1.0 - es[:, None]) * c2
        local = ls[:, None] * c3 + (1.0 - ls[:, None]) * c4
        scores = (
            wd[:, 0, None] * eff
            + wd[:, 1, None] * local
            + wd[:, 2, None] * c5
            + wd[:, 3, None] * c6
        )

        order = np.argsort(-scores, axis=1, kind="stable")
        strict_scores = scores.copy()
        strict_scores[:, one_mm] = -np.inf
        strict_order = np.argsort(-strict_scores, axis=1, kind="stable")

        for rank_index in range(n):
            np.add.at(rank_counts[:, rank_index], order[:, rank_index], 1)
            np.add.at(
                strict_counts[:, rank_index],
                strict_order[:, rank_index],
                1,
            )

        score_sum += scores.sum(axis=0)
        score_sq_sum += np.square(scores).sum(axis=0)

        criterion_weights = np.column_stack([
            wd[:, 0] * es,
            wd[:, 0] * (1.0 - es),
            wd[:, 1] * ls,
            wd[:, 1] * (1.0 - ls),
            wd[:, 2],
            wd[:, 3],
        ])
        winners = order[:, 0]
        np.add.at(winner_count, winners, 1)
        for j in range(6):
            np.add.at(central_sum[:, j], winners, criterion_weights[:, j])

        done += batch

    mean_score = score_sum / N_ITER
    score_sd = np.sqrt(
        np.maximum(score_sq_sum / N_ITER - np.square(mean_score), 0.0)
    )
    central = np.full((n, 6), np.nan)
    mask = winner_count > 0
    central[mask] = central_sum[mask] / winner_count[mask, None]
    return rank_counts, strict_counts, mean_score, score_sd, winner_count, central


def main():
    rows = read_tsv(INPUT)
    if len(rows) != 86:
        raise RuntimeError("Expected 86 candidates.")

    ids = np.array([r["candidate_id"] for r in rows], dtype=object)
    d = np.array([float(r["midpoint_distance_to_tss_bp"]) for r in rows])
    atac = np.array([float(r["C2_A549_ATAC_support"]) for r in rows])
    talam = 1.0 - np.array([float(r["C3_TALAM1_safety"]) for r in rows])
    prom = 1.0 - np.array(
        [float(r["C4_nearby_lncRNA_promoter_safety"]) for r in rows]
    )
    seq = np.array([float(r["C5_sequence_feasibility"]) for r in rows])
    jost = np.array(
        [float(r["C6_Jost_CRISPRi_specificity"]) for r in rows]
    )
    one_mm_count = np.array([int(r["one_mismatch_hits"]) for r in rows])
    two_mm_count = np.array([int(r["two_mismatch_hits"]) for r in rows])
    one_mm = one_mm_count > 0

    profile_rows = []
    rank_rows = []
    central_rows = []
    by_candidate = {candidate_id: {} for candidate_id in ids}

    for i, (profile, anchor) in enumerate(PROFILES.items()):
        counts, strict, mean_score, score_sd, wins, central = simulate(
            anchor, BASE_SEED + i, ids, d, atac, talam, prom, seq, jost, one_mm
        )
        for idx, candidate_id in enumerate(ids):
            p1 = counts[idx, 0] / N_ITER
            p5 = counts[idx, :5].sum() / N_ITER
            mean_rank = (
                counts[idx] * np.arange(1, len(ids) + 1)
            ).sum() / N_ITER
            strict_p1 = strict[idx, 0] / N_ITER
            strict_p5 = strict[idx, :5].sum() / N_ITER
            strict_mean = (
                strict[idx] * np.arange(1, len(ids) + 1)
            ).sum() / N_ITER

            row = OrderedDict([
                ("profile", profile),
                ("candidate_id", candidate_id),
                ("p_rank_1", round(float(p1), 8)),
                ("p_top_5", round(float(p5), 8)),
                ("mean_rank", round(float(mean_rank), 6)),
                ("median_rank", median_rank(counts[idx])),
                ("mean_utility", round(float(mean_score[idx]), 8)),
                ("utility_sd", round(float(score_sd[idx]), 8)),
                ("strict_no_1mm_p_rank_1", round(float(strict_p1), 8)),
                ("strict_no_1mm_p_top_5", round(float(strict_p5), 8)),
                ("strict_no_1mm_mean_rank", round(float(strict_mean), 6)),
                ("one_mismatch_hits", int(one_mm_count[idx])),
                ("two_mismatch_hits", int(two_mm_count[idx])),
                ("sequence_feasibility_pass", int(seq[idx])),
            ])
            profile_rows.append(row)
            by_candidate[candidate_id][profile] = row

            for rank_index in range(len(ids)):
                rank_rows.append(OrderedDict([
                    ("profile", profile),
                    ("candidate_id", candidate_id),
                    ("rank", rank_index + 1),
                    (
                        "rank_acceptability",
                        round(counts[idx, rank_index] / N_ITER, 8),
                    ),
                ]))

            cw = central[idx]
            central_row = OrderedDict([
                ("profile", profile),
                ("candidate_id", candidate_id),
                ("first_rank_count", int(wins[idx])),
                ("p_rank_1", round(float(p1), 8)),
            ])
            for j, name in enumerate(CRITERIA):
                central_row["central_weight_" + name] = (
                    "" if np.isnan(cw[j]) else round(float(cw[j]), 8)
                )
            central_rows.append(central_row)

    overall_rows = []
    for candidate_id in ids:
        p = by_candidate[candidate_id]
        p1 = np.array([p[name]["p_rank_1"] for name in PROFILES])
        p5 = np.array([p[name]["p_top_5"] for name in PROFILES])
        mr = np.array([p[name]["mean_rank"] for name in PROFILES])
        idx = int(np.flatnonzero(ids == candidate_id)[0])
        overall_rows.append(OrderedDict([
            ("candidate_id", candidate_id),
            ("mean_p_rank_1_across_profiles", round(float(p1.mean()), 8)),
            ("max_p_rank_1_across_profiles", round(float(p1.max()), 8)),
            ("mean_p_top_5_across_profiles", round(float(p5.mean()), 8)),
            ("minimum_p_top_5_across_profiles", round(float(p5.min()), 8)),
            ("mean_rank_across_profiles", round(float(mr.mean()), 6)),
            ("worst_mean_rank_across_profiles", round(float(mr.max()), 6)),
            ("balanced_p_top_5", p["balanced"]["p_top_5"]),
            ("efficacy_p_top_5", p["efficacy_focused"]["p_top_5"]),
            ("safety_p_top_5", p["safety_focused"]["p_top_5"]),
            ("one_mismatch_hits", int(one_mm_count[idx])),
            ("two_mismatch_hits", int(two_mm_count[idx])),
            ("sequence_feasibility_pass", int(seq[idx])),
            ("jost_crispri_specificity", round(float(jost[idx]), 8)),
            ("midpoint_distance_to_tss_bp", float(d[idx])),
            ("ATAC_support", int(atac[idx])),
            ("TALAM1_safety", int(1 - talam[idx])),
            ("nearby_promoter_safety", int(1 - prom[idx])),
        ]))

    overall_rows.sort(
        key=lambda r: (
            -float(r["minimum_p_top_5_across_profiles"]),
            -float(r["mean_p_top_5_across_profiles"]),
            float(r["mean_rank_across_profiles"]),
        )
    )
    for i, row in enumerate(overall_rows, 1):
        row["robustness_order"] = i
        row.move_to_end("robustness_order", last=False)

    write_tsv(OUT / "malat1_a549_smaa_profile_summary_v0_1.tsv", profile_rows)
    write_tsv(
        OUT / "malat1_a549_smaa_rank_acceptability_v0_1.tsv", rank_rows
    )
    write_tsv(
        OUT / "malat1_a549_smaa_central_weights_v0_1.tsv", central_rows
    )
    write_tsv(
        OUT / "malat1_a549_smaa_overall_robustness_v0_1.tsv",
        overall_rows,
    )

    manifest = {
        "candidate_count": len(ids),
        "iterations_per_profile": N_ITER,
        "profiles": list(PROFILES.keys()),
        "seed_base": BASE_SEED,
        "maximum_binomial_monte_carlo_se": math.sqrt(0.25 / N_ITER),
        "status": "complete",
    }
    (OUT / "malat1_a549_smaa_run_manifest_v0_1.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print("SMAA outputs written to {}".format(OUT))


if __name__ == "__main__":
    raise SystemExit(main())
