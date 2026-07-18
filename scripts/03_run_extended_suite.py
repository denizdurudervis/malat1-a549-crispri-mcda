# Extended reproducibility and sensitivity analyses.

from pathlib import Path
from collections import OrderedDict, defaultdict
import csv
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE / "generated" / "extended"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
DOCS = ROOT / "docs"
SCRIPTS = ROOT / "scripts"
for p in [ROOT, TABLES, FIGURES, DOCS, SCRIPTS]:
    p.mkdir(parents=True, exist_ok=True)

RAW_PATH = BASE / "inputs" / "malat1_a549_raw_decision_matrix_v0_2.tsv"
NORM_PATH = BASE / "generated" / "normalized" / "malat1_a549_normalized_decision_matrix_v0_3.tsv"
SMAA_OVERALL_PATH = BASE / "generated" / "profile_smaa" / "malat1_a549_smaa_overall_robustness_v0_1.tsv"
SMAA_PROFILE_PATH = BASE / "generated" / "profile_smaa" / "malat1_a549_smaa_profile_summary_v0_1.tsv"
PANEL_PATH = BASE / "inputs" / "malat1_a549_provisional_validation_panel_v0_1.tsv"
STOJIC_PATH = BASE / "inputs" / "stojic2018_guides_grch38_corrected.tsv"
CRINCL_PATH = BASE / "inputs" / "crincl_malat1_unique_guides_grch38_v2.tsv"

SEED = 20260715
rng_master = np.random.default_rng(SEED)

PROFILES = OrderedDict([
    ("balanced", np.array([0.30, 0.30, 0.15, 0.25], float)),
    ("efficacy_focused", np.array([0.45, 0.20, 0.10, 0.25], float)),
    ("safety_focused", np.array([0.20, 0.35, 0.15, 0.30], float)),
])


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)


def save_png(fig, path, dpi=300):
    """Write a PNG and fail immediately if the render is incomplete."""
    path.unlink(missing_ok=True)
    fig.savefig(path, dpi=dpi)
    if path.stat().st_size <= 100 or path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"Invalid PNG render: {path}")


def as_bool(v):
    return str(v).strip().lower() in {"true", "1", "yes"}


def clamp01(x):
    return np.clip(x, 0.0, 1.0)


def ranks_desc_1d(values):
    values = np.asarray(values, float)
    order = np.argsort(-values, kind="stable")
    ranks = np.empty(len(values), dtype=int)
    ranks[order] = np.arange(1, len(values) + 1)
    return ranks


def average_ranks(values):
    values = np.asarray(values, float)
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg
        i = j
    return ranks


def spearman(a, b):
    ra, rb = average_ranks(a), average_ranks(b)
    if np.std(ra) == 0 or np.std(rb) == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else 1.0


def percentile_rank(values):
    r = average_ranks(values)
    return (r - 1.0) / max(len(values) - 1, 1)


raw = read_tsv(RAW_PATH)
norm = read_tsv(NORM_PATH)
smaa_overall = read_tsv(SMAA_OVERALL_PATH)
smaa_profile = read_tsv(SMAA_PROFILE_PATH)
panel = read_tsv(PANEL_PATH)
stojic = read_tsv(STOJIC_PATH)
crincl = read_tsv(CRINCL_PATH)

if len(raw) != 86 or len(norm) != 86:
    raise RuntimeError("Expected 86 candidates in raw and normalized matrices.")

ids = np.array([r["candidate_id"] for r in norm], dtype=object)
n = len(ids)
id_to_i = {cid: i for i, cid in enumerate(ids)}
raw_by_id = {r["candidate_id"]: r for r in raw}
norm_by_id = {r["candidate_id"]: r for r in norm}
panel_ids = [r["candidate_id"] for r in panel]

d = np.array([float(r["midpoint_distance_to_tss_bp"]) for r in norm])
atac = np.array([float(r["C2_A549_ATAC_support"]) for r in norm])
talam_safety = np.array([float(r["C3_TALAM1_safety"]) for r in norm])
prom_safety = np.array([float(r["C4_nearby_lncRNA_promoter_safety"]) for r in norm])
seq_locked = np.array([float(r["C5_sequence_feasibility"]) for r in norm])
jost = np.array([float(r["C6_Jost_CRISPRi_specificity"]) for r in norm])
cfd = np.array([float(r["doench_cfd_specificity"]) for r in norm])
hsu = np.array([float(r["hsu2013_specificity"]) for r in norm])
one_mm = np.array([int(float(r["one_mismatch_hits"])) for r in norm])
two_mm = np.array([int(float(r["two_mismatch_hits"])) for r in norm])
gc = np.array([float(r["gc_percent"]) for r in norm])
poly_t = np.array([as_bool(r["has_poly_t_4"]) for r in norm])
homopoly = np.array([int(float(r["max_homopolymer_run"])) for r in norm])

c1_locked = np.array([float(r["C1_position_soft_support"]) for r in norm])
eff_locked = (c1_locked + atac) / 2
local_locked = (talam_safety + prom_safety) / 2


# ---------------------------------------------------------------------
# 1. Position model-form sensitivity
# ---------------------------------------------------------------------
def trapezoid_score(x, core_left, core_right, zero_left=-500.0, zero_right=500.0):
    x = np.asarray(x, float)
    score = np.ones_like(x)
    left = x < core_left
    right = x > core_right
    score[left] = (x[left] - zero_left) / (core_left - zero_left)
    score[right] = (zero_right - x[right]) / (zero_right - core_right)
    return clamp01(score)


def triangular_score(x, peak=100.0, left=-500.0, right=500.0):
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    m1 = (x >= left) & (x <= peak)
    m2 = (x > peak) & (x <= right)
    out[m1] = (x[m1] - left) / (peak - left)
    out[m2] = (right - x[m2]) / (right - peak)
    return clamp01(out)


def gaussian_score(x, center=100.0, sigma=200.0):
    x = np.asarray(x, float)
    return np.exp(-0.5 * ((x - center) / sigma) ** 2)


position_models = OrderedDict([
    ("baseline_trapezoid", c1_locked),
    ("narrow_trapezoid", trapezoid_score(d, -25, 200, -300, 400)),
    ("broad_trapezoid", trapezoid_score(d, -100, 350, -500, 500)),
    ("triangular_peak_plus100", triangular_score(d, 100, -500, 500)),
    ("gaussian_center_plus100_sigma200", gaussian_score(d, 100, 200)),
])

position_rows = []
position_rank_map = {}
for model_name, c1 in position_models.items():
    eff = (c1 + atac) / 2
    for profile_name, w in PROFILES.items():
        utility = w[0] * eff + w[1] * local_locked + w[2] * seq_locked + w[3] * jost
        ranks = ranks_desc_1d(utility)
        position_rank_map[(model_name, profile_name)] = ranks
        for i, cid in enumerate(ids):
            position_rows.append(OrderedDict([
                ("position_model", model_name),
                ("profile", profile_name),
                ("candidate_id", cid),
                ("position_score", round(float(c1[i]), 8)),
                ("utility", round(float(utility[i]), 8)),
                ("rank", int(ranks[i])),
                ("panel_flag", int(cid in panel_ids)),
            ]))
write_tsv(TABLES / "position_model_sensitivity.tsv", position_rows)

baseline_pos_ranks = position_rank_map[("baseline_trapezoid", "balanced")]
baseline_pos_top10 = ids[np.argsort(baseline_pos_ranks)[:10]].tolist()
position_summary_rows = []
for model_name in position_models:
    for profile_name in PROFILES:
        ranks = position_rank_map[(model_name, profile_name)]
        top10 = ids[np.argsort(ranks)[:10]].tolist()
        position_summary_rows.append(OrderedDict([
            ("position_model", model_name),
            ("profile", profile_name),
            ("spearman_vs_baseline_balanced", round(spearman(-baseline_pos_ranks, -ranks), 6)),
            ("top10_jaccard_vs_baseline_balanced", round(jaccard(baseline_pos_top10, top10), 6)),
            ("top_5_candidates", ";".join(ids[np.argsort(ranks)[:5]].tolist())),
        ]))
write_tsv(TABLES / "position_model_summary.tsv", position_summary_rows)

position_candidate_summary = []
for i, cid in enumerate(ids):
    rr = np.array([position_rank_map[k][i] for k in position_rank_map], int)
    position_candidate_summary.append(OrderedDict([
        ("candidate_id", cid),
        ("minimum_rank_across_position_models_profiles", int(rr.min())),
        ("maximum_rank_across_position_models_profiles", int(rr.max())),
        ("mean_rank_across_position_models_profiles", round(float(rr.mean()), 6)),
        ("rank_range", int(rr.max() - rr.min())),
        ("panel_flag", int(cid in panel_ids)),
    ]))
position_candidate_summary.sort(key=lambda r: (r["mean_rank_across_position_models_profiles"], r["maximum_rank_across_position_models_profiles"]))
write_tsv(TABLES / "position_candidate_rank_stability.tsv", position_candidate_summary)


# ---------------------------------------------------------------------
# 2. Sequence-QC sensitivity
# ---------------------------------------------------------------------
moderate_seq = (
    (gc >= 25) & (gc <= 75) & (~poly_t) & (homopoly <= 5)
).astype(float)
poly_t_only_seq = (~poly_t).astype(float)
continuous_gc = np.ones(n, float)
continuous_gc[gc < 30] = np.clip((gc[gc < 30] - 20) / 10, 0, 1)
continuous_gc[gc > 70] = np.clip((80 - gc[gc > 70]) / 10, 0, 1)
continuous_hom = np.where(homopoly <= 4, 1.0, np.where(homopoly == 5, 0.5, 0.0))
continuous_seq = continuous_gc * (~poly_t).astype(float) * continuous_hom

sequence_scenarios = OrderedDict([
    ("locked_strict_composite", seq_locked),
    ("moderate_25_75_no_polyT_homopoly_le5", moderate_seq),
    ("polyT_only", poly_t_only_seq),
    ("continuous_penalty", continuous_seq),
    ("no_sequence_domain", np.full(n, np.nan)),
    ("hard_locked_eligibility", seq_locked),
])

sequence_rows = []
sequence_rank_map = {}
for scenario_name, seq_value in sequence_scenarios.items():
    for profile_name, w_orig in PROFILES.items():
        w = w_orig.copy()
        if scenario_name == "no_sequence_domain":
            w = w.copy()
            kept = np.array([w[0], w[1], w[3]])
            kept = kept / kept.sum()
            utility = kept[0] * eff_locked + kept[1] * local_locked + kept[2] * jost
        else:
            utility = w[0] * eff_locked + w[1] * local_locked + w[2] * seq_value + w[3] * jost
        if scenario_name == "hard_locked_eligibility":
            utility = utility.copy()
            utility[seq_locked < 0.5] = -np.inf
        ranks = ranks_desc_1d(utility)
        sequence_rank_map[(scenario_name, profile_name)] = ranks
        for i, cid in enumerate(ids):
            sequence_rows.append(OrderedDict([
                ("sequence_scenario", scenario_name),
                ("profile", profile_name),
                ("candidate_id", cid),
                ("sequence_value", "" if np.isnan(seq_value[i]) else round(float(seq_value[i]), 8)),
                ("eligible", int(not (scenario_name == "hard_locked_eligibility" and seq_locked[i] < 0.5))),
                ("utility", "" if not np.isfinite(utility[i]) else round(float(utility[i]), 8)),
                ("rank", int(ranks[i])),
                ("panel_flag", int(cid in panel_ids)),
            ]))
write_tsv(TABLES / "sequence_qc_sensitivity.tsv", sequence_rows)

baseline_seq_ranks = sequence_rank_map[("locked_strict_composite", "balanced")]
baseline_seq_top10 = ids[np.argsort(baseline_seq_ranks)[:10]].tolist()
sequence_summary_rows = []
for scenario_name in sequence_scenarios:
    for profile_name in PROFILES:
        ranks = sequence_rank_map[(scenario_name, profile_name)]
        top10 = ids[np.argsort(ranks)[:10]].tolist()
        sequence_summary_rows.append(OrderedDict([
            ("sequence_scenario", scenario_name),
            ("profile", profile_name),
            ("spearman_vs_locked_balanced", round(spearman(-baseline_seq_ranks, -ranks), 6)),
            ("top10_jaccard_vs_locked_balanced", round(jaccard(baseline_seq_top10, top10), 6)),
            ("top_5_candidates", ";".join(ids[np.argsort(ranks)[:5]].tolist())),
        ]))
write_tsv(TABLES / "sequence_qc_summary.tsv", sequence_summary_rows)

sequence_candidate_summary = []
for i, cid in enumerate(ids):
    rr = np.array([sequence_rank_map[k][i] for k in sequence_rank_map], int)
    sequence_candidate_summary.append(OrderedDict([
        ("candidate_id", cid),
        ("minimum_rank_across_sequence_scenarios_profiles", int(rr.min())),
        ("maximum_rank_across_sequence_scenarios_profiles", int(rr.max())),
        ("mean_rank_across_sequence_scenarios_profiles", round(float(rr.mean()), 6)),
        ("rank_range", int(rr.max() - rr.min())),
        ("locked_sequence_pass", int(seq_locked[i])),
        ("panel_flag", int(cid in panel_ids)),
    ]))
sequence_candidate_summary.sort(key=lambda r: (r["mean_rank_across_sequence_scenarios_profiles"], r["maximum_rank_across_sequence_scenarios_profiles"]))
write_tsv(TABLES / "sequence_candidate_rank_stability.tsv", sequence_candidate_summary)


# ---------------------------------------------------------------------
# Shared uncertainty simulator
# ---------------------------------------------------------------------
def simulate_uncertain(
    n_iter,
    seed,
    weight_mode="anchor",
    anchor=None,
    vary_weights=True,
    vary_position=True,
    vary_atac=True,
    vary_local=True,
    vary_jost=True,
    strict_no_1mm=False,
    checkpoints=None,
):
    rng = np.random.default_rng(seed)
    rank_counts = np.zeros((n, n), dtype=np.int32)
    score_sum = np.zeros(n, float)
    score_sq = np.zeros(n, float)
    checkpoint_data = {}
    completed = 0
    batch_size = 2000
    checkpoints = sorted(checkpoints or [])
    check_set = set(checkpoints)

    # Central values for fixed-source runs.
    c1_center = c1_locked
    c2_center = atac * 0.95 + (1 - atac) * 0.05
    c3_center = 1.0 - 0.30 * (1.0 - talam_safety)
    c4_center = 1.0 - 0.70 * (1.0 - prom_safety)
    c6_center = jost

    while completed < n_iter:
        next_checkpoint = min([c for c in checkpoints if c > completed], default=n_iter)
        batch = min(batch_size, n_iter - completed, next_checkpoint - completed)

        if vary_weights:
            if weight_mode == "uniform":
                wd = rng.dirichlet(np.ones(4), size=batch)
                es = rng.beta(1.0, 1.0, size=batch)
                ls = rng.beta(1.0, 1.0, size=batch)
            else:
                if anchor is None:
                    anchor = PROFILES["balanced"]
                wd = rng.dirichlet(anchor * 60.0, size=batch)
                es = rng.beta(20.0, 20.0, size=batch)
                ls = rng.beta(20.0, 20.0, size=batch)
        else:
            if anchor is None:
                anchor = PROFILES["balanced"]
            wd = np.broadcast_to(anchor[None, :], (batch, 4))
            es = np.full(batch, 0.5)
            ls = np.full(batch, 0.5)

        if vary_position:
            left_core = rng.triangular(-75.0, -50.0, -25.0, size=batch)
            right_core = rng.triangular(250.0, 300.0, 350.0, size=batch)
            tss_shift = np.clip(rng.normal(0.0, 2.0, size=batch), -6, 6)
            adj = d[None, :] - tss_shift[:, None]
            c1 = np.ones((batch, n))
            left = adj < left_core[:, None]
            right = adj > right_core[:, None]
            lv = (adj + 500) / (left_core[:, None] + 500)
            rv = (500 - adj) / (500 - right_core[:, None])
            c1[left] = lv[left]
            c1[right] = rv[right]
            c1 = clamp01(c1)
        else:
            c1 = np.broadcast_to(c1_center[None, :], (batch, n))

        if vary_atac:
            rel = rng.beta(19.0, 1.0, size=batch)
            c2 = np.where(atac[None, :] > 0.5, rel[:, None], 1 - rel[:, None])
        else:
            c2 = np.broadcast_to(c2_center[None, :], (batch, n))

        if vary_local:
            tp = rng.triangular(0.10, 0.30, 0.60, size=batch)
            pp = rng.triangular(0.40, 0.70, 0.95, size=batch)
            c3 = 1 - tp[:, None] * (1 - talam_safety)[None, :]
            c4 = 1 - pp[:, None] * (1 - prom_safety)[None, :]
        else:
            c3 = np.broadcast_to(c3_center[None, :], (batch, n))
            c4 = np.broadcast_to(c4_center[None, :], (batch, n))

        c5 = np.broadcast_to(seq_locked[None, :], (batch, n))

        if vary_jost:
            a = np.maximum(jost * 40.0, 0.01)
            b = np.maximum((1 - jost) * 40.0, 0.01)
            c6 = rng.beta(a[None, :], b[None, :], size=(batch, n))
        else:
            c6 = np.broadcast_to(c6_center[None, :], (batch, n))

        eff = es[:, None] * c1 + (1 - es[:, None]) * c2
        local = ls[:, None] * c3 + (1 - ls[:, None]) * c4
        score = (
            wd[:, 0, None] * eff
            + wd[:, 1, None] * local
            + wd[:, 2, None] * c5
            + wd[:, 3, None] * c6
        )
        if strict_no_1mm:
            score[:, one_mm > 0] = -np.inf

        order = np.argsort(-score, axis=1, kind="stable")
        for rank_idx in range(n):
            np.add.at(rank_counts[:, rank_idx], order[:, rank_idx], 1)
        finite_score = np.where(np.isfinite(score), score, 0.0)
        score_sum += finite_score.sum(axis=0)
        score_sq += np.square(finite_score).sum(axis=0)

        completed += batch
        if completed in check_set:
            checkpoint_data[completed] = rank_counts.copy()

    mean_score = score_sum / n_iter
    sd_score = np.sqrt(np.maximum(score_sq / n_iter - mean_score**2, 0))
    return rank_counts, mean_score, sd_score, checkpoint_data


# ---------------------------------------------------------------------
# 3. Global/uninformed SMAA
# ---------------------------------------------------------------------
GLOBAL_ITER = 50000
global_counts, global_mean, global_sd, _ = simulate_uncertain(
    GLOBAL_ITER,
    SEED + 100,
    weight_mode="uniform",
    vary_weights=True,
    vary_position=True,
    vary_atac=True,
    vary_local=True,
    vary_jost=True,
)
global_strict_counts, _, _, _ = simulate_uncertain(
    GLOBAL_ITER,
    SEED + 101,
    weight_mode="uniform",
    vary_weights=True,
    vary_position=True,
    vary_atac=True,
    vary_local=True,
    vary_jost=True,
    strict_no_1mm=True,
)

global_rows = []
for i, cid in enumerate(ids):
    ranks_axis = np.arange(1, n + 1)
    p1 = global_counts[i, 0] / GLOBAL_ITER
    p5 = global_counts[i, :5].sum() / GLOBAL_ITER
    mean_rank = float((global_counts[i] * ranks_axis).sum() / GLOBAL_ITER)
    strict_p5 = global_strict_counts[i, :5].sum() / GLOBAL_ITER
    global_rows.append(OrderedDict([
        ("candidate_id", cid),
        ("p_rank_1", round(float(p1), 8)),
        ("p_top_5", round(float(p5), 8)),
        ("mean_rank", round(mean_rank, 6)),
        ("mean_utility", round(float(global_mean[i]), 8)),
        ("utility_sd", round(float(global_sd[i]), 8)),
        ("strict_no_1mm_p_top_5", round(float(strict_p5), 8)),
        ("one_mismatch_hits", int(one_mm[i])),
        ("panel_flag", int(cid in panel_ids)),
    ]))
global_rows.sort(key=lambda r: (-r["p_top_5"], r["mean_rank"]))
for k, r in enumerate(global_rows, 1):
    r["global_robustness_order"] = k
    r.move_to_end("global_robustness_order", last=False)
write_tsv(TABLES / "global_weight_smaa_summary.tsv", global_rows)


# Pairwise win probability among global top 10.
global_top10_ids = [r["candidate_id"] for r in global_rows[:10]]
top_idx = np.array([id_to_i[c] for c in global_top10_ids])
PAIR_ITER = 30000
rng = np.random.default_rng(SEED + 102)
pair_wins = np.zeros((10, 10), float)
done = 0
while done < PAIR_ITER:
    batch = min(2000, PAIR_ITER - done)
    wd = rng.dirichlet(np.ones(4), size=batch)
    es = rng.beta(1, 1, size=batch)
    ls = rng.beta(1, 1, size=batch)
    left_core = rng.triangular(-75, -50, -25, size=batch)
    right_core = rng.triangular(250, 300, 350, size=batch)
    shift = np.clip(rng.normal(0, 2, size=batch), -6, 6)
    adj = d[None, :] - shift[:, None]
    c1 = np.ones((batch, n))
    left = adj < left_core[:, None]
    right = adj > right_core[:, None]
    lv = (adj + 500) / (left_core[:, None] + 500)
    rv = (500 - adj) / (500 - right_core[:, None])
    c1[left] = lv[left]
    c1[right] = rv[right]
    c1 = clamp01(c1)
    rel = rng.beta(19, 1, size=batch)
    c2 = np.where(atac[None, :] > .5, rel[:, None], 1-rel[:, None])
    tp = rng.triangular(.1, .3, .6, size=batch)
    pp = rng.triangular(.4, .7, .95, size=batch)
    c3 = 1 - tp[:, None] * (1-talam_safety)[None, :]
    c4 = 1 - pp[:, None] * (1-prom_safety)[None, :]
    a = np.maximum(jost*40, .01)
    b = np.maximum((1-jost)*40, .01)
    c6 = rng.beta(a[None, :], b[None, :], size=(batch, n))
    eff = es[:, None]*c1 + (1-es[:, None])*c2
    loc = ls[:, None]*c3 + (1-ls[:, None])*c4
    score = wd[:,0,None]*eff + wd[:,1,None]*loc + wd[:,2,None]*seq_locked + wd[:,3,None]*c6
    s = score[:, top_idx]
    for i in range(10):
        for j in range(10):
            if i != j:
                pair_wins[i,j] += np.mean(s[:,i] > s[:,j]) * batch
    done += batch
pair_wins /= PAIR_ITER
pair_rows = []
for i, cid_i in enumerate(global_top10_ids):
    row = OrderedDict([("candidate_id", cid_i)])
    for j, cid_j in enumerate(global_top10_ids):
        row[cid_j] = "" if i == j else round(float(pair_wins[i,j]), 6)
    pair_rows.append(row)
write_tsv(TABLES / "global_top10_pairwise_win_probabilities.tsv", pair_rows)


# ---------------------------------------------------------------------
# 4. Pareto front and layers
# ---------------------------------------------------------------------
domains = np.column_stack([eff_locked, local_locked, seq_locked, jost])


def nondominated_layers(values):
    remaining = set(range(values.shape[0]))
    layers = np.zeros(values.shape[0], dtype=int)
    layer = 1
    while remaining:
        front = []
        rem_list = list(remaining)
        for i in rem_list:
            dominated = False
            vi = values[i]
            for j in rem_list:
                if i == j:
                    continue
                vj = values[j]
                if np.all(vj >= vi) and np.any(vj > vi):
                    dominated = True
                    break
            if not dominated:
                front.append(i)
        for i in front:
            layers[i] = layer
            remaining.remove(i)
        layer += 1
    return layers


pareto_layer = nondominated_layers(domains)
strict_mask = one_mm == 0
strict_layers = np.full(n, np.nan)
strict_layers[strict_mask] = nondominated_layers(domains[strict_mask])

pareto_rows = []
for i, cid in enumerate(ids):
    pareto_rows.append(OrderedDict([
        ("candidate_id", cid),
        ("pareto_layer_4_domains", int(pareto_layer[i])),
        ("strict_no_1mm_pareto_layer", "" if np.isnan(strict_layers[i]) else int(strict_layers[i])),
        ("efficacy_domain", round(float(eff_locked[i]), 8)),
        ("local_safety_domain", round(float(local_locked[i]), 8)),
        ("sequence_domain", round(float(seq_locked[i]), 8)),
        ("offtarget_domain_jost", round(float(jost[i]), 8)),
        ("one_mismatch_hits", int(one_mm[i])),
        ("panel_flag", int(cid in panel_ids)),
    ]))
pareto_rows.sort(key=lambda r: (r["pareto_layer_4_domains"], -r["efficacy_domain"], -r["local_safety_domain"], -r["offtarget_domain_jost"]))
write_tsv(TABLES / "pareto_front_and_layers.tsv", pareto_rows)


# ---------------------------------------------------------------------
# 5. Leave-one-domain-out
# ---------------------------------------------------------------------
baseline_utility = (
    PROFILES["balanced"][0]*eff_locked
    + PROFILES["balanced"][1]*local_locked
    + PROFILES["balanced"][2]*seq_locked
    + PROFILES["balanced"][3]*jost
)
baseline_rank = ranks_desc_1d(baseline_utility)
baseline_top10 = ids[np.argsort(baseline_rank)[:10]].tolist()

loo_defs = OrderedDict([
    ("baseline_all_domains", [0,1,2,3]),
    ("without_efficacy", [1,2,3]),
    ("without_local_safety", [0,2,3]),
    ("without_sequence", [0,1,3]),
    ("without_offtarget", [0,1,2]),
])
base_w = PROFILES["balanced"]
domain_matrix = np.column_stack([eff_locked, local_locked, seq_locked, jost])
loo_rows = []
loo_summary = []
for name, keep in loo_defs.items():
    w = base_w[keep]
    w = w / w.sum()
    utility = domain_matrix[:,keep] @ w
    ranks = ranks_desc_1d(utility)
    top10 = ids[np.argsort(ranks)[:10]].tolist()
    loo_summary.append(OrderedDict([
        ("scenario", name),
        ("spearman_vs_baseline", round(spearman(-baseline_rank, -ranks), 6)),
        ("top10_jaccard_vs_baseline", round(jaccard(baseline_top10, top10), 6)),
        ("top_5_candidates", ";".join(ids[np.argsort(ranks)[:5]].tolist())),
    ]))
    for i,cid in enumerate(ids):
        loo_rows.append(OrderedDict([
            ("scenario", name),
            ("candidate_id", cid),
            ("utility", round(float(utility[i]), 8)),
            ("rank", int(ranks[i])),
            ("rank_change_vs_baseline", int(ranks[i] - baseline_rank[i])),
            ("panel_flag", int(cid in panel_ids)),
        ]))
write_tsv(TABLES / "leave_one_domain_out.tsv", loo_rows)
write_tsv(TABLES / "leave_one_domain_out_summary.tsv", loo_summary)


# ---------------------------------------------------------------------
# 6. Off-target model substitution and ensemble
# ---------------------------------------------------------------------
off_models = OrderedDict([
    ("Jost_CRISPRi_percentile", percentile_rank(jost)),
    ("Hsu2013_percentile", percentile_rank(hsu)),
    ("Doench_CFD_percentile", percentile_rank(cfd)),
])
off_models["three_model_mean_percentile"] = np.mean(np.column_stack(list(off_models.values())), axis=1)

off_rows = []
off_summary = []
off_rank_map = {}
for model_name, off_value in off_models.items():
    for profile_name,w in PROFILES.items():
        utility = w[0]*eff_locked + w[1]*local_locked + w[2]*seq_locked + w[3]*off_value
        ranks = ranks_desc_1d(utility)
        off_rank_map[(model_name, profile_name)] = ranks
        top10 = ids[np.argsort(ranks)[:10]].tolist()
        reference = off_rank_map.get(("Jost_CRISPRi_percentile", profile_name), ranks)
        ref_top10 = ids[np.argsort(reference)[:10]].tolist()
        off_summary.append(OrderedDict([
            ("offtarget_model", model_name),
            ("profile", profile_name),
            ("spearman_vs_Jost_same_profile", round(spearman(-reference, -ranks), 6)),
            ("top10_jaccard_vs_Jost_same_profile", round(jaccard(ref_top10, top10), 6)),
            ("top_5_candidates", ";".join(ids[np.argsort(ranks)[:5]].tolist())),
        ]))
        for i,cid in enumerate(ids):
            off_rows.append(OrderedDict([
                ("offtarget_model", model_name),
                ("profile", profile_name),
                ("candidate_id", cid),
                ("offtarget_value_percentile", round(float(off_value[i]), 8)),
                ("utility", round(float(utility[i]), 8)),
                ("rank", int(ranks[i])),
                ("panel_flag", int(cid in panel_ids)),
            ]))
write_tsv(TABLES / "offtarget_model_sensitivity.tsv", off_rows)
write_tsv(TABLES / "offtarget_model_summary.tsv", off_summary)

off_candidate = []
for i,cid in enumerate(ids):
    rr = np.array([off_rank_map[k][i] for k in off_rank_map])
    off_candidate.append(OrderedDict([
        ("candidate_id", cid),
        ("minimum_rank_across_models_profiles", int(rr.min())),
        ("maximum_rank_across_models_profiles", int(rr.max())),
        ("mean_rank_across_models_profiles", round(float(rr.mean()), 6)),
        ("rank_range", int(rr.max()-rr.min())),
        ("panel_flag", int(cid in panel_ids)),
    ]))
off_candidate.sort(key=lambda r: (r["mean_rank_across_models_profiles"], r["maximum_rank_across_models_profiles"]))
write_tsv(TABLES / "offtarget_candidate_rank_stability.tsv", off_candidate)


# ---------------------------------------------------------------------
# 7. Classical MDS
# ---------------------------------------------------------------------
X = np.column_stack([c1_locked, atac, talam_safety, prom_safety, seq_locked, jost])
D = np.sqrt(((X[:,None,:] - X[None,:,:])**2).sum(axis=2))
J = np.eye(n) - np.ones((n,n))/n
B = -0.5 * J @ (D**2) @ J
eigvals, eigvecs = np.linalg.eigh(B)
idx = np.argsort(eigvals)[::-1]
eigvals = eigvals[idx]
eigvecs = eigvecs[:,idx]
positive = eigvals > 1e-12
coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2],0))[None,:]
D2 = np.sqrt(((coords[:,None,:]-coords[None,:,:])**2).sum(axis=2))
upper = np.triu_indices(n,1)
stress1 = math.sqrt(np.sum((D[upper]-D2[upper])**2) / np.sum(D[upper]**2))
distance_r2 = float(np.corrcoef(D[upper], D2[upper])[0,1]**2)
positive_sum = eigvals[positive].sum()
explained = eigvals[:2].sum()/positive_sum if positive_sum>0 else 0

mds_rows = []
robust_by_id = {r["candidate_id"]: r for r in smaa_overall}
for i,cid in enumerate(ids):
    mds_rows.append(OrderedDict([
        ("candidate_id", cid),
        ("MDS1", round(float(coords[i,0]), 10)),
        ("MDS2", round(float(coords[i,1]), 10)),
        ("pareto_layer", int(pareto_layer[i])),
        ("mean_p_top5_across_profiles", robust_by_id[cid]["mean_p_top_5_across_profiles"]),
        ("panel_flag", int(cid in panel_ids)),
    ]))
write_tsv(TABLES / "mds_coordinates.tsv", mds_rows)
eig_rows = []
for i,val in enumerate(eigvals[:10],1):
    eig_rows.append(OrderedDict([
        ("dimension", i),
        ("eigenvalue", round(float(val), 10)),
        ("positive_variance_fraction", round(float(val/positive_sum), 10) if val>0 and positive_sum>0 else 0),
    ]))
write_tsv(TABLES / "mds_eigenvalues.tsv", eig_rows)
mds_fit = {
    "criteria": ["C1_position","C2_ATAC","C3_TALAM1_safety","C4_promoter_safety","C5_sequence","C6_Jost"],
    "stress_1": stress1,
    "distance_R2": distance_r2,
    "positive_eigenvalue_variance_explained_by_2D": explained,
}
(DOCS / "mds_fit_summary.json").write_text(json.dumps(mds_fit, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------
# 8. One-at-a-time uncertainty-source audit
# ---------------------------------------------------------------------
OAT_ITER = 20000
source_configs = OrderedDict([
    ("weights_only", dict(vary_weights=True, vary_position=False, vary_atac=False, vary_local=False, vary_jost=False)),
    ("position_only", dict(vary_weights=False, vary_position=True, vary_atac=False, vary_local=False, vary_jost=False)),
    ("ATAC_only", dict(vary_weights=False, vary_position=False, vary_atac=True, vary_local=False, vary_jost=False)),
    ("local_context_only", dict(vary_weights=False, vary_position=False, vary_atac=False, vary_local=True, vary_jost=False)),
    ("Jost_only", dict(vary_weights=False, vary_position=False, vary_atac=False, vary_local=False, vary_jost=True)),
    ("all_sources_combined", dict(vary_weights=True, vary_position=True, vary_atac=True, vary_local=True, vary_jost=True)),
])

# Baseline deterministic central score/rank matching simulator centers.
c2_center = atac*.95 + (1-atac)*.05
c3_center = 1 - .30*(1-talam_safety)
c4_center = 1 - .70*(1-prom_safety)
eff_center = .5*c1_locked + .5*c2_center
local_center = .5*c3_center + .5*c4_center
baseline_center_score = .30*eff_center + .30*local_center + .15*seq_locked + .25*jost
baseline_center_rank = ranks_desc_1d(baseline_center_score)
baseline_center_top5 = ids[np.argsort(baseline_center_rank)[:5]].tolist()

oat_candidate_rows = []
oat_summary_rows = []
for k,(source, cfg) in enumerate(source_configs.items()):
    counts, mean_s, sd_s, _ = simulate_uncertain(
        OAT_ITER,
        SEED + 200 + k,
        anchor=PROFILES["balanced"],
        weight_mode="anchor",
        **cfg,
    )
    mean_ranks = (counts * np.arange(1,n+1)[None,:]).sum(axis=1)/OAT_ITER
    ptop5 = counts[:,:5].sum(axis=1)/OAT_ITER
    prank1 = counts[:,0]/OAT_ITER
    winner_prob = prank1[prank1>0]
    entropy = -float(np.sum(winner_prob*np.log(winner_prob))) if len(winner_prob) else 0
    expected_jaccard = 0.0
    # Approximate expected overlap of random top-5 with deterministic top-5:
    expected_overlap = float(np.sum(ptop5[[id_to_i[c] for c in baseline_center_top5]]))
    expected_jaccard = expected_overlap / (10 - expected_overlap) if expected_overlap < 10 else 1.0
    oat_summary_rows.append(OrderedDict([
        ("uncertainty_source", source),
        ("mean_absolute_rank_shift_vs_central", round(float(np.mean(np.abs(mean_ranks-baseline_center_rank))), 6)),
        ("maximum_absolute_rank_shift_vs_central", round(float(np.max(np.abs(mean_ranks-baseline_center_rank))), 6)),
        ("mean_score_sd", round(float(np.mean(sd_s)), 8)),
        ("maximum_score_sd", round(float(np.max(sd_s)), 8)),
        ("winner_entropy_nats", round(entropy, 8)),
        ("candidates_with_p_top5_ge_0_10", int(np.sum(ptop5>=.10))),
        ("expected_top5_jaccard_vs_central_approx", round(expected_jaccard, 6)),
    ]))
    for i,cid in enumerate(ids):
        oat_candidate_rows.append(OrderedDict([
            ("uncertainty_source", source),
            ("candidate_id", cid),
            ("p_rank_1", round(float(prank1[i]), 8)),
            ("p_top_5", round(float(ptop5[i]), 8)),
            ("mean_rank", round(float(mean_ranks[i]), 6)),
            ("mean_rank_shift_vs_central", round(float(mean_ranks[i]-baseline_center_rank[i]), 6)),
            ("utility_sd", round(float(sd_s[i]), 8)),
            ("panel_flag", int(cid in panel_ids)),
        ]))
write_tsv(TABLES / "uncertainty_source_oat_candidate.tsv", oat_candidate_rows)
write_tsv(TABLES / "uncertainty_source_oat_summary.tsv", oat_summary_rows)


# ---------------------------------------------------------------------
# 9. Monte Carlo convergence and seed stability
# ---------------------------------------------------------------------
checkpoints = [1000, 5000, 10000, 20000, 30000]
counts, _, _, cp_data = simulate_uncertain(
    30000,
    SEED + 300,
    anchor=PROFILES["balanced"],
    vary_weights=True,
    vary_position=True,
    vary_atac=True,
    vary_local=True,
    vary_jost=True,
    checkpoints=checkpoints,
)
final_p5 = counts[:,:5].sum(axis=1)/30000
final_p1 = counts[:,0]/30000
final_top10 = ids[np.argsort(-final_p5)[:10]].tolist()
conv_rows = []
for cp in checkpoints:
    cc = cp_data[cp]
    p5 = cc[:,:5].sum(axis=1)/cp
    p1 = cc[:,0]/cp
    top10 = ids[np.argsort(-p5)[:10]].tolist()
    conv_rows.append(OrderedDict([
        ("iterations", cp),
        ("max_abs_p_top5_difference_vs_30000_all_candidates", round(float(np.max(np.abs(p5-final_p5))), 8)),
        ("mean_abs_p_top5_difference_vs_30000_all_candidates", round(float(np.mean(np.abs(p5-final_p5))), 8)),
        ("max_abs_p_rank1_difference_vs_30000_all_candidates", round(float(np.max(np.abs(p1-final_p1))), 8)),
        ("top10_jaccard_vs_30000", round(jaccard(top10, final_top10), 6)),
        ("top_5_candidates", ";".join(ids[np.argsort(-p5)[:5]].tolist())),
    ]))
write_tsv(TABLES / "monte_carlo_convergence.tsv", conv_rows)

seed_rows = []
seed_result_p5 = []
SEED_N = 5
for s in range(SEED_N):
    cc, _, _, _ = simulate_uncertain(
        10000,
        SEED + 400 + s,
        anchor=PROFILES["balanced"],
        vary_weights=True,
        vary_position=True,
        vary_atac=True,
        vary_local=True,
        vary_jost=True,
    )
    p5 = cc[:,:5].sum(axis=1)/10000
    seed_result_p5.append(p5)
seed_result_p5 = np.array(seed_result_p5)
seed_mean = seed_result_p5.mean(axis=0)
seed_sd = seed_result_p5.std(axis=0, ddof=1)
for i,cid in enumerate(ids):
    seed_rows.append(OrderedDict([
        ("candidate_id", cid),
        ("mean_p_top5_across_5_seeds", round(float(seed_mean[i]), 8)),
        ("sd_p_top5_across_5_seeds", round(float(seed_sd[i]), 8)),
        ("min_p_top5_across_5_seeds", round(float(seed_result_p5[:,i].min()), 8)),
        ("max_p_top5_across_5_seeds", round(float(seed_result_p5[:,i].max()), 8)),
        ("panel_flag", int(cid in panel_ids)),
    ]))
seed_rows.sort(key=lambda r: -r["mean_p_top5_across_5_seeds"])
write_tsv(TABLES / "random_seed_stability.tsv", seed_rows)

seed_summary_rows = []
ref = seed_result_p5[0]
ref_top10 = ids[np.argsort(-ref)[:10]].tolist()
for s in range(SEED_N):
    p5 = seed_result_p5[s]
    top10 = ids[np.argsort(-p5)[:10]].tolist()
    seed_summary_rows.append(OrderedDict([
        ("seed_index", s+1),
        ("spearman_p_top5_vs_seed1", round(spearman(ref,p5), 6)),
        ("top10_jaccard_vs_seed1", round(jaccard(ref_top10,top10), 6)),
        ("top_5_candidates", ";".join(ids[np.argsort(-p5)[:5]].tolist())),
    ]))
write_tsv(TABLES / "random_seed_summary.tsv", seed_summary_rows)


# ---------------------------------------------------------------------
# 10. Historical locus-context comparison
# ---------------------------------------------------------------------
historical_rows = []
for cid in panel_ids:
    r = raw_by_id[cid]
    historical_rows.append(OrderedDict([
        ("source_group", "Current_A549_panel"),
        ("guide_id", cid),
        ("guide_sequence", r["protospacer_sequence_20nt"]),
        ("cell_or_library_context", "A549;ZIM3-KRAB-dCas9;in_silico"),
        ("tss_distance_bp", r["midpoint_distance_to_tss_bp"]),
        ("in_minus50_plus300", int(as_bool(r["midpoint_in_minus50_plus300"]))),
        ("ATAC_overlap", int(as_bool(r["overlaps_atac_peak"]))),
        ("TALAM1_overlap", int(as_bool(r["overlaps_talam1_annotation"]))),
        ("nearby_promoter_overlap", int(as_bool(r["overlaps_nearby_lncRNA_promoter_window"]))),
        ("evidence_class", "current_in_silico_candidate"),
        ("known_knockdown_result", ""),
    ]))
for r in stojic:
    historical_rows.append(OrderedDict([
        ("source_group", "Stojic2018_HeLa"),
        ("guide_id", r["guide_label"]),
        ("guide_sequence", r["published_sequence"]),
        ("cell_or_library_context", f'{r["cell_line"]};{r["crispri_effector"]}'),
        ("tss_distance_bp", r["selected_midpoint_distance_to_anchor_bp"]),
        ("in_minus50_plus300", int(as_bool(r["midpoint_in_crispri_minus50_plus300"]))),
        ("ATAC_overlap", int(as_bool(r["overlaps_atac_peak"]))),
        ("TALAM1_overlap", int(as_bool(r["overlaps_talam1_annotation"]))),
        ("nearby_promoter_overlap", int(as_bool(r["overlaps_nearby_lncRNA_promoter_window"]))),
        ("evidence_class", r["evidence_class"]),
        ("known_knockdown_result", r["knockdown_result"]),
    ]))
for idx,r in enumerate(crincl,1):
    historical_rows.append(OrderedDict([
        ("source_group", "CRiNCL2017_design"),
        ("guide_id", f"CRiNCL_{idx:02d}"),
        ("guide_sequence", r["library_spacer_20nt"]),
        ("cell_or_library_context", r["source_sublibraries"]),
        ("tss_distance_bp", r["selected_midpoint_distance_to_anchor_bp"]),
        ("in_minus50_plus300", int(as_bool(r["midpoint_in_crispri_minus50_plus300"]))),
        ("ATAC_overlap", int(as_bool(r["overlaps_atac_peak"]))),
        ("TALAM1_overlap", int(as_bool(r["overlaps_talam1_annotation"]))),
        ("nearby_promoter_overlap", int(as_bool(r["overlaps_nearby_lncRNA_promoter_window"]))),
        ("evidence_class", r["evidence_class"]),
        ("known_knockdown_result", ""),
    ]))
write_tsv(TABLES / "historical_locus_context_comparison.tsv", historical_rows)

hist_summary = []
for group in ["Current_A549_panel","Stojic2018_HeLa","CRiNCL2017_design"]:
    group_rows = [r for r in historical_rows if r["source_group"]==group]
    distances = np.array([float(r["tss_distance_bp"]) for r in group_rows])
    hist_summary.append(OrderedDict([
        ("source_group", group),
        ("guide_count", len(group_rows)),
        ("median_tss_distance_bp", round(float(np.median(distances)), 3)),
        ("minimum_absolute_tss_distance_bp", round(float(np.min(np.abs(distances))), 3)),
        ("in_minus50_plus300_count", sum(int(r["in_minus50_plus300"]) for r in group_rows)),
        ("ATAC_overlap_count", sum(int(r["ATAC_overlap"]) for r in group_rows)),
        ("TALAM1_overlap_count", sum(int(r["TALAM1_overlap"]) for r in group_rows)),
        ("nearby_promoter_overlap_count", sum(int(r["nearby_promoter_overlap"]) for r in group_rows)),
    ]))
write_tsv(TABLES / "historical_locus_context_summary.tsv", hist_summary)


# ---------------------------------------------------------------------
# 11. Integrated candidate technical status
# ---------------------------------------------------------------------
orig_robust = {r["candidate_id"]: r for r in smaa_overall}
global_by_id = {r["candidate_id"]: r for r in global_rows}
pos_stab = {r["candidate_id"]: r for r in position_candidate_summary}
seq_stab = {r["candidate_id"]: r for r in sequence_candidate_summary}
off_stab = {r["candidate_id"]: r for r in off_candidate}

integrated_rows = []
for i,cid in enumerate(ids):
    o = orig_robust[cid]
    g = global_by_id[cid]
    p = pos_stab[cid]
    s = seq_stab[cid]
    m = off_stab[cid]
    # Transparent evidence flags, not another weighted score.
    robust_three_profiles = float(o["minimum_p_top_5_across_profiles"]) >= 0.50
    robust_global = float(g["p_top_5"]) >= 0.50
    position_stable = int(p["maximum_rank_across_position_models_profiles"]) <= 15
    sequence_stable = int(s["maximum_rank_across_sequence_scenarios_profiles"]) <= 15
    model_stable = int(m["maximum_rank_across_models_profiles"]) <= 15
    pareto_front = int(pareto_layer[i]) == 1
    no_one_mm = one_mm[i] == 0
    flag_count = sum([
        robust_three_profiles,
        robust_global,
        position_stable,
        sequence_stable,
        model_stable,
        pareto_front,
        no_one_mm,
    ])
    if flag_count >= 6 and robust_three_profiles:
        tier = "Tier_A_highly_robust"
    elif flag_count >= 5:
        tier = "Tier_B_robust_tradeoff"
    elif flag_count >= 3:
        tier = "Tier_C_context_dependent"
    else:
        tier = "Tier_D_low_priority"
    integrated_rows.append(OrderedDict([
        ("candidate_id", cid),
        ("technical_tier", tier),
        ("evidence_flag_count_0_to_7", flag_count),
        ("three_profile_min_p_top5", o["minimum_p_top_5_across_profiles"]),
        ("three_profile_mean_p_top5", o["mean_p_top_5_across_profiles"]),
        ("global_weight_p_top5", g["p_top_5"]),
        ("position_max_rank", p["maximum_rank_across_position_models_profiles"]),
        ("sequence_max_rank", s["maximum_rank_across_sequence_scenarios_profiles"]),
        ("offtarget_model_max_rank", m["maximum_rank_across_models_profiles"]),
        ("pareto_layer", int(pareto_layer[i])),
        ("one_mismatch_hits", int(one_mm[i])),
        ("panel_flag", int(cid in panel_ids)),
        ("robust_three_profiles_flag", int(robust_three_profiles)),
        ("robust_global_flag", int(robust_global)),
        ("position_stable_top15_flag", position_stable),
        ("sequence_stable_top15_flag", sequence_stable),
        ("offtarget_model_stable_top15_flag", model_stable),
        ("pareto_front_flag", int(pareto_front)),
        ("no_one_mismatch_flag", int(no_one_mm)),
    ]))
integrated_rows.sort(key=lambda r: (-r["evidence_flag_count_0_to_7"], -float(r["three_profile_mean_p_top5"]), -float(r["global_weight_p_top5"])))
for k,r in enumerate(integrated_rows,1):
    r["integrated_order"] = k
    r.move_to_end("integrated_order", last=False)
write_tsv(TABLES / "integrated_candidate_technical_status.tsv", integrated_rows)


# ---------------------------------------------------------------------
# 12. Figures
# ---------------------------------------------------------------------
# Position rank ranges for top integrated 12.
top12 = integrated_rows[:12]
top12_ids = [r["candidate_id"] for r in top12]
fig, ax = plt.subplots(figsize=(11,7))
for y,cid in enumerate(top12_ids):
    rr = [position_rank_map[k][id_to_i[cid]] for k in position_rank_map]
    ax.plot([min(rr), max(rr)], [y,y], marker="|")
    ax.scatter(np.mean(rr), y, s=50)
ax.set_yticks(range(len(top12_ids)))
ax.set_yticklabels([c.replace("MALAT1_A549_NGG_","") for c in top12_ids])
ax.invert_yaxis()
ax.set_xlabel("Rank across position models and profiles")
ax.set_ylabel("Candidate ID suffix")
ax.set_title("Position-model sensitivity of leading candidates")
ax.grid(axis="x", alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_7_position_model_rank_ranges.png")
plt.close(fig)

# Sequence rank ranges.
fig, ax = plt.subplots(figsize=(11,7))
for y,cid in enumerate(top12_ids):
    rr = [sequence_rank_map[k][id_to_i[cid]] for k in sequence_rank_map]
    ax.plot([min(rr), max(rr)], [y,y], marker="|")
    ax.scatter(np.mean(rr), y, s=50)
ax.set_yticks(range(len(top12_ids)))
ax.set_yticklabels([c.replace("MALAT1_A549_NGG_","") for c in top12_ids])
ax.invert_yaxis()
ax.set_xlabel("Rank across sequence-QC scenarios and profiles")
ax.set_ylabel("Candidate ID suffix")
ax.set_title("Sequence-QC sensitivity of leading candidates")
ax.grid(axis="x", alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_8_sequence_qc_rank_ranges.png")
plt.close(fig)

# Global SMAA top 10.
gtop = global_rows[:10][::-1]
fig, ax = plt.subplots(figsize=(10,7))
ax.barh(
    range(len(gtop)),
    [float(r["p_top_5"]) for r in gtop],
)
ax.set_yticks(range(len(gtop)))
ax.set_yticklabels([r["candidate_id"].replace("MALAT1_A549_NGG_","") for r in gtop])
ax.set_xlabel("Top-5 probability under globally sampled weights")
ax.set_ylabel("Candidate ID suffix")
ax.set_title("Global-weight SMAA robustness")
ax.set_xlim(0,1)
ax.grid(axis="x",alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_9_global_weight_smaa_top10.png")
plt.close(fig)

# Pareto trade-off.
fig, ax = plt.subplots(figsize=(10,7))
size = 70 + 500*jost
ax.scatter(eff_locked, local_locked, s=size, alpha=.55)
front_idx = np.where(pareto_layer==1)[0]
ax.scatter(eff_locked[front_idx], local_locked[front_idx], s=size[front_idx], facecolors="none", linewidths=1.5)
for cid in panel_ids:
    i=id_to_i[cid]
    ax.annotate(cid.replace("MALAT1_A549_NGG_",""),(eff_locked[i],local_locked[i]),xytext=(5,5),textcoords="offset points",fontsize=8)
ax.set_xlabel("Efficacy-context domain")
ax.set_ylabel("Local-safety domain")
ax.set_title("Pareto structure of efficacy and local safety")
ax.grid(alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_10_pareto_efficacy_local_safety.png")
plt.close(fig)

# MDS.
fig, ax = plt.subplots(figsize=(10,7))
ax.scatter(coords[:,0],coords[:,1],s=60+400*np.array([float(orig_robust[c]["mean_p_top_5_across_profiles"]) for c in ids]),alpha=.55)
for cid in panel_ids:
    i=id_to_i[cid]
    ax.annotate(cid.replace("MALAT1_A549_NGG_",""),(coords[i,0],coords[i,1]),xytext=(5,5),textcoords="offset points",fontsize=8)
ax.set_xlabel("MDS dimension 1")
ax.set_ylabel("MDS dimension 2")
ax.set_title("Two-dimensional criterion-space representation")
ax.grid(alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_11_mds_criterion_space.png")
plt.close(fig)

# Uncertainty source impact.
fig, ax = plt.subplots(figsize=(10,6))
labels=[r["uncertainty_source"] for r in oat_summary_rows]
vals=[float(r["mean_absolute_rank_shift_vs_central"]) for r in oat_summary_rows]
ax.bar(range(len(labels)),vals)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels,rotation=35,ha="right")
ax.set_ylabel("Mean absolute rank shift")
ax.set_title("One-at-a-time uncertainty-source impact")
ax.grid(axis="y",alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_12_uncertainty_source_rank_impact.png")
plt.close(fig)

# Convergence.
fig, ax = plt.subplots(figsize=(9,6))
ax.plot([r["iterations"] for r in conv_rows],[r["max_abs_p_top5_difference_vs_30000_all_candidates"] for r in conv_rows],marker="o")
ax.set_xlabel("Monte Carlo iterations")
ax.set_ylabel("Maximum |P(top 5) - final P(top 5)|")
ax.set_title("Monte Carlo convergence")
ax.grid(alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_13_monte_carlo_convergence.png")
plt.close(fig)

# Historical TSS distance.
group_order=["Current_A549_panel","Stojic2018_HeLa","CRiNCL2017_design"]
fig, ax = plt.subplots(figsize=(10,6))
for x,group in enumerate(group_order):
    vals=[float(r["tss_distance_bp"]) for r in historical_rows if r["source_group"]==group]
    jitter=np.linspace(-.15,.15,len(vals)) if len(vals)>1 else [0]
    ax.scatter(np.full(len(vals),x)+jitter,vals,alpha=.65)
ax.axhspan(-50,300,alpha=.12)
ax.axhline(0,linewidth=1)
ax.set_xticks(range(len(group_order)))
ax.set_xticklabels(group_order,rotation=20,ha="right")
ax.set_ylabel("Midpoint distance to A549 active TSS (bp)")
ax.set_title("Current and historical MALAT1 guide locations")
ax.grid(axis="y",alpha=.25)
fig.tight_layout()
save_png(fig, FIGURES/"Figure_14_historical_guide_tss_distance.png")
plt.close(fig)


# ---------------------------------------------------------------------
# 13. Documentation
# ---------------------------------------------------------------------
technical_report = f"""# MALAT1 A549 CRISPRi — Final Technical Completion Report v1

## Scope

This package closes the computational and decision-analytic work for the
single-locus MALAT1/A549/ZIM3-KRAB-dCas9 case study. It extends the locked
86-candidate matrix and the original three-profile SMAA analysis with all
pre-specified and reasonable optional stress tests that can be performed using
the frozen project data.

## Completed analyses

1. Alternative position-response functions:
   - baseline trapezoid;
   - narrow trapezoid;
   - broad trapezoid;
   - triangular peak at +100 bp;
   - Gaussian center at +100 bp.

2. Sequence-QC sensitivity:
   - locked strict composite;
   - moderate 25–75% GC / no poly-T / homopolymer ≤5;
   - poly-T-only screen;
   - continuous feasibility penalty;
   - removal of the sequence domain;
   - hard strict eligibility.

3. Global-weight SMAA:
   - {GLOBAL_ITER:,} simulations with Dirichlet(1,1,1,1) domain weights;
   - uniform within-domain splits;
   - the same position, ATAC, local-context and Jost uncertainty families;
   - strict no-1-mismatch companion analysis.

4. Pareto non-dominance over four domains.

5. Leave-one-domain-out stress testing.

6. Off-target model substitution:
   - Jost CRISPRi;
   - Hsu2013;
   - Doench CFD;
   - three-model percentile ensemble.

7. Classical MDS on the six normalized criteria.

8. One-at-a-time uncertainty-source audit.

9. Monte Carlo convergence at 1k, 5k, 10k, 20k and 30k iterations.

10. Five-seed stability audit.

11. Locus-context comparison with mapped Stojic 2018 and CRiNCL historical
    guides.

12. Integrated technical-status table using transparent binary evidence flags,
    not a new compensatory score.

## MDS fit

- Kruskal stress-1: {stress1:.4f}
- Distance R²: {distance_r2:.4f}
- Positive-eigenvalue variance represented in 2D: {explained:.4f}

The MDS plot is descriptive. It is not used to rank candidates.

## Interpretation rule

No output in this package is experimental validation. Rank probabilities are
conditional on the frozen locus, A549 evidence, candidate set, criteria,
normalization functions and uncertainty distributions.

## Data-acquisition scope freeze

No additional cCRE criterion, second cell line or additional guide library is
added. Those would create a new biological study scope rather than complete the
current calculations. ENCODE peak-set duplication and cCRE evidence reuse are
therefore not treated as independent criteria.

## Remaining non-computational work

- external peer review and any venue-specific manuscript revision;
- experimental repression testing;
- measurement of unintended local transcriptional effects;
- empirical genome-wide off-target validation.

These cannot be completed by further calculation on the present data.
"""
(DOCS/"TECHNICAL_COMPLETION_REPORT.md").write_text(technical_report,encoding="utf-8")

# Data dictionary for new outputs.
dictionary_rows = [
    OrderedDict([("file","position_model_sensitivity.tsv"),("purpose","Candidate ranks under five alternative position transforms and three preference profiles."),("key_interpretation","Model-form stress test; none is claimed to be a validated efficacy predictor.")]),
    OrderedDict([("file","sequence_qc_sensitivity.tsv"),("purpose","Candidate ranks under six sequence-QC treatments."),("key_interpretation","Tests dependence on QC thresholds and compensatory versus hard-screen handling.")]),
    OrderedDict([("file","global_weight_smaa_summary.tsv"),("purpose","Rank acceptability with globally sampled domain weights."),("key_interpretation","Preference-agnostic robustness stress test.")]),
    OrderedDict([("file","pareto_front_and_layers.tsv"),("purpose","Non-dominated sorting over four decision domains."),("key_interpretation","Shows trade-off-efficient candidates without weights.")]),
    OrderedDict([("file","leave_one_domain_out.tsv"),("purpose","Ranks after deleting one decision domain at a time."),("key_interpretation","Detects dependence on a single domain.")]),
    OrderedDict([("file","offtarget_model_sensitivity.tsv"),("purpose","Rankings with Jost, Hsu, CFD and an ensemble."),("key_interpretation","Cross-model robustness; scores remain predictions.")]),
    OrderedDict([("file","mds_coordinates.tsv"),("purpose","Two-dimensional classical MDS coordinates."),("key_interpretation","Descriptive visualization only.")]),
    OrderedDict([("file","uncertainty_source_oat_summary.tsv"),("purpose","One-at-a-time uncertainty impact summary."),("key_interpretation","Not a formal Sobol decomposition.")]),
    OrderedDict([("file","monte_carlo_convergence.tsv"),("purpose","Convergence against the 30,000-iteration estimate."),("key_interpretation","Numerical stability check.")]),
    OrderedDict([("file","random_seed_stability.tsv"),("purpose","Top-5 probability stability across five seeds."),("key_interpretation","Numerical reproducibility check.")]),
    OrderedDict([("file","historical_locus_context_comparison.tsv"),("purpose","Current panel versus mapped historical guides."),("key_interpretation","Locus-context comparison, not head-to-head efficacy testing.")]),
    OrderedDict([("file","integrated_candidate_technical_status.tsv"),("purpose","Transparent evidence-flag synthesis."),("key_interpretation","Tiering is based on explicit flags, not a new hidden weighted score.")]),
]
write_tsv(TABLES/"DATA_DICTIONARY.tsv",dictionary_rows)

summary = {
    "status": "technical_calculations_complete",
    "candidate_count": n,
    "global_smaa_iterations": GLOBAL_ITER,
    "oat_iterations_per_source": OAT_ITER,
    "convergence_final_iterations": 30000,
    "seed_stability_runs": SEED_N,
    "position_models": list(position_models.keys()),
    "sequence_scenarios": list(sequence_scenarios.keys()),
    "offtarget_models": list(off_models.keys()),
    "pareto_front_candidate_count": int(np.sum(pareto_layer==1)),
    "mds": mds_fit,
    "integrated_top_10": [r["candidate_id"] for r in integrated_rows[:10]],
    "locked_panel": panel_ids,
    "interpretation_guard": "All candidate performance is in-silico and conditional on the frozen study design.",
}
(DOCS/"technical_completion_manifest.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")

print("Technical analyses completed.")
print("Top integrated candidates:")
for r in integrated_rows[:10]:
    print(r["integrated_order"], r["candidate_id"], r["technical_tier"], r["evidence_flag_count_0_to_7"])
print("Created root:", ROOT)
