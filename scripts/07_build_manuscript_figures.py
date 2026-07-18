#!/usr/bin/env python3
"""Generate publication-oriented main and supplementary figures from frozen TSV outputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.colors import (
    BoundaryNorm,
    LinearSegmentedColormap,
    ListedColormap,
    Normalize,
)
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"
OUT = GENERATED / "manuscript_figures"
MAIN = OUT / "main_figures"
SUPP = OUT / "supplementary_figures"
for path in (MAIN, SUPP):
    path.mkdir(parents=True, exist_ok=True)


COLORS = {
    "navy": "#164B73",
    "blue": "#2F7EAA",
    "teal": "#2A9D8F",
    "orange": "#E28E2C",
    "red": "#C44E52",
    "purple": "#7A5195",
    "gray": "#8A9197",
    "light_gray": "#E8ECEF",
    "dark": "#20262E",
    "green": "#4E9A51",
}

PROFILE_COLORS = {
    "balanced": COLORS["navy"],
    "efficacy_focused": COLORS["orange"],
    "safety_focused": COLORS["teal"],
}

PANEL_IDS = [
    "MALAT1_A549_NGG_0046",
    "MALAT1_A549_NGG_0047",
    "MALAT1_A549_NGG_0048",
    "MALAT1_A549_NGG_0045",
    "MALAT1_A549_NGG_0076",
    "MALAT1_A549_NGG_0071",
]

PANEL_LABELS = {
    "MALAT1_A549_NGG_0046": "0046 · core",
    "MALAT1_A549_NGG_0047": "0047 · core",
    "MALAT1_A549_NGG_0048": "0048 · core",
    "MALAT1_A549_NGG_0045": "0045 · conditional",
    "MALAT1_A549_NGG_0076": "0076 · specificity contrast",
    "MALAT1_A549_NGG_0071": "0071 · local-context contrast",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def suffix(candidate_id: str) -> str:
    return candidate_id.rsplit("_", 1)[-1]


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def style_ax(ax, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(colors=COLORS["dark"], labelsize=9)
    ax.xaxis.label.set_color(COLORS["dark"])
    ax.yaxis.label.set_color(COLORS["dark"])
    if grid_axis:
        ax.grid(axis=grid_axis, color="#C7CDD3", linewidth=0.6, alpha=0.55, zorder=0)


def save(fig, path: Path) -> None:
    png = path.with_suffix(".png")
    pdf = path.with_suffix(".pdf")
    png.unlink(missing_ok=True)
    pdf.unlink(missing_ok=True)
    fig.savefig(
        png,
        dpi=320,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
        metadata={"Software": "MALAT1 A549 reproducibility pipeline"},
    )
    fig.savefig(
        pdf,
        bbox_inches="tight",
        facecolor=fig.get_facecolor(),
        metadata={
            "Creator": "MALAT1 A549 reproducibility pipeline",
            "Producer": "Matplotlib",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    if png.stat().st_size <= 100 or png.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"Invalid PNG render: {png}")
    if pdf.stat().st_size <= 100 or pdf.read_bytes()[:4] != b"%PDF":
        raise RuntimeError(f"Invalid PDF render: {pdf}")
    plt.close(fig)


def panel_letter(ax, letter: str) -> None:
    ax.text(
        -0.08,
        1.04,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        color=COLORS["dark"],
    )


def build_figure_1(raw: list[dict[str, str]]) -> None:
    anchor = 65_499_045
    window_left, window_right = anchor - 500, anchor + 500
    atac_left, atac_right, summit = 65_498_986, 65_499_152, 65_499_071
    nearby_tss = anchor - 141
    candidate_midpoints = np.array([float(r["protospacer_midpoint_1based"]) for r in raw])
    candidate_colors = [COLORS["navy"] if r["guide_orientation"] == "+" else COLORS["orange"] for r in raw]

    fig = plt.figure(figsize=(7.15, 5.7))
    grid = fig.add_gridspec(2, 1, height_ratios=(1.55, 1.15), hspace=0.42)
    ax = fig.add_subplot(grid[0])
    ax.axvspan(window_left, window_right, color=COLORS["light_gray"], alpha=0.9, zorder=0)
    ax.axvspan(atac_left, atac_right, color="#8FD0C7", alpha=0.8, zorder=1)
    ax.axvline(anchor, color=COLORS["red"], linewidth=2.0, zorder=4)
    ax.axvline(summit, color=COLORS["teal"], linewidth=1.3, linestyle="--", zorder=4)
    ax.axvline(nearby_tss, color=COLORS["purple"], linewidth=1.3, linestyle=":", zorder=4)
    for x, color in zip(candidate_midpoints, candidate_colors):
        ax.vlines(x, 0.04, 0.29, color=color, linewidth=0.8, alpha=0.75, zorder=3)
    ax.hlines(0.165, window_left, window_right, color="#AEB5BC", linewidth=0.7, zorder=2)
    ax.text(anchor, 0.89, "A549 active-TSS anchor\nchr11:65,499,045", ha="center", va="top", color=COLORS["red"], fontsize=9, fontweight="bold")
    ax.text(atac_left + 10, 0.58, "A549 ATAC peak", ha="left", va="bottom", fontsize=8.5, color="#125E56")
    ax.text(summit + 3, 0.31, "summit", rotation=90, va="bottom", ha="left", fontsize=8, color="#125E56")
    ax.text(nearby_tss - 3, 0.38, "nearby antisense\nTSS (−141 bp)", rotation=90, va="bottom", ha="right", fontsize=8, color=COLORS["purple"])
    ax.text(window_left + 8, 0.305, "86 NGG candidates", fontsize=8.5, color=COLORS["dark"], va="bottom")
    ax.set_xlim(window_left - 35, window_right + 35)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("GRCh38 coordinate on chromosome 11")
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda value, _: f"{int(value):,}"))
    style_ax(ax)
    panel_letter(ax, "A")

    ax2 = fig.add_subplot(grid[1])
    evidence = [
        ("GENCODE 50", "version-pinned transcript annotation", COLORS["gray"]),
        ("FANTOM5 CAGE", "A549 transcription-initiation evidence", COLORS["red"]),
        ("ENCODE ATAC-seq", "A549 accessibility evidence", COLORS["teal"]),
        ("ENCODE poly(A)− RNA-seq", "A549 expression evidence", COLORS["blue"]),
    ]
    for y, (source, role, color) in enumerate(evidence[::-1]):
        ax2.scatter(0.06, y, s=95, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        ax2.hlines(y, 0.09, 0.96, color="#CDD2D7", linewidth=1.0, zorder=1)
        ax2.text(0.12, y + 0.14, source, fontsize=9, fontweight="bold", color=COLORS["dark"], va="center")
        ax2.text(0.12, y - 0.16, role, fontsize=8.3, color="#4C5661", va="center")
        ax2.scatter(0.96, y, marker="D", s=54, color=COLORS["navy"], edgecolor="white", linewidth=0.7, zorder=3)
    ax2.text(0.96, 3.48, "convergent support for\nthe reporting anchor", ha="center", va="bottom", fontsize=8.6, color=COLORS["navy"], fontweight="bold")
    ax2.set_xlim(0, 1.05)
    ax2.set_ylim(-0.55, 3.75)
    ax2.axis("off")
    panel_letter(ax2, "B")
    fig.suptitle("Cell-context evidence and candidate window at the MALAT1 locus", fontsize=13, fontweight="bold", color=COLORS["dark"], y=0.98)
    save(fig, MAIN / "Figure_1_A549_TSS_evidence_and_candidate_window")


def build_figure_2(norm_by_id: dict[str, dict[str, str]]) -> None:
    columns = [
        ("C1_position_soft_support", "Position\nsupport"),
        ("C2_A549_ATAC_support", "A549 ATAC\nsupport"),
        ("C3_TALAM1_safety", "TALAM1\nnon-overlap"),
        ("C4_nearby_lncRNA_promoter_safety", "Nearby TSS\nnon-overlap"),
        ("C5_sequence_feasibility", "Sequence\nfeasibility"),
        ("C6_Jost_CRISPRi_specificity", "Jost CRISPRi\nspecificity"),
    ]
    values = np.array([[float(norm_by_id[cid][key]) for key, _ in columns] for cid in PANEL_IDS])
    palette = LinearSegmentedColormap.from_list(
        "khaki_sage",
        ["#354032", "#596650", "#78886A", "#98A98A", "#BAC7AA", "#DCE3CF"],
        N=256,
    )
    value_norm = Normalize(vmin=0, vmax=1)
    fig, ax = plt.subplots(figsize=(7.15, 4.6), facecolor="#FBFCF8")
    ax.set_facecolor("#FBFCF8")
    image = ax.imshow(
        values,
        cmap=palette,
        norm=value_norm,
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(-0.5, values.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, values.shape[0], 1), minor=True)
    ax.grid(which="minor", color="#F4F6EF", linewidth=1.15)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            red, green, blue, _ = palette(value_norm(value))
            def linearise(channel):
                return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4
            luminance = (
                0.2126 * linearise(red)
                + 0.7152 * linearise(green)
                + 0.0722 * linearise(blue)
            )
            color = "#172019" if luminance > 0.39 else "#FFFFFF"
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8.5, color=color, fontweight="bold")
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([label for _, label in columns], fontsize=8.5, color="#263025")
    ax.set_yticks(range(len(PANEL_IDS)))
    ax.set_yticklabels([PANEL_LABELS[cid] for cid in PANEL_IDS], fontsize=8.7, color="#263025")
    ax.set_xlabel("Normalized criterion value (higher is more favourable)", labelpad=12, color="#263025")
    ax.set_title("Role-structured six-guide panel across the six decision criteria", fontsize=12.5, fontweight="bold", pad=12, color="#263025")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Normalized value", fontsize=9, color="#263025")
    cbar.ax.tick_params(labelsize=8, colors="#263025")
    cbar.outline.set_edgecolor("#596650")
    cbar.outline.set_linewidth(0.7)
    fig.tight_layout()
    save(fig, MAIN / "Figure_2_role_structured_panel_criteria")


def build_figure_3(profile_rows, overall_by_id, global_by_id) -> None:
    profile_by = {(r["profile"], r["candidate_id"]): r for r in profile_rows}
    fig, ax = plt.subplots(figsize=(7.15, 4.8))
    y = np.arange(len(PANEL_IDS))
    offsets = {"balanced": -0.16, "efficacy_focused": 0.0, "safety_focused": 0.16}
    for idx, cid in enumerate(PANEL_IDS):
        vals = [float(profile_by[(profile, cid)]["p_top_5"]) for profile in offsets]
        ax.hlines(idx, min(vals), max(vals), color="#B8C0C7", linewidth=2.0, zorder=1)
        for profile, offset in offsets.items():
            ax.scatter(
                float(profile_by[(profile, cid)]["p_top_5"]),
                idx + offset,
                s=42,
                color=PROFILE_COLORS[profile],
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
        ax.scatter(
            float(global_by_id[cid]["p_top_5"]),
            idx,
            marker="D",
            s=42,
            facecolor="white",
            edgecolor=COLORS["purple"],
            linewidth=1.4,
            zorder=4,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([PANEL_LABELS[cid] for cid in PANEL_IDS])
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Probability of ranking in the top five")
    ax.set_title("Profile robustness and global-weight stress test of the panel", fontsize=12.5, fontweight="bold", pad=12)
    style_ax(ax, "x")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PROFILE_COLORS["balanced"], markeredgecolor="white", markersize=7, label="Balanced"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PROFILE_COLORS["efficacy_focused"], markeredgecolor="white", markersize=7, label="Efficacy-focused"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=PROFILE_COLORS["safety_focused"], markeredgecolor="white", markersize=7, label="Safety-focused"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white", markeredgecolor=COLORS["purple"], markersize=6.5, label="Global weights"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8.5, ncol=2)
    fig.tight_layout()
    save(fig, MAIN / "Figure_3_panel_profile_and_global_robustness")


def build_figure_4(norm, overall_by_id, pareto_by_id) -> None:
    ids = [r["candidate_id"] for r in norm]
    x = np.array([float(r["efficacy_domain_score"]) for r in norm])
    y = np.array([float(r["whole_genome_safety_domain_score"]) for r in norm])
    local = np.array([float(r["local_safety_domain_score"]) for r in norm])
    robust = np.array([float(overall_by_id[cid]["mean_p_top_5_across_profiles"]) for cid in ids])
    layers = np.array([int(pareto_by_id[cid]["pareto_layer_4_domains"]) for cid in ids])
    sizes = 22 + 240 * np.sqrt(robust)
    cmap = ListedColormap(["#C44E52", "#E2B13C", "#2A9D8F"])
    norm_color = BoundaryNorm([-0.1, 0.25, 0.75, 1.1], cmap.N)
    fig, ax = plt.subplots(figsize=(7.15, 5.3))
    sc = ax.scatter(x, y, c=local, cmap=cmap, norm=norm_color, s=sizes, alpha=0.72, edgecolor="white", linewidth=0.65, zorder=2)
    front = layers == 1
    ax.scatter(x[front], y[front], s=sizes[front] + 38, facecolors="none", edgecolors=COLORS["dark"], linewidths=1.15, zorder=3)
    offsets = {
        "MALAT1_A549_NGG_0046": (-50, 22),
        "MALAT1_A549_NGG_0047": (-58, 35),
        "MALAT1_A549_NGG_0048": (-58, -34),
        "MALAT1_A549_NGG_0045": (-58, 9),
        "MALAT1_A549_NGG_0076": (8, 8),
        "MALAT1_A549_NGG_0071": (8, -18),
    }
    for cid in PANEL_IDS:
        i = ids.index(cid)
        ax.annotate(
            suffix(cid),
            (x[i], y[i]),
            xytext=offsets[cid],
            textcoords="offset points",
            fontsize=8.3,
            fontweight="bold",
            color=COLORS["dark"],
            arrowprops=dict(arrowstyle="-", color="#808890", linewidth=0.7),
            path_effects=[pe.withStroke(linewidth=2.5, foreground="white")],
            zorder=5,
        )
    ax.set_xlabel("Efficacy-context domain score")
    ax.set_ylabel("Predicted Jost CRISPRi specificity")
    ax.set_xlim(-0.02, 1.04)
    ax.set_ylim(-0.03, 1.04)
    ax.set_title("Trade-off structure across efficacy, local context, and predicted specificity", fontsize=12.1, fontweight="bold", pad=12)
    style_ax(ax, "both")
    color_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=cmap(0), markersize=7, label="Local-context score 0.0"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=cmap(1), markersize=7, label="Local-context score 0.5"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=cmap(2), markersize=7, label="Local-context score 1.0"),
        Line2D([0], [0], marker="o", color=COLORS["dark"], markerfacecolor="none", markersize=8, label="Four-domain Pareto front"),
    ]
    ax.legend(handles=color_handles, loc="lower left", frameon=True, framealpha=0.94, fontsize=7.8, ncol=2)
    ax.text(0.99, 0.03, "Point area ∝ mean P(top five) across profiles", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.7, color="#58616A")
    fig.tight_layout()
    save(fig, MAIN / "Figure_4_tradeoff_and_four_domain_pareto_projection")


def build_figure_5(oat_rows) -> None:
    mapping = {
        "weights_only": "Preference\nweights",
        "position_only": "Position\nmodel",
        "ATAC_only": "ATAC\nreliability",
        "local_context_only": "Local-context\npenalties",
        "Jost_only": "Jost\nspecificity",
        "all_sources_combined": "All sources\ncombined",
    }
    order = ["weights_only", "Jost_only", "local_context_only", "ATAC_only", "position_only", "all_sources_combined"]
    by_source = {r["uncertainty_source"]: r for r in oat_rows}
    values = [float(by_source[key]["mean_absolute_rank_shift_vs_central"]) for key in order]
    colors = [COLORS["navy"], COLORS["purple"], COLORS["orange"], COLORS["teal"], COLORS["green"], COLORS["red"]]
    fig, ax = plt.subplots(figsize=(7.15, 4.45))
    x = np.arange(len(order))
    bars = ax.bar(x, values, color=colors, width=0.72, zorder=2)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.045, f"{value:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=COLORS["dark"])
    ax.axvline(4.5, color="#A7ADB3", linewidth=1.0, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels([mapping[key] for key in order])
    ax.set_ylabel("Mean absolute rank shift from modal central configuration")
    ax.set_ylim(0, max(values) + 0.35)
    ax.set_title("One-at-a-time audit identifies preference weights as the leading single source", fontsize=12.0, fontweight="bold", pad=12)
    style_ax(ax, "y")
    ax.text(0.985, 0.96, "Combined scenario is not an OAT component", transform=ax.transAxes, ha="right", va="top", fontsize=7.8, color="#5D6670")
    fig.tight_layout()
    save(fig, MAIN / "Figure_5_uncertainty_source_rank_impact")


def build_figure_6(hist_rows) -> None:
    group_order = ["Current_A549_panel", "Stojic2018_HeLa", "CRiNCL2017_design"]
    labels = ["Current A549\npanel", "Stojic 2018\nHeLa", "CRiNCL 2017\ndesign"]
    group_colors = [COLORS["navy"], COLORS["orange"], COLORS["gray"]]
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 4.7), gridspec_kw={"width_ratios": [1.05, 1.25]}, sharex=True)
    rng = np.random.default_rng(20260718)
    for ax_idx, ax in enumerate(axes):
        for x_idx, group in enumerate(group_order):
            vals = np.array([float(r["tss_distance_bp"]) for r in hist_rows if r["source_group"] == group])
            jitter = rng.uniform(-0.14, 0.14, size=len(vals))
            ax.scatter(np.full(len(vals), x_idx) + jitter, vals, s=28 if group != "Current_A549_panel" else 42, color=group_colors[x_idx], edgecolor="white", linewidth=0.5, alpha=0.78, zorder=3)
        ax.axhspan(-50, 300, color="#B9D8F0", alpha=0.42, zorder=0)
        ax.axhline(0, color=COLORS["navy"], linewidth=1.0, zorder=1)
        ax.set_xticks(range(3))
        ax.set_xticklabels(labels, fontsize=8.2)
        style_ax(ax, "y")
    axes[0].set_ylim(-1800, 12500)
    axes[0].set_ylabel("Guide midpoint distance to the A549 active-TSS anchor (bp)")
    axes[0].set_title("Full range", fontsize=10.5, fontweight="bold")
    axes[1].set_ylim(-1650, 500)
    axes[1].set_title("TSS-proximal view", fontsize=10.5, fontweight="bold")
    panel_letter(axes[0], "A")
    panel_letter(axes[1], "B")
    axes[1].text(0.98, 0.965, "Shaded region: −50 to +300 bp", transform=axes[1].transAxes, ha="right", va="top", fontsize=7.7, color=COLORS["navy"])
    fig.suptitle("Historical MALAT1 guides require re-evaluation against the A549 active TSS", fontsize=12.3, fontweight="bold", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, MAIN / "Figure_6_historical_guide_distance_full_and_zoom")


def build_supp_rank_heatmap(rank_rows, overall_rows) -> None:
    top = [r["candidate_id"] for r in sorted(overall_rows, key=lambda r: int(r["robustness_order"]))[:8]]
    lookup = {(r["profile"], r["candidate_id"], int(r["rank"])): float(r["rank_acceptability"]) for r in rank_rows}
    values = np.array([[lookup[("balanced", cid, rank)] for rank in range(1, 16)] for cid in top])
    fig, ax = plt.subplots(figsize=(7.15, 4.1))
    im = ax.imshow(values, cmap="magma", vmin=0, vmax=max(0.75, values.max()), aspect="auto")
    ax.set_xticks(range(15))
    ax.set_xticklabels(range(1, 16))
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([suffix(cid) for cid in top])
    ax.set_xlabel("Rank (display truncated at rank 15)")
    ax.set_ylabel("Candidate suffix")
    ax.set_title("Balanced-profile rank acceptability for the top eight robustness-ordered candidates", fontsize=11.4, fontweight="bold", pad=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Rank acceptability")
    ax.tick_params(length=0)
    fig.tight_layout()
    save(fig, SUPP / "Figure_S1_balanced_rank_acceptability")


def build_supp_global(global_rows) -> None:
    top = sorted(global_rows, key=lambda r: int(r["global_robustness_order"]))[:10][::-1]
    fig, ax = plt.subplots(figsize=(7.15, 4.4))
    vals = [float(r["p_top_5"]) for r in top]
    bars = ax.barh(range(10), vals, color=[COLORS["navy"] if as_bool(r["panel_flag"]) else COLORS["gray"] for r in top], zorder=2)
    for bar, value in zip(bars, vals):
        ax.text(value + 0.012, bar.get_y() + bar.get_height() / 2, f"{value:.3f}", va="center", fontsize=8)
    ax.set_yticks(range(10))
    ax.set_yticklabels([suffix(r["candidate_id"]) for r in top])
    ax.set_xlim(0, 1)
    ax.set_xlabel("P(top five) under globally sampled weights")
    ax.set_ylabel("Candidate suffix")
    ax.set_title("Global-weight SMAA robustness", fontsize=11.5, fontweight="bold")
    style_ax(ax, "x")
    ax.legend(handles=[Line2D([0], [0], color=COLORS["navy"], lw=6, label="Panel candidate"), Line2D([0], [0], color=COLORS["gray"], lw=6, label="Other candidate")], frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    save(fig, SUPP / "Figure_S2_global_weight_top10")


def range_panel(rows, min_col, max_col, mean_col, title, xlabel, outfile) -> None:
    by_id = {r["candidate_id"]: r for r in rows}
    fig, ax = plt.subplots(figsize=(7.15, 3.8))
    for y, cid in enumerate(PANEL_IDS):
        row = by_id[cid]
        lo, hi, mean = float(row[min_col]), float(row[max_col]), float(row[mean_col])
        ax.hlines(y, lo, hi, color=COLORS["gray"], linewidth=2.4, zorder=1)
        ax.scatter(mean, y, s=48, color=COLORS["navy"] if y < 3 else COLORS["orange"], edgecolor="white", linewidth=0.6, zorder=3)
        ax.scatter([lo, hi], [y, y], marker="|", s=95, color=COLORS["gray"], zorder=2)
    ax.set_yticks(range(len(PANEL_IDS)))
    ax.set_yticklabels([PANEL_LABELS[cid] for cid in PANEL_IDS])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=11.4, fontweight="bold")
    style_ax(ax, "x")
    fig.tight_layout()
    save(fig, SUPP / outfile)


def build_supp_seed_convergence(seed_rows, convergence_rows) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.8))
    seeds = [int(r["seed_index"]) for r in seed_rows]
    rho = [float(r["spearman_p_top5_vs_seed1"]) for r in seed_rows]
    jac = [float(r["top10_jaccard_vs_seed1"]) for r in seed_rows]
    axes[0].plot(seeds, rho, marker="o", color=COLORS["navy"], label="Spearman ρ")
    axes[0].plot(seeds, jac, marker="s", color=COLORS["orange"], label="Top-10 Jaccard")
    axes[0].set_ylim(0.75, 1.02)
    axes[0].set_xticks(seeds)
    axes[0].set_xlabel("Seed run")
    axes[0].set_ylabel("Agreement with seed 1")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title("Seed stability", fontsize=10.5, fontweight="bold")
    style_ax(axes[0], "y")
    panel_letter(axes[0], "A")
    iters = [int(r["iterations"]) for r in convergence_rows]
    diff = [float(r["max_abs_p_top5_difference_vs_30000_all_candidates"]) for r in convergence_rows]
    axes[1].plot(iters, diff, marker="o", color=COLORS["teal"])
    axes[1].set_xlabel("Iterations")
    axes[1].set_ylabel("Maximum |P(top five) − 30k estimate|")
    axes[1].set_title("Monte Carlo convergence", fontsize=10.5, fontweight="bold")
    style_ax(axes[1], "both")
    panel_letter(axes[1], "B")
    fig.tight_layout()
    save(fig, SUPP / "Figure_S6_seed_stability_and_convergence")


def build_supp_mds(mds_rows) -> None:
    ids = [r["candidate_id"] for r in mds_rows]
    x = np.array([float(r["MDS1"]) for r in mds_rows])
    y = np.array([float(r["MDS2"]) for r in mds_rows])
    robust = np.array([float(r["mean_p_top5_across_profiles"]) for r in mds_rows])
    fig, ax = plt.subplots(figsize=(7.15, 4.5))
    ax.scatter(x, y, s=28 + 230 * np.sqrt(robust), color=COLORS["blue"], alpha=0.6, edgecolor="white", linewidth=0.5)
    offsets = [(6, 7), (6, -13), (6, 7), (-28, 8), (6, 7), (6, 7)]
    for cid, offset in zip(PANEL_IDS, offsets):
        i = ids.index(cid)
        ax.annotate(suffix(cid), (x[i], y[i]), xytext=offset, textcoords="offset points", fontsize=8, fontweight="bold", path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    ax.set_xlabel("MDS dimension 1")
    ax.set_ylabel("MDS dimension 2")
    ax.set_title("Descriptive two-dimensional representation of six-criterion space", fontsize=11.2, fontweight="bold")
    ax.text(0.99, 0.02, "Stress-1 = 0.302; distance R² = 0.782; 2D positive-eigenvalue fraction = 0.662", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="#555E67")
    style_ax(ax, "both")
    fig.tight_layout()
    save(fig, SUPP / "Figure_S7_MDS_criterion_space")


def build_supp_loo(loo_rows) -> None:
    display = {
        "without_efficacy": "Without efficacy",
        "without_local_safety": "Without local context",
        "without_sequence": "Without sequence",
        "without_offtarget": "Without specificity",
    }
    rows = [r for r in loo_rows if r["scenario"] != "baseline_all_domains"]
    fig, ax = plt.subplots(figsize=(7.15, 3.8))
    x = np.arange(len(rows))
    jac = [float(r["top10_jaccard_vs_baseline"]) for r in rows]
    rho = [float(r["spearman_vs_baseline"]) for r in rows]
    width = 0.36
    ax.bar(x - width / 2, jac, width=width, color=COLORS["orange"], label="Top-10 Jaccard")
    ax.bar(x + width / 2, rho, width=width, color=COLORS["navy"], label="Spearman ρ")
    ax.set_xticks(x)
    ax.set_xticklabels([display[r["scenario"]] for r in rows])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Agreement with balanced all-domain ranking")
    ax.set_title("Leave-one-domain-out sensitivity", fontsize=11.4, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5, ncol=2, loc="lower right")
    style_ax(ax, "y")
    fig.tight_layout()
    save(fig, SUPP / "Figure_S8_leave_one_domain_out")


def main() -> None:
    for directory in (MAIN, SUPP):
        for pattern in ("*.tmp.png", "*.tmp.pdf"):
            for stale in directory.glob(pattern):
                stale.unlink()
    raw = read_tsv(ROOT / "inputs" / "malat1_a549_raw_decision_matrix_v0_2.tsv")
    norm = read_tsv(GENERATED / "normalized" / "malat1_a549_normalized_decision_matrix_v0_3.tsv")
    profile = read_tsv(GENERATED / "profile_smaa" / "malat1_a549_smaa_profile_summary_v0_1.tsv")
    overall = read_tsv(GENERATED / "profile_smaa" / "malat1_a549_smaa_overall_robustness_v0_1.tsv")
    rank_accept = read_tsv(GENERATED / "profile_smaa" / "malat1_a549_smaa_rank_acceptability_v0_1.tsv")
    global_rows = read_tsv(GENERATED / "extended" / "tables" / "global_weight_smaa_summary.tsv")
    pareto = read_tsv(GENERATED / "extended" / "tables" / "pareto_front_and_layers.tsv")
    oat = read_tsv(GENERATED / "extended" / "tables" / "uncertainty_source_oat_summary.tsv")
    hist = read_tsv(GENERATED / "extended" / "tables" / "historical_locus_context_comparison.tsv")
    position = read_tsv(GENERATED / "extended" / "tables" / "position_candidate_rank_stability.tsv")
    sequence = read_tsv(GENERATED / "extended" / "tables" / "sequence_candidate_rank_stability.tsv")
    off_target = read_tsv(GENERATED / "extended" / "tables" / "offtarget_candidate_rank_stability.tsv")
    seed = read_tsv(GENERATED / "extended" / "tables" / "random_seed_summary.tsv")
    convergence = read_tsv(GENERATED / "extended" / "tables" / "monte_carlo_convergence.tsv")
    mds = read_tsv(GENERATED / "extended" / "tables" / "mds_coordinates.tsv")
    loo = read_tsv(GENERATED / "extended" / "tables" / "leave_one_domain_out_summary.tsv")

    norm_by_id = {r["candidate_id"]: r for r in norm}
    overall_by_id = {r["candidate_id"]: r for r in overall}
    global_by_id = {r["candidate_id"]: r for r in global_rows}
    pareto_by_id = {r["candidate_id"]: r for r in pareto}

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9.3,
        "axes.titleweight": "bold",
        "figure.dpi": 140,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    build_figure_1(raw)
    build_figure_2(norm_by_id)
    build_figure_3(profile, overall_by_id, global_by_id)
    build_figure_4(norm, overall_by_id, pareto_by_id)
    build_figure_5(oat)
    build_figure_6(hist)
    build_supp_rank_heatmap(rank_accept, overall)
    build_supp_global(global_rows)
    range_panel(
        position,
        "minimum_rank_across_position_models_profiles",
        "maximum_rank_across_position_models_profiles",
        "mean_rank_across_position_models_profiles",
        "Position-model sensitivity of the role-structured panel",
        "Rank across five position models and three profiles",
        "Figure_S3_position_model_rank_ranges",
    )
    range_panel(
        sequence,
        "minimum_rank_across_sequence_scenarios_profiles",
        "maximum_rank_across_sequence_scenarios_profiles",
        "mean_rank_across_sequence_scenarios_profiles",
        "Sequence-QC sensitivity of the role-structured panel",
        "Rank across six sequence-QC scenarios and three profiles",
        "Figure_S4_sequence_QC_rank_ranges",
    )
    range_panel(
        off_target,
        "minimum_rank_across_models_profiles",
        "maximum_rank_across_models_profiles",
        "mean_rank_across_models_profiles",
        "Off-target-model sensitivity of the role-structured panel",
        "Rank across four specificity models and three profiles",
        "Figure_S5_offtarget_model_rank_ranges",
    )
    build_supp_seed_convergence(seed, convergence)
    build_supp_mds(mds)
    build_supp_loo(loo)

    manifest = {
        "main_figures": sorted(p.name for p in MAIN.glob("*.png") if ".tmp." not in p.name),
        "supplementary_figures": sorted(p.name for p in SUPP.glob("*.png") if ".tmp." not in p.name),
        "source": "Frozen repository TSVs regenerated by the validated post-FlashFry pipeline.",
        "figure_1_note": "Schematic evidence summary; not a raw signal track.",
    }
    (OUT / "figure_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
