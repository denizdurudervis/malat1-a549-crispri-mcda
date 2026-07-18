# Reproducibility validation report

**Overall: PASS**

Tolerance: `5e-08`

## Invariants
- PASS — raw_candidate_count_86
- PASS — normalized_candidate_count_86
- PASS — unique_candidate_ids
- PASS — all_jost_present
- PASS — one_mismatch_candidate_count_6
- PASS — no_extra_exact_offtarget
- PASS — core_integrated_top3
- PASS — profile_robustness_top1
- PASS — all_22_expected_png_figures_valid
- PASS — all_14_expected_pdf_figures_valid

## Output comparisons
- PASS — `generated/normalized/malat1_a549_normalized_decision_matrix_v0_3.tsv`; max difference `0.0`
- PASS — `generated/profile_smaa/malat1_a549_smaa_central_weights_v0_1.tsv`; max difference `0.0`
- PASS — `generated/profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv`; max difference `0.0`
- PASS — `generated/profile_smaa/malat1_a549_smaa_profile_summary_v0_1.tsv`; max difference `0.0`
- PASS — `generated/profile_smaa/malat1_a549_smaa_rank_acceptability_v0_1.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/global_weight_smaa_summary.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/integrated_candidate_technical_status.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/monte_carlo_convergence.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/offtarget_model_summary.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/pareto_front_and_layers.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/position_model_summary.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/random_seed_summary.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/sequence_qc_summary.tsv`; max difference `0.0`

Validated locked scientific TSV outputs: `13`.

All expected figure files were checked for valid nonempty PNG or PDF artifacts. Figure byte equality is not used as a scientific equality test because rendering and font metadata can differ across systems.