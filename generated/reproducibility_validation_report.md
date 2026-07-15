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

## Output comparisons
- PASS — `generated/normalized/malat1_a549_normalized_decision_matrix_v0_3.tsv`; max difference `0.0`
- PASS — `generated/profile_smaa/malat1_a549_smaa_profile_summary_v0_1.tsv`; max difference `0.0`
- PASS — `generated/profile_smaa/malat1_a549_smaa_overall_robustness_v0_1.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/global_weight_smaa_summary.tsv`; max difference `0.0`
- PASS — `generated/extended/tables/integrated_candidate_technical_status.tsv`; max difference `0.0`

PNG and XLSX binary hashes are not used as scientific equality tests because rendering and archive metadata can differ across systems.