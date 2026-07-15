# Claude Technical Handoff — Read This First

You are receiving the complete technical package for a computational CRISPRi
guide-prioritization case study.

## Locked study identity

- Target: MALAT1
- Cell context: A549
- Effector: ZIM3-KRAB-dCas9 CRISPRi
- Genome: GRCh38.p14
- Annotation: GENCODE 50 / Ensembl 116
- Active A549 TSS anchor: chr11:65,499,045
- Candidate set: 86 SpCas9-NGG guides in TSS ±500 bp
- Primary genome-wide off-target score: Jost–Santos CRISPRi specificity
- Decision framework: hierarchical multi-criteria analysis with SMAA-like
  Monte Carlo rank acceptability

## Files to start with

1. `tables/integrated_candidate_technical_status.tsv`
2. `tables/global_weight_smaa_summary.tsv`
3. `tables/pareto_front_and_layers.tsv`
4. `tables/position_candidate_rank_stability.tsv`
5. `tables/sequence_candidate_rank_stability.tsv`
6. `tables/offtarget_candidate_rank_stability.tsv`
7. `docs/TECHNICAL_COMPLETION_REPORT.md`
8. the original `malat1_a549_raw_decision_matrix_v0_2.tsv`
9. the original `malat1_a549_normalized_decision_matrix_v0_3.tsv`

## Writing constraints

- Describe results as in-silico prioritization, not validated efficacy.
- Do not call a guide clinically safe or biologically validated.
- Do not describe CFD/Hsu/Jost predictions as experimental off-target data.
- Explain that CAGE defines the TSS and is not separately weighted.
- Explain that ATAC peak files from the same ENCODE experiment are processed
  alternatives, not independent biological replicates.
- Treat local TALAM1 and nearby antisense-promoter effects as graded risks.
- State that historical guides were designed around an older/inactive TSS
  relative to the A549 evidence.
- Report uncertainty distributions as modelling assumptions.
- Do not infer unpublished knockdown percentages from bar heights.

## Preferred manuscript logic

Problem → cell-context-specific TSS audit → candidate generation → local and
whole-genome safety audit → hierarchical criteria → profile SMAA → global-weight
stress test → sensitivity/Pareto/MDS → provisional validation panel →
limitations.
