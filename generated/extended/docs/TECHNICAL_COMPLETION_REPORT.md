# MALAT1 A549 CRISPRi — Final Technical Completion Report v1

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
   - 50,000 simulations with Dirichlet(1,1,1,1) domain weights;
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

- Kruskal stress-1: 0.3016
- Distance R²: 0.7818
- Positive-eigenvalue variance represented in 2D: 0.6617

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
