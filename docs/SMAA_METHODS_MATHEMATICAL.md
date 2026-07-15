# Mathematical specification — customized SMAA-2-type framework

## Alternatives and criteria

Let alternative guide `i = 1,...,86` have stochastic criterion vector

`x_i^(s) = (x_i1^(s), ..., x_i6^(s))`

for simulation draw `s`. The six criteria are position support, A549 ATAC support, TALAM1 safety, nearby-promoter safety, sequence feasibility, and Jost CRISPRi specificity.

## Additive utility

For criterion weights `w_k^(s) >= 0` and `sum_k w_k^(s)=1`:

`U_i^(s) = sum_k w_k^(s) x_ik^(s)`.

The hierarchy is represented by four domain weights—efficacy, local safety, sequence feasibility, and whole-genome safety—and stochastic within-domain splits for the two paired domains.

## Weight distributions

Profile runs use Dirichlet distributions centred on:

- balanced: `(0.30, 0.30, 0.15, 0.25)`;
- efficacy-focused: `(0.45, 0.20, 0.10, 0.25)`;
- safety-focused: `(0.20, 0.35, 0.15, 0.30)`.

The Dirichlet concentration is 60. Efficacy and local-safety within-domain splits use `Beta(20,20)`. The global stress test uses `Dirichlet(1,1,1,1)` and uniform within-domain splits.

## Criterion uncertainty

- TSS shift: clipped Normal with SD 2 bp and limits ±6 bp.
- Full-support position boundaries: triangular `(-75,-50,-25)` upstream and `(250,300,350)` downstream.
- ATAC classification reliability: `Beta(19,1)`.
- TALAM1 overlap penalty: triangular `(0.10,0.30,0.60)`.
- Nearby-promoter overlap penalty: triangular `(0.40,0.70,0.95)`.
- Jost specificity: beta distribution centred on each observed value, concentration 40.
- Sequence feasibility: fixed under the primary run; varied separately in threshold sensitivity analyses.

These distributions are transparent modelling assumptions, not experimentally estimated error distributions.

## Rank acceptability

For rank `r`, the estimated rank-acceptability index is

`b_i^r = (1/S) sum_s I(rank_i^(s)=r)`.

The top-five probability is

`P_i(Top5) = sum_(r=1)^5 b_i^r`.

## Central weights and confidence-factor analogue

For alternatives that ranked first at least once, the central criterion-weight vector is the mean criterion-weight vector conditional on the alternative ranking first. The close-out pipeline fixes that vector, resamples criterion uncertainty, and estimates the probability that the same candidate ranks first. This is reported as a profile-specific SMAA-2-type confidence-factor analogue.

## Naming discipline

The method should be described as a **customized SMAA-2-type stochastic rank-acceptability framework** unless the manuscript adopts and documents the complete standard SMAA-2 output definitions. It is not a new SMAA theory and does not experimentally validate guide efficacy or safety.
