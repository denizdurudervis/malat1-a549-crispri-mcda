# Clean-run completion status

**Result: PASS**

The standalone post-FlashFry pipeline was executed in a newly created isolated Python virtual environment using the pinned package versions. The generated scientific TSV outputs matched the locked reference outputs within tolerance `5e-08`.

Validated items include:

- 86 unique candidate guides;
- six candidates with a 1-mismatch hit;
- no extra exact genomic target;
- profile-based stochastic rankings;
- global-weight stochastic rankings;
- extended sensitivity analyses;
- integrated core ordering `0046`, `0047`, `0048`;
- profile-specific central-weight confidence-factor analogues.

This clean run starts from the frozen post-FlashFry raw decision matrix. Whole-genome index construction remains documented through the original manifests and hashes rather than repeated here.
