# Reproducibility close-out protocol

## What this validates

The pipeline verifies that the frozen post-FlashFry input regenerates:

- the 86-row normalized matrix;
- the three-profile rank-acceptability results;
- global-weight SMAA summaries;
- sensitivity and integrated-status tables;
- the locked core ordering `0046`, `0047`, `0048`.

## What this does not validate

It does not repeat public-data downloads, GENCODE parsing, whole-genome index construction, or FlashFry discovery. Those stages are preserved through frozen inputs, commands, versions, hashes and manifests. It also does not provide wet-lab validation.

## Pass criteria

- all comparison keys match;
- selected numeric results differ by no more than `5e-8`;
- 86 unique candidates remain;
- six candidates have a one-mismatch hit;
- no candidate has an extra exact genomic match;
- the primary robust candidate is `NGG_0046` and the integrated top three are `0046–0048`.
