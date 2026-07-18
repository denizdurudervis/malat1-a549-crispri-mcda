# Uncertainty-aware CRISPRi guide prioritization for the MALAT1 locus in A549

This repository contains the reproducible, post-FlashFry decision-analytic
pipeline for a single-locus CRISPR interference case study.

The study prioritizes SpCas9-NGG guides for **MALAT1 repression in A549 cells
under the ZIM3-KRAB–dCas9 modality**. It combines cell-context-specific
transcription evidence, local locus risks, sequence feasibility and
whole-genome off-target predictions in a customized **SMAA-2-type stochastic
multicriteria ranking framework**.

> **Scope guard:** The output is an in-silico prioritization panel for
> experimental validation. It does not establish repression efficiency,
> biological safety, clinical suitability or experimentally measured
> off-target activity.

## Study snapshot

- Reference genome: GRCh38.p14
- Annotation: GENCODE 50 / Ensembl 116
- A549 active MALAT1 TSS anchor: chr11:65,499,045
- Candidate set: 86 SpCas9-NGG guides within TSS ±500 bp
- Primary off-target criterion: Jost–Santos CRISPRi specificity
- Alternative off-target sensitivity models: Hsu2013 and Doench CFD
- Main uncertainty analysis:
  - 30,000 simulations per preference profile
  - 50,000 globally sampled-weight simulations
  - position-model, sequence-QC and off-target-model sensitivity
  - Pareto, leave-one-domain-out, convergence and seed-stability analyses

## Main technical interpretation

The strongest core candidates within the prespecified model family are:

1. `MALAT1_A549_NGG_0046`
2. `MALAT1_A549_NGG_0047`
3. `MALAT1_A549_NGG_0048`

`MALAT1_A549_NGG_0045` is a conditional fourth core candidate.
`MALAT1_A549_NGG_0076` and `MALAT1_A549_NGG_0071` are deliberate trade-off
comparators rather than equal-confidence winners.

These statements are conditional on the frozen A549 context, candidate set,
criteria, normalization functions and uncertainty distributions.

## Quick start

### Windows PowerShell

Python 3.13 (64-bit) is required.

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\clean_run.ps1
```

A successful run ends with:

```text
Validation PASS
POST-FLASHFRY REPRODUCIBILITY PIPELINE COMPLETE
```

### Linux or macOS

```bash
bash clean_run.sh
```

### Existing Python 3.13 environment

```bash
python -m pip install -r requirements-lock.txt
python scripts/06_validate_file_manifest.py
python scripts/run_pipeline.py
```

The final validation report is written to:

```text
generated/reproducibility_validation_report.md
```

## Colab

Open `notebooks/MALAT1_A549_post_FlashFry_reproducibility.ipynb` in Google
Colab. The public repository URL is preconfigured; the notebook also supports
running from an extracted repository ZIP.

The notebook deliberately starts from the frozen post-FlashFry input matrix.
It does not rebuild the large GRCh38 whole-genome FlashFry index.

## Pipeline

1. Rebuild the normalized decision matrix.
2. Re-run three-profile stochastic rank-acceptability analysis.
3. Compute central-weight confidence-factor analogues.
4. Re-run global-weight and sensitivity analyses.
5. Rebuild the six manuscript and eight supplementary publication figures.
6. Compare generated scientific TSV outputs against locked references.
7. Confirm that all 22 PNG and 14 PDF figure artifacts are valid.

## Repository structure

```text
.
├── .github/workflows/       # continuous reproducibility check
├── docs/                    # methods, scope and reproducibility notes
├── generated/               # reproduced outputs from the latest run
├── inputs/                  # frozen post-FlashFry inputs and provenance
├── notebooks/               # Colab-ready reproducibility notebook
├── reference/               # locked scientific reference outputs
├── scripts/                 # executable analysis pipeline
├── clean_run.ps1            # Windows clean run
├── clean_run.sh             # Linux/macOS clean run
├── requirements-lock.txt
├── environment.yml
├── software_versions.json
├── CITATION.cff
└── LICENSE
```

## Reproducibility status

The post-FlashFry pipeline has passed clean-run validation in:

- an isolated Linux/Python 3.13 environment;
- a local Windows/Python 3.13 environment.

All columns in 13 locked scientific TSV outputs are compared by keys and numeric
tolerance. The pipeline also checks 22 expected PNG and 14 expected PDF figure
artifacts. Figure binary hashes are intentionally not treated as scientific
equality tests because rendering, compression and font metadata can differ
without changing the numerical results. `FILE_MANIFEST_SHA256.tsv` separately
records the byte-level integrity of the deposited release package.

## Important methodological notes

- CAGE was used to define the active A549 TSS and was not counted again as an
  independent weighted criterion.
- ATAC peak sets processed from the same experiment are not treated as
  independent biological replicates.
- TALAM1 and the nearby antisense-promoter context are modeled as graded risks,
  not experimentally proven adverse effects.
- Jost, Hsu and CFD are predictive scoring models, not empirical off-target
  measurements.
- Historical MALAT1 guides are interpreted as having limited transferability to
  the current A549 TSS context, not as inherently incorrect designs.
- The 2D MDS representation is supplementary and descriptive because its stress
  value indicates limited low-dimensional fidelity.

## Citation

Citation metadata and the repository URL are provided in `CITATION.cff`. The
archived release DOI and associated-manuscript record should be added after
those persistent identifiers exist.

## License

Code is released under the MIT License. Data files retain the terms of their
original public sources; consult the corresponding source records before reuse.
