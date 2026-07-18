# Repository package validation

The GitHub-ready package was validated before delivery.

## Checks performed

- `python scripts/run_pipeline.py`
  - Result: PASS
  - Every column in all 13 locked scientific TSV outputs matched the reference
    tables within the configured tolerance.
  - All 22 expected generated PNG and 14 expected generated PDF files were
    valid nonempty artifacts.

- `python scripts/06_validate_file_manifest.py`
  - Result: PASS
  - Deposited files matched the release-level size and SHA-256 manifest.

- `notebooks/MALAT1_A549_post_FlashFry_reproducibility.ipynb`
  - Result: PASS
  - Every notebook cell executed successfully in repository-root context.
  - The notebook reproduced the pipeline and confirmed the validation report.

## Interpretation

This validates the deposited post-FlashFry computational workflow. It does not
provide experimental validation of guide repression or safety.
