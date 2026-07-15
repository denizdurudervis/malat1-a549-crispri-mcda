# Windows local clean run

1. Extract this folder somewhere outside the main project, for example the Desktop.
2. Open the folder in File Explorer.
3. Click the address bar, type `powershell`, and press Enter.
4. Run:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\clean_run.ps1
```

The script creates a separate `.venv-repro`, installs the pinned packages and runs every post-FlashFry analysis. Success is shown by:

```text
Validation PASS
POST-FLASHFRY REPRODUCIBILITY PIPELINE COMPLETE
```

Then open:

```text
generated\reproducibility_validation_report.md
```

This local run is an independent Windows confirmation. It does not modify the main project folder.
