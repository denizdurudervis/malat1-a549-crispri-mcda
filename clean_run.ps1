$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host "Checking for Python 3.13..."

& py -3.13 -c "import sys; print(sys.version)"
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "Python 3.13 is not installed." -ForegroundColor Red
    Write-Host "Install a 64-bit Python 3.13 release, close PowerShell, reopen it, and run this script again."
    Write-Host "Installed versions can be checked with: py -0p"
    exit 1
}

$venv = ".venv-repro"

if (Test-Path $venv) {
    Write-Host "Removing previous virtual environment..."
    Remove-Item -Recurse -Force $venv
}

Write-Host "Creating clean Python 3.13 virtual environment..."
& py -3.13 -m venv $venv
if ($LASTEXITCODE -ne 0 -or -not (Test-Path "$venv\Scripts\python.exe")) {
    throw "Virtual environment creation failed."
}

$python = Join-Path $venv "Scripts\python.exe"

Write-Host "Installing locked dependencies..."
& $python -m pip install --upgrade pip
if ($LASTEXITCODE -ne 0) { throw "pip upgrade failed." }

& $python -m pip install -r requirements-lock.txt
if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }

Write-Host "Running the reproducibility pipeline..."
& $python scripts\run_pipeline.py
if ($LASTEXITCODE -ne 0) { throw "Pipeline failed." }

Write-Host ""
Write-Host "POST-FLASHFRY REPRODUCIBILITY PIPELINE COMPLETE" -ForegroundColor Green
Write-Host "Open generated\reproducibility_validation_report.md"
