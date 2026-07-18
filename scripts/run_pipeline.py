#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys, shutil, argparse

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

def run(name):
    cmd=[sys.executable,str(SCRIPTS/name)]
    print("\n>>> " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=ROOT)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--keep-generated',action='store_true',help='do not delete prior generated outputs')
    args=p.parse_args()
    generated=ROOT/'generated'
    if generated.exists() and not args.keep_generated:
        shutil.rmtree(generated)
    generated.mkdir(parents=True,exist_ok=True)
    for script in [
        '01_build_normalized_matrix.py',
        '02_run_profile_smaa.py',
        '04_compute_confidence_factors.py',
        '03_run_extended_suite.py',
        '07_build_manuscript_figures.py',
        '05_validate_outputs.py',
    ]:
        run(script)
    print("\nPOST-FLASHFRY REPRODUCIBILITY PIPELINE COMPLETE")

if __name__=='__main__': main()
