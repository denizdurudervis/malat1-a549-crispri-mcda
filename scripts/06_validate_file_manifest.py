#!/usr/bin/env python3
"""Build or validate the byte-level repository release manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "FILE_MANIFEST_SHA256.tsv"
EXCLUDED_DIRS = {".git", ".venv-repro", "__pycache__"}
EXCLUDED_SUFFIXES = {".pyc", ".log"}
EXCLUDED_NAMES = {".DS_Store", MANIFEST.name}


def included_files():
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if any(part in EXCLUDED_DIRS for part in relative.parts):
            continue
        if path.name in EXCLUDED_NAMES or path.suffix in EXCLUDED_SUFFIXES:
            continue
        yield path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def current_rows():
    return [
        {
            "relative_path": path.relative_to(ROOT).as_posix(),
            "size_bytes": str(path.stat().st_size),
            "sha256": sha256(path),
        }
        for path in included_files()
    ]


def write_manifest(rows):
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "size_bytes", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} entries -> {MANIFEST.relative_to(ROOT)}")


def validate_manifest(rows):
    if not MANIFEST.exists():
        raise SystemExit("Manifest is missing. Run with --write before release.")
    with MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
        deposited = list(csv.DictReader(handle, delimiter="\t"))
    current = {row["relative_path"]: row for row in rows}
    expected = {row["relative_path"]: row for row in deposited}
    issues = []
    for relative in sorted(set(current) | set(expected)):
        if relative not in expected:
            issues.append(f"unlisted:{relative}")
        elif relative not in current:
            issues.append(f"missing:{relative}")
        elif current[relative] != expected[relative]:
            issues.append(f"content_mismatch:{relative}")
    if issues:
        print("Repository manifest validation FAIL")
        for issue in issues[:30]:
            print(f"- {issue}")
        if len(issues) > 30:
            print(f"- ... and {len(issues) - 30} more")
        raise SystemExit(1)
    print(f"Repository manifest validation PASS ({len(rows)} files)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write",
        action="store_true",
        help="refresh FILE_MANIFEST_SHA256.tsv after an intentional release edit",
    )
    args = parser.parse_args()
    rows = current_rows()
    if args.write:
        write_manifest(rows)
    else:
        validate_manifest(rows)


if __name__ == "__main__":
    main()
