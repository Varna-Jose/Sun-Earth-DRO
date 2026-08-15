"""Shared repository paths for the DRO campaign.

Source code stays under ``DRO_injection``. Runtime artifacts are kept out of
source directories:

- Python products: ``DRO_injection/output/python``
- GMAT reports: ``DRO_injection/output/gmat``
- Generated GMAT scripts: ``DRO_injection/generated/gmat``
- Shareable result package: ``dist/A_within_one_year``
"""
from pathlib import Path

DRO_DIR = Path(__file__).resolve().parent
ROOT = DRO_DIR.parent
SCRIPTS_DIR = ROOT / "scripts"
PYTHON_OUTPUT = DRO_DIR / "output" / "python"
GMAT_OUTPUT = DRO_DIR / "output" / "gmat"
GENERATED_GMAT = DRO_DIR / "generated" / "gmat"
DIST_ROOT = ROOT / "dist"


def ensure_runtime_dirs() -> None:
    """Create runtime directories used by the campaign."""
    for path in (PYTHON_OUTPUT, GMAT_OUTPUT, GENERATED_GMAT, DIST_ROOT):
        path.mkdir(parents=True, exist_ok=True)
