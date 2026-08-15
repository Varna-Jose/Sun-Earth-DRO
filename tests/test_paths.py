from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "DRO_injection"))

from paths import DRO_DIR, GENERATED_GMAT, GMAT_OUTPUT, PYTHON_OUTPUT, ROOT as PROJECT_ROOT


def test_paths_stay_inside_repo():
    assert PROJECT_ROOT == ROOT
    for path in (DRO_DIR, GENERATED_GMAT, GMAT_OUTPUT, PYTHON_OUTPUT):
        assert path.is_relative_to(ROOT)
