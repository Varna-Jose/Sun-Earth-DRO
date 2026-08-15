from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "DRO_injection"))

from gmat.runner import GmatNotFoundError, find_gmat_executable, run_gmat


def test_explicit_gmat_executable(tmp_path):
    exe = tmp_path / "GmatConsole"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    assert find_gmat_executable(exe) == exe.resolve()


def test_missing_gmat(monkeypatch):
    monkeypatch.delenv("GMAT_BIN", raising=False)
    monkeypatch.setenv("PATH", "")
    with pytest.raises(GmatNotFoundError):
        find_gmat_executable()


def test_runner_passes_script_to_executable(tmp_path):
    exe = tmp_path / "GmatConsole"
    exe.write_text("#!/bin/sh\necho SCRIPT=$1\nexit 0\n")
    exe.chmod(0o755)
    script = tmp_path / "case.script"
    script.write_text("% test\n")

    result = run_gmat(script, executable=exe)

    assert result.success
    assert str(script.resolve()) in result.stdout
