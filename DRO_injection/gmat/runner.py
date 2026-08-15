# DRO_injection/gmat/runner.py

from __future__ import annotations

from dataclasses import dataclass

import os
import shutil
import subprocess
from pathlib import Path


class GmatNotFoundError(RuntimeError):
    pass


def find_gmat_executable(explicit: str | Path | None = None) -> Path:
    """Resolve the GMAT console executable.

    Resolution order:
    1. Explicit argument
    2. GMAT_BIN environment variable
    3. Executable found on PATH
    """

    candidates: list[str | Path] = []

    if explicit is not None:
        candidates.append(explicit)

    env_path = os.environ.get("GMAT_BIN")
    if env_path:
        candidates.append(env_path)

    for candidate in candidates:
        path = Path(candidate).expanduser().resolve()
        if path.is_file():
            return path

    for executable_name in (
        "GmatConsole",
        "GMATConsole",
    ):
        discovered = shutil.which(executable_name)
        if discovered:
            return Path(discovered).resolve()

    raise GmatNotFoundError(
        "GMAT console executable not found. "
        "Set GMAT_BIN or add GmatConsole to PATH."
    )


@dataclass(frozen=True)
class GmatResult:
    script: Path
    returncode: int
    stdout: str
    stderr: str

    @property
    def success(self) -> bool:
        return self.returncode == 0


def run_gmat(
    script: str | Path,
    *,
    executable: str | Path | None = None,
    cwd: str | Path | None = None,
    timeout: float | None = None,
) -> GmatResult:

    script_path = Path(script).resolve()

    if not script_path.is_file():
        raise FileNotFoundError(f"GMAT script does not exist: {script_path}")

    gmat = find_gmat_executable(executable)

    working_directory = (
        Path(cwd).resolve()
        if cwd is not None
        else script_path.parent
    )

    completed = subprocess.run(
        [str(gmat), str(script_path)],
        cwd=working_directory,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    return GmatResult(
        script=script_path,
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
    )