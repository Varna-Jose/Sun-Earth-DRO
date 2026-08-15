from pathlib import Path
import subprocess
import sys
import time

from gmat.runner import run_gmat
from paths import GENERATED_GMAT, ensure_runtime_dirs


HERE = Path(__file__).resolve().parent
PYTHON = sys.executable

PYTHON_STEPS = [
    ("Python optimisation", "phase9_single_launch.py"),
    ("Departure design", "departure.py"),
    ("Figure pack", "phase9_figures.py"),
    ("Snapshot figure", "phase9_snapshot.py"),
    ("GMAT validation generation", "gmat/gen_phase9_gmat.py"),
    ("Replay generator", "gmat/gen_all9_replay.py"),
    ("RESULTS package", "build_results_package.py"),
]


def run_python_step(name: str, relative_path: str) -> None:
    path = HERE / relative_path

    print(f"\n=== {name} ===", flush=True)

    result = subprocess.run(
        [PYTHON, str(path)],
        cwd=path.parent,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        print("\n".join(result.stdout.strip().splitlines()[-12:]))

    if result.returncode != 0:
        if result.stderr:
            print(result.stderr[-2000:])
        raise RuntimeError(f"STEP FAILED: {name}")


def main() -> None:
    ensure_runtime_dirs()
    started = time.time()

    for name, path in PYTHON_STEPS[:-1]:
        run_python_step(name, path)

    print("\n=== Replay GMAT run ===")
    result = run_gmat(
        GENERATED_GMAT / "g03_all9_replay.script",
        timeout=1800,
    )

    if not result.success:
        raise RuntimeError(
            f"GMAT replay failed:\n{result.stderr}"
        )

    run_python_step(*PYTHON_STEPS[-1])

    minutes = (time.time() - started) / 60
    print(f"\nCampaign complete in {minutes:.1f} min.")


if __name__ == "__main__":
    main()