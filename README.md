# Sun-Earth DRO

Python campaign-design and GMAT-validation workflow for a nine-spacecraft Sun-Earth distant retrograde orbit (DRO) phasing study.

## Repository layout

```text
Sun-Earth-DRO/
├── DRO_injection/
│   ├── campaign_config.py          # launch epoch / campaign configuration
│   ├── paths.py                    # one source of truth for repo/runtime paths
│   ├── scan_launch_dates.py
│   ├── phase9_single_launch.py     # Python optimisation
│   ├── departure.py                # grouped departure design
│   ├── phase9_figures.py
│   ├── phase9_snapshot.py
│   ├── run_campaign.py             # end-to-end orchestrator
│   ├── build_results_package.py    # builds dist/A_within_one_year/
│   ├── gmat/
│   │   ├── runner.py               # portable GmatConsole discovery/execution
│   │   ├── gen_phase9_gmat.py      # 9 slot + 9 transfer validation scripts
│   │   └── gen_all9_replay.py      # all-spacecraft replay generator
│   ├── generated/gmat/             # generated .script files (gitignored)
│   └── output/
│       ├── python/                  # Python CSV/JSON/figures (gitignored)
│       └── gmat/                    # GMAT report files (gitignored)
├── scripts/
│   └── lambert.py
├── tests/
├── dist/                            # exported result packages (gitignored)
├── pyproject.toml
└── .env.example
```

Source code is kept separate from generated scripts and runtime results. The campaign never writes into the GMAT installation directory.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

GMAT is an external application. Either place `GmatConsole` on `PATH` or set:

```bash
export GMAT_BIN="/path/to/GmatConsole"
```

Do not commit your machine-specific GMAT path.

## Run

From the repository root:

```bash
python DRO_injection/run_campaign.py
```

The pipeline performs the Python optimisation and departure design, creates figures, generates and executes the per-spacecraft GMAT validation scripts, generates and executes the fleet replay, and finally builds a shareable package under `dist/A_within_one_year/`.

## Runtime artifacts

The three runtime areas are intentionally separate:

- `DRO_injection/output/python/`: Python-generated campaign products.
- `DRO_injection/generated/gmat/`: generated GMAT input scripts.
- `DRO_injection/output/gmat/`: GMAT-generated report files.

All three are reproducible and ignored by Git except for `.gitkeep` placeholders.
