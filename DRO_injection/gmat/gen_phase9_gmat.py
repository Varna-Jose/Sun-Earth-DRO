"""
GMAT-fidelity rerun of the 9-spacecraft single-launch phasing solution.

Single launch, but no longer a single shared asymptote: Ariane 6's
restartable Vinci + APU support multiple distinct targeted separation
events on one flight, so the 9 spacecraft release in 3 sequential groups
(G1 = S1-S2, G2 = S3-S5, G3 = S6-S9), each with its own optimised C3/angle
(see departure.py). All 9 still leave Earth on the same calendar launch
day - only the post-Ariane escape state differs by group.

Pipeline:
  A. write + run g01_slot<k>.script x9
  B. write + run g02_sc<k>.script x9
  C. parse everything and print the GMAT-vs-Python delta-v table

Generated GMAT scripts:
    DRO_injection/gmat/

Generated GMAT reports:
    DRO_injection/gmat/output/

GMAT executable:
    configured via GMAT_BIN environment variable
"""

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
DRO_DIR = HERE.parent
ROOT = DRO_DIR.parent

sys.path.insert(0, str(DRO_DIR))
from paths import GMAT_OUTPUT, GENERATED_GMAT, PYTHON_OUTPUT, ensure_runtime_dirs
ensure_runtime_dirs()

sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(DRO_DIR))

from lambert import lambert, MU, AU, DAY
from campaign_config import A1_LAUNCH, LAUNCH


GMAT_EXE = os.environ.get("GMAT_BIN")

GMAT_OUT = GMAT_OUTPUT


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

YEAR_D = 365.25636
YEAR = YEAR_D * DAY
N_E = 2 * np.pi / YEAR

AMP = 22122321.422754
ECC = AMP / AU

NSAT = 9

C3 = 4.0
ANG = np.deg2rad(269.0)

R_SEP = 1.5e6

VY_GUESS = 1.026 * 2 * AMP * N_E


# ---------------------------------------------------------------------------
# analytic helpers
# ---------------------------------------------------------------------------

def kepler_E(M, e):
    E = float(M)

    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))

    return E


def earth_state(t):
    th = N_E * t

    return (
        AU * np.array([np.cos(th), np.sin(th), 0.0]),
        np.sqrt(MU / AU)
        * np.array([-np.sin(th), np.cos(th), 0.0]),
    )


def slot_state(k, t):
    tp = k * YEAR / NSAT
    w = N_E * tp

    M = np.mod(N_E * (t - tp), 2 * np.pi)
    E = kepler_E(M, ECC)

    nu = 2 * np.arctan2(
        np.sqrt(1 + ECC) * np.sin(E / 2),
        np.sqrt(1 - ECC) * np.cos(E / 2),
    )

    r = AU * (1 - ECC * np.cos(E))
    p = AU * (1 - ECC**2)

    th = w + nu

    pos = r * np.array(
        [np.cos(th), np.sin(th), 0.0]
    )

    vr = np.sqrt(MU / p) * ECC * np.sin(nu)
    vt = np.sqrt(MU / p) * (1 + ECC * np.cos(nu))

    vel = (
        vr * np.array([np.cos(th), np.sin(th), 0.0])
        + vt * np.array([-np.sin(th), np.cos(th), 0.0])
    )

    return pos, vel


# ---------------------------------------------------------------------------
# departure design
# ---------------------------------------------------------------------------

DEPARTURE_FILE = PYTHON_OUTPUT / "campaign_departure.json"

with DEPARTURE_FILE.open("r", encoding="utf-8") as f:
    DEP = json.load(f)


sc = [
    dict(
        day=m["day"],
        group=m["group"],
        vinf_vec=np.array(m["vinf_vec"]),
        dv1v=np.array(m["dv1v"]),
        dv1=m["dv1"],
        dv2=m["dv2"],
    )
    for m in DEP["members"]
]


def departure_state(vinf_vec):
    """Departure state in SunEarthRot frame."""

    u = vinf_vec / np.linalg.norm(vinf_vec)

    r_rot = R_SEP * u

    w_cross_r = N_E * np.array(
        [-r_rot[1], r_rot[0], 0.0]
    )

    v_rot = vinf_vec - w_cross_r

    return r_rot, v_rot


# ---------------------------------------------------------------------------
# common GMAT resources
# ---------------------------------------------------------------------------

COMMON = """
Create CoordinateSystem SunEcl;
GMAT SunEcl.Origin = Sun;
GMAT SunEcl.Axes   = MJ2000Ec;

Create CoordinateSystem SunEarthRot;
GMAT SunEarthRot.Origin    = Earth;
GMAT SunEarthRot.Axes      = ObjectReferenced;
GMAT SunEarthRot.XAxis     = R;
GMAT SunEarthRot.ZAxis     = N;
GMAT SunEarthRot.Primary   = Sun;
GMAT SunEarthRot.Secondary = Earth;

Create ForceModel HelioFM;
GMAT HelioFM.CentralBody  = Sun;
GMAT HelioFM.PointMasses  = {Sun, Earth, Luna, Venus, Mars, Jupiter};
GMAT HelioFM.Drag         = None;
GMAT HelioFM.SRP          = On;
GMAT HelioFM.SRP.Flux     = 1367;
GMAT HelioFM.SRP.SRPModel = Spherical;
GMAT HelioFM.ErrorControl = RSSStep;

Create Propagator Helio;
GMAT Helio.FM              = HelioFM;
GMAT Helio.Type            = PrinceDormand78;
GMAT Helio.InitialStepSize = 3600;
GMAT Helio.Accuracy        = 1e-11;
GMAT Helio.MinStep         = 0;
GMAT Helio.MaxStep         = 86400;
GMAT Helio.MaxStepAttempts = 50;

Create DifferentialCorrector DC;
GMAT DC.ShowProgress      = true;
GMAT DC.MaximumIterations = 40;
GMAT DC.DerivativeMethod  = ForwardDifference;
GMAT DC.Algorithm         = NewtonRaphson;
"""


def sat_block(name, a1_epoch, state6):
    x, y, z, vx, vy, vz = state6

    return f"""
Create Spacecraft {name};
GMAT {name}.DateFormat       = A1ModJulian;
GMAT {name}.Epoch            = '{a1_epoch:.9f}';
GMAT {name}.CoordinateSystem = SunEarthRot;
GMAT {name}.DisplayStateType = Cartesian;
GMAT {name}.X  = {x:.6f};
GMAT {name}.Y  = {y:.6f};
GMAT {name}.Z  = {z:.6f};
GMAT {name}.VX = {vx:.9f};
GMAT {name}.VY = {vy:.9f};
GMAT {name}.VZ = {vz:.9f};
GMAT {name}.DryMass = 265;
GMAT {name}.Cr      = 1.3;
GMAT {name}.SRPArea = 4.0;
"""


# ---------------------------------------------------------------------------
# phase A
# ---------------------------------------------------------------------------

def write_g01(k):
    a1_apse = A1_LAUNCH + (k * YEAR_D / NSAT) - YEAR_D
    a1_ins = A1_LAUNCH + sc[k]["day"]

    report_path = (
        GMAT_OUT / f"g01_slot{k}.txt"
    ).resolve()

    report_path.unlink(missing_ok=True)

    s = (
        "% auto-generated by gen_phase9_gmat.py "
        "- slot DRO convergence\n"
    )

    s += COMMON

    s += sat_block(
        "Slot",
        a1_apse,
        (-AMP, 0, 0, 0, VY_GUESS, 0),
    )

    s += f"""
Create Variable vy_solved;

Create ReportFile rpt;
GMAT rpt.Filename         = '{report_path.as_posix()}';
GMAT rpt.WriteHeaders     = false;
GMAT rpt.SolverIterations = None;
GMAT rpt.Precision        = 16;

BeginMissionSequence;

Target DC {{SolveMode = Solve, ExitMode = SaveAndContinue}};
   Vary DC(Slot.SunEarthRot.VY = {VY_GUESS:.6f}, {{Perturbation = 1e-6, Lower = 6, Upper = 12, MaxStep = 0.3}});
   vy_solved = Slot.SunEarthRot.VY;
   Propagate Helio(Slot) {{Slot.ElapsedDays = 60}};
   Propagate Helio(Slot) {{Slot.SunEarthRot.Y = 0}};
   Achieve DC(Slot.SunEarthRot.VX = 0, {{Tolerance = 1e-5}});
EndTarget;

Propagate Helio(Slot) {{Slot.A1ModJulian = {a1_ins:.9f}}};
Report rpt vy_solved Slot.A1ModJulian Slot.SunEarthRot.X Slot.SunEarthRot.Y Slot.SunEarthRot.Z Slot.SunEarthRot.VX Slot.SunEarthRot.VY Slot.SunEarthRot.VZ Slot.Sun.RMAG;
"""

    p = GENERATED_GMAT / f"g01_slot{k}.script"
    p.write_text(s, encoding="utf-8")

    return p


# ---------------------------------------------------------------------------
# phase B
# ---------------------------------------------------------------------------

def write_g02(k, tgt):
    a1_ins = A1_LAUNCH + sc[k]["day"]

    g = sc[k]["dv1v"]

    r_rot, v_rot = departure_state(
        sc[k]["vinf_vec"]
    )

    summary_path = (
        GMAT_OUT / f"g02_sc{k}_summary.txt"
    ).resolve()

    traj_path = (
        GMAT_OUT / f"g02_sc{k}_traj.txt"
    ).resolve()

    summary_path.unlink(missing_ok=True)
    traj_path.unlink(missing_ok=True)

    s = (
        "% auto-generated by gen_phase9_gmat.py "
        f"- spacecraft transfer + insertion "
        f"(release group G{sc[k]['group'] + 1})\n"
    )

    s += COMMON

    s += sat_block(
        "Sat",
        A1_LAUNCH,
        (
            r_rot[0],
            r_rot[1],
            0,
            v_rot[0],
            v_rot[1],
            0,
        ),
    )

    s += f"""
Create ImpulsiveBurn B1;
GMAT B1.CoordinateSystem = SunEarthRot;
GMAT B1.Element1 = {g[0]:.9f};
GMAT B1.Element2 = {g[1]:.9f};
GMAT B1.Element3 = 0.0;

Create ImpulsiveBurn B2;
GMAT B2.CoordinateSystem = SunEarthRot;

Create Variable dv1 dv2 rx ry rz;

Create ReportFile rpt;
GMAT rpt.Filename         = '{summary_path.as_posix()}';
GMAT rpt.WriteHeaders     = false;
GMAT rpt.SolverIterations = None;
GMAT rpt.Precision        = 16;

Create ReportFile traj;
GMAT traj.Filename         = '{traj_path.as_posix()}';
GMAT traj.WriteHeaders     = false;
GMAT traj.SolverIterations = None;
GMAT traj.Precision        = 12;
GMAT traj.WriteReport      = false;
GMAT traj.Add = {{Sat.A1ModJulian, Sat.SunEarthRot.X, Sat.SunEarthRot.Y, Sat.Sun.RMAG, Sat.Earth.RMAG}};

Create OrbitView vw;
GMAT vw.Add                = {{Sat, Earth}};
GMAT vw.CoordinateSystem   = SunEarthRot;
GMAT vw.ViewPointReference = Earth;
GMAT vw.ViewPointVector    = [0 0 150000000];
GMAT vw.ViewDirection      = Earth;
GMAT vw.ViewUpAxis         = Y;
GMAT vw.ShowPlot           = true;
GMAT vw.SolverIterations   = None;

Create Variable q;

BeginMissionSequence;

Target DC {{SolveMode = Solve, ExitMode = SaveAndContinue}};
   Vary DC(B1.Element1 = {g[0]:.9f}, {{Perturbation = 1e-7, MaxStep = 0.2}});
   Vary DC(B1.Element2 = {g[1]:.9f}, {{Perturbation = 1e-7, MaxStep = 0.2}});
   Vary DC(B1.Element3 = 0.0, {{Perturbation = 1e-7, MaxStep = 0.05}});
   Maneuver B1(Sat);
   Propagate Helio(Sat) {{Sat.A1ModJulian = {a1_ins:.9f}}};
   Achieve DC(Sat.SunEarthRot.X = {tgt[0]:.6f}, {{Tolerance = 500}});
   Achieve DC(Sat.SunEarthRot.Y = {tgt[1]:.6f}, {{Tolerance = 500}});
   Achieve DC(Sat.SunEarthRot.Z = {tgt[2]:.6f}, {{Tolerance = 500}});
EndTarget;

dv1 = sqrt(B1.Element1^2 + B1.Element2^2 + B1.Element3^2);

rx = Sat.SunEarthRot.X - {tgt[0]:.6f};
ry = Sat.SunEarthRot.Y - {tgt[1]:.6f};
rz = Sat.SunEarthRot.Z - {tgt[2]:.6f};

B2.Element1 = {tgt[3]:.9f} - Sat.SunEarthRot.VX;
B2.Element2 = {tgt[4]:.9f} - Sat.SunEarthRot.VY;
B2.Element3 = {tgt[5]:.9f} - Sat.SunEarthRot.VZ;

Maneuver B2(Sat);

dv2 = sqrt(B2.Element1^2 + B2.Element2^2 + B2.Element3^2);

Report rpt dv1 dv2 B1.Element1 B1.Element2 B1.Element3 rx ry rz Sat.A1ModJulian;

GMAT traj.WriteReport = true;

For q = 1:1:400;
   Propagate Helio(Sat) {{Sat.ElapsedDays = 1}};
EndFor;
"""

    p = GENERATED_GMAT / f"g02_sc{k}.script"
    p.write_text(s, encoding="utf-8")

    return p


# ---------------------------------------------------------------------------
# GMAT runner
# ---------------------------------------------------------------------------

def run(script):
    if not GMAT_EXE:
        raise RuntimeError(
            "GMAT_BIN is not set. "
            "Set it to the full path of GmatConsole."
        )

    gmat_exe = Path(
        GMAT_EXE
    ).expanduser().resolve()

    if not gmat_exe.is_file():
        raise RuntimeError(
            "GMAT_BIN does not point to an existing file:\n"
            f"  {gmat_exe}"
        )

    script = Path(
        script
    ).resolve()

    r = subprocess.run(
        [
            str(gmat_exe),
            str(script),
        ],
        cwd=script.parent,
        capture_output=True,
        text=True,
        timeout=900,
    )

    ok = (
        r.returncode == 0
        and "Mission run completed" in r.stdout
    )

    if not ok:
        print(
            f"\nGMAT FAILED: "
            f"{script.name}"
        )

        print(
            f"return code: "
            f"{r.returncode}"
        )

        if r.stdout:
            print(
                "\n--- GMAT stdout ---"
            )
            print(
                r.stdout[-8000:]
            )

        if r.stderr:
            print(
                "\n--- GMAT stderr ---"
            )
            print(
                r.stderr[-8000:]
            )

    return ok, r.stdout


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print(
        f"launch A1MJD "
        f"{A1_LAUNCH:.6f} "
        f"({LAUNCH} 00:00 UTC)"
    )

    print(
        "Python guesses:",
        ", ".join(
            f"S{k + 1}"
            f"(G{s['group'] + 1}) "
            f"d{s['day']:.0f} "
            f"{1000 * (s['dv1'] + s['dv2']):.0f}"
            for k, s in enumerate(sc)
        ),
        "m/s",
    )

    # -----------------------------------------------------------------------
    # phase A
    # -----------------------------------------------------------------------

    targets = {}

    for k in range(NSAT):
        p = write_g01(k)

        ok, out = run(p)

        if not ok:
            raise RuntimeError(
                f"GMAT failed during DRO-slot validation "
                f"for slot {k}.\n"
                f"Script:\n  {p}"
            )

        report_path = (
            GMAT_OUT
            / f"g01_slot{k}.txt"
        )

        if not report_path.is_file():
            raise RuntimeError(
                "GMAT completed but did not produce "
                "the expected report:\n"
                f"  {report_path}"
            )

        row = (
            report_path
            .read_text(
                encoding="utf-8"
            )
            .split()
        )

        if len(row) < 9:
            raise RuntimeError(
                f"Unexpected report format in "
                f"{report_path}: "
                f"expected at least 9 values, "
                f"got {len(row)}"
            )

        (
            vy,
            a1,
            X,
            Y,
            Z,
            VX,
            VY,
            VZ,
            rsun,
        ) = map(
            float,
            row[:9],
        )

        targets[k] = (
            X,
            Y,
            Z,
            VX,
            VY,
            VZ,
        )

        print(
            f"g01 slot{k}: "
            f"OK  "
            f"VY* {vy:.4f}  "
            f"state at ins "
            f"({X / 1e6:.2f}, "
            f"{Y / 1e6:.2f}) Mkm"
        )

    # -----------------------------------------------------------------------
    # phase B
    # -----------------------------------------------------------------------

    print()

    results = {}

    for k in range(NSAT):
        p = write_g02(
            k,
            targets[k],
        )

        ok, out = run(p)

        if not ok:
            raise RuntimeError(
                f"GMAT failed during spacecraft "
                f"validation for spacecraft {k}.\n"
                f"Script:\n  {p}"
            )

        summary_path = (
            GMAT_OUT
            / f"g02_sc{k}_summary.txt"
        )

        if not summary_path.is_file():
            raise RuntimeError(
                "GMAT completed but did not produce "
                "the expected summary:\n"
                f"  {summary_path}"
            )

        row = (
            summary_path
            .read_text(
                encoding="utf-8"
            )
            .split()
        )

        if len(row) < 8:
            raise RuntimeError(
                f"Unexpected report format in "
                f"{summary_path}: "
                f"expected at least 8 values, "
                f"got {len(row)}"
            )

        dv1 = float(
            row[0]
        )

        dv2 = float(
            row[1]
        )

        res = np.sqrt(
            float(row[5]) ** 2
            + float(row[6]) ** 2
            + float(row[7]) ** 2
        )

        results[k] = (
            dv1,
            dv2,
            res,
            ok,
        )

        print(
            f"g02 sc{k}: "
            f"OK  "
            f"dv1 {dv1 * 1000:7.1f}  "
            f"dv2 {dv2 * 1000:7.1f}  "
            f"pos residual {res:.0f} km"
        )

    # -----------------------------------------------------------------------
    # phase C
    # -----------------------------------------------------------------------

    print(
        "\n=== GMAT (full ephemeris) "
        "vs Python (2-body) ==="
    )

    print(
        f"{'s/c':>4} "
        f"{'day':>4} | "
        f"{'dv1 py':>7} "
        f"{'dv1 gmat':>8} | "
        f"{'dv2 py':>7} "
        f"{'dv2 gmat':>8} | "
        f"{'tot py':>7} "
        f"{'tot gmat':>8} "
        f"{'diff':>6}"
    )

    for k in range(NSAT):
        d1g, d2g, res, ok = results[k]

        d1p = sc[k]["dv1"]
        d2p = sc[k]["dv2"]

        print(
            f"S{k + 1:<3} "
            f"{sc[k]['day']:>4.0f} | "
            f"{d1p * 1000:7.0f} "
            f"{d1g * 1000:8.0f} | "
            f"{d2p * 1000:7.0f} "
            f"{d2g * 1000:8.0f} | "
            f"{(d1p + d2p) * 1000:7.0f} "
            f"{(d1g + d2g) * 1000:8.0f} "
            f"{((d1g + d2g) - (d1p + d2p)) * 1000:+6.0f}"
        )


if __name__ == "__main__":
    main()