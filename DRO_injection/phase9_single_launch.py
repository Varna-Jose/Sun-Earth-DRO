"""9-spacecraft single-launch phasing into the JPL Sun-Earth DRO ID 365.

Question answered
    All 9 spacecraft leave on ONE Ariane 64 (shared escape asymptote, C3 <= 20).
    Each then makes two impulsive burns of its own:
        burn 1  departure trim, just outside the SOI  (re-shapes its heliocentric orbit)
        burn 2  DRO insertion, anywhere along its own slot's track
    All 9 insertions must happen within DEADLINE days of launch.
    What is the per-spacecraft delta-v, and does it fit the 2700 m/s budget?

Model
    Planar heliocentric two-body, Earth circular at 1 AU.  The DRO is the
    equivalent heliocentric ellipse (a = 1 AU, e = amplitude/AU) locked 1:1 to
    Earth -- the approximation validated against GMAT/FD to 0.73 % on flux
    (orbit_gmat/README.md).  Slot k has its perihelion (sunward apse) epoch
    offset by k*T/9.  Burn-1 is costed with NO Oberth discount (applied at the
    SOI edge) -- conservative; a perigee burn would be cheaper.

Method
    For each slot k and each candidate insertion date t (60..DEADLINE days,
    1-day grid), solve Lambert from Earth@t0 to slot-k state@t.  This gives
    dv2(k,t) and the required heliocentric departure velocity.  The shared
    launcher asymptote v_inf is then chosen (2-D grid, |v_inf| <= sqrt(20))
    to minimise the WORST member's total dv1+dv2, where each member takes its
    own best insertion date.

Outputs (out/)
    phase9_summary.csv    per-member: slot, insertion day, dv1, dv2, total
    phase9_rotating.png   rotating-frame picture of the 9 transfers
    printed summary       worst / mean member, chosen C3, sensitivity to deadline
"""
import os, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from paths import PYTHON_OUTPUT, ensure_runtime_dirs
ensure_runtime_dirs()
from lambert import lambert, MU, AU, DAY

YEAR   = 365.25636 * DAY          # sidereal year, s
N_E    = 2 * np.pi / YEAR         # Earth mean motion
AMP    = 22122321.422754          # km, DRO ID 365 max sunward excursion (Varna slide)
ECC    = AMP / AU                 # equivalent heliocentric eccentricity
NSAT   = 9
C3MAX  = 20.0                     # km2/s2, their stated launcher cap
DEADLINE = 365.0                  # days, all insertions inside this
TOF_MIN  = 60.0                   # days, commissioning before burn 2
OUT = str(PYTHON_OUTPUT)


def kepler_E(M, e):
    E = M.copy() if isinstance(M, np.ndarray) else float(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def earth_state(t):
    th = N_E * t
    r = AU * np.array([np.cos(th), np.sin(th), 0.0])
    v = np.sqrt(MU / AU) * np.array([-np.sin(th), np.cos(th), 0.0])
    return r, v


def slot_state(k, t):
    """Heliocentric state of slot k's DRO track at time t (s)."""
    tp = k * YEAR / NSAT                    # perihelion epoch of slot k
    w  = N_E * tp                           # perihelion longitude (mean-longitude lock)
    M  = N_E * (t - tp)
    E  = kepler_E(np.mod(M, 2 * np.pi), ECC)
    a  = AU
    nu = 2 * np.arctan2(np.sqrt(1 + ECC) * np.sin(E / 2), np.sqrt(1 - ECC) * np.cos(E / 2))
    r  = a * (1 - ECC * np.cos(E))
    p  = a * (1 - ECC**2)
    th = w + nu
    pos = r * np.array([np.cos(th), np.sin(th), 0.0])
    vr  = np.sqrt(MU / p) * ECC * np.sin(nu)
    vt  = np.sqrt(MU / p) * (1 + ECC * np.cos(nu))
    vel = vr * np.array([np.cos(th), np.sin(th), 0.0]) + vt * np.array([-np.sin(th), np.cos(th), 0.0])
    return pos, vel


def rotating(pos, t):
    """Heliocentric inertial -> Earth-centred Sun-Earth rotating frame (Mkm)."""
    th = N_E * t
    c, s = np.cos(-th), np.sin(-th)
    x = c * pos[0] - s * pos[1] - AU
    y = s * pos[0] + c * pos[1]
    return x / 1e6, y / 1e6


# ---------------------------------------------------------------- 1. Lambert pre-scan
t0 = 0.0
rE, vE = earth_state(t0)
days = np.arange(TOF_MIN, DEADLINE + 0.5, 1.0)

vdep = np.full((NSAT, len(days), 3), np.nan)   # required heliocentric departure velocity
dv2  = np.full((NSAT, len(days)), np.nan)      # insertion burn

for k in range(NSAT):
    for j, d in enumerate(days):
        t = t0 + d * DAY
        rT, vT = slot_state(k, t)
        sol = lambert(rE, rT, d * DAY)
        if sol is None:
            continue
        v1, v2 = sol
        vdep[k, j] = v1
        dv2[k, j] = np.linalg.norm(vT - v2)

# ---------------------------------------------------------------- 2. shared-asymptote search
best = None
for c3 in np.arange(0.0, C3MAX + 0.01, 0.25):  # physical floor is 0
    vinf = np.sqrt(c3)
    for ang in np.deg2rad(np.arange(0, 360, 1.0)):
        vstack = vE + vinf * np.array([np.cos(ang), np.sin(ang), 0.0])
        dv1 = np.linalg.norm(vdep - vstack, axis=2)          # (NSAT, ndays)
        tot = dv1 + dv2
        pick = np.nanargmin(tot, axis=1)                     # best insertion day per member
        per = tot[np.arange(NSAT), pick]
        worst = per.max()
        if best is None or worst < best["worst"]:
            best = dict(worst=worst, c3=c3, ang=np.rad2deg(ang), pick=pick, per=per,
                        dv1=dv1[np.arange(NSAT), pick], dv2=dv2[np.arange(NSAT), pick])

pick, per = best["pick"], best["per"]
print(f"shared asymptote: C3 = {best['c3']:.1f} km2/s2, in-plane angle {best['ang']:.0f} deg "
      f"(vs Earth velocity direction)")
print(f"worst member total = {best['worst']*1000:.0f} m/s   mean = {per.mean()*1000:.0f} m/s\n")
print(f"{'slot':>4} {'insert day':>10} {'dv1 m/s':>8} {'dv2 m/s':>8} {'total m/s':>9}")
rows = []
for k in range(NSAT):
    d = days[pick[k]]
    print(f"{k:>4} {d:>10.0f} {best['dv1'][k]*1000:>8.0f} {best['dv2'][k]*1000:>8.0f} {per[k]*1000:>9.0f}")
    rows.append((k, d, best['dv1'][k] * 1000, best['dv2'][k] * 1000, per[k] * 1000))

os.makedirs(OUT, exist_ok=True)
import csv
with open(os.path.join(OUT, "phase9_summary.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["slot", "insertion_day", "dv1_ms", "dv2_ms", "total_ms",
                f"C3_{best['c3']:.1f}", f"asymptote_deg_{best['ang']:.0f}", f"deadline_d_{DEADLINE:.0f}"])
    w.writerows(rows)

# ---------------------------------------------------------------- 3. checks
# sanity: cheapest single-member insertion should reproduce the known one-off
# transfer (~C3 12.5, dv2 ~0.5-0.7 km/s at a perpendicular crossing)
for k in range(NSAT):
    j = np.nanargmin(dv2[k])
    c3_solo = np.linalg.norm(vdep[k, j] - vE) ** 2
    print(f"solo slot {k}: best dv2 {dv2[k, j]*1000:5.0f} m/s at day {days[j]:3.0f}, "
          f"needs C3 {c3_solo:5.1f} km2/s2")

# ---------------------------------------------------------------- 4. figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(9, 9))
tt = np.linspace(0, YEAR, 720)
dx, dy = zip(*[rotating(slot_state(0, t)[0], t) for t in tt])
ax.plot(dx, dy, color="0.65", lw=1.2, label="DRO ID 365 (rotating frame)")
cmap = plt.get_cmap("viridis")
for k in range(NSAT):
    d = days[pick[k]]
    t_ins = t0 + d * DAY
    sol = lambert(rE, slot_state(k, t_ins)[0], d * DAY)
    v1 = sol[0]
    # propagate the transfer arc by Kepler stepping (elements from r,v)
    from numpy.linalg import norm
    r0, vv = rE.copy(), v1.copy()
    npts = 300
    arc = []
    # universal-variable propagation via repeated small Lambert-free steps: use f&g series with eccentric orbit elements
    hvec = np.cross(r0, vv); h = norm(hvec)
    evec = np.cross(vv, hvec) / MU - r0 / norm(r0); e = norm(evec)
    a = 1 / (2 / norm(r0) - norm(vv)**2 / MU)
    n = np.sqrt(MU / a**3)
    nu0 = np.arctan2(np.dot(np.cross(evec, r0), hvec) / (h), np.dot(evec, r0))
    E0 = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu0 / 2), np.sqrt(1 + e) * np.cos(nu0 / 2))
    M0 = E0 - e * np.sin(E0)
    w_t = np.arctan2(evec[1], evec[0])
    for tq in np.linspace(0, d * DAY, npts):
        Mq = M0 + n * tq
        Eq = kepler_E(np.mod(Mq, 2 * np.pi), e)
        nuq = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Eq / 2), np.sqrt(1 - e) * np.cos(Eq / 2))
        rq = a * (1 - e * np.cos(Eq))
        posq = rq * np.array([np.cos(w_t + nuq), np.sin(w_t + nuq), 0.0])
        arc.append(rotating(posq, t0 + tq))
    ax.plot(*zip(*arc), color=cmap(k / NSAT), lw=1.0)
    ax.plot(*rotating(slot_state(k, t_ins)[0], t_ins), "o", color=cmap(k / NSAT), ms=6,
            label=f"S{k+1} day {d:.0f}, {per[k]*1000:.0f} m/s")
ax.plot(0, 0, "o", color="tab:blue", ms=8)
ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(8, 6))
ax.arrow(-5, 52, -12, 0, head_width=2, color="orange")
ax.annotate("to Sun", (-20, 55), color="orange")
ax.set_xlabel("rotating x, Mkm (sunward negative)")
ax.set_ylabel("rotating y, Mkm")
ax.set_title(f"9 s/c, one launch, all inserted < {DEADLINE:.0f} d - worst {best['worst']*1000:.0f} m/s, "
             f"C3 {best['c3']:.1f}")
ax.legend(fontsize=7, loc="lower right")
ax.set_aspect("equal")
fig.savefig(os.path.join(OUT, "phase9_rotating.png"), dpi=150, bbox_inches="tight")
print("\nwrote out/phase9_summary.csv, out/phase9_rotating.png")
