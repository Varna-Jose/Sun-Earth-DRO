"""Two-panel deployment snapshot for the 9-s/c single-launch phasing study.

Left panel : the day the FIRST spacecraft inserts (rest still in transfer)
Right panel: the day the LAST spacecraft inserts (deployment complete)

Filled diamonds = in the DRO.  Open circles = still cruising (thin line =
transfer path flown so far).  Same solution as phase9_single_launch.py
(shared asymptote C3 4.0 at 269 deg, per-member best insertion day).
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from paths import PYTHON_OUTPUT, ensure_runtime_dirs
ensure_runtime_dirs()
from lambert import lambert, MU, AU, DAY

YEAR = 365.25636 * DAY
N_E  = 2 * np.pi / YEAR
AMP  = 22122321.422754
ECC  = AMP / AU
NSAT = 9
C3, ANG = 4.0, np.deg2rad(269.0)
CMAP = plt.get_cmap("viridis")


def kepler_E(M, e):
    E = float(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def earth_state(t):
    th = N_E * t
    return (AU * np.array([np.cos(th), np.sin(th), 0.0]),
            np.sqrt(MU / AU) * np.array([-np.sin(th), np.cos(th), 0.0]))


def slot_state(k, t):
    tp = k * YEAR / NSAT
    w = N_E * tp
    M = np.mod(N_E * (t - tp), 2 * np.pi)
    E = kepler_E(M, ECC)
    nu = 2 * np.arctan2(np.sqrt(1 + ECC) * np.sin(E / 2), np.sqrt(1 - ECC) * np.cos(E / 2))
    r = AU * (1 - ECC * np.cos(E))
    p = AU * (1 - ECC**2)
    th = w + nu
    pos = r * np.array([np.cos(th), np.sin(th), 0.0])
    vr = np.sqrt(MU / p) * ECC * np.sin(nu)
    vt = np.sqrt(MU / p) * (1 + ECC * np.cos(nu))
    vel = vr * np.array([np.cos(th), np.sin(th), 0.0]) + vt * np.array([-np.sin(th), np.cos(th), 0.0])
    return pos, vel


def rotating(pos, t):
    th = N_E * t
    c, s = np.cos(-th), np.sin(-th)
    return (c * pos[0] - s * pos[1] - AU) / 1e6, (s * pos[0] + c * pos[1]) / 1e6


def elements(r0, v0):
    h = np.cross(r0, v0)
    evec = np.cross(v0, h) / MU - r0 / np.linalg.norm(r0)
    e = np.linalg.norm(evec)
    a = 1 / (2 / np.linalg.norm(r0) - np.linalg.norm(v0)**2 / MU)
    nu0 = np.arctan2(np.dot(np.cross(evec, r0), h) / np.linalg.norm(h), np.dot(evec, r0))
    E0 = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu0 / 2), np.sqrt(1 + e) * np.cos(nu0 / 2))
    return a, e, np.arctan2(evec[1], evec[0]), E0 - e * np.sin(E0)


def arc_pos(el, dt):
    a, e, w, M0 = el
    E = kepler_E(np.mod(M0 + np.sqrt(MU / a**3) * dt, 2 * np.pi), e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e * np.cos(E))
    return r * np.array([np.cos(w + nu), np.sin(w + nu), 0.0])


# ----- re-derive the chosen transfers (fast: one Lambert per day per slot)
t0 = 0.0
rE, vE = earth_state(t0)
vstack_default = vE + np.sqrt(C3) * np.array([np.cos(ANG), np.sin(ANG), 0.0])

import json
_dep_path = PYTHON_OUTPUT / "campaign_departure.json"
DEP = json.loads(_dep_path.read_text()) if _dep_path.exists() else None


def member_vstack(k):
    if not DEP:
        return vstack_default
    m = DEP["members"][k]
    if "vinf_vec" in m:                 # grouped format: per-member asymptote
        return vE + np.array(m["vinf_vec"])
    return vE + np.array(DEP["vinf_vec"])   # old single-asymptote format


sc = []
for k in range(NSAT):
    vstack = member_vstack(k)
    best = None
    for d in ([DEP["members"][k]["day"]] if DEP else np.arange(60.0, 365.5, 1.0)):
        rT, vT = slot_state(k, t0 + d * DAY)
        sol = lambert(rE, rT, d * DAY)
        if sol is None:
            continue
        cost = np.linalg.norm(sol[0] - vstack) + np.linalg.norm(vT - sol[1])
        if best is None or cost < best[0]:
            best = (cost, d, sol[0])
    sc.append(dict(day=best[1], el=elements(rE, best[2])))

first = min(s["day"] for s in sc)
last  = max(s["day"] for s in sc)
loop = [rotating(slot_state(0, t)[0], t) for t in np.linspace(0, YEAR, 720)]

fig, axes = plt.subplots(1, 2, figsize=(13, 7.6), sharex=True, sharey=True)
for ax, snap in zip(axes, [first, last]):
    n_in = sum(1 for s in sc if s["day"] <= snap)
    ax.plot(*zip(*loop), color="0.75", lw=1.2, zorder=1)
    for k, s in enumerate(sc):
        col = CMAP(k / NSAT)
        if s["day"] <= snap:
            x, y = rotating(slot_state(k, t0 + snap * DAY)[0], t0 + snap * DAY)
            ax.plot(x, y, "D", color=col, ms=9, zorder=4)
        else:
            qs = np.linspace(0, snap * DAY, 200)
            pts = [rotating(arc_pos(s["el"], q), t0 + q) for q in qs]
            ax.plot(*zip(*pts), color=col, lw=0.8, alpha=0.45, zorder=2)
            x, y = pts[-1]
            ax.plot(x, y, "o", mfc="none", mec=col, ms=8, mew=1.6, zorder=4)
        ax.annotate(f"S{k+1}", (x, y), textcoords="offset points", xytext=(7, 5), fontsize=8.5)
    ax.plot(0, 0, "o", color="tab:blue", ms=7, zorder=5)
    ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(7, 5), fontsize=9)
    ax.set_title(f"Day {snap:.0f} - {n_in}/9 inserted" +
                 ("  (first insertion)" if snap == first else "  (deployment complete)"))
    ax.set_xlabel("rotating x, Mkm (sunward negative)")
    ax.set_aspect("equal")
axes[0].set_ylabel("rotating y, Mkm")
axes[0].arrow(8, 52, -14, 0, head_width=2, color="orange")
axes[0].annotate("to Sun", (-14, 55), color="orange")
h = [plt.Line2D([], [], marker="D", color="0.3", ls="", ms=8, label="in the DRO"),
     plt.Line2D([], [], marker="o", mfc="none", mec="0.3", ls="", ms=8, label="still in transfer")]
axes[0].legend(handles=h, loc="lower right", fontsize=9)
fig.suptitle("Deployment snapshots - one launch, nine spacecraft, DRO ID 365 (Sun-Earth rotating frame)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = PYTHON_OUTPUT / "phase9_snapshot.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out, f"(panels at day {first:.0f} and day {last:.0f})")
