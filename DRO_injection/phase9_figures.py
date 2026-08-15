"""Figure pack for the 9-spacecraft single-launch phasing study (DRO ID 365).

Re-derives the phase9_single_launch.py / departure.py solution (per-member
asymptote if campaign_departure.json has release groups, else the fixed
C3 = 4.0 / 269 deg fallback below; per-member best insertion day), then draws:

    out/phase9_grid.png       3x3 small multiples, one transfer per spacecraft
    out/phase9_combined.png   all transfers, time-graded colour, cleaned up
    out/phase9_timeline.png   deployment timeline (cruise bars + insertion)
    out/phase9_sundist.png    Sun distance vs time per spacecraft (thermal)
    out/phase9_anim.gif       rotating-frame animation of the deployment

All in the planar heliocentric two-body model (see phase9_single_launch.py).
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS_DIR = (HERE / ".." / "scripts").resolve()

sys.path.insert(0, str(SCRIPTS_DIR))

from lambert import lambert, MU, AU, DAY

YEAR = 365.25636 * DAY
N_E  = 2 * np.pi / YEAR
AMP  = 22122321.422754
ECC  = AMP / AU
NSAT = 9
C3   = 4.0
ANG  = np.deg2rad(269.0)
TOF_MIN, DEADLINE = 60.0, 365.0
from paths import PYTHON_OUTPUT, ensure_runtime_dirs
ensure_runtime_dirs()
OUT = PYTHON_OUTPUT
CMAP = plt.get_cmap("viridis")


def kepler_E(M, e):
    E = float(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def earth_state(t):
    th = N_E * t
    r = AU * np.array([np.cos(th), np.sin(th), 0.0])
    v = np.sqrt(MU / AU) * np.array([-np.sin(th), np.cos(th), 0.0])
    return r, v


def slot_state(k, t):
    tp = k * YEAR / NSAT
    w  = N_E * tp
    M  = np.mod(N_E * (t - tp), 2 * np.pi)
    E  = kepler_E(M, ECC)
    nu = 2 * np.arctan2(np.sqrt(1 + ECC) * np.sin(E / 2), np.sqrt(1 - ECC) * np.cos(E / 2))
    r  = AU * (1 - ECC * np.cos(E))
    p  = AU * (1 - ECC**2)
    th = w + nu
    pos = r * np.array([np.cos(th), np.sin(th), 0.0])
    vr  = np.sqrt(MU / p) * ECC * np.sin(nu)
    vt  = np.sqrt(MU / p) * (1 + ECC * np.cos(nu))
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
    M0 = E0 - e * np.sin(E0)
    return a, e, np.arctan2(evec[1], evec[0]), M0


def arc_pos(a, e, w, M0, dt):
    """Inertial position on the transfer ellipse dt seconds after departure."""
    n = np.sqrt(MU / a**3)
    E = kepler_E(np.mod(M0 + n * dt, 2 * np.pi), e)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e * np.cos(E))
    return r * np.array([np.cos(w + nu), np.sin(w + nu), 0.0])


# ------------------------------------------------------------- solve the 9 transfers
t0 = 0.0
rE, vE = earth_state(t0)
vstack_default = vE + np.sqrt(C3) * np.array([np.cos(ANG), np.sin(ANG), 0.0])
days = np.arange(TOF_MIN, DEADLINE + 0.5, 1.0)

# if the campaign departure design exists, use its asymptote(s) + insertion days
# so the figures match the GMAT/spreadsheet solution exactly. Each spacecraft
# may belong to its own release group with its own asymptote (Ariane 6
# supports multiple targeted separations); fall back to a single global
# vinf_vec for older single-asymptote campaign_departure.json files.
import json
_dep_path = PYTHON_OUTPUT / "campaign_departure.json"
DEP = json.load(open(_dep_path)) if os.path.exists(_dep_path) else None


def member_vstack(k):
    if not DEP:
        return vstack_default
    m = DEP["members"][k]
    if "vinf_vec" in m:                 # grouped format: per-member asymptote
        return vE + np.array(m["vinf_vec"])
    return vE + np.array(DEP["vinf_vec"])   # old single-asymptote format


sc = []          # per s/c dict: day, dv1, dv2, elements, insertion state
for k in range(NSAT):
    vstack = member_vstack(k)
    best = None
    for d in ([DEP["members"][k]["day"]] if DEP else days):
        t = t0 + d * DAY
        rT, vT = slot_state(k, t)
        sol = lambert(rE, rT, d * DAY)
        if sol is None:
            continue
        v1, v2 = sol
        dv1 = np.linalg.norm(v1 - vstack)
        dv2 = np.linalg.norm(vT - v2)
        if best is None or dv1 + dv2 < best["dv1"] + best["dv2"]:
            best = dict(day=d, dv1=dv1, dv2=dv2, v1=v1)
    a, e, w, M0 = elements(rE, best["v1"])
    best.update(el=(a, e, w, M0))
    sc.append(best)
    print(f"S{k+1}: day {best['day']:.0f}  dv1 {best['dv1']*1000:.0f}  dv2 {best['dv2']*1000:.0f}  "
          f"tot {(best['dv1']+best['dv2'])*1000:.0f} m/s")

dro_loop = [rotating(slot_state(0, t)[0], t) for t in np.linspace(0, YEAR, 720)]


def draw_base(ax, loop_lw=1.0):
    ax.plot(*zip(*dro_loop), color="0.75", lw=loop_lw, zorder=1)
    ax.plot(0, 0, "o", color="tab:blue", ms=6, zorder=5)
    ax.set_aspect("equal")


# ------------------------------------------------------------- 1. small multiples
fig, axes = plt.subplots(3, 3, figsize=(12, 12.5), sharex=True, sharey=True)
for k, ax in enumerate(axes.flat):
    s = sc[k]
    a, e, w, M0 = s["el"]
    pts = [rotating(arc_pos(a, e, w, M0, q), t0 + q)
           for q in np.linspace(0, s["day"] * DAY, 300)]
    draw_base(ax)
    ax.plot(*zip(*pts), color=CMAP(k / NSAT), lw=1.6, zorder=3)
    t_ins = t0 + s["day"] * DAY
    ax.plot(*rotating(slot_state(k, t_ins)[0], t_ins), "s", color=CMAP(k / NSAT), ms=7, zorder=4)
    ax.set_title(f"S{k+1}  insert day {s['day']:.0f}  "
                 f"{(s['dv1']+s['dv2'])*1000:.0f} m/s", fontsize=10)
    if k == 0:
        ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(6, 5), fontsize=8)
for ax in axes[-1]:
    ax.set_xlabel("rotating x, Mkm")
for ax in axes[:, 0]:
    ax.set_ylabel("rotating y, Mkm")
fig.suptitle("One launch, nine transfers - Sun-Earth rotating frame "
             "(square marker = DRO insertion)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(os.path.join(OUT, "phase9_grid.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------- 2. combined, time-graded
fig, ax = plt.subplots(figsize=(9, 9.5))
draw_base(ax, loop_lw=1.4)
for k, s in enumerate(sc):
    a, e, w, M0 = s["el"]
    qs = np.linspace(0, s["day"] * DAY, 250)
    pts = np.array([rotating(arc_pos(a, e, w, M0, q), t0 + q) for q in qs])
    ax.scatter(pts[:, 0], pts[:, 1], c=qs / DAY, cmap="plasma", s=2, vmin=0, vmax=365, zorder=3)
    t_ins = t0 + s["day"] * DAY
    x, y = rotating(slot_state(k, t_ins)[0], t_ins)
    ax.plot(x, y, "o", mfc="none", mec="k", ms=9, zorder=4)
    ax.annotate(f"S{k+1}", (x, y), textcoords="offset points", xytext=(7, 5), fontsize=9)
sm = plt.cm.ScalarMappable(cmap="plasma", norm=plt.Normalize(0, 365))
fig.colorbar(sm, ax=ax, label="days after launch", shrink=0.75)
ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(7, 5))
ax.arrow(5, 52, -14, 0, head_width=2, color="orange")
ax.annotate("to Sun", (-16, 55), color="orange")
ax.set_xlabel("rotating x, Mkm (sunward negative)")
ax.set_ylabel("rotating y, Mkm")
ax.set_title("Nine transfers from one launch, coloured by time since launch\n"
             "(circles = each spacecraft's DRO insertion point)")
fig.savefig(os.path.join(OUT, "phase9_combined.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------- 3. timeline
fig, ax = plt.subplots(figsize=(9, 4.8))
order = np.argsort([s["day"] for s in sc])
for row, k in enumerate(order):
    s = sc[k]
    ax.barh(row, s["day"], left=0, height=0.55, color=CMAP(k / NSAT), alpha=0.55)
    ax.plot(s["day"], row, "s", color=CMAP(k / NSAT), ms=8)
    ax.annotate(f"{(s['dv1']+s['dv2'])*1000:.0f} m/s",
                (s["day"] + 6, row), va="center", fontsize=8.5)
ax.set_yticks(range(NSAT))
ax.set_yticklabels([f"S{k+1}" for k in order])
ax.axvline(365, color="crimson", lw=1.2, ls="--")
ax.annotate("1 year", (365, NSAT - 0.3), color="crimson", ha="right", fontsize=9,
            textcoords="offset points", xytext=(-4, 0))
ax.set_xlabel("days after launch (bar = cruise, square = DRO insertion)")
ax.set_title("Deployment timeline - all nine inserted inside one year")
ax.set_xlim(0, 400)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase9_timeline.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------- 4. Sun distance vs time
fig, ax = plt.subplots(figsize=(9, 5))
tt = np.arange(0, 500.0, 1.0)
for k, s in enumerate(sc):
    a, e, w, M0 = s["el"]
    r = []
    for d in tt:
        if d <= s["day"]:
            r.append(np.linalg.norm(arc_pos(a, e, w, M0, d * DAY)) / AU)
        else:
            r.append(np.linalg.norm(slot_state(k, t0 + d * DAY)[0]) / AU)
    ax.plot(tt, r, color=CMAP(k / NSAT), lw=1.1, label=f"S{k+1}")
    ax.plot(s["day"], r[int(s["day"])], "s", color=CMAP(k / NSAT), ms=5)
ax.axhline(1 - ECC, color="0.4", lw=0.8, ls=":")
ax.axhline(1 + ECC, color="0.4", lw=0.8, ls=":")
ax.annotate("DRO perihelion 0.852 au", (2, 1 - ECC - 0.006), fontsize=8, color="0.35")
ax.annotate("DRO aphelion 1.148 au", (2, 1 + ECC + 0.002), fontsize=8, color="0.35")
ax.set_xlabel("days after launch")
ax.set_ylabel("Sun distance, au")
ax.set_title("Sun distance per spacecraft (squares = DRO insertion) - cruise thermal envelope")
ax.legend(fontsize=8, ncol=3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase9_sundist.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------- 5. animation
from matplotlib.animation import FuncAnimation, PillowWriter

fig, ax = plt.subplots(figsize=(6.4, 7))
draw_base(ax, loop_lw=1.2)
ax.annotate("Earth", (0, 0), textcoords="offset points", xytext=(6, 5), fontsize=9)
ax.set_xlabel("rotating x, Mkm")
ax.set_ylabel("rotating y, Mkm")
ax.set_xlim(-30, 30)
ax.set_ylim(-52, 52)
dots = [ax.plot([], [], "o", color=CMAP(k / NSAT), ms=7)[0] for k in range(NSAT)]
trails = [ax.plot([], [], "-", color=CMAP(k / NSAT), lw=0.8, alpha=0.6)[0] for k in range(NSAT)]
title = ax.set_title("")
STEP = 2.0
hist = [[] for _ in range(NSAT)]


def pos_at(k, d):
    s = sc[k]
    if d <= s["day"]:
        a, e, w, M0 = s["el"]
        return rotating(arc_pos(a, e, w, M0, d * DAY), t0 + d * DAY)
    return rotating(slot_state(k, t0 + d * DAY)[0], t0 + d * DAY)


def frame(i):
    d = i * STEP
    for k in range(NSAT):
        x, y = pos_at(k, d)
        dots[k].set_data([x], [y])
        hist[k].append((x, y))
        trails[k].set_data(*zip(*hist[k][-90:]))
        dots[k].set_marker("o" if d <= sc[k]["day"] else "D")
    title.set_text(f"day {d:.0f}   ({sum(1 for s in sc if s['day'] <= d)}/9 inserted)")
    return dots + trails + [title]


anim = FuncAnimation(fig, frame, frames=int(500 / STEP), blit=False)
anim.save(os.path.join(OUT, "phase9_anim.gif"), writer=PillowWriter(fps=18), dpi=90)
plt.close(fig)
print("wrote grid / combined / timeline / sundist / anim to out/")
