"""Departure design for the configured launch date (campaign_config.LAUNCH).

Ariane 6's restartable Vinci + APU support multiple, distinct targeted
separation events on one flight, so instead of one shared post-Ariane
asymptote for all 9 spacecraft, this optimises an INDEPENDENT shared
asymptote per release GROUP. Groups are contiguous in slot-phase order
(natural for a sequential dispenser release) and were chosen by brute-force
search over all contiguous partitions of the 9 slots, minimising the worst
member's total delta-v (see notebook history - the 3-way split below beat
1-way by ~2,200 m/s and 2-way by ~280 m/s on the worst member; a 4th group
was not worth the extra release event).

Also fixes a bug in the old single-asymptote version: the C3 search floor
was hardcoded at 4.0 km2/s2, and the found optimum sat exactly on that
boundary. C3 = 0 (no launcher-imparted excess velocity beyond Earth's own)
is the true physical floor and is included in the grid below.

Writes out/campaign_departure.json for the rest of the chain (figure
scripts, GMAT generators, results package). Each member now carries its
OWN vinf_vec (its group's asymptote), not a single global one.
"""
import os, sys, json
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
sys.path.insert(0, str(HERE))
from paths import PYTHON_OUTPUT, ensure_runtime_dirs
ensure_runtime_dirs()
from lambert import lambert, MU, AU, DAY
from campaign_config import LAUNCH
from scan_launch_dates import earth_departure_velocity, slot_state, days, vdep, dv2, NSAT

C3_LO, C3_HI, C3_STEP = 0.0, 20.0, 0.25     # C3=0 is the physical floor, not 4.0
ANG_STEP_DEG = 1.0

# Contiguous 3-way split in slot-phase order (0-indexed slots = S1..S9):
# G1: S1-S2, G2: S3-S5, G3: S6-S9 - found by brute-force search over all
# contiguous partitions, refined on the fine C3/angle grid.
GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7, 8]]

vE = earth_departure_velocity(LAUNCH)


def optimize_group(members):
    """Best shared asymptote (C3, angle) for just this subset of slots."""
    best = None
    for c3 in np.arange(C3_LO, C3_HI + 0.01, C3_STEP):
        vinf = np.sqrt(c3)
        for ang in np.deg2rad(np.arange(0, 360, ANG_STEP_DEG)):
            vstack = vE + vinf * np.array([np.cos(ang), np.sin(ang), 0.0])
            dv1 = np.linalg.norm(vdep[members] - vstack, axis=2)
            tot = dv1 + dv2[members]
            pick = np.nanargmin(tot, axis=1)
            per = tot[np.arange(len(members)), pick]
            worst = per.max()
            if best is None or worst < best["worst"]:
                best = dict(worst=float(worst), mean=float(per.mean()),
                            c3=float(c3), ang=float(np.rad2deg(ang)),
                            pick=pick, per=per, vstack=vstack)
    return best


members = [None] * NSAT
group_summaries = []
for gi, group in enumerate(GROUPS):
    b = optimize_group(group)
    vinf_vec = b["vstack"] - vE
    for i, k in enumerate(group):
        j = b["pick"][i]
        dv1v = vdep[k, j] - b["vstack"]
        members[k] = dict(slot=k, group=gi, day=float(days[j]),
                          vinf_vec=[float(x) for x in vinf_vec],
                          dv1v=[float(x) for x in dv1v],
                          dv1=float(np.linalg.norm(dv1v)), dv2=float(dv2[k, j]),
                          total=float(np.linalg.norm(dv1v) + dv2[k, j]))
    group_summaries.append(dict(group=gi, slots=[k + 1 for k in group],
                                C3=b["c3"], ang_deg=b["ang"],
                                vinf_vec=[float(x) for x in vinf_vec],
                                worst_ms=b["worst"] * 1000, mean_ms=b["mean"] * 1000))

worst_ms = max(m["total"] for m in members) * 1000
mean_ms = np.mean([m["total"] for m in members]) * 1000

out = dict(launch=str(LAUNCH), vE_kms=[float(x) for x in vE],
           n_groups=len(GROUPS), groups=group_summaries,
           worst_ms=float(worst_ms), mean_ms=float(mean_ms), members=members)
PYTHON_OUTPUT.mkdir(parents=True, exist_ok=True)
with (PYTHON_OUTPUT / "campaign_departure.json").open("w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)

print(f"launch {LAUNCH}: |vE| = {np.linalg.norm(vE):.3f} km/s (circular mean 29.785)")
print(f"{len(GROUPS)} release groups (sequential Ariane 6 APU separations):")
for g in group_summaries:
    print(f"  G{g['group']+1} slots {g['slots']}: C3 = {g['C3']:.2f} km2/s2, "
          f"angle {g['ang_deg']:.0f} deg  (worst {g['worst_ms']:.0f}, "
          f"mean {g['mean_ms']:.0f} m/s)")
print()
for m in members:
    print(f"  S{m['slot']+1} (G{m['group']+1}): day {m['day']:.0f}  dv1 {m['dv1']*1000:5.0f}  "
          f"dv2 {m['dv2']*1000:5.0f}  total {m['total']*1000:5.0f} m/s")
print(f"\nworst member {worst_ms:.0f} m/s, fleet mean {mean_ms:.0f} m/s -> output/python/campaign_departure.json")
