"""Rank candidate launch dates for the 9-s/c single-launch phasing deployment.

A Sun-Earth DRO has NO launch window - C3 and the transfer geometry repeat
daily. The only date-dependent term at this fidelity is Earth's DEPARTURE
VELOCITY: orbital eccentricity modulates it by +/-0.5 km/s through the year
(perihelion ~3 Jan: fastest; aphelion ~4 Jul: slowest), which shifts the
shared-asymptote / burn-1 costs by a few hundred m/s.

Model: the circular-Earth phasing machinery of phase9_single_launch.py, with
the departure velocity replaced by Earth's true velocity vector on the
candidate date (radial + transverse components from the eccentric orbit).
Slots and positions stay ideal; GMAT provides truth at the chosen date.

Output: a ranked table (1st + 15th of each 2045 month) and the recommended
date to put into campaign_config.py.
"""
import sys, os
from datetime import date, timedelta
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
from lambert import lambert, MU, AU, DAY

YEAR = 365.25636 * DAY
N_E = 2 * np.pi / YEAR
AMP = 22122321.422754
ECC_DRO = AMP / AU
NSAT = 9
C3MAX = 20.0
E_E = 0.0167086                      # Earth orbital eccentricity
PERIHELION = date(2045, 1, 3)        # Earth perihelion passage (approx)


def kepler_E(M, e):
    E = float(M)
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E


def earth_departure_velocity(d):
    """Earth heliocentric velocity vector on date d, expressed in the scan
    frame (Earth on +x axis, transverse = +y). Radial/transverse components
    from the eccentric orbit; magnitude 29.29 (Jul) - 30.29 (Jan) km/s."""
    M = 2 * np.pi * ((d - PERIHELION).days % 365.25636) / 365.25636
    E = kepler_E(M, E_E)
    nu = 2 * np.arctan2(np.sqrt(1 + E_E) * np.sin(E / 2), np.sqrt(1 - E_E) * np.cos(E / 2))
    p = AU * (1 - E_E**2)
    v_r = np.sqrt(MU / p) * E_E * np.sin(nu)
    v_t = np.sqrt(MU / p) * (1 + E_E * np.cos(nu))
    return np.array([v_r, v_t, 0.0])


def slot_state(k, t):
    tp = k * YEAR / NSAT
    w = N_E * tp
    M = np.mod(N_E * (t - tp), 2 * np.pi)
    E = kepler_E(M, ECC_DRO)
    nu = 2 * np.arctan2(np.sqrt(1 + ECC_DRO) * np.sin(E / 2), np.sqrt(1 - ECC_DRO) * np.cos(E / 2))
    r = AU * (1 - ECC_DRO * np.cos(E))
    p = AU * (1 - ECC_DRO**2)
    th = w + nu
    pos = r * np.array([np.cos(th), np.sin(th), 0.0])
    vr = np.sqrt(MU / p) * ECC_DRO * np.sin(nu)
    vt = np.sqrt(MU / p) * (1 + ECC_DRO * np.cos(nu))
    vel = vr * np.array([np.cos(th), np.sin(th), 0.0]) + vt * np.array([-np.sin(th), np.cos(th), 0.0])
    return pos, vel


# Lambert pre-scan is date-independent in this frame (positions ideal):
rE = AU * np.array([1.0, 0.0, 0.0])
days = np.arange(60.0, 365.5, 1.0)
vdep = np.full((NSAT, len(days), 3), np.nan)
dv2 = np.full((NSAT, len(days)), np.nan)
for k in range(NSAT):
    for j, d in enumerate(days):
        rT, vT = slot_state(k, d * DAY)
        sol = lambert(rE, rT, d * DAY)
        if sol is None:
            continue
        vdep[k, j] = sol[0]
        dv2[k, j] = np.linalg.norm(vT - sol[1])

def _scan():
    results = []
    cand = [date(2045, m, dd) for m in range(1, 13) for dd in (1, 15)]
    for D in cand:
        vE = earth_departure_velocity(D)
        best = None
        for c3 in np.arange(4.0, C3MAX + 0.01, 1.0):
            vinf = np.sqrt(c3)
            for ang in np.deg2rad(np.arange(0, 360, 2.0)):
                vstack = vE + vinf * np.array([np.cos(ang), np.sin(ang), 0.0])
                dv1 = np.linalg.norm(vdep - vstack, axis=2)
                tot = dv1 + dv2
                worst = np.nanmin(tot, axis=1).max()
                mean = np.nanmin(tot, axis=1).mean()
                if best is None or worst < best[0]:
                    best = (worst, mean, c3, np.rad2deg(ang))
        results.append((D, *best))
        print(f"{D}  worst {best[0]*1000:5.0f} m/s   mean {best[1]*1000:5.0f} m/s   "
              f"C3 {best[2]:4.1f}  ang {best[3]:3.0f}")

    results.sort(key=lambda r: r[1])
    D, worst, mean, c3, ang = results[0]
    print(f"\nBEST 2045 EPOCH: {D}  (worst member {worst*1000:.0f} m/s, mean {mean*1000:.0f} m/s)")
    print(f"WORST 2045 EPOCH: {results[-1][0]}  (worst member {results[-1][1]*1000:.0f} m/s)")
    print(f"spread best-to-worst across the year: {(results[-1][1]-worst)*1000:.0f} m/s")
    print("\n-> set campaign_config.LAUNCH to the best epoch and run run_campaign.py")


if __name__ == "__main__":
    _scan()
