"""Universal-variable Lambert solver (Vallado) and the DRO/Earth state loader.

Shared by transfer_scan.py and pick_arc.py so the solver exists in one place.
"""
import os
from pathlib import Path
import numpy as np

MU = 1.32712440018e11          # km3/s2, Sun
AU = 149597870.691
DAY = 86400.0

# Optional legacy state table. Keep machine-specific locations outside source code.
# Set DRO_STATES_FILE when load_states() is used without an explicit path.
STATES = os.environ.get("DRO_STATES_FILE")


def _c2c3(psi):
    if psi > 1e-6:
        s = np.sqrt(psi)
        return (1 - np.cos(s)) / psi, (s - np.sin(s)) / s**3
    if psi < -1e-6:
        s = np.sqrt(-psi)
        return (1 - np.cosh(s)) / psi, (np.sinh(s) - s) / s**3
    return 0.5, 1.0 / 6.0


def lambert(r1, r2, tof, prograde=True):
    """Single-revolution Lambert. r1, r2 in km, tof in s. Returns (v1, v2) or None."""
    R1, R2 = np.linalg.norm(r1), np.linalg.norm(r2)
    cosdnu = np.clip(np.dot(r1, r2) / (R1 * R2), -1.0, 1.0)
    dnu = np.arccos(cosdnu)
    # the arc must run the same way round the Sun as the Earth does
    if (np.cross(r1, r2)[2] < 0) == prograde:
        dnu = 2 * np.pi - dnu
    dm = 1.0 if dnu < np.pi else -1.0
    A = dm * np.sqrt(R1 * R2 * (1 + cosdnu))
    if abs(A) < 1e-9:
        return None

    psi, lo, hi = 0.0, -4 * np.pi, 4 * np.pi**2
    dt = np.inf
    for _ in range(200):
        c2, c3 = _c2c3(psi)
        y = R1 + R2 + A * (psi * c3 - 1) / np.sqrt(c2)
        if A > 0 and y < 0:
            for _ in range(100):
                psi += 0.1
                c2, c3 = _c2c3(psi)
                y = R1 + R2 + A * (psi * c3 - 1) / np.sqrt(c2)
                if y > 0:
                    break
            lo = psi
        if y < 0:
            return None
        chi = np.sqrt(y / c2)
        dt = (chi**3 * c3 + A * np.sqrt(y)) / np.sqrt(MU)
        lo, hi = (psi, hi) if dt <= tof else (lo, psi)
        new = (hi + lo) / 2
        if abs(new - psi) < 1e-12:
            psi = new
            break
        psi = new

    c2, c3 = _c2c3(psi)
    y = R1 + R2 + A * (psi * c3 - 1) / np.sqrt(c2)
    if y <= 0:
        return None
    f = 1 - y / R1
    g = A * np.sqrt(y / MU)
    gdot = 1 - y / R2
    if abs(g) < 1e-12 or abs(dt - tof) / tof > 1e-4:
        return None
    return (r2 - f * r1) / g, (gdot * r2 - r1) / g


def load_states(path=None):
    """Return t, DRO state, and Earth state from a legacy state table.

    Pass ``path`` explicitly or set the ``DRO_STATES_FILE`` environment variable.
    """
    source = path or STATES
    if not source:
        raise ValueError("State table not configured; pass path or set DRO_STATES_FILE")
    D = np.genfromtxt(Path(source), skip_header=1)
    D = D[np.argsort(D[:, 1])]
    a1mjd, t = D[:, 0], D[:, 1]
    sc_r, sc_v = D[:, 2:5], D[:, 5:8]          # spacecraft from Sun
    rel_r, rel_v = D[:, 8:11], D[:, 11:14]     # spacecraft from Earth, same axes
    return a1mjd, t, sc_r, sc_v, sc_r - rel_r, sc_v - rel_v
