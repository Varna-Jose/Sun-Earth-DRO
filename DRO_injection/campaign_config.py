"""Single source of truth for the phasing campaign: the launch epoch.

Change LAUNCH and re-run run_campaign.py - everything (GMAT scripts, replay,
figures, RESULTS package, spreadsheet dates) rebuilds for the new date.
"""
from datetime import date

LAUNCH = date(2045, 1, 1)           # best 2045 epoch per scan_launch_dates.py
                                    # (worst member 5161 m/s; year spread only 314 m/s)

A1_ANCHOR = 31407.00042863868       # A1MJD of 01 Jan 2027 12:00:00 UTC
def a1mjd(d):
    """A1 modified Julian date of 00:00 UTC on calendar date d."""
    return A1_ANCHOR + (d - date(2027, 1, 1)).days - 0.5

A1_LAUNCH = a1mjd(LAUNCH)
