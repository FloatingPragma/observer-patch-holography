# Recording an unresolved clock in the donor release

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The sealed donor release replays the deposited sixteen signed clock settings and the zero-clock control. Its driver aborts when a well-formed clock bank fails to identify a magnitude on a nonzero setting. The [operating protocol](DONOR_RELEASE.md) instead calls that result inconclusive. The [checked acquisition](donor_release_checked.py) records the bank and its disposition before any release. It reuses the original clock decoder, native preparation, detector and flight laws. Original receipts and manifest files retain their checksums.

A valid unresolved bank gives `inconclusive_clock`, regardless of the simulator setting. A malformed count bank gives `invalid_clock_bank`. An absent bank gives `missing_clock_bank`. A resolved interval that misses the known simulated magnitude gives `clock_rejected` through a separate simulator check. Only calibrated banks run the original pulse programs. A failed quantitative release or changed-law gate gives `scientific_rejection` with the associated outcomes. Success of the executable means records were generated or verified. The scientific dispositions are in the records.

Run from the repository root to verify the complete resolved results and zero-clock disposition against the sealed archive.

~~~sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.donor_release_checked
~~~

To examine another clock acquisition on the same seventeen scheduled settings, retain the slot order from `clock.json` and mark an absent individual bank with an empty object in the corresponding `runs` entry. The input uses the same declared number of trials and time settings. Give a separate, unoccupied file path for the output.

~~~sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.donor_release_checked \
    --clock-receipt path/to/clock-counts.json --record path/to/checked-release.json
~~~

The focused regression uses a valid all-zero bank at strength $+1/2$ and an out-of-range count at strength $-1/2$. Both original banks and their separate dispositions appear in the revised evidence. The next calibrated slot still releases both donor columns. The checked acquisition is a local research entry point. A prospective filing requires authorization, an immutable nomination and a registered replay on the platform. [Regression](test_donor_release_checked.py)
