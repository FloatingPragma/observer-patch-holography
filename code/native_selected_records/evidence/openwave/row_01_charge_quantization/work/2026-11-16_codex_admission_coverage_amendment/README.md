# Admission coverage checker amendment

The independent checker now enforces the exact catalog and registered finite
input domain. The original source, controller, receipt and evidence remain
unchanged under seal
`ba4fc9e36bf41592b272feee134c4bdd54d67c51416c8f992416535d536c9f2f`.
This is verification-only work, with no new native source runs or awards.

The [old-checker probe](GAP_PROBE.json) retains the actual gap: dropping a
consumer column, duplicating a column in place, duplicating one public-stream
history or replacing a different-total input with a duplicate all exited zero.
A simple uncompensated missing row already rejected. The original strict
producer replay continued to compare the full saved evidence; this gap was in
the independently implemented semantic reader.

The [amended checker](check.py) is a complete 148-line reader with no producer
imports. In addition to its existing current, source-law, extraction and
complete finite-prefix checks, it requires:

- All 30 unique `(position, edge, coin)` columns for each of 864 successful
  native histories, giving 25,920 checked columns.
- Exactly one of each registered B/M/decoded-input/public-stream/law case,
  giving 1,752 histories. Initial states and totals come from census currents;
  stream identity comes from the entire observed coin tape.
- The registered vacuum/marker preparation and exact twelve signed prefix
  keys, so an omitted case cannot be concealed by repeating a valid case.

The stream comparison is an audit of the filed finite experiment, separate
from the operational controller. Both same public streams occur for every
input. The checker never uses seed, total, filename or order labels to infer
a case. A shuffled file with all those labels false remains valid. It still
checks that the independent full-state audit agrees with the current-derived
state; a physical audit vector is not an arbitrary metadata label.

Run the read-only headline from this directory:

```sh
python3 -B check.py
```

It reads the unchanged sibling admission packet and local public graph file,
requires Python 3 and NumPy, writes nothing and exits nonzero on failure.
`--evidence PATH` is an explicit semantic-check interface for disposable
counterexamples; the default additionally validates the filed receipt and
evidence digest. All scientific checks remain active under `python -O`.

[VALIDATION.json](VALIDATION.json) contains clean normal/optimized checks,
ten intended-predicate rejections (five malformed records in both modes),
and an optimized shuffled/false-label positive control. Each subprocess has
a 60-second limit. Saved source bytes, sizes, mtimes and inodes are checked
before and after. `python3 -B validate.py` repeats these checks and compares
the saved validation; only `--record-to NEW_FILE` permits exclusive creation.
The old-gap replay is separately available as `python3 -B probe_gaps.py`.

[DESIGN.md](DESIGN.md), [GAP_PREREG.json](GAP_PREREG.json) and
[AMENDMENT_PREREG.json](AMENDMENT_PREREG.json) preserve the order: root's
static finding, frozen old-reader probes, observed old-reader failures,
frozen repair and then amended verification. This is known-result verification,
not group-blind discovery. The native law, controller, admitted source family,
12 executed prefix scope and physical hypotheses are unchanged. The amendment
supports local E1/E2 and D2 coverage only; it does not close external G6,
public-source, final-outgoing or release gates.
