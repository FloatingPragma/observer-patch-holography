# Native count-clock correction

The original separated signed-pair source now has a source-derived additive
generator under an explicitly implemented Poisson event-count clock. All
four separated compact products were checked on their own native source:
generator value12 equals the two single values6+6. The old fixed-attempt
logarithmic excess remains correct for its earlier clock.

All9 prospective finite gates pass. The3.352-second run executes242,175
native calls:100,800 source-mode columns,131,070 coin-prefix branches,
9,800 auxiliary events and505 clocked body events. Eight64-tick histories
complete. The exact16-bit prefix instrument retains failure mass1/8192;
the128-bit/count32 sampler has a rigorous failure bound3.185e−37 per count.
No failed or unresolved mass is normalized away.

The [proof](PROOF.md) distinguishes the exact ideal Poisson clock from its
complete finite flagged implementation. The finite rare-event alternative
retains positive excess at every nonzero step and bounds its convergence
to zero. [Receipt](RECEIPT.json), [records](RECORDS.json.gz),
[freeze](PREREG.json), [independent check](INDEPENDENT_CHECK.json).

Run offline using only Python's standard library:

```
PYTHONDONTWRITEBYTECODE=1 python3 headline.py
PYTHONDONTWRITEBYTECODE=1 python3 verify.py
```

The independent implementation reconstructs every event using sparse
floor/ceiling updates and verifies the exponential bounds using a separate
positive-series reciprocal enclosure. It passes in1.305seconds and imports
no producer modules. Both implementations have the same Codex author; this
is not other-group review.

This source-count generator does not itself establish physical rest energy,
unitary evolution or unread quantum persistence. Its coherent U_L reader
is a separate declared mode being implemented in the energy packet. That
consumer must retain the clock, preparation and readout errors. The public
row2 decision remains PR723's🚧; this local result changes no pin or award.
