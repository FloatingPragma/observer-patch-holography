# Native Maxwell and Klein–Gordon row evidence

This package preserves the full-row decisions and independent evidence for
OpenWave rows29 and30 on their declared source domains. Each decision binds
two independent full-criterion reviews. Four isolated modes reproduce the
complete additional Maxwell review, its size-three frame/native instrument,
and both Klein–Gordon reviews. The original scientific files are unchanged.

Use Python3.13 and the exact numerical dependencies in `requirements.txt`:

```sh
python3 code/openwave_wave_rows/replay.py --mode all --output /tmp/wave-replay.json
```

Each mode copies the frozen snapshot into a temporary directory, deletes its
scientific outputs, and regenerates them with site startup disabled. The
controller permits the listed numerical runtime and denies outside project
reads, subprocesses and network access. Input hashes and output bytes are
checked before and after execution. The second Klein–Gordon review's private
input layout is materialized solely from the same flat pinned snapshot.

Row29 concerns the free noncompact radiative source, exact oriented curl,
transverse waves and native first-moment/clock instruments. Its coefficients
are declared and chosen for unit speed. Row30 concerns the scalar mass-coin
walk, the full slow subspace, native preparation/readout and a controlled
massive scalar limit. Its carrier-coordinate speed is1/sqrt(48). Coordinate
rescaling does not prove a shared physical speed or an interacting Dirac
transfer. Architecture, mass, preparation, coherent control and clock inputs
remain explicit. Finite enormous source words are certified compositionally;
they are not described as literally executed histories.

These are internal full-criterion decisions. No official OpenWave status is
changed, and the package does not claim all31 rows or one interacting
continuum construction. `VALIDATION.json` gives the precise review evidence;
`PLANS.json` lists executable replay modes and outputs.
