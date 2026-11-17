# Native charge readout with record integrity

A source emits native current records, a thirty-event contrast sweep checks
the reconstructed charge, and a second native carrier reads its finite
counter through a winding. The pointer reports the charge after the sweep
and retains timeouts. The result is conditional on a known preparation,
conserved load and faithful closing records, for one or two loads. The
charge unit and observer apparatus are declared construction choices;
their physical identification requires a separate argument.

Run from any directory with Python 3 and its standard library:

```sh
python3 /path/to/openwave_row1_integrity/run.py
```

The command writes `RESULTS.json` beside the code. It verifies frozen bytes,
actual native pointer tapes, complete bounded probability laws, failed
tape controls and adversarial-channel counterexamples. The executable
uses no workspace imports or external packages. [PROOF.md](PROOF.md)
states the exact theorem, clock accounting and boundaries.
