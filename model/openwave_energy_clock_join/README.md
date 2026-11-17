# Static transfer energy and reflection clocks

One calibration of the compact signed-source mode fails to predict a second
full-graph reflection frequency under the declared linear energy conversion.
The same coordinate apparatus resolves the discrepancy, with preparation
failures retained. The [proof and apparatus contract](PROOF.md) specify the
two modes, exact domain and certified error margins.

This is the diagnostic `static-two-mode-arc-diagnostic-v1`, an ideal apparatus
extension of the public compact arc clock. Its rank-one preparation filters
and full coordinate detector are declared services. Physical energy, a native
duration per reflection word and transfer to the official OpenWave candidate
are outside its result. Review status is provisional: the second method is
an internal implementation crosscheck, with no independent reviewer.

From the repository root, using the installed NumPy and SymPy dependencies:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B model/openwave_energy_clock_join/check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B model/openwave_energy_clock_join/crosscheck.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B -m pytest -q -W error model/openwave_energy_clock_join
```

Both calculators read only the sealed graph. They import no producer,
published spectrum or success receipt. Their shared conceptual dependencies
are the native endpoint formula, the graph and the stated reflection
architecture. Their stdout is regenerated evidence; no stored result or
global manifest is required.
