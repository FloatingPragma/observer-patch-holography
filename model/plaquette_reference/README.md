# Full plaquette source reference

This closes a finite preparation gap in the same-action regional modular
calculation: all 256 states of the four-site, eight-mode Weyl/contact
plaquette participate in an executed source preparation. The original action
has `c=1/8`, `U=1`, and the declared inverse temperature is `beta=1/32`.
No particle-number sector is truncated. The reference is a controlled
approximation to the Gibbs state, with an exact rational error bound.

Write `Hc=H-I`, `Q8=sum_{k=0}^8 (-Hc/64)^k/k!`. The exact onsite frame
already established in `field_symmetry` makes `Hc` real. Integer matrix
powers compute every coefficient of `Q8`, without an eigenbasis oracle.
For column `j`, set `p_j=||Q8[:,j]||^2/Tr(Q8^2)` and
`v_j=Q8[:,j]/||Q8[:,j]||`. The source executes:

1. Prepare the eight-bit probability amplitudes `sqrt(p_j)` with the complete
   native amplitude-tree schedule.
2. Read all eight bits through the existing complete two-dimensional-code
   dephasing export, retaining all 256 records and every abort outcome.
3. For each record, reversibly map `|j>` to blank using native `X` pulses
   and the compensating phase; source-prepare the real column `v_j`.
4. Execute the inverse onsite frame with native Pauli words.

The literal tree is derived directly from the exact integer polynomial
columns and squared weights. Common integer factors are cancelled before
angle evaluation in a private 192-bit `mpmath` context; Walsh sums of the
represented angles are exact rational sums before conversion to binary64.
Thus BLAS norm reductions and host `atan2` implementations cannot alter
the retained word hashes or the number of distinct pulse controls.
All zero-angle requests remain. This is a reproducible angle construction,
not a claim of correct rounding: the independent `certified_reference`
interval checks still bound every actual emitted angle and integrated
source control against the exact integer target.

The unconditioned field density is `Q8^2/Tr(Q8^2)` in the original Weyl
frame. Its extension retains a **classical column record**. It is not a
pure thermofield purification, and a selected record's conditional field
is not claimed to be Gibbs. The native word is run for every one of the
256 branches; failures remain separate outputs and are never normalized
away. The full source snapshots and lazy failure vectors are available
from `prepare()`. A trajectory uses at most 10,855 native pulses, eight
readouts, eight initialized field bits and eight retained record bits.
Abort labels additionally retain stage, input-record prefix, first failing
event and Kraus outcome; the receipt counts these separately, including
the surviving spectator quantum register.
Enumeration of all branches is a host audit, not a requirement that one
physical trajectory execute all conditional words.

The exact row-sum proof gives `||Hc||<=3`. With `q=3/64`,
`r=q^9/(9! (1-q))` bounds `||Q8-exp(-Hc/64)||`. Normalized vectorization
then gives a density trace error `delta<=4r/(1-q)`. The signed CAR partial
trace is CPTP, so the same error controls the regional state. A rational
common lower eigenvalue bound `m>1/20` gives
`||log rho_A(Q8)-log rho_A(Gibbs)||<=delta/m<3e-16`.
This proof uses the full original action and an explicit finite circuit.
Ideal analog angles are supplied source controls; floating source
execution has a separate `2e-11` complete-output admission and is not
presented as an interval proof of hardware precision.

The receipt compares the **actually executed, unnormalized success field**
with the old regional expansion, while retaining its abort complement.
For the normalized Gibbs modular generator, the regional expansion is
`P8+log Tr_A(E exp(-beta Hc)) I`; its exact scalar is retained in the
absolute operator comparison. The scalar drops only on trace-zero entropy
tangents. The rational expansion remainder and the preparation bridge are
reported separately. The original source action is also continued at
modular parameter `s=1/2` on every record, retaining every export failure.
Since `Q8` is a polynomial in the original `H`, its ideal field marginal
commutes with that action exactly.

A separate portable bound is conditional on the declared complete-output
admission `eps=2e-11`: set `d_exec=(2+eps)eps` and add
`(delta+d_exec)/(m-d_exec)` to the exact regional expansion remainder.
The resulting rational bound is below `2.71e-9`. This does not assert
hardware precision. The tighter bridge calculated from the measured
floating residual is labeled a numerical diagnostic, not used to choose
the declared tolerance.

The temperature, action, coherent load/export, initialization, complete
readout, conditional scheduler, memory and analog controls are declared
simulator/observer services. They are not omitted assumptions or claimed
emergent physical calibrations. This construction establishes no
thermodynamic or interacting continuum limit, physical clock, cap geometry,
Newton normalization, or Einstein dynamics.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/plaquette_reference/test_preparation.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.plaquette_reference --check model/plaquette_reference/results.json
```
