# DYN02: complete instruments for the proposed common LOCAL12 walk

This package attaches position and relative-phase instruments to DYN01's one
candidate, `native-local12-current-v0-proposed`. Its sole evolution input is the
frozen `native_walk_bridge/common_walk.py` / `compile_walk.py` source word for
`W = S R`. It supplies no Hamiltonian, mass parameter, eigenbasis preparation,
time-to-second conversion or empirical calibration. The earlier 192-CAR bridge
is a separately frozen comparison and is not imported here. All cases below
are development evidence; no official row status is changed.

## Source and clock

The source is a padded 128-position by 16-event memory. There are 96 physical
vertices and 12 active event slots. Each of six matched directions has a hold
and a move slot. `R = 2|u12><u12| - I16` acts on every position row; `S` is the
owner's native-event/record-recoding/geometric-move/uncomputation word. The
original native labels remain recoverable by the owner's inverse recoder.
Repeated coherent steps do not generate independent fresh uniform native
attempts. Global attempts, local attempts, source pulses and coherent W steps
remain distinct clocks.

The common preparation has amplitude 1/12 on the first cell's 12 ports and
12 active event slots, and zero elsewhere. The owner also retains the law,
endpoint, record and scratch memory, with native reset, source export and
padding records. Our reduced 128 by 16 numerical representation is justified
by the owner's reversible clean-workspace contract, rather than by deleting
failed or dirty outputs.

## Position instruments and their disturbance

`native_phase_probe.actual_native_pointer` executes the six-pulse source
pointer followed by the literal native C or E instrument. It checks the
complete one-bit map on every matrix unit, retaining coin, current, read bit
and source-export complements. Seven such readouts give the instrument

`rho -> P_x rho P_x`, for every `0 <= x < 128`,

with identity on the event register, every retained current coordinate and
any external reference. This is a sharp intervention, not a terminal hidden
position lookup. Two readouts require 84 source pulses and 14 fresh ready
auxiliary carriers. They give the complete unnormalized joint outputs

`P_y W^n P_x W^t rho (W*)^t P_x (W*)^n P_y`.

`SAME_WALK_RECORDS.json` contains every 128 by 128 joint probability, including
zero and padding outcomes, for `t = 0, 1` and `n = 0, 1, 2, 4, 8`. Integer
numerator recurrences prove normalization and preserve the conditional
poststate needed for every later W. Measurement disturbance is explicitly
nonzero. No projection onto the smaller cyclic span of `(A, SA)` is used after
a readout: the full flat and dispersive sectors continue under the same W.

## Relative-phase instrument

Each primitive source pulse is controlled on an added retained probe. The
construction covers all six primitive kinds I, X, Y, Z, H and J, including
their scalar phases. Complete local source execution checks every input
column and every export complement. The exact word identity

`C(U_N ... U_1) = C(U_N) ... C(U_1)`

then controls the owner's unchanged W word. The probe is prepared in plus,
followed by controlled W^n, an X or Y analysis word, and the actual native
pointer. Each quadrature is a declared setting with probability 1/2; C and E
are separately configured native laws. No undocumented prior is applied to
the step count. The full native fine branches remain part of the output.

The resulting ideal reference contrasts are the real or imaginary part
of `C_n = <psi|W^n|psi>`. The exact values for the common preparation are:

| W steps n | C_n |
| --- | --- |
| 0 | 1 |
| 1 | 11/12 |
| 2 | 25/36 |
| 4 | 5/54 |
| 8 | -13087/26244 |

The imaginary contrasts vanish for this real source and real preparation.
X outcomes have probabilities `(1 +/- C_n)/2`; Y outcomes each have 1/2.
Dephasing the probe before analysis destroys the contrast. Replacing W by -W
changes the odd-step contrast. A separate actual scalar-I control detects
the same convention error: source-global phase becomes relative probe phase.

The spectral interpretation uses the proven compression
`T = A* S A = I - L/12` and `A* W^n A = T_n(T)` for the unmeasured common
preparation. Its principal quasiphases are
`+/- arccos(1 - lambda(L)/12)` modulo 2 pi per coherent W step. Eigenvectors
are used only for a diagnostic comparison; they are not execution inputs.
This interpretation is not reapplied to a projected cyclic sector after a
position measurement. It establishes neither a physical rest mass nor an
absolute energy or time unit.

## Retained current and external references

DYN03's proposed finite current attachment uses `G = delta(vertex0)` and
`F = P_x + C z` in `(Z6)^288`, with 193 independent cycle coordinates. Its
normalized isometry is

`J_G |x,e> = |x,e> 6^(-193/2) sum_z |P_x + C z>`.

Here x ranges over valid source vertices; invalid position outputs are still
retained by the full padded instrument. The detector transfer requires
DYN03's exact, separately owned identities
`W_F J_G = J_G W` and `P_x J_G = J_G P_x`, tensor identity on an arbitrary
external reference. Controlled W obeys the same identity with identity on
the probe. Consequently all present phase and position instruments transfer
to that coherent correlated preparation while retaining the full current
fiber. A position intervention can leave span(A,SA) while remaining in the
full J_G image; these are different subspaces. This package does not trace a
definite current tape and assert preservation of the bare coherent density.
The current attachment is conditional on DYN03's final frozen proof/receipt;
it is not certified by the bare-walk numerical comparison alone. Nonzero
cycle characters and definite cycle coordinates are different preparations.

## Complete failures and quantitative bounds

`INSTRUMENT_CONTRACT.json` composes the source-control words with outward
integer-interval bounds from the existing certified-reference machinery.
One owner W has 432,425 source pulses. Its explicitly controllable word has
6,948,708 source pulses. The additional controlled-word implementation error
is at most `9.375998013913197e-10` per W. This bound compares the complete
controlled source output with the controlled product of the owner's ideal
proper-code gates. It must be added to the owner's separate error relative
to the intended geometric W; it does not replace that error.

Added apparatus resets have the declared initial carrier in the receipt and
use the existing complete native reset. The two-position protocol uses 14
ready carriers, with timeout upper bound `7.374596908588228e-8`. The phase
protocol adds one blank probe and one ready pointer carrier. Fresh independent
apparatus is reset before common source preparation. Each first-reset failure
has its exact product of preceding success probabilities; every fine native
history and timeout remains. Proper-code load/export and its internal
receiver initialization remain admitted coherent services. These added
carrier resets do not derive those services.

The complete bound recipes are:

* Position: common preparation error + `(t+n)` common-W error + 14 pointer
  errors + added reset timeout.
* Phase: common preparation error + `n` times (common-W error + controlled
  wrapper error) + probe preparation/analysis/pointer error + added reset
  timeout.

All source first-export failures retain their source/environment/reference
outputs. None of these bounds is obtained by normalizing successful runs.
The common provider's current receipt deliberately does not infer an
accumulated physical whole-W bound from floating local residuals. Therefore
this package does not report a completed numerical end-to-end physical error
budget until that separate bound is supplied. The exact logical records,
native local instrument checks, wrapper certificate and error recipes are
the bounded present result.

## Reproduction and custody

From this directory run, in order:

```
python native_phase_probe.py
python same_walk_records.py
python instrument_contract.py
```

The three receipts pin their sources. `MANIFEST.json` additionally binds the
frozen owner program/receipt and the independent kinematic review. Exact
records are rational; the reported maximum comparison error against the
owner's full padded source callback is `2.9249267195438145e-14`. Local source
float checks remain diagnostics alongside, not substitutes for, the stated
complete-error certificates.
