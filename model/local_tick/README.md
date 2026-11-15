# N4 memory survives the observer's local tick clock

This construction derives the occupation-memory test when the observer counts
its own incident attempts, while remote repairs continue.
**The N4 distinction remains strictly positive at the two certified couplings.**
Neither a preferred global physical clock nor a full-network matrix is needed
for this calculation. The event schedule, preparation and readout are supplied.

This is an extension of the [existing N4 finite result](../../code/native_response/README.md)
and uses the [agreement research foundations](../../docs/AGREEMENT_RESEARCH_FOUNDATIONS.md)
and [OPH interpretation](../../docs/INTERPRETATION_AND_TERMS.md). The mathematical
tools are the first-hit geometric series and Schur complement; the contribution
is their mapping to the declared source and a certified observable separation.
The candidate law and its state factors are unchanged.

## Result and operational contract

Keep the binary total-two shell, all twelve six-level private factors, full
normalized trace `tau`, and IID uniform choices among the thirty seams in the
source chart. A tick is **any attempted incident seam**, including an event
that leaves the observer's load unchanged. Its probability per scheduler
attempt is `p=1/6`. The waiting count is unbounded but almost surely finite.
This is an observer-local event count, not calibrated proper duration.

For each independent trial, measure only `N_v`, allow the unmodified source to
run until one or two local ticks have occurred, then measure only `N_v` again.
No intermediate measurement or private reset occurs. Both endpoint outcomes
are retained. The initial local Lüders measurement preserves the tracial
reference. Its joint law therefore supplies

```text
f = N_v - 1/6,     ||f||_tau² = 5/36,
c_k = E[(X_initial-1/6)(X_after k ticks-1/6)] / (5/36).
```

Let `r=1-c_1`. The corrected statistic is

```text
b_tick = c_2 - c_1²
         - r²(209r+180)/(1331r²+3300r+1800).
```

The native and one-term exchange controls give exactly

```text
c_1 = 5/11,       c_2 = 93/407,       b_tick = 0.
```

For both N4 completions C and E, the existing interval Gram margin gives:

| Coupling | Certified lower bound on `b_tick` | Decimal rendering of that lower bound |
| --- | --- | --- |
| `1/16` | `139996257/12500000000000` | `0.00001119970056` |
| `1/8` | `11335290321/50000000000000` | `0.00022670580642` |

These are inherited rigorous lower bounds, not estimates of the actual N4
statistic. The full N4 `c_2` is not numerically evaluated. C and E share `c_1`
and the lower-bound proof, but their remote resolvents may differ; equality
of their one/two-global-attempt correlations does not imply equality here.
Their local-tick values satisfy `c_2,C >= c_2,E` and `b_tick,C >= b_tick,E`;
strict inequality is not proved.

The closed comparator below permits the N4 compression factor `alpha` to be
matched using `c_1`. Thus positivity excludes even that refitted closed native
mean description under this same stopping instrument. It does not uniquely
identify N4, quantumness, a learning law or the physical Universe.

## First-hit channel and its ordering

Work in normalized Hilbert–Schmidt operator space. Each seam channel `Phi_e`
is self-adjoint and positive semidefinite as an operator on this space. For
native C its spectrum is contained in `{0,1/2,1}`; for the conditional
expectation E it is in `{0,1}`. The N4 unitary conjugations preserve these
facts. This is Hilbert-space operator positivity, in addition to complete
positivity of the channel; the latter alone would not justify this step.

Split the average into remote and incident parts:

```text
A = sum(remote e) Phi_e/30,       B = sum(incident e) Phi_e/30,
Phi = A+B,                       D = I-A.
```

Then `0<=A<=(1-p)I` and `pI<=D<=I`. A word of n remote events followed by one
incident event contributes `B A^n` in Schrödinger convention. Thus

```text
K = B D^(-1),                    K* = D^(-1) B.
```

The series converges in operator norm. Each stopped map is trace preserving
and unital, since every event probability is state independent and its total
success probability is one. The tracial state remains stationary at ticks.
The endpoint correlation is `c_k=<f,(K*)^k f>/||f||²`. `K*` is generally not
self-adjoint; treating it as such would give the wrong memory formula.

Every remote seam fixes `f`, so `Af=(1-p)f` and `D^(-1)f=f/p`. Put
`d=(I-Phi)f`, hence `Bf=pf-d`, and let `r_global=<f,d>/||f||²=alpha/11`.
Using the adjoint on the left factor gives exactly

```text
c_1 = 1 - r_global/p = 1 - 6 alpha/11,
c_2 = 1 - 2 r_global/p
      + <d,D^(-1)d>/(p ||f||²).
```

Remote events vanish from the first scalar moment by locality and stationarity.
They remain inside the second moment's resolvent. They cannot be deleted,
replaced by a mean waiting count of six, or treated as a deterministic time
rescaling.

## Closed comparator and the exact correction

Here `alpha` is the occupation Dirichlet compression factor and `P` below
is an operator-space projection. They are distinct from the electromagnetic
coupling `alpha_em` and geometric pixel parameter used by the
[binary detector](../observer/p_response/README.md).

Let `P` project onto centered linear occupations, and let `L_R` be the carrier
Laplacian containing only remote seams. The proved seam compression
gives, in coefficient coordinates on the centered subspace,

```text
D_00 = P D P = p I + alpha L_R/60,
d_0 = P d = (alpha/60) L e_v.
```

For the closed comparator, every seam acts on means by
`M_e=I-(alpha/2) w_e w_e^T`, where `w_e` is its oriented incidence vector.
For `0<=alpha<=1` this is a valid lazy-swap mean law, with the same constant
occupation sector and stationary covariance. Its local-tick map follows from
the same first-hit formula. This comparator is specified separately from N4;
compression alone is not claimed to make N4 means close.

The distance shells have sizes `1,5,5,1`. On shell-constant coefficient
vectors, the remote Laplacian is

```text
L_R,shell = [[0, 0,  0,  0],
             [0, 2, -2,  0],
             [0,-2,  3, -1],
             [0, 0, -5,  5]].
```

The last two unequal off-diagonal entries reflect shell multiplicities;
this displayed coordinate matrix is self-adjoint in the shell-size weighted
inner product. Direct inversion gives

```text
(L e_v)^T (10I+alpha L_R)^(-1) L e_v
 = 5(2alpha+3)(3alpha+10)/(11alpha²+50alpha+50).

c_2,closed - c_1²
 = 6alpha²(19alpha+30)/[121(11alpha²+50alpha+50)].
```

Substituting `alpha=11(1-c_1)/6` produces the correction in `b_tick`. Exact
12-port calculations agree for all observer ports at four rational control
couplings, and a symbolic identity covers the general parameter. An independent
floating evaluation of the complete 66-state native chain agrees with the exact
values `5/11` and `93/407` to the declared numerical tolerance.

The global-attempt correction `(1-c_1)²/10` does not hold at local ticks.
It assigns `-16/2035` to the native null; the test suite retains that negative
control. It is a null-model error, not a positive false detection of memory.

## Why the N4 excess stays strictly positive

Decompose `D` on `P` and `Q=I-P`, with blocks `D_ij`. The Schur complement
`S=D_11-D_10 D_00^(-1) D_01` is positive and obeys `S<=I`. With

```text
g = D_00^(-1)d_0,
z = Qd - D_10 g,
```

the standard block inverse identity gives

```text
b_tick = <z,S^(-1)z>/(p ||f||²) >= ||z||²/(p ||f||²) >= 0.
```

This also explains why the scalar correction is a valid null despite the
stopped channel's lack of self-adjointness. The inequality `S<=I` follows
from `D<=I` and subtraction of a positive Schur term.

The existing N4 construction defines odd-block edge residuals
`S_e=Phi_e(Z_e)-(1-alpha)Z_e` (extended by zero outside odd occupation), with
full-reference Gram matrix `G_ef=<S_e,S_f>`. The interval certificate proves
`G>=m I` at `theta=1/16` and `1/8`. To avoid confusing the two objects, this
packet calls the edge Gram `G`; the global-attempt packet calls it `B`.

For an occupation coefficient vector `u`, the hidden part of the all-seam
defect is `-(1/60) sum_e (partial u)_e S_e`. The hidden part of the remote
defect has the same expression restricted to remote seams. Consequently

```text
z = -(1/60) sum_e t_e S_e,
t_e = (partial e_v)_e                 for incident e,
t_e = -(partial g)_e                  for remote e.
```

Here `g` in the coefficient formula denotes the centered coefficient vector
representing the operator above. The five incident coefficients remain `+1`
or `-1`, irrespective of the remote resolvent, so `sum_e t_e²>=5`. Therefore

```text
||z||² >= 5m/3600,
b_tick >= (5m/3600)/[(1/6)(5/36)] = 3m/50.
```

This is six times the global-attempt lower bound `m/100`. It compares
two differently normalized event-count statistics; it is not a claim of six
times greater physical power or a sample-size improvement. The proof uses
the full Gram bound, not a reset, small-coupling truncation or mean-field
approximation. Both certified positive rational margins are recomputed from
the unchanged interval producer and propagated exactly.

There is a further non-strict completion comparison. On each native seam,
`C-E` is positive semidefinite: it is one half of the projection onto C's
`1/2` eigenspace. Unitary conjugation preserves this order, hence
`A_C>=A_E`, `D_C<=D_E` and `D_C^(-1)>=D_E^(-1)`. The vector
`d=(I-Phi)f` and `c_1` are identical for C and E. Their exact second-moment
formula therefore implies `c_2,C>=c_2,E`, and subtracting the same scalar null
gives `b_tick,C>=b_tick,E`. This ordering does not establish strictness or
numerically determine either moment.

For native and one-term exchange controls the occupations are invariant under
each seam's mean map, so all hidden residuals vanish and the exact null follows.
The one-term controls have normal exchange amplitude and `alpha=1` in the
existing source proof. The present packet does not transfer any higher-moment
closure from that fact to N4.

## What is supplied and what is open

The exact local-tick mathematical interface is constructed. Detecting a
tick, preparing the full trace, counting to one or two and retaining endpoint
records require an internal apparatus if the model is to realize its
own observer. The branch's dressed current label is not relabeled a bare
occupation measurement. The stopping proof supplies neither a calibrated
physical duration nor an autonomous schedule.

This experiment integrates an unbounded number of remote events. A finite
global-attempt cap H has missed-tick probability `(5/6)^H` for one tick;
two-tick completion needs its own retained failure accounting. Such a cap is
an auxiliary scheduler resource bound, not an observer-local timeout. No
finite-horizon postselection is silently substituted for the stated protocol.

A measured value needs a prospective error/confidence rule for the entire
nonlinear statistic, including correlated estimates, preparation/readout
errors and scheduler deviations. No sampling run or finite-sample detection
claim is made. Another enlarged hidden-state law may also give a positive
value. These finite source results do not select the physical law or complete
OPH's algebra-to-observer construction.

## Reproduction and independent controls

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/local_tick
python3 -m model.local_tick --check model/local_tick/results.json
```

The calculation requires the existing NumPy, SciPy, SymPy, mpmath and pytest
environment. `--write NEW_PATH` uses exclusive creation; source byte hashes
are checked before any retained calculation. `results.json` is deterministic
and contains exact comparator values plus the recomputed interval enclosures.

Tests compare the 12-port reduction against the full 66-state native stopping
kernel, preserve finite-prefix failure mass, and reject reversed factor order
using arbitrary preparations. A separate three-port one-particle quantum
model with a private qubit at the middle port constructs full Kraus maps and
checks the resolvent, non-self-adjoint stopped channel and Schur identity in
its entire operator space. That dimension-six model is an independent check
of the argument, not substituted N4 evidence. The N4 positivity comes from
the source-bound interval certificate and the analytic transfer above.
