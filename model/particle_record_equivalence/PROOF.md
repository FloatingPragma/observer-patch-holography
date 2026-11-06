# Exact repair and the recorded spectrum

The native matching code supplies an exact recorded recovery identity. Its
family memory also has a smaller encoding with the same complete observations.
These are positive finite results for the declared source process. Physical
mass identification, laboratory clocks and capacity-independent physical
preparation are separate requirements.

## Recorded correction preserves the logical process

Let `E` encode a six-dimensional logical space in the 512-dimensional native
carrier. Each column is a uniform sum of six distinct occupation addresses.
For a fixed seam, placement coin and raw record `y`, the native branch maps
that column to `c_y` distinct addresses with amplitude `1/sqrt(12)`. The factor
12 includes the six code terms and the fair placement coin.

The executable checks both matching codes, both C/E laws and all thirty native
seams. Every one of the 1,440 raw branches has the same `c_y` for all six
logical columns, with disjoint image supports between columns. Thus

```
(K_y E)* K_y E = w_y I6,          w_y = c_y/12,
sum_y w_y = 1.
```

For `w_y>0`, put `V_y=K_y E/sqrt(w_y)`. This is an isometry. On its image the
ideal recorded recovery is `R_y=E V_y*`, hence

```
R_y K_y E = sqrt(w_y) E.
```

An isometry between these equal-dimensional subspaces extends to a unitary
on the complete carrier. The source recovery uses the corresponding branch
permutation, inverse uniform mixing and code encoder. Zero branches retain
their outcome slots. This exact algebra concerns the ideal logical operations;
the numerical source compiler and physical services retain their own errors.

For every logical density matrix, including entanglement with a reference,
the corrected complete channel is

```
rho -> sum_y w_y |y><y| tensor E rho E*.
```

The raw record is independent of the logical input. Inserting any specified
logical evolution between such corrected repairs therefore preserves that
evolution, conditioned on every raw record. The correction cannot change its
frequencies. A conditional logical feedback operation changes this statement
unless it is included in the evolution being compared.

More generally, `T E=E B` and `L E=I` imply `L T^n E=B^n`. Consequently

```
L (I-z T)^(-1) E = (I-z B)^(-1)
```

on the common resolvent domain. The logical poles agree exactly; inaccessible
carrier modes can be absent from the compressed response. A different state
or lossy detector can make a pole's residue vanish. A faithful re-encoding
cannot shift the logical poles, while adding memory, feedback or a different
time rule changes the dynamics and requires its own derivation. No unmodelled
mass correction follows merely from observers writing or comparing records.

## Complete quantum histories in four memory addresses

The native family instrument has 128 raw C/E labels and twelve old-memory
inputs. Its complete branch is

```
K_(r,m) = sqrt(w_r) |F(m,r)><m| tensor U_F(m,r)(h),
U_j(h) = exp(-ih X tensor P_j).
```

Both the old memory and the raw labels are retained as outcome/environment
indices. The selected reset is address 2. The stationary law is supported on
`S={2,3,9,10}`. Every nonzero raw branch maps `S` into `S`; the executable
checks the source branch table directly, before aggregating its weights.

Let `E_S` send four encoded memory basis vectors to these four addresses.
Define the smaller branch using exactly the same raw and original old-memory
labels, the corresponding encoded destination, the same weight and the same
quantum pulse. Labels with an old memory outside `S` have zero branch maps.
Then, for every retained label and every real `h`,

```
K_(r,m) (E_S tensor I) = (E_S tensor I) B_(r,m).
```

The identity holds on the entire encoded memory/work space and remains valid
with an arbitrary reference. Multiplying branch identities proves equality of
all finite joint record probabilities and terminal work/reference states.
No stationary-product assumption is needed: every initially encoded state is
admitted, including memory/work correlations. Reset2 and stationary memory
satisfy this preparation condition.

Initial memory address 0 is outside that direct encoding, but every one of
its nonzero first-event branches lands at 2 or 9. Retaining the actual first
instrument and its old0 record gives an exact four-address continuation. Its
first-event marginal is `(71/72,0,1/72,0)`; replacing this with stationary
memory or reset2 would change the process. The complete one-step startup
domain is `{0,1,2,3,9,10}`, with arbitrary quantum correlations retained under
the actual first instrument. The other initial addresses need their transient
process or a separate comparison.

The exact four-address transition in the ordered basis `(2,3,9,10)` is

```
[67/72  71/72  67/72  67/72]
[ 5/72      0      0      0]
[    0      0      0   5/72]
[    0   1/72   5/72      0].
```

This is a complete instrument equivalence for an encoding of the same law.
It does not identify arbitrary initial states in the other eight addresses,
another feedback rule, a physical clock or a larger physical universe. The
proof explicitly supplies the outcome-wise maps that determinism alone lacks.

## Pulse duty is an action choice

With stationary memory `pi`, the first-order average for an admitted pulse
policy `D` is

```
Z_D = sum_(r,m) pi_m w_r chi_D(r) P_F(m,r).
```

All records remain present when a pulse is disabled; its quantum operation
is identity. Source enumeration gives traces `1,p,q,pq` for the four policies
always, C, E and CE, with `p=1/12`, `q=1/6`. In particular,

```
Z_C = Y-(1-p)P2,                Tr(Z_C)=p.
```

The exact normalized quadratic traces are different:

```
Tr(Y^2)                    = 3004995/3333362,
Tr(Z_C^2)/Tr(Z_C)^2         = 1200289/1666681.
```

Therefore `Y` and `Z_C/p` are not isospectral and cannot be related by a
faithful coordinate change or a common clock scaling. The diagnostic ordered
ratios are `1:195.11703757:3579.82624571` for always-duty and
`1:53.58800864:271.40307014` for C duty. Dividing by the C event probability
changes a common scale and leaves those ratios unchanged. The E and CE
policies have exact rank two; floating residual eigenvalues cannot be used as
denominators for a three-mass ratio.

Selecting a policy changes the coherent operation applied in the source.
It is a declared dynamical choice, with a different generator, rather than a
correction of the same measured quantity. These stationary averages are
first-order generators under the stated averaging limit, not exact
finite-step particle poles.

## Finite preparation and the two-clock inverse

The complete fixed-basis observer readout uses eight burn-in transitions and
one final recorded transition. From reset2 the exact memory error is

```
||T^9 e2-pi||_1 = 9765625/2001949865847816192
                <= 2(1/12)^9 = 1/2579890176.
```

Each `P_j` is a positive rank-one projector of norm one. The same upper bound
therefore controls the difference of finite and stationary averaged kernels.
For either coherent step `h`, the specified contrast has the exact form
`A_h=Y_finite+(cos(h)-1)Z_finite`. The two-clock weights cancel the second term
exactly. The replay retains all 82,944 scalar pointer probabilities and
recovers `Y_finite` with floating operator discrepancy below `3e-15` in this
execution. That number is a numerical diagnostic, not a physical error bound.

The producer reuses the finite flavor package's rational interval LDL inertia
certificate to enclose all three eigenvalues of `Y`. Its 55-step bisections
use exact rational arithmetic and certified interval pivot signs. Weyl's
inequality and the rational mixing bound then give exact rational second/first
and third/first ratio intervals for the finite preparation. Their approximate
decimal displays are

```
[195.116749, 195.117326],
[3579.820990, 3579.831501].
```

The exact rational endpoints are retained in the receipt; its floating endpoint
copies are diagnostic approximations. The finite preparation and finite pulse
inverse in this declared protocol
cannot explain a percent-scale hierarchy discrepancy. Unknown physical
service errors, finite-shot statistics or a changed protocol are not assigned
zero. No empirical mass table enters these calculations.
