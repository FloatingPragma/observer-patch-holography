# A faithful reference that survives the original field and complete repair

The same literal source pulse, repeated once on each occupation mode, now
prepares a faithful nontracial **stationary field weight**. Its stationarity
does not require perfect control angles. The actual represented controls
determine two positive numbers `a,b`, and the unnormalized successful field
block is exactly

```
sigma_M = diag(a,b) tensor M = a^M exp(gamma N),  gamma=log(b/a)>0.
```

Here `N` is the original total CAR number. Both the unchanged interacting
Weyl/contact action and every unchanged complete charged-block repair
preserve it. All local export aborts and all later repair records remain in
the output. This is a number-biased reference, not a Gibbs state for the
interacting Hamiltonian. Its modular flow is trivial on fixed-number
observables, including the previous fixed-N particle sector.

## Actual source preparation and complete outputs

For each mode, initialize the existing six-level source in `|0>`, apply its
native `Y(-pi/3)` request, and use the existing complete dephasing export.
The processor is reset by that specified export service before reuse.
Interpret the twelve emitted binary64 integrated controls as exact rational
coefficients of the retained `Q(sqrt(5))` response matrices. If
`U=exp(sum c_p D_p)`, then `a=|U00|^2`, `b=|U10|^2`. The other four amplitudes
go to explicitly labeled abort records. No floating matrix exponential is
used to prove these definitions or their bounds.

Identical initialized inputs, pulse controls and export maps give identical
local instruments. Their product retains all `6^M` records. The proper field
block therefore has trace `(a+b)^M`, and total abort mass is
`1-(a+b)^M`. Neither quantity is silently replaced by one or zero. Product
record access in `preparation.numerical_branch` retains even zero-weight
slots. The total source cost is M pulses and M complete exports; the
initialized modes, receiver memory and reset/readout service are explicit
supplied apparatus. Preparing the field requires no eigensolver, imported
thermal state, outcome rejection or conditional preparation tree.

Independent outward dyadic arithmetic bounds the represented angle and all
twelve represented controls on the complete six-level processor. Duhamel's
inequality for anti-Hermitian generators and the exact complete export give

```
epsilon < 1.894e-16                 one-mode complete isometry error
|a-1/4|, |b-3/4| <= delta=(2+epsilon)epsilon
one-mode abort probability <= epsilon^2 < 3.587e-32
|gamma-log(3)| < 2.021e-15.
```

Thus `a,b>0` and `b>a` without trusting numerical cancellation or an
assumed execution tolerance. The ideal reference would have occupation
probability3/4. The actual exact mathematical source generally prepares a
nearby number-biased weight; **that nearby weight itself is stationary**.
The displayed floating execution is a diagnostic of the implementation,
not a guarantee of physical hardware accuracy or correctly rounded `expm`.

## Continuing process and scheduler premise

Every term `c_i^* c_j` and `n_i n_j` in the original action commutes with N.
The existing dimer ground has number2, the hole/particle doublets have
number1/3, and graded whole-block swaps preserve total number. Thus every
native projector, partial swap and outside projector commutes with N.
The integer native routing checks both completeness and unitality on every
existing seam. Adding `Q=I-P` gives unitality on the complete field algebra.
For each such instrument,

```
sum_r L_r sigma_M L_r^* = sigma_M sum_r L_r L_r^* = sigma_M.
```

This identity holds for arbitrary finite geometry of the original field,
and any placement of the existing finite native instruments on its modes;
it does not require the aligned spin frame or a larger native port graph.
Sequential reuse of the specified finite processor and tensor-identity
extension are the stated placement services. No arbitrary full-Fock matrix
is introduced as a primitive.

After preparing every mode, retain and bypass any register sector containing
an export-abort flag. Apply the field/repair schedule on the all-proper
sector. Both the proper block and untouched abort blocks are then invariant.
Any finite **record-independent** schedule of Hamiltonian durations and
these complete repairs preserves the field-plus-abort marginal. The same
claim covers an independent classical random choice of schedules.
This marginal traces all successful occupation-readout records as well as
the later repair records; those records remain retained in the full output.

It does not say that the growing classical-history state is stationary.
Individual records normally condition the field. For the ideal eight-mode
number reference, a two-dimer native hold has probability `9/32768`; its
normalized conditioned state is `P/4`, at trace distance `16375/16384` from
the unconditioned reference. That normalized state describes the record's
mathematical conditional state, not postselection in this protocol.

Record-dependent feedback also needs a separate theorem. An exact native
counterexample starts in the uniform three-position highest-spin corner,
applies E on seam(0,1), and applies E on seam(0,2) only after
`native(coin0,raw5)`. All33 resulting records, including outside/zero slots,
are retained. The final distribution is `(5/12,1/3,1/4)`, at trace distance
`1/12` from uniform. The composite is trace preserving but not unital,
although both available repair channels are unital. No feedback claim is
hidden in the finite-schedule result.

`Lean/QFT/ContinuingReference.lean` proves the generic complete-sum
stationarity identity, its power specialization and finite fixed-history
composition. It does not encode the concrete source matrices, positivity,
fixed-point classification or physical clock identification as Lean values.

## Local modular bounds without total-volume amplification

First trace the **complete** output outside a k-mode region, then compress
that region to its proper field code. Its weight is `diag(a,b) tensor k`,
independently of total M. Relative to `diag(1/4,3/4) tensor k`,

```
local complete-output trace error <= 2 k epsilon
minimum proper regional weight >= (1/4-delta)^k
||-log sigma_k + log rho_k,ideal|| <= k delta/(1/4-delta).
```

The last inequality is the scalar log Lipschitz bound applied separately
to each tensor factor. It does not divide a global error by a vanishing
global minimum eigenvalue. These bounds remain unchanged under the ideal
continuing process above, because its complete field-plus-abort marginal
is exactly fixed. Projecting **globally** onto success before tracing would
instead multiply the regional weight by `(a+b)^(M-k)` and add an explicit
scalar to its logarithm; that different operation is not used in the local
bound. No global success conditioning or missing scalar is involved.

The modular automorphism is conjugation by `exp(it gamma N)`, with the vacuum
phase chosen as zero. The absolute scalar `a^(itM)` cancels only in that
automorphism. A native Z/I pair per mode implements a controlled-precision
approximation, including all clock-export aborts and arbitrary input/reference
columns. Its bound includes the exact interval for gamma, the log3 series,
angle rounding and port rounding. This is an accessible charge clock;
coherent charge references are the specified source encoding, not a derived
physical superselection resource. Imperfect compiled clock/field words add
their own errors per use. Exact stationarity of the mathematical field and
repair laws does not erase those implementation errors.

## Which stationary clocks survive inside the particle sector

For a finite complete unital channel and Hermitian X,

```
sum_r ||[L_r,X]||_HS^2
 = 2 (||X||_HS^2 - Re Tr(X* Phi(X))).
```

Consequently every fixed density commutes with all Kraus operators.
Requiring stationarity for every Hamiltonian time also requires commutation
with H. The already source-bound joint algebra on six-site `N=4,S=1`
therefore classifies **all** common stationary references there:

```
I90 tensor A  direct_sum  B0 direct_sum B1 ... direct_sum B12,
A>0 in M3; B0>0 in M9; Bj>0 in M3;
90 Tr(A)+sum_j Tr(Bj)=1.
```

There are198 complex commutant dimensions. The active factor is scalar on
its90 energy multiplicities. The dark energy at K=-1 has spin-inclusive
dimension9; each other distinct dark energy has dimension3. Positive
spin weights or distinct dark energy weights give nontracial common
references. This classification is for every constituent map and every
Hamiltonian time, not fixed points of one arbitrarily chosen stroboscopic
composite. It is not a claim to have source-prepared every matrix in the
classification.

The new exact projector
`P=(4I-Splus* Splus-eta_minus* eta_minus)/4` on the120-dimensional highest
weight space gives `Tr(PK)=-135`, `Tr(PK^2)=1815/8`, hence active energy
variance `13/48`. Popoviciu's spectral variance bound gives distance at
least `|beta|sqrt(13/48)>|beta|/2` between any active stationary modular
generator and `beta K`, even allowing an arbitrary scalar and protected
spin term. Thus the original nonzero-beta Gibbs state is not a common
stationary reference. This is a quantified obstruction to identifying
these clocks, not to the positive continuing reference just constructed.

A separate exact target on aligned chains is
`3^(2S_z_hidden)/(16/3)^L`; it has local eigenvalues
`(3,9,1,3)/16` and nontrivial modular action within fixed N. Its invariance
follows from the previously proved common spin. Its alternating mode
preparation is not asserted to have the identical-pulse exact-stationarity
property, and unaligned three-dimensional swaps remain outside its domain.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.continuing_reference --check model/continuing_reference/results.json
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/continuing_reference
```
