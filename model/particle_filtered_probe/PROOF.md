# A complete finite apparatus for a time-filtered physical probe

## 1. Source and reused constructions

A member of the [primary contract](../particle_common_contract/PROOF.md)
fixes the finite Wilson Hamiltonian H, its original energy zero, the
source-Gauss polynomial density rho and its canonical conjugate reference.
The density commutes with H. The reference remains idle throughout this
experiment and is retained in the final output; its joint purification is
not assumed stationary. No other vacuum, group factor or Hamiltonian is
substituted.

The [common physical field algebra](../particle_primary_field_algebra/PROOF.md)
also supplies authenticated charge projections and normalized Hermitian
Wilson plaquette contractions at its declared centered positions. These
are even and commute with all source Gauss constraints. Their exact
descriptor, filter width, member and clock are checked by `verify_probe`;
an inactive early-level probe remains zero. The construction below covers
this entire bounded menu as well as the explicit parity example. A parity
variance floor is not transferred to another probe.

The [physical boundary apparatus](../particle_physical_boundary_state/PROOF.md)
supplies complete spectral projectors E_q for integer-normalized endpoint
hypercharge Q6 on the origin-owned first link. Choose

```
A = sum_(q in S_ell) (-1)^q E_q.
```

This is an even, gauge-invariant Hermitian unitary of norm one on the
entire finite source. Its action is identity on the other field factors.
A coherent evaluation of the existing complete charge menu followed by
a parity phase and uncomputation implements controlled A. Measurement
of the charge label before uncomputation would be a different instrument.

The construction reuses the complete linear-combination interference
argument in [the spin packet proof](../spin_packet_readout/README.md),
the all-reference instrument compiler in
[the finite instrument proof](../spin_operational_position/PROOF.md),
and `spectral_window_access.cayley_service_bound`. The inherited finite
LCU matrix control and compiler planner are executed in the receipt.
Their spin or dimer states are not used as the primary preparation.
The [positive-energy result](../particle_joint_reconstruction/POSITIVE_ENERGY.md)
gives the ideal triangle multiplier and spectral-tail theorem used below.

## 2. Positive triangular quadrature

Fix a source-time half-width T>0 and an even integer M. Divide [-T,T]
into M cells of length h=2T/M, with midpoints t_j. Since zero is a cell
boundary, the triangular density is linear on every cell. Its exact
cell integral is

```
f_T(t)=(1-|t|/T)_+/T,
w_j=integral_cell f_T(t)dt = (2/M)(1-|t_j|/T),
w_j>0,       sum_j w_j=1.
```

Define alpha_t(A)=exp(itH) A exp(-itH), the continuous filtered probe A_f
and its quadrature B_M by

```
A_f=integral f_T(t) alpha_t(A)dt,
B_M=sum_j w_j alpha_tj(A).
```

If [L,U] bounds the original H, its commutator satisfies
`||[H,A]||<=U-L` for this unit probe. Subtracting the midpoint only in
this norm estimate does not alter any executed clock. Therefore

```
||B_M-A_f|| <= delta_q = min(2,(U-L)T/M).               (1)
```

More generally a proved bound J on `||[H,A]||` replaces U-L. The bound
uses exact cell integrals, so there is no separate coefficient or tail
error. The finite source norm estimate can grow rapidly with refinement;
no uniform bounded implementation cost is inferred.

## 3. An exact complete finite instrument

Use a neutral M-state selector and let v_j=sqrt(w_j). A real Householder
R maps its blank e_0 to v. Every entry is algebraic. Set

```
SELECT=sum_j |j><j| tensor alpha_tj(A),
W=(R† tensor I) SELECT (R tensor I),
K_r=<r|W|0>.
```

All blocks in SELECT are unitaries. Thus W is unitary,
`sum_r K_r†K_r=I`, and K_0=B_M. The apparatus retains every r, including
zero-weight outcomes, with the source and its untouched reference.
It does not repeat until success or divide by the probability of r=0.
The original success probability is

```
p_0=Tr rho B_M†B_M.
```

The other Kraus maps depend on the specified Householder, so the complete
instrument is fixed as well as its success effect. A random classical
choice of t_j gives `sum_j w_j alpha_tj(A) rho alpha_tj(A)†`; it has no
coherent cross terms and does not implement this filter.

For a Hermitian contraction A, introduce a neutral two-state probe ancilla
and the exact algebraic unitary

```
U_A = [[A,sqrt(I-A²)],[sqrt(I-A²),-A]].
```

Functional calculus makes the off-diagonal products commute, so
`U_A† U_A=I`. The source clock on this carrier is `I_probe tensor H`.
Apply the same selector to the conjugates of U_A, then retain both the
selector r and probe-ancilla outcome p. There are 2M fine Kraus operators
`K_(r,p)=<r,p|W|0,0>` and their squared effects sum to identity. The
double-blank block is exactly B_M for the original A. In particular A=0
has zero blank success, with all weight retained in other outcomes. The
dilation and its source conjugations preserve the physical Gauss carrier;
the added ancilla has trivial gauge action.

## 4. Original Gram correlations from complete Ramsey records

To measure the ordered correlation at source time tau, coherently
choose between W and its source conjugate
`(I tensor exp(itau H))W(I tensor exp(-itau H))`. A neutral Ramsey
qubit is prepared in |+>, with analyzer phase phi=0 or -pi/2. Its
complete fine Kraus menu is

```
L_(r,b) = [K_r+(-1)^b exp(i phi) alpha_tau(K_r)]/2,
r=0,...,M-1,        b=0,1.                              (2)
```

Summing b cancels cross terms, and summing r gives identity. On the
stationary source the r=0 probabilities obey

```
P(0,b)=[p_0+(-1)^b Re(exp(i phi)F_B(tau))]/2,
F_B(tau)=Tr rho B_M† alpha_tau(B_M).
```

Give score (-1)^b to r=0 and score zero to all r!=0. Its unconditional
mean is Re F_B for phi=0 and Im F_B for phi=-pi/2. The sum over the two
r=0 records is p_0. All 2M records are retained separately in each phase
setting. These are distinct experiments on fresh copies of the same
preparation, not sequential measurements on an undisturbed source.
For a nonstationary input the effects remain complete, but this displayed
identification of their sum with the stationary p_0 requires its premise.

The canonical reference is never counter-evolved. The two source echoes
already implement the Heisenberg conjugations. Discarding the reference
would weaken the observation domain, and replacing an echo by a one-sided
unitary gives a different Gram response.

For the contraction dilation replace r in (2) by the complete pair (r,p).
There are 4M records per Ramsey phase setting; only the double-blank pair
has a nonzero signed score. The stationary correlation identity and all
norm bounds remain valid. Sequential composition retains every tuple of
fine labels, with Kraus operator `K_last ... K_first`; its blank block is
the corresponding ordered product. This supplies finite-word original
Gram readouts without conditioning any intermediate record.

## 5. Finite algebraic clocks and complete implementation errors

For a signed time t use the exact algebraic unitary

```
C_r(t)=[(I-itH/(2r))(I+itH/(2r))^(-1)]^r.
```

For R_H=max(|L|,|U|), the inherited Cayley bound gives

```
||C_r(t)-exp(-itH)|| <= d(t)
 <= min(2,|t|^3 R_H^3/(12r²)).                          (3)
```

`cayley_service_bound(|t|,|t|,R_H,r)` computes this bound. The denominator
is invertible for every finite self-adjoint H. The inverse at -t is the
adjoint, and each transform commutes with H and every Gauss constraint.
Every negative-time use is an explicitly controlled reversible operation.
It is not an assertion about autonomous backwards time. Scalar shifts
cancel from the ideal Heisenberg echo. A finite Cayley function does not
have that exact shift covariance, so the source energy zero is retained
in both the matrix and error calculation.

Replacing every node clock by C_r changes each SELECT block by at most
2d(T). The finite implemented success operator B_tilde therefore obeys

```
||B_tilde-A_f|| <= e_B=delta_q+2d(T),                   (4)
```

with the trivial cap two. Comparing the complete finite Ramsey instrument
with its exact-node counterpart gives half-trace error at most
`2d(T)+2d(tau)` on every input and every untouched reference. Add a
compiler half-diamond allowance sigma and preparation half-trace allowance
zeta by composition. This comparison retains all fine failure maps.
Quadrature error (1) compares success operators and their original
responses. It does not identify their different failure maps with a
continuum ancilla instrument.

Let p_tilde and F_tilde denote the measured original success weight and
signed Ramsey scores. The unitary and contraction norm bounds give

```
|p_tilde-F_f(0)| <= 2e_B+sigma+zeta,
|Re F_tilde-Re F_f(tau)|,
|Im F_tilde-Im F_f(tau)| <= 2e_B+2d(tau)+2(sigma+zeta).  (5)
```

Probability errors may be capped by one and signed-score errors by two.
Subtract both allowances when forming the lower contrast
`F_f(0)-Re F_f(tau)`. Finite-shot sampling uncertainty, if present, is
additional; the receipt executes exact finite probabilities and generates
no empirical counts.

The exact apparatus route declares controlled algebraic Cayley gates,
neutral Householder preparation and the coherent boundary parity phase,
or the authenticated Hermitian-contraction dilation, as
program services on the selected finite carrier. This extends the
apparatus menu without changing H, rho or the reference. It is consistent
with the primary's declared exact algebraic preparation service. The
rational and algebraic matrices make the complete finite operation
computable.

An alternative compilation route uses the inherited addressed coherent
load/export service and finite register encoding. That service must be
admitted for the chosen source code; the spin-code proof alone does not
supply such a physical embedding. Its universal instrument formula is
`20D²+4D+10q` pulses for q source bits and padded source-plus-record
dimension D. `factored_compiler_plan` records this formula without
allocating D, and specializes its angle and port precision estimates.
Unused code inputs act by a declared unitary extension. Approximate
compilation can leave the physical code; its original code-complement
and export-failure outputs are retained and covered by sigma. Gatewise
Gauss preservation and a literal enormous compiled spatial execution are
not claimed. The exact primitive route preserves the physical carrier.

`precision_parameters` chooses finite M and r for any positive requested
response tolerance, using only source bounds and times. It gives
`delta_q<=epsilon/16`, `d(T),d(tau)<=epsilon/64` and
`sigma=epsilon/32`, so the total contrast error is at most epsilon/2
when zeta=0. Counts can be enormous. This is a complete finite
construction, without an efficiency or full spatial execution claim.

## 6. Complete Lüders histories for the common algebra

The coherent filter's blank effect is `B_tilde† B_tilde`. It is not the
effect used by the common algebra's binary Lüders measurements. For its
charge menu set `E=A_f`; for its Hermitian plaquette menu set
`E=(I+A_f)/2`. Positivity of the triangular density gives `0<=E<=I`.
The finite Cayley quadrature gives respectively

```
E_M=B_tilde,           E_M=(I+B_tilde)/2.
```

Each node conjugation is an exact unitary and every weight is positive,
so `0<=E_M<=I` holds before any approximation comparison. Let
`delta=min(1,e_B)` for a charge projection and `delta=min(1,e_B/2)`
for a plaquette. Then `||E_M-E||<=delta`. Define the distinct complete
instrument

```
J_0=sqrt(E_M),         J_1=sqrt(I-E_M),
V_M=[[J_0,-J_1],[J_1,J_0]].
```

The two roots commute and `J_0²+J_1²=I`, hence V_M is unitary. Both
binary records and the quantum reference are retained. These positive
algebraic roots and their neutral dilation are declared finite apparatus
services; the alternative inherited arbitrary-instrument compiler has
two fine outcomes and the same explicit source-code interface premise
as Section 5. Functional calculus preserves every gauge commutant.

For positive operators the square-root norm bound is
`||sqrt(X)-sqrt(Y)||<=sqrt(||X-Y||)`. One direct justification is
`X<=Y+delta I`, operator monotonicity of the square root, and
`sqrt(Y+delta I)<=sqrt(Y)+sqrt(delta)I`; reversing X,Y bounds the
self-adjoint root difference on both sides. Thus the two-outcome
Stinespring columns obey

```
||[sqrt(E_M);sqrt(I-E_M)]-[sqrt(E);sqrt(I-E)]||
 <= sqrt(2 delta).
```

Pure-state trace distance is at most the distance between their chosen
unit vectors. Purifying an arbitrary input and then writing the record
therefore gives complete half-trace distance at most
`min(1,sqrt(2 delta))`, uniformly over all untouched quantum references.
This compares the actual complete Lüders maps, not just their effects.

For a history of N probes with a common upper bound delta, compiler
allowance sigma per probe, initial preparation allowance zeta and total
interleaved source-clock operator error d_history, channel telescoping
and trace-norm contraction give

```
complete history half-trace error
 <= min(1,N*(sqrt(2 delta)+sigma)+zeta+d_history).       (6)
```

The clock allowance is a sum over the actual interleaved evolutions;
zero means their exact selected source clock. Conditional controls on
retained earlier records satisfy the same bound with their largest
branch allowance. Every binary history, including zero-probability
records, remains in the comparison domain. `lueders_error_budget`
returns a rational upper bound for (6) without a fixed absolute
square-root precision floor.

At each finite source and fixed finite history the schedule in Section 5
can make delta and sigma arbitrarily small. For example choosing its
response target at most `eta²/(16N²)` makes the instrument part of (6)
less than eta for `0<eta<1`; preparation and interleaved clocks require
their own allowances. This closes the finite apparatus premise for the
common algebra's complete effect histories. It does not assert uniform
cost in source size, a nonzero limiting response or particle content.

## 7. Spectral weights, finite survival and continuum boundary

For the original stationary source measure mu_A,

```
dmu_(A_f)(E)=sinc(ET/2)^4 dmu_A(E),
F_f(0)=mu_(A_f)(R),
F_f(0)-Re F_f(tau)=integral [1-cos(tau E)] dmu_(A_f)(E).
```

The multiplier is a probe operation. No weight is replaced by a
conditional success probability. For finite Cayley clocks the response
is compared by (4)-(5); their rational phases are not silently identified
with a different exact Hamiltonian clock.

The actual boundary-state proof supplies both even and odd endpoint
charge sectors with probability at least c_n>0. Consequently the
unfiltered parity probe has `Var_rho(A)>=4c_n²`. Stationarity preserves
the mean under time filtering, and

```
||A_f-A|| <= (U-L)T/3,
sqrt(Var_rho(A_f)) >= max(0,sqrt(Var_rho(A))-||A_f-A||).
```

This gives a finite survival condition for a sufficiently short filter.
The certified c_n can tend to zero, and connected variance can include
zero-frequency conserved fluctuations. It supplies neither a uniform
positive-energy window nor continuum noncollapse. A lower original
Ramsey contrast with controlled first moment and negative tail can test
such a window separately. The primary contrast is not computed here.

The [positive-window bound](../particle_positive_window/PROOF.md) uses the
same original measure: for a certified contrast D, mass bound M_0,
first absolute moment M_1 and negative tail eta_delta,

```
mu_(A_f)([delta,R])
 >= max(0,(D-M_0*(delta*tau)²/2-2*eta_delta-2*M_1/R)/2).
```

The measured contrast is reduced by both allowances in (5) before this
bound is applied. A positive value excludes an explanation entirely by
zero energy or weight above R. It is not an isolated pole or continuum
mass assertion, and no primary contrast is supplied by an apparatus plan.

## 8. Executed controls

The independent three-level source has eigenvalues 0,1,3 in a rotated
coordinate basis and stationary mean ground excess 1/1024. Its probe
swaps the first two energy directions and acts by -1 on the third.
The code executes all four selector outcomes and both Ramsey bits for
each phase, with T=1/16, tau=1 and sixteen Cayley subdivisions. It checks
normalization, the original success weight, both original Gram scores
and a strictly positive contrast allowance against the continuous
triangle. The reference dimension is three. This source is separate
from the spatial Wilson family.

For this control, the triangle moment bound is 32 and the original
negative tail below -1/4 is at most 1/256. With delta=1/4 and R=1024,
the exact recorded contrast and all errors above certify original
positive-band weight greater than 1/16. An independent spectral sum
checks it against the actual continuous filtered measure. This is a
finite integration check on the existing three-level control.

A separate exact two-level control uses H=diag(0,4), a rank-one
projection, T=1/2, two nodes and one Cayley subdivision. Its filtered
effect has diagonal 1/2 and off-diagonal 3/10. The contraction dilation
retains four outputs; the Lüders instrument retains two. Their original
blank probabilities on the stated diagonal state are respectively
17/50 and 1/2. The code executes both instruments and the four records
of a two-measurement Lüders history. Entangled-reference and perturbed
rank-deficient-effect controls check completeness and (6).

Tests independently compare the whole selector/Ramsey circuit, retained
reference marginal, exact spectral projector sums and numerical time
integrals. Countercontrols retain the error from classical selector
mixing, one-sided clocks, success normalization and a changed primary
member. No measured mass or learned parameter enters either construction.
