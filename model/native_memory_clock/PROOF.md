# A finite autonomous rail for native event records

## Construction and ancestry

The construction has a finite event horizon `T >= 1`. Its moving source is
the twelve directed bits on the tetrahedron, with declared background
`eta=(0,0,0,0)` or `eta=(1,0,0,0)`.
The directed-bit update and the address/overwritten-bit record are those of
[`native_join.py`](https://github.com/muellerberndt/openwave/blob/0c738545b10b655c72cf2f877ec55b983592f535/experiments/2026-11-17_mass_free_link_search/native_join.py).
The tetrahedron and incidence convention are specified in
[`source.py`](https://github.com/muellerberndt/openwave/blob/0c738545b10b655c72cf2f877ec55b983592f535/experiments/2026-11-17_mass_free_link_search/source.py).
The native crossing function is
[`Law.crossing`](https://github.com/muellerberndt/openwave/blob/0c738545b10b655c72cf2f877ec55b983592f535/candidates/A-address-bit-v1/a_address_bit_v1.py).
These are source mappings for the finite event rule. The rail, reversible
extension, finite memory allocation and clock Hamiltonian are supplied
construction choices.

For an edge `e=(a,b)` with `a<b`, incidence column `B_e` has entries `+1`
at `a` and `-1` at `b`. A source word `x` has flux
`z_e=x_(2e+1)-x_(2e)` and load vector `q=Bz+eta`. Address `2e` has ordered
endpoints `(a,b)` and address `2e+1` has `(b,a)`. For an address `k` with
ordered endpoints `(a,b)`, define

```text
j_k(x) =  1  if x_k=0 and q_a>q_b,
         -1  if x_k=1 and q_b>q_a,
          0  otherwise.
```

The map `f_k` flips bit `k` precisely when `j_k` is nonzero. All other bits
are retained. Each address has probability `1/12`, including addresses
whose update is a no-operation. This fixes the source event and its native
classical interpretation without an energy or particle dictionary.

## A reversible event on every memory state

Write `x=(s,b)`, where `b=x_k` and `s` contains the other eleven bits, and
write `g_k(s,b)=(f_k(x))_k`. The memory bit is `r`. Define

```text
W_k : (s,b,r) -> (s, g_k(s,b) XOR r, b).
```

This is a permutation of the entire basis, including dirty memory. For
output `(s,c,d)`, its unique inverse is

```text
b=d,  r=c XOR g_k(s,d).
```

Consequently `W_k` extends to a unitary. It has the explicit reversible
circuit: swap `b` and `r`, then XOR `g_k(s,r)` into the source bit. The
Boolean function reads only the bits incident to the two selected endpoints
and the selected bit. Its bounded support on this fixed tetrahedron does
not require a global state lookup as a physical primitive. A gate-synthesis
cost for a particular elementary gate set is a separate resource question.

Let `A` be a twelve-state address register and set
`W=sum_k |k><k|_A tensor W_k`. On address state
`|u>=sum_k |k>/sqrt(12)` and blank memory, this gives

```text
W (|psi> |u> |0>) = (1/sqrt(12)) sum_(x,k) psi_x |f_k(x)> |k> |x_k>.
```

This is precisely the source-bound one-event isometry. The permutation
proof applies to arbitrary superpositions and to an arbitrary untouched
quantum reference. It does not assume that the input is a classical bit
distribution.

Fresh memory is a hypothesis of the native-event interpretation. At zero
background, source integer `1`, address `0` and memory `0` give source `0`
and memory `1`. Reusing that memory at the same address restores source
`1`; two native events with fresh memory give source `0`. The dirty-memory
operation is unitary and fails the claimed native repeated-event law.
Reusing one address register also changes the law: it gives identical
successive addresses, whereas two independent uniform registers give all
`12^2` ordered address pairs.

## Complete carrier and local Hamiltonian

Cells `0,...,T` lie on a declared rail. Cell `j` carries a twelve-bit source
buffer `S_j`; for `j<T` it also carries address `A_j` and memory bit `R_j`.
The clock carrier is exactly the one-excitation space with orthonormal
basis `|0>,...,|T>`. No zero-excitation or multiple-excitation clock states
are part of this declared carrier.

The shared data dimension and total dimension are

```text
d_data = 2^(12(T+1)) 24^T,
d_total = (T+1) d_data.
```

Storage consists of `12(T+1)` source bits, `T` memory bits, `T`
twelve-state address registers and a clock of dimension `T+1`. A unary
embedding of the clock uses `T+1` two-state sites restricted to one
excitation. These are finite resources growing linearly in the horizon
when counted as fixed-size local registers. No memory reset occurs during
free evolution.

`W_j` applies the controlled reversible event to `S_j,A_j,R_j`. Set

```text
U_j = SWAP(S_j,S_(j+1)) W_j,
g_j = sqrt((j+1)(T-j))/T,
H = sum_(j=0)^(T-1) g_j [ |j+1><j| tensor U_j
                         + |j><j+1| tensor U_j^* ].
```

Each summand acts on two adjacent rail cells. Its data support is two
twelve-bit buffers, one address register and one memory bit; its clock
support is the corresponding neighboring clock sites. Thus the interaction
support is bounded independently of `T`, with a fixed finite Hilbert space
per rail cell. The summands are supplied interactions of whole patches.
This statement does not identify the rail with the tetrahedron's original
seam geometry, establish an elementary two-qubit Hamiltonian compilation,
or derive a physical spatial metric. The coupling profile depends on the
declared horizon and is not translation invariant.

Initial preparation places an arbitrary source/reference state on `S_0`,
sets every other buffer and every `R_j` to zero, prepares each `A_j` in
`|u>`, and puts the clock at `0`. Uniform addresses, blank storage and the
localized clock are supplied preparation resources. The free law is the
single time-independent Hamiltonian `H`.

## Full-carrier spectral identity

Define `P_0=I`, `P_n=U_(n-1)...U_0` on the entire shared data carrier, and

```text
D = sum_(n=0)^T |n><n| tensor P_n,
A = sum_(j=0)^(T-1) g_j (|j+1><j| + |j><j+1|).
```

Every `U_j` is unitary even on dirty buffers and memory. Therefore `D` is
unitary on the full carrier, and direct multiplication gives

```text
H = D (A tensor I_data) D^*.
```

This identity is stronger than a compression to the prepared history.
There are no omitted data sectors in its spectral conclusion. With
spin `J=T/2`, the matrix `A` is `(2/T) J_x`. One direct verification uses
the permutation-symmetric subspace of `T` spins: the normalized state with
`n` excitations has matrix element `sqrt((n+1)(T-n))` under the sum of
single-spin raising operators. Thus the complete spectrum is

```text
lambda_k = -1 + 2k/T,  k=0,...,T,
```

and every listed eigenvalue has multiplicity `d_data`. All negative,
zero and positive modes, with their multiplicities, belong to the law.
In particular `||H||=1`. The spacing `2/T` is a clock-carrier diagnostic;
it has no identified particle-mass interpretation.

## Binary source and numerical error budget

The spectral and transfer identities above concern the ideal real
couplings. The executable stores positive binary floating-point
couplings `x_j`. Their exact rational values are obtained with
`Fraction.from_float`. For the transfer profile let
`r_j=(j+1)(T-j)/T^2`, so the ideal coupling is `sqrt(r_j)`. The exact
rational number

```text
e_j = |x_j^2-r_j|/x_j
```

is an upper bound on `|x_j-sqrt(r_j)|`, since
`|x-sqrt(r)|=|x^2-r|/(x+sqrt(r)) <= |x^2-r|/x` for `x>0`.
For the uniform control use `r_j=1/4`; its binary coupling `1/2` is exact.
Set the nonexistent endpoint errors to zero and define

```text
epsilon = max_(n=0,...,T) (e_(n-1)+e_n).
```

The absolute row sums of the symmetric clock perturbation are bounded
by `epsilon`, which bounds its operator norm. Both clocks use the same
exact data gates and the same unitary `D`. Thus, writing `H_binary` for
the Hamiltonian with the stored couplings,
`||H_binary-H|| <= epsilon` on the full carrier, including every dirty
data sector and every untouched reference. The executable records this
bound as an exact rational. Its binary clock retains full-carrier gauge
conjugacy. Every eigenvalue of its clock matrix has multiplicity
`d_data`, while its ordered eigenvalues differ from the ideal values by
at most `epsilon`. Positive neighboring couplings make that clock matrix
an irreducible Jacobi matrix, whose eigenvalues are simple.

At a fixed real time `tau`, Duhamel's formula gives

```text
||exp(-i tau H_binary)-exp(-i tau H)|| <= min(2, |tau| epsilon).
```

For the same normalized initial vector, this bounds state-vector
distance. For any complete measurement, the corresponding record-law
`L1` distance is at most twice that state distance, and is also at most
two. Purification and trace-norm contraction give the same record bound
for mixed inputs with arbitrary quantum references. These bounds are
unconditional record bounds. Renormalizing a rare clock outcome requires
a separate conditional-error bound.

The recorded runtime `tau` is a binary number interpreted as an exact
real argument in this comparison. Exact arrival at `pi T/2` and exact
return at `pi T` are statements about the ideal real law and those exact
times. Comparing a binary runtime `tau_binary` to an ideal target `tau_*`
adds at most `|tau_binary-tau_*| ||H||` to the state bound. An initial
preparation-vector error adds its norm as a further term.

The rational coupling bound does not certify roundoff in initial square
roots, `pi`, trigonometric evaluation, probability accumulation, or
SciPy's numerical exponential. The fixed numerical norm, complete-law
and spectral discrepancy gates are implementation checks; a small norm
residual alone cannot bound a coherent phase error. The exact theorems,
the binary-coupling bound and numerical replay tolerances have separate
roles and are not physical uncertainty estimates.

## Clock amplitudes and native event instruments

In units with `hbar=1`, let `tau` be the external argument of
`exp(-i tau H)`. For initial data vector `|Psi>`,

```text
exp(-i tau H) |0>|Psi> = sum_(n=0)^T a_n(tau) |n> P_n|Psi>,
a_n(tau) = sqrt(binomial(T,n)) (-i sin(tau/T))^n cos(tau/T)^(T-n).
```

The formula follows by applying the independent single-spin rotation
`exp(-i (tau/T) X)` to the symmetric spin realization of `A`. Hence

```text
p_n(tau) = binomial(T,n) sin(tau/T)^(2n) cos(tau/T)^(2(T-n)).
```

The clock law is independent of the source and of its quantum reference.
At `tau=pi T/2` the clock reaches `T` with probability one. At
`tau=pi T`, `exp(-i tau H)=(-1)^T I` on the complete carrier: the source,
addresses, memory and buffers return together. The same law that writes
the finite archive also uncomputes it. Couplings and observation arguments
are supplied dimensionless choices, with no selection from empirical
masses, couplings or durations.

For address `a` and old bit `b`, define the source operator

```text
F_(a,b) = sum_(x:x_a=b) |f_a(x)><x|.
```

The map `f_a` is injective on the set with fixed `x_a=b`, since its other
bits retain the input. Thus `F_(a,b)^* F_(a,b)` is the projector onto that
set, and `sum_(a,b) F_(a,b)^* F_(a,b)=12 I`. For history
`h=((a_0,b_0),...,(a_(n-1),b_(n-1)))`, put

```text
K_h = 12^(-n/2) F_(a_(n-1),b_(n-1)) ... F_(a_0,b_0).
```

The empty history has `K_empty=I`. Induction gives
`sum_h K_h^* K_h=I`. Conditional on clock position `n`, the source is in
`S_n`, earlier buffers are blank, and the first `n` record cells contain
the orthogonal history labels of this instrument. Later address registers
retain their supplied `|u>` state; they do not denote occurred events.
Later memories and buffers are blank. These unused quantum resources are
retained in the carrier.

For arbitrary density operator `rho_SR`, including source/reference
entanglement, the unnormalized state for terminal outcome `(n,h)` is

```text
p_n(tau) (K_h tensor I_R) rho_SR (K_h^* tensor I_R).
```

All positions, histories and no-operation records are part of the complete
instrument. A further source-basis readout multiplies its Kraus operator
by the corresponding basis projector. For diagonal inputs this is the
exact native `n`-event classical process. For coherent inputs it is the
specified coherent completion, whose phases are additional architectural
data. Neither diagonal agreement nor the one-event source theorem proves
that this completion is unique.

The endpoint word and retained history also determine the primary
candidate's complete `(seam,coin,current)` records without an initial-word
side channel. In reverse event order, read the successor's addressed bit
`c`, recover the directed current `j=c-b` from the retained old bit `b`,
and restore that bit to obtain the predecessor. Verify that its actual
native update gives the successor. The public seam is `floor(a/2)`, the
coin is `a mod 2`, and the seam-oriented current is `(-1)^(a mod 2) j`.
Induction reconstructs the initial word and every event. `decode_history`
implements this map and rejects incompatible records. This native-history
interpretation requires the blank-memory preparation; dirty complements
retain the total unitary law without acquiring that interpretation.

Clock position supplies an internal finite ordinal record. It does not
identify the external parameter `tau` with an observer's proper time.
Reading only the arrival outcome and renormalizing would omit the other
clock outcomes. Persistent intermediate observers require their apparatus
to be included in the dynamics; the terminal instrument alone provides no
such construction.

The numerical prepared span contains every address word and every clock
position for each supplied source-basis component. Invertibility makes
these basis states distinct, giving dimension
`(T+1) 12^T M` for `M` supplied source labels. Evaluating every literal
Hamiltonian column verifies that no complementary output escapes that
span. It is therefore a reducing span for that preparation. The spectral
theorem uses the full-carrier gauge identity, rather than treating this
prepared span as the complete data carrier.

## Measurement and energy

The autonomous energy is `H`, conserved during the free interval. Native
disagreement `Q(x)=sum_v q_v^2` obeys the distinct event identity

```text
Q(x)-Q(f_k(x)) = 2 j_k(x)(q_a-q_b)-2 j_k(x)^2 >= 0.
```

This is a directed event-order identity. It is not conservation of the
autonomous rail energy, nor a thermodynamic cost assigned to the record.
At general source time the source is correlated with multiple clock
positions; the law includes reverse propagation and record erasure.

Let `Phi` be nonselective dephasing in clock position, with identity on
the data. Every term of `H` is off-diagonal in that basis, so
`Phi(H)=0` and `Tr(H Phi(rho))=0` for every `rho`. An initial position
state has mean energy zero, and its autonomous evolution retains that
mean. Clock readout therefore has zero mean energy change on this
preparation. For a general coherent clock state the mean change is
`-Tr(H rho)` and need not vanish. Its absolute value is at most
`||H||=1`; the general difference-of-expectations bound is `2||H||=2`.
For the binary-coupling source, replace the norm-one bound by
`||H_binary|| <= 1+epsilon`.
Mean energy accounting alone supplies no apparatus Hamiltonian or
zero-cost measurement implementation.

The dephasing intervention does not commute with the free dynamics. In
particular, dephasing at `tau=pi T/2` would act trivially on the occupied
endpoint; a nontrivial coherence control must choose an interior clock
distribution. Dephasing at `tau=pi T/4` and continuing to `tau=pi T`
gives clock-return probability

```text
sum_(n=0)^T [binomial(T,n)/2^T]^2 = binomial(2T,T)/4^T < 1.
```

Without that intervention the return probability is one. This analytic
control concerns a declared measurement intervention, with all clock
outcomes retained; it is separate from uninterrupted autonomous evolution.

## Finite-resource obstructions

Perfect retention of all address strings for `n` events requires at least
`12^n` distinguishable record states, since every string has nonzero
probability. Thus a record register decoding all strings without error
needs dimension at least `12^n`. The explicit address/old-bit allocation
has dimension `24^n`; no minimality claim is made for that allocation.
This bound is about complete distinguishable histories, not a general
prohibition on indefinitely operating finite reversible systems.

A second obstruction concerns the native marginal without requiring its
complete history. Suppose one finite closed unitary `V` and a fixed finite
environment reproduced the native classical process at every integer
event number. Choose a source input with a strict expected drop of `Q`
after one event. Source integer `1` on the zero-background tetrahedron
has such a drop. The native rule then gives `E[Q_n] <= E[Q_1] < Q_0`
for all `n>=1`. Finite-dimensional unitary recurrence supplies a sequence
`n_l -> infinity` with `V^(n_l)` tending to the identity, and hence
`E[Q_(n_l)] -> Q_0`, a contradiction. Recurrence follows by simultaneous
approximation of the finitely many eigenphases on their compact torus.
The argument applies to fixed-time samples of a finite autonomous
Hamiltonian as well. Infinite resources, external reset, approximation
over a stated horizon, and changed marginal laws lie outside its premises.

If a finite observable's expectation is nondecreasing under a fixed
Hamiltonian for every initial state and all times, that observable is
conserved. The premise requires `i[H,N]` to be positive semidefinite.
Its trace is zero, so it must vanish. For a discrete unitary the same proof uses
`V^* N V-N`. Monotone mean progress for one prepared finite interval is
compatible with this obstruction; a permanent all-state clock is not.

## Sampling and complete record laws

Each sampled table describes independent repetitions of preparation,
free evolution and one terminal instrument. It is not a continuously
measured trajectory. The headline protocol has twelve tables, each with
`N=256` outcomes drawn from a `K`-label numerical law. For independent
ideal draws with exact categorical probabilities `q_i`, Hoeffding's
inequality and a union bound give the per-label radius

```text
r_K = sqrt(log(2 K 12 / 0.01)/(2 N)).
```

For each table, the probability that any label exceeds this radius is
at most `0.01/12`; the bound over all twelve tables is `0.01`. Different
tables may share random numbers. Independence between tables is not
required by the union bound. Independence within a table is a hypothesis
of the probability statement. A fixed-seed pseudorandom replay does not
prove that hypothesis or supply independent physical measurement data.

The implemented inverse-CDF sampler uses stored nonnegative binary
probabilities `p_i`, floating cumulative sums `C_i`, and a uniform grid
with spacing `2^-53`. Let `S_i=sum_(j<=i) p_j` denote exact rational sums,
`delta=max_i |C_i-S_i|`, and `eta=|S_K-1|`. Clipping interval endpoints
to the unit interval and counting uniform grid points gives the
conservative per-label discrepancy

```text
|q_i-p_i| <= beta = 2 delta + eta + 2^-52.
```

The executable converts each stored probability and actual floating
cumulative sum to its exact `Fraction` value, computes `delta` and `eta`
with rational arithmetic, and records the exact rational `beta`. Its
floating bias is rounded upward with `nextafter`. The frequency gate
uses `r_K+beta` when comparing with the stored table, separately from
error between that numerical table and the ideal quantum law.

If the floating CDF ends below one, the sampler has a missing-tail
category. Its probability is at most
`m=max(0,1-C_K)+2^-53`. The program rejects any sampled index outside the
table; it does not discard or resample that draw. The ideal-IID failure
probability for the combined frequency and no-overflow gates is bounded
by `0.01 + sum_tables N m`, using the recorded missing-tail bounds.
The executable rounds each floating tail bound and the combined bound
upward; the frequency radius retains its separate numerical evaluation.
An observed absence of overflow does not establish zero tail probability.
Extended checks have no shots: their exact categorical descriptions and
numerical complete-law comparisons carry no sampling confidence claim.

## Controls, source boundaries and possible falsification

The same reversible construction accepts a competing nonzero update law
that reverses the native inequalities and drives load disagreement uphill.
It uses the same rail, address preparation, storage and terminal apparatus.
Unitarity and the clock spectrum cannot distinguish those two laws; their
source/history statistics must do so. This is a control against treating
installed unitary structure as evidence for the native repair law.

A clock control retains the same `U_j` and apparatus with uniform nonzero
nearest-neighbor couplings. The gauge identity and conditional native
instrument survive, while the engineered binomial clock and perfect
endpoint-transfer formula are not transferred to that changed source.
Its full spectrum and clock distribution require their own calculation.
Dirty-memory reuse and interior clock dephasing supply the explicit
failure mechanisms above.

The phase control multiplies each forward gate by `(-1)^b`, where `b`
is its overwritten bit, and uses the adjoint on reverse propagation.
For each resolved history this changes `K_h` by the scalar
`(-1)^(sum_j b_j)`. Consequently every classical archive instrument map
is unchanged, including on source/reference inputs. With `T=1`, prepare
`(|0>+|1>)/sqrt(2)` and apply the same inverse native data gate to both
completions at the endpoint, followed by a complete `X` measurement on
source bit zero. The native completion gives minus probability zero.
Only address zero makes the phase control act as `Z` on that input, so
its minus probability is `1/12`. The address and clock complements are
retained by the instrument.

The classical archive readout and the inverse-decoder `X` readout are
alternative terminal instruments on separate preparations. The inverse
decoder uncomputes the overwritten-bit archive. Its `X` result therefore
does not coexist with a committed classical value of that overwritten
bit in this experiment. Copying or reading that archive first changes the
coherent instrument and invalidates the stated `X` probabilities. No
joint committed-old-bit/`X` interpretation follows from comparing the
two terminal settings.

The central claim is falsified by a failure of any total-basis inverse,
an incorrect source-bound event, a disagreement with the complete clock
formula, a reference-dependent clock probability, or a conditional record
instrument differing from the literal native product. Controls that fail
to distinguish the native and competing law on a declared informative
preparation invalidate that discrimination claim. Numerical verification
requires source-bound gates and complete outcome comparisons; numerical
residuals alone do not prove spectral multiplicities.

The finite Hermitian and complete-instrument hypotheses described in
[CCF-LT-2's finite bridges](../ccf_lt2/bridges/PROOF.md) are proved directly
for this carrier above. The matter algebra, calibrated Hamiltonian,
capacity response, field preparations and source-specific conclusions of
CCF-LT-2 are not used. No map is asserted from the rail to the independent
OPH ratios `1/6:1/12:1/12`, `sqrt(2)`, `2/9`, the `Z6` class map, `3/8`,
flux gaps, `6`, or `5-sqrt(5)`. The rail supplies no relativistic field
limit, particle mass, gravitational coupling, rest clock or repair of the
separated-pair energy anomaly. It validates no complete OpenWave physical
row. Empirical comparison values do not enter this construction.
