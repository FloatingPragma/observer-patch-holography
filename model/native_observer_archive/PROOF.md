# Native records, local observers and a finite readiness interval

## Declared source and carrier

The native event is the twelve-directed-bit tetrahedron rule in
[`native_memory_clock/native.py`](../native_memory_clock/native.py), on either
declared background `(0,0,0,0)` or `(1,0,0,0)`. Its primary source mapping and
total reversible extension are proved in
[`native_memory_clock/PROOF.md`](../native_memory_clock/PROOF.md).
For an address `k`, write `f_k(x)` for the native successor and `b=x_k`.
The event gate acts on a source word and a memory bit by

```text
W_k(s,b,r) = (s, (f_k(s,b))_k XOR r, b).
```

For output `(s,c,d)`, the inverse recovers `b=d` and
`r=c XOR (f_k(s,d))_k`. This is a permutation on the full dirty-memory
domain. Blank memory supplies its native-event interpretation. Source
events retain no-operations and use the strict native load inequalities.

The clock has `N` sites `0,...,N-1`, restricted to exactly one excitation.
Four consecutive edges starting at `p` carry the data program below. Other
edges act as identity on the data. Five adjacent program cells carry five
twelve-bit source buffers, four twelve-state first-address registers, one
twelve-state action-address register, three two-bit couriers, two two-bit
observer pointers, and two native old-bit memories. A pointer or courier
has a value bit and a validity bit. The complete data and total dimensions
are

```text
D = 2^(5*12 + 3*2 + 2*2 + 2) 12^5 = 2^72 12^5,
dim(H_total) = N D.
```

Padding sites carry only the clock. A unary embedding uses `N` two-state
clock sites restricted to one excitation. All resources are finite; the
clock allocation grows with the declared rail length. No memory reset,
outside address draw, or measurement intervention occurs during free
evolution. Preparation of the packet, blank storage and address state is
a supplied resource.

## A total local observer program

Let `A_0,...,A_3` be the first-address registers, `B` the action address,
`R_1,R_2` the old-bit memories, `C_1,C_2,C_3` the couriers, and `P_A,P_B`
the observer pointers. Pairwise XOR acts on both bits of a pointer or
courier. The four forward data gates, in order, are:

1. Apply `W_(A_0)` to `S_0,R_1`; swap `S_0,S_1` and `A_0,A_1`.
2. XOR the bit `(S_1)_(A_1)` into `P_A.value`, toggle `P_A.valid`, and
   XOR `P_A` into `C_1`. Swap `S_1,S_2`, `A_1,A_2`, and `C_1,C_2`.
3. Swap `S_2,S_3`, `A_2,A_3`, and `C_2,C_3`.
4. XOR `C_3` into `P_B`. Add `6 P_B.value P_B.valid` to `B` modulo twelve,
   apply `W_B` to `S_3,R_2`, and swap `S_3,S_4`.

Each instruction is a permutation on every data basis state. In particular,
pointer publication is XOR into the existing pointer, rather than an
overwrite, and the controlled address change is a modular translation.
Reversing the instruction order and inverting each instruction therefore
defines every adjoint on dirty buffers, couriers, pointers and addresses.
The copy operations copy orthogonal computational records. They do not
clone arbitrary quantum states.

Gate `j` acts on program cells `j,j+1`. A whole tetrahedral buffer is a
bounded-size local register, so each gate and its associated clock hop
have bounded support independent of `N`. These are supplied patch
interactions; no elementary two-qubit Hamiltonian compilation, original
tetrahedron seam geometry, or physical spatial metric is inferred.

The preparation places an arbitrary source/reference density on `S_0`,
sets all other data registers to zero, and prepares
`A_0 = sum_(a=0)^11 |a>/sqrt(12)`. In particular `B=0`; the second address
is a declared adaptive intervention, not a second independent uniform
draw. If the initial word is `x` and the first address is `a`, put

```text
y = f_a(x),       b = x_a,       r = y_a,
d = 6r,          z = f_d(y),    b' = y_d.
```

After four gates, `S_4=z`, `A_3=a`, `B=d`, `R_1=b`, `R_2=b'`, and
`P_A=P_B=C_3=(r,1)`. The other transported buffers, first-address registers
and couriers are blank. Both observers thus hold the same actual updated
bit, and the second native event uses the address selected by that record.
This is an adaptive two-event native instrument with a supplied observation
and feedback rule.

## Complete Hamiltonian and spectrum

Let `U_n` be the relevant program gate on edges `p,...,p+3` and identity
on every other edge. Define

```text
H = (1/2) sum_(n=0)^(N-2)
    [ |n+1><n| tensor U_n + |n><n+1| tensor U_n^* ],
P_0 = I,    P_n = U_(n-1) ... U_0,
G = sum_(n=0)^(N-1) |n><n| tensor P_n.
```

All `U_n` are total unitaries, including their dirty-data sectors. Hence
`G` is unitary and direct multiplication gives

```text
H = G (A_N tensor I_data) G^*,
A_N = (1/2) sum_n (|n+1><n| + |n><n+1|).
```

The sine vectors `sin((n+1)k pi/(N+1))`, for `k=1,...,N`, diagonalize
`A_N`. Every eigenvalue

```text
lambda_k = cos(k pi/(N+1))
```

has multiplicity exactly `D` in the full Hamiltonian. Its norm is
`cos(pi/(N+1)) < 1`. No low-energy or prepared-history projection is used
in this conclusion. The coupling `1/2` is exactly representable in binary
and does not depend on the horizon or program. The spectrum describes
this finite clock carrier, with no particle-mass interpretation.

## Packet moments on the infinite line

Choose integers `m>=1` and `a>=0`, with `a+m<=p`. The normalized initial
clock packet is supported entirely before the first program edge:

```text
phi_(a+j) = (-i)^j binomial(m,j) / sqrt(binomial(2m,m)),  0<=j<=m.
```

Thus `P_n=I` throughout its support, and the prepared clock and data are
uncorrelated. The normalization follows from Vandermonde's identity.
Extend the same packet by zero to the infinite integer lattice. On that
lattice write `S|n>=|n+1>`, `A=(S+S^*)/2`, and `X|n>=n|n>`. The bounded
velocity is

```text
V = i[A,X] = i(S^*-S)/2,     [A,V]=0,
X(t) = exp(i t A) X exp(-i t A) = X+tV.
```

These identities apply to the finite-support packet and its finite
position moments. Write `f_j=binomial(m,j)/sqrt(binomial(2m,m))`, with
`f_j=0` outside `0,...,m`. Vandermonde's identity gives

```text
sum_j f_j f_(j+1) = m/(m+1),
sum_j f_j f_(j+2) = m(m-1)/((m+1)(m+2)).
```

The phase `(-i)^j` makes the velocity expectation positive. The position
probabilities are the hypergeometric weights
`binomial(m,j)^2/binomial(2m,m)`. Their mean and variance, and the velocity
moments, are

```text
mu_0 = a+m/2,                 v = E[V] = m/(m+1),
s_0 = Var(X) = m^2/(4(2m-1)),
E[V^2] = (1 + m(m-1)/((m+1)(m+2)))/2,
b = Var(V) = (2m+1)/((m+1)^2(m+2)).
```

The symmetrized covariance of `X,V` vanishes. To see this directly, each
neighboring pair contributes `(j+1/2) f_j f_(j+1)` to `E[(X-a)V]`;
reflection about the packet center gives `(m/2) E[V]`. Consequently the
infinite-line evolution has the exact moments

```text
mu(t) = a+m/2+vt,             s(t) = s_0+b t^2.
```

Let `c=p+4` be the first completed clock position. For `mu(t)>c-1`, the
one-sided variance inequality gives

```text
Pr_infinite[X<c] <= q(t)
  = s(t) / (s(t) + (mu(t)-(c-1))^2).
```

This bound retains every unfinished position. It does not divide by an
arrival probability. If `B=c-1-mu_0>0`, then `mu(t)-(c-1)=vt-B` and

```text
d/dt [ (s_0+b t^2)/(vt-B)^2 ]
  = -2(bBt+v s_0)/(vt-B)^3 < 0       for t>B/v.
```

Hence `q(t)` decreases on that whole interval. Its value at the left
endpoint certifies a continuous readiness window; a time grid is not
needed for this step.

## Finite boundary and a uniform readiness bound

Embed the finite clock in the infinite line by cutting the bonds
`(-1,0)` and `(N-1,N)`. The resulting operator `A_cut` is the direct sum
of `A_N` and two disconnected exterior half-lines. Both `A` and `A_cut`
have norm at most one. Set

```text
d_boundary = min(a+1, N-a-m).
```

A walk from the initial support needs at least `d_boundary` hops to leave
the finite rail. Thus `A^k phi=A_cut^k phi` for every `k<d_boundary`.
Expanding the two exponentials and using their norm bounds gives, uniformly
for `0<=t<=t_max`,

```text
eta = 2 sum_(k=d_boundary)^infinity t_max^k/k!,
||exp(-it A)phi-exp(-it A_cut)phi|| <= eta.
```

For `t_max<d_boundary+1`, the ratios of successive tail terms are at most
`t_max/(d_boundary+1)`. The factorial bound
`d! >= (d/e)^d > (d/3)^d` therefore supplies the rational enclosure

```text
eta <= 2 (3 t_max/d_boundary)^d_boundary
       / (1-t_max/(d_boundary+1)).
```

Pure-state trace distance is at most vector distance. Every event
probability, including the entire unfinished event, consequently changes
by at most this `eta`. Combining the two bounds gives

```text
Pr_finite[clock<c] <= q(t_min)+eta
```

throughout any window `[t_min,t_max]` satisfying the stated conditions.
Gauge conjugacy makes this probability independent of the source and of
every untouched quantum reference.

The declared parameters are

```text
N=4096,   a=2048,   m=64,   p=2112,   c=2116,
t_min=256,   t_max=512.
```

Here `d_boundary=1984`, `v=64/65`, `s_0=1024/127`, and `b=43/92950`.
At the left endpoint,

```text
mu(256)-(c-1) = 14109/65,
s(256) = 226536448/5902325,
q(256) = 226536448/278318778205.
```

The finite-boundary enclosure is

```text
eta <= (3970/1473) (24/31)^1984 < 2^(-659).
```

The last inequality follows from `(24/31)^3<1/2` and `3970/1473<3`.
Exact rational comparison then yields

```text
Pr_finite[clock>=2116] > 999/1000
                   for every real t in [256,512].
```

The interval has duration `256` in the supplied dimensionless Hamiltonian
parameter. It is a uniform single-snapshot statement. It supplies no
joint probability for an unmeasured trajectory, repeated outside
measurements, or observer proper time.

## Coherent adaptive instrument and unfinished outcomes

For first address `a` and old bit `b`, let
`F_(a,b)=sum_(x:x_a=b) |f_a(x)><x|`. Let `Q_(a,r)` project onto source
words with addressed bit `r`. The coordinates of the completed adaptive
isometry are

```text
L_(a,b,r,b') = (1/sqrt(12)) F_(6r,b') Q_(a,r) F_(a,b).
```

The corresponding retained data label contains `A_3=a`, `R_1=b`,
`R_2=b'`, `C_3=(r,1)`, `B=6r`, and `P_A=P_B=(r,1)`, together with all
blank complementary registers and the final source in `S_4`. Denote this
complete isometry by `E=sum_h |label_h> tensor L_h`, with the source
factor understood in `S_4`. The labels are orthogonal. For a fixed
address, `sum_b F_(a,b)^* F_(a,b)=I`; summing the second old bit, observed
bit, first old bit and first address in that order proves `E^*E=I`.
Its basis action is precisely the four-gate action above.

The public instrument measures only the actual pointer pair. For the
completed outcome `P_A=P_B=(r,1)`, its operator is

```text
E_r = sum_(a,b,b') |label_(a,b,r,b')> tensor L_(a,b,r,b').
```

For an arbitrary source/reference density `rho_SR`, its unnormalized
output is `(E_r tensor I_R) rho_SR (E_r^* tensor I_R)`. The sum over private
histories occurs before forming this density: address and old-bit archive
coherences within a public outcome are retained. Replacing it by a sum
of separate history densities would specify a different instrument.
The source and the reference also remain quantum.

Let `E_j` be the prepared data isometry after `j` program gates and
`j(n)=min(4,max(0,n-p))`. The exact joint evolved state has clock
amplitudes `alpha_n(t)` and data isometries `E_(j(n))`. Tracing the clock
and applying the pointer instrument therefore gives the complete channel

```text
Phi_t = sum_(n=0)^(N-1) |alpha_n(t)|^2 Phi_(j(n)),
```

where each `Phi_j` retains the private data and untouched reference.
For the declared blank preparation and communicated program,
`P_B.valid=1` exactly for `j=4`. Thus the pointer register itself supplies
an unfinished flag; it is not manufactured from a clock postselection.
Writing `p_t=Pr[clock>=c]`,

```text
Phi_t = p_t Phi_4 + (1-p_t) Phi_unfinished.
```

All unfinished outcomes are retained. The two terms have orthogonal
public validity flags, and their weights are independent of the input.
The half-diamond distance from the ideal completed instrument is exactly
`1-p_t`: its upper bound follows from the channel mixture, and the
validity-flag measurement attains the matching lower bound. This equality
allows arbitrary untouched references. In particular its value is below
`1/1000` throughout the declared window. Conditioning on a positive-weight
completed outcome gives the exact adaptive native instrument, while the
unconditional statement keeps its failure probability.

A separate computational diagnostic measures the clock and every data
basis label. For a source-basis preparation `x`, its full law has weight
`|alpha_n(t)|^2/12` on each of the twelve data words obtained by applying
prefix `j(n)` to `(x,a)`. For a fixed prefix these words are distinct,
because the total data gate is a permutation. The full `N`-position clock
law and the five twelve-word prefix tables thus describe every supported
diagnostic outcome without repeating the same data table for each clock
site. Positions of zero weight are retained. This diagnostic is a
different instrument from pointer-only readout; its classical tables do
not certify private or source/reference coherences.

## Measurements, controls and finite recurrence

The Hamiltonian is conserved during free evolution. The prepared packet
has zero mean energy because adjacent amplitudes have imaginary products.
For a complete computational-basis diagnostic, every Hamiltonian term is
off-diagonal in the measured clock basis, so nonselective readout also has
zero mean energy. This does not specify an apparatus Hamiltonian or a
zero-cost implementation.

Pointer-only dephasing can change energy for a general clock/data state.
For this preparation it has zero mean energy change as well. Indeed the
real path Hamiltonian preserves the form
`alpha_n(t)=(-i)^(n-a) r_n(t)` with real `r_n(t)`. If `Q_l` are the pointer
projectors, the data coefficient on a dephased forward bond is the
expectation of

```text
P_n^* [sum_l U_n^* Q_l U_n Q_l] P_n.
```

Each `U_n` and `P_n` is a basis permutation, so the bracket and its
conjugate are diagonal self-adjoint operators. Their expectations are
real even for an arbitrary source/reference density. Multiplication by
`alpha_(n+1)^* alpha_n`, which is imaginary, gives zero real contribution
after adding the adjoint bond. This proves the stated mean-energy result
for the actual pointer instrument, without transferring it to arbitrary
preparations.

The uphill, no-communication, no-feedback and no-publication programs
retain their own specified gates and complete outputs. Their clock
spectrum and clock-sector probability are unchanged whenever their gates
remain unitary. Their pointer semantics and adaptive native channel need
not agree with the communicated program. In particular a completed clock
sector does not turn an invalid pointer into a valid record under a
communication or publication control.

Finite-dimensional unitary recurrence implies arbitrarily close returns
to the full initial state along an unbounded sequence of external times.
This includes the pointer blanks, transported buffers and native archive.
No exact `pi` return formula is asserted for the uniform rail. The
readiness interval is therefore a finite persistence guarantee for
single terminal snapshots, not a permanent archive or indefinitely
operating dissipative native process.

The proof concerns the ideal normalized packet and exact finite unitary
law. The implementation bounds the binary packet vector error by an exact
rational `epsilon`. If `v_hat` is that represented vector and `v` the ideal
unit vector, then `||v_hat/||v_hat||-v|| <= 2 epsilon`. Unitary propagation
preserves this bound, which adds at most `2 epsilon` to every event
probability and the half-diamond bound for the supplied product preparation.
The receipt retains this preparation allowance separately and verifies that
the combined bound is below `1/1000`. Matrix exponentiation and probability
evaluation have numerical residual checks without a certified SciPy
roundoff claim. A time-grid
check is an implementation check, not the certificate for the whole
interval. No source extraction certificate, physical clock calibration,
continuum limit, particle dictionary, or transfer to the official
CCF-LT-2 Hamiltonian follows from this construction.
