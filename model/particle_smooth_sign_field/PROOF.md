# A smooth field action with controlled local quantum records

A positive smoothing parameter gives a declared gauge-coordinate action
whose full interaction is exponentially summable in spatial support.
Its polynomial truncations approximate complete local detector records
uniformly in exterior volume, on one state and clock. This result treats
the quantum electric field and every Fock occupation; it requires no
Wilson gap or link-derivative hypothesis. Spacing refinement requires
a separate argument. Zero smoothing is required only if the chosen
target is the discontinuous-overlap action.

## 1. Source conventions and the changed regulator

Use the actual M15 representation, sixty Dirac/CAR modes per cell,
coordinate-link multiplication, positive electric Casimirs and periodic
Wilson kernel of [the field source](../particle_observer_outputs/FIELD_ROUTE.md).
Its Hermitian kernel `H_W=beta(aD_W-I)` obeys `||H_W||<=5<8` for every
unitary coordinate configuration, without a gap assumption. Choose
`delta>0`, `a>0` and `kappa>=0`, and define

```
s_delta(H)=H(H²+delta² I)^(-1/2),    c_delta=sqrt(1+delta²),
h_delta=(beta+c_delta s_delta(H_W))/a,
H_delta=H_E+dGamma(h_delta)-kappa sum_p Re Tr_M15(U_p)/15. (1)
```

Equation(1) is a different declared action from the discontinuous
overlap sign and from its retained Legendre polynomial. It does not
establish equality of their records. The coefficient c_delta follows
from the flat-kernel value `H_W(0)=-beta`: it makes `h_delta(0)=0` and
fixes the leading free Dirac coefficient to one. It contains no measured
mass or fitted coupling. Every comparison in this proof targets this
same smooth action, with its own explicit polynomial family below.

The source norm satisfies
`||h_delta||<=(1+sqrt(1+delta²))/a<=(2+delta)/a`. This gives
full-Fock norm density `b_delta=60(2+delta)/a`. It replaces the kinetic
allowance of any different regulator when applying energy conservation.
The coordinate completion and all apparatus services are admitted
design choices. The proof makes no claim that a narrower native repair
alphabet implements them without the existing compiler contract.

## 2. A uniformly convergent finite-path expansion

Set

```
D=64+delta²,    q=64/D<1,    u=1-q,
Z=(64I-H_W²)/D,    c_n=binomial(2n,n)/4^n.
```

The source bound implies `0<=Z<=qI`. The ordinary binomial series gives

```
c_delta s_delta(H_W)
 =sum_(n>=0) F_n,
F_n=(c_delta/sqrt(D)) c_n H_W Z^n.                     (2)
```

It converges in norm on every coordinate configuration, uniformly over
those configurations and finite volumes. It defines the same bounded
coordinate multiplication operator on the quantum link Hilbert space.
The truncation `p_M=sum_(n=0)^M F_n` has coordinate-word degree and
path range at most `2M+1`. Its coefficients are those of (2), not the
old Legendre approximation. Because c_n>=0 and the scalar series has
one sign at each x, `|p_M(x)|<=c_delta` throughout `|x|<=8`.
Thus every truncated one-body Hamiltonian has the same bound
`||(beta+p_M(H_W))/a||<=(2+delta)/a`.

At order n, `c_n<=1` and `8/sqrt(D)<=1`, so

```
||F_n/a|| <=(1+delta)q^n/a.                            (3)
```

At infinite volume, every finite coefficient is defined directly by
the same finite Wilson paths. Their local interaction estimates below
construct the limit without positing an infinite-volume matrix square
root as extra source data.

## 3. Exponentially summable full-Fock interactions

Let P_x project onto all sixty one-body modes at cell x. As in the
[existing locality proof](../particle_field_causality/PROOF.md), set

```
Phi_(x,n)=dGamma((P_x F_n+F_n P_x)/(2a)).
```

The trace-norm estimate on the rank-sixty projector gives
`||Phi_(x,n)||<=C q^n`, with `C=60(1+delta)/a`, on the entire Fock
space. These even interactions are supported in the ball of radius
`d_n=2n+1`, including every link appearing in its paths. Their sum over
x is exactly the nth source term. A configuration-dependent coefficient
is treated as a bounded operator on its actual quantum link factors;
it is not replaced by a classical background.

Choose the rational spatial weight

```
lambda=1+u/8>1,    r=q lambda^4<1.
```

Expanding `(1-u)(1+u/8)^4` gives constant1, linear coefficient−1/2
and nonpositive higher coefficients for0<=u<=1. Hence
`r<=1-u/2`. The number N_d of cubic cells in the radius-d ball is
at most `(2d+1)^3`. A cell belongs to at most N_d rooted stars and
each star contains at most N_d cells. The weighted interaction sum
therefore obeys

```
sup_y sum_(Z contains y) |Z| ||Phi_Z|| lambda^diam(Z)
 <= C lambda² sum_(n>=0) N_(2n+1)² r^n
 <= C lambda² 4096 sum_(n>=0)(n+1)^6 r^n.
```

The generating series of `(n+1)^6` has numerator
`1+57r+302r²+302r³+57r⁴+r⁵` and denominator `(1-r)^7`.
Its numerator is at most720. Adding the twelve incident plaquettes,
each of support size four and diameter two, gives the finite bound

```
W_delta=C lambda² 4096*720/(1-r)^7+48 kappa lambda².     (4)
```

The coordinate-independent beta term is onsite. Electric Casimirs and
that term are removed in the onsite interaction picture. Their unitary
conjugation preserves the support, parity and norm of every Phi_(x,n).
No bound on derivatives of a sign function with respect to quantum
links is used. The control comes from convergence of the full local
interaction expansion, including its link dependence.

For unbounded electric generators use their bounded spectral clippings,
exactly as in the owning locality proof. Finite-volume evolutions
converge strongly on vectors as the clipping ceiling grows. All local
norm bounds are uniform in that ceiling, so they pass to the limit.
The finite source is self-adjoint on the electric domain because its
remaining interaction is bounded at finite volume. On the full local
`B(L²(G))` algebra this construction does not assert point-norm
continuity in time.

## 4. A complete local-record perturbation bound

Let `N=M+1` be the first omitted term in (2). Its tail interaction has
the same weighted norm bounded by

```
J_tail <= A (N+1)^6 r^N,
A=C lambda² 4096*720/(1-r)^7.                          (5)
```

Indeed `(N+j+1)^6<=(N+1)^6(j+1)^6`, so summing the shifted tail
reduces to the same generating-series bound. Both the full smooth
interaction and its truncation have weighted norm at most W_delta.

The graded commutator-chain proof of the owning field theorem works
with lambda in place of two. Including the zero-length chain for
possibly overlapping supports gives

```
||[alpha_t(A_X),B_Z]_gr||
 <=2||A_X||||B_Z|| exp(2 W_delta |t|)
       sum_(x in X,y in Z) lambda^(-dist(x,y)).         (6)
```

The lattice weight sum is finite and explicit:
`S_lambda=sum_(z in Z³) lambda^(-|z|_1)=((lambda+1)/(lambda-1))³`.
Summing (6) over the omitted interaction terms and applying Duhamel's
identity proves, for any bounded local operator,

```
||alpha_T^smooth(A_X)-alpha_T^M(A_X)||
 <=2|X| ||A_X|| T exp(2 W_delta T) S_lambda J_tail.     (7)
```

The bound is independent of total volume, state, fermion number and
untouched reference dimension. It acts on the whole quantum field,
so local electric or link observables are covered as bounded operators.
It is not an inference from an expected one-body Hilbert–Schmidt norm.

For k specified intermediate local instrument dilations and a final
local output, use the same full input state, policy and pointer services.
Defer the recorded outcomes coherently and telescope the instruments
and the final output. The existing complete-history argument gives

```
half_trace(output_smooth,output_M)
 <=min(1,(2k+1)2|X| T exp(2 W_delta T) S_lambda J_tail). (8)
```

This keeps all outcomes and untouched references. Apparatus supports
must lie in the specified observation region and use the even/graded
algebra or explicit compensating references required by the source.
The output does not include the whole unobserved dynamical exterior.
No stationarity, thermal state, translation symmetry or low-mode
probability is required for (7)–(8).

## 5. Fixed-time planning without giant exponentials

For rational delta,a,kappa,T, all majorants in (4)–(8) are rational.
Since `r<=1-u/2`,
`r^N<=2^(-floor(uN/2))`. Also
`exp(2 W_delta T)<=2^ceil(4 W_delta T)` because e<4.
Consequently a sufficient logarithmic error exponent is

```
ceil(log2[(2k+1)2|X|T S_lambda A])
 +ceil(4 W_delta T)+6ceil(log2(N+1))-floor(uN/2).        (9)
```

Integer doubling of N makes (9) smaller than every prescribed negative
bit target. The code calculates exact integer/rational inequalities
and stores the final polynomial degree in factored resource metadata.
It does not materialize that enormous polynomial or a giant exponential.
T is fixed throughout this search.

The mixed [local gauge state](../particle_local_gauge_state/PROOF.md)
has exactly matching source marginals across volumes. Therefore its
truncated polynomial records have the
[fixed-time local limit](../particle_local_history_family/PROOF.md).
Taking M large using (8) gives a norm-Cauchy local smooth-field limit.
The state and clock remain those of (1) and its approximants. This
does not identify canonical purifiers from different volumes or the
global exterior state; the reference restrictions of the owning local
history proof remain in force.

For a concrete finite simulator, the code also chooses an open box
for that same binomial polynomial. The observation cells are contained
in a declared radius-r cube, and the boundary separation is
`D_boundary=R-r-2max(d,2)>0`. Its one-body norm is bounded by
`(2+delta)/a`, its range is `d=2M+1`, and its weighted locality rate
is bounded by W_delta uniformly in M. With the boundary constants
of the owning volume proof, use lambda^(-D_boundary) in place of
2^(-D_boundary). Since
`log(lambda)>=u/(8+u)>=u/9`, its tail is bounded by
`2^(-floor(u D_boundary/9))`. A sufficiently large integer radius
exists at every fixed T. The box input is the actual full local
marginal with outgoing electric spectators and retained boundary
charge references.

Finally apply the existing representation-cutoff theorem to the
binomial polynomial, with word step `max(4,2M+1)`, full-Fock norm density
`b_delta`, and initial degree `3|box|L`. The generic integer cutoff
planner consumes these quantities only; its legacy degree-key name
does not import Legendre coefficients. Complete local occupation
pointers preserve coordinate degree. The whole-history Dyson count
includes all intervals before imposing the cutoff, with no restart
of initial support after an instrument.

The default certificate assigns20 error bits to each of smooth-to-
polynomial, polynomial-to-open-volume and representation compression.
Their sum is at most `3/2^20` at T=1. The action, source state, positive
electric generator, kappa and clock are shared within these comparisons.
Actual preparation, compiled apparatus and numerical evolution errors
are additional. This is a finite resource certificate, not execution
of the enormous state or pulse tape.

## 6. Source controls and the free spatial symbol

The exact native sixty-column control from the local Wilson calculation
has H_W² eigenvalues0,1,3 with ranks16,28,16. At delta=1 the normalized
smooth term has squared eigenvalues0,1,3/2. Its sixteen zero-mode
columns are retained and mapped continuously to zero. The code checks
these complete projectors and an explicit small binomial approximation
without a floating spectral admission. This coordinate control is
not called a normalizable gauge-field preparation.

For free coordinate links put `k_i=a p_i`,
`w=sum_i(1-cos k_i)-1`, `s²=sum_i sin² k_i` and
`R²=w²+s²=1+2sum_(i<j)(1-cos k_i)(1-cos k_j)`. The actual kernel is
`H_W=w beta+sum_i alpha_i sin k_i`, with `alpha_i=i beta gamma_i`.
Then

```
h_delta² = a^(-2)[(1+c_delta w/sqrt(R²+delta²))²
                  +c_delta² s²/(R²+delta²)] I.        (10)
```

Both summands are nonnegative. A zero requires every sine to vanish
and w<0. Among all eight corners only k=0 has w=−1, and normalization
by c_delta makes that corner exactly zero. At every other corner
`w=2j-1>=1`; `c_delta w/sqrt(w²+delta²)>=1`, so its positive energy
is at least2/a. Both full Dirac branches and every internal component
are retained; this is a vectorlike census.

At fixed delta and fixed continuum momentum,
`R²=1+O(a⁴|p|⁴)`. Expanding sine and cosine in (10)'s Hermitian
matrix symbol gives

```
h_delta(p)=sum_i alpha_i p_i + (a|p|²/2) beta + O(a²) . (11)
```

Thus the leading free Dirac velocity is one in the declared spacing
and time convention. A finite polynomial truncation approximates this
normalized action with its own error; it need not have an exactly
zero corner at finite order. No interacting pole, physical mass,
spin assignment or chiral removal follows from this free calculation.

## 7. State, energy and refinement boundaries

The default uses the fixed-L referenced mixed preparation and kappa
of its degree-zero density planner. It changes the regulator action
explicitly to delta=1. The certificate recomputes the kinetic allowance
as `b_delta=180` at a=1 before quoting any same-action roughness bound.
It does not reuse the old regulator's kinetic constant silently.
The record comparison itself is state independent and uses no roughness
bound. The maximally mixed full Fock input is not a Lorentz vacuum,
and its conjugate charged reference remains part of the preparation.

The positive delta action avoids the discontinuity by an explicit
change of architecture. Exact Ginsparg–Wilson or chiral-projector
identities of the discontinuous sign cannot be transferred to it:
`s_delta(H_W)²` is generally not identity. A chiral interacting theory
requires its own regulator and matter construction.

As delta tends to zero, the spatial weight lambda tends to one and
the displayed bounds diverge. No estimate uniform in that limit is
claimed. Recovering the discontinuous-sign architecture would require
that extra limit; it is not necessary for every permissible field
architecture. A fixed positive delta can instead be retained during
spacing refinement, provided its own interacting and chiral claims
are proved.

There is also an exact simpler endpoint of the declared design family:
as delta tends to infinity, `c_delta s_delta(H_W)` tends in operator
norm to H_W on its bounded spectrum. Hence
`h_delta -> (beta+H_W)/a`, the ordinary finite-range Wilson
Hamiltonian. The convergence is uniform for `||H_W||<=8`; it follows
directly from `sqrt(1+delta²)/sqrt(x²+delta²)->1`. The finite-range
Wilson action may be chosen directly. It retains the full vectorlike
Dirac space and is not a chiral Standard Model construction. It shows
that an exact-sign comparison is not an obligatory intermediate step
when the architecture itself is an allowed design choice.

Spacing refinement, interacting state selection, surviving
connected responses and a common stress/geometry map are also separate
requirements. This fixed-delta thermodynamic construction supplies no
measured Standard Model masses, scattering match or Einstein limit.
