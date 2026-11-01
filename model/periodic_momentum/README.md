# Conserved periodic momentum, winding and the locality boundary

The same local Weyl/contact action has **three exactly conserved global
lattice momenta** on a declared periodic box. They are Hermitian spectral
logarithms of the graded translations, generally nonquadratic on Fock space.
This supplies an interacting conserved quantity beyond the quadratic
classification in [lattice_momentum_balance](../lattice_momentum_balance/README.md).
It does not supply local conserved momentum columns or an Einstein stress
tensor. Periodic boundary conditions are an explicit input; the previous
open-box action does not acquire translation symmetry.

There is also a constructive local approximation on an initially smooth,
bounded-number sector, with a finite-time error including interaction
leakage. No sector is silently made invariant.

## Same action and complete graded translations

On a Cartesian torus, with two CAR modes at each site, retain

```
H = sum_(x,i) [ (ic/2) c_x^dag sigma_i c_(x+ei) + h.c. ]
      + U sum_x n_x0 n_x1 .
```

Each active period is at least three, so a directed wrapping bond is not a
duplicate short-ring edge. A side of length one denotes an inactive axis.
The local coefficients and contact are unchanged. The four-site z-ring
is exactly the existing [particle_response](../particle_response/dispersion.py)
action, including its wrap orientation. The cubic convention agrees with
[bound_pair](../bound_pair/protection.py). The algebraic results hold for
all finite real c,U; the executable source domain retains nonnegative
exact rational c,U and at most 64 sites, with c=1/8,U=1 in the receipt.

Let t_i shift one-particle modes by one site. Its complete exterior lift
T_i=Gamma(t_i) fixes the vacuum and multiplies an occupation column by the
sign of the induced permutation of its occupied modes. Consequently

```
T_i c_xs^dag T_i^dag = c_(x+ei),s^dag,
T_i^Li = I,       [T_i,T_j]=0,       [H,T_i]=0.
```

These identities hold on every number sector, its coherences and arbitrary
untouched references. The last equality follows by permuting the entire
bond and onsite sums. The exact Clifford certificate checks H, H0 and V
separately under all three shifts on the complete 54-mode `3x3x3` algebra:
406 nonzero Majorana coefficients, nine covariances, all nine commuting
shift pairs and their orders. This is symbolic coverage of all `2^54`
states, not execution of those amplitudes. Exterior lifting is a group
representation, so no fermion-parity factor remains after a full period.

## A conserved Hermitian logarithm, with its branch fixed

Write zeta=exp(2 pi i/L) and let b(j) be the representative of j modulo L
in `(-L/2,L/2]`. Define

```
Pi_j = (1/L) sum_(r=0)^(L-1) zeta^(-jr) T^r,
P    = (2 pi/L) sum_(j=0)^(L-1) b(j) Pi_j.
```

The finite geometric sum gives orthogonal Hermitian projectors summing to
I, with `T Pi_j=zeta^j Pi_j`. Thus P is Hermitian, `exp(iP)=T`,
`||P||<=pi`, and `[H,P]=0`. The three Cartesian P_i commute. Here “principal”
means Arg in **(-pi,pi]**. At an even-period eigenvalue -1 we retain +pi;
the alternate -pi convention is the distinct conserved logarithm
`P-2 pi Pi_(L/2)`. No numerical logarithm or eigenvalue cut defines either.

The complete 256-state ring certificate uses Gaussian rational matrices
for P/pi and all four projectors, each of exact rank 64. It checks their
products, H commutation and both -1 conventions. Three elementary Lean
theorems in `Lean/QFT/PeriodicMomentum.lean` verify conservation and mutual
commutation of the finite scalar polynomials and the integer branch
coefficient contradiction below. The source CAR and spectral/locality
attachments are the separate exact certificate and proofs here.

## Additive momentum fails by an explicit integer winding

Let p be the same branch logarithm of the **one-particle** t, and put
Q=dGamma(p). On a momentum occupation basis, Q adds individual momenta,
whereas P reduces that sum to its principal representative. Exactly,

```
Q = P + 2 pi W,
W = (sum occupied b(j) - b(sum occupied j))/L in Z.
```

Hence W is a well-defined Hermitian integer-spectrum operator. It is not
asserted to be conserved. On the retained ring its multiplicities at
`W=-1,0,1` are `4,152,100`, and `W^3-W=0` as an exact full matrix.
The vacuum fixes the additive scalar in any proposed quadratic formula
for P to zero; the one-particle sector fixes that quadratic to Q. But on
the completely filled ring P=0 and Q=2 pi. Therefore P is not a scalar
plus a homogeneous number-preserving quadratic.

The contact vertex gives a more informative two-particle witness. Use
normalized opposite-spin Fourier pairs with the `exp(-ikx)` convention.
The incoming labels `(1,1)` and outgoing `(-1,-1)` on L=4 have the same
T eigenvalue -1 and P=pi. Their Q values are pi and -pi respectively.
The exact matrix element is

```
<out|V|in> = U/4,
<out|[H,Q]|in> = pi U/2.
```

The contact changes additive momentum by 2 pi while preserving crystal
momentum. This is an explicit Umklapp transition, not a failure of global
conservation. Full matrices verify `[H,Q]=2 pi [H,W]` at c=1/8,U=1.

## Every branch has a finite-range obstruction on odd-prime rings

This statement concerns **any number-preserving Hermitian** Q with
`exp(iQ)=T`, not just a polynomial chosen above and not just quadratic Q.
Number preservation makes its one-particle restriction q a logarithm of
t. Every matrix commutes with its own exponential, so q commutes with t.
At momentum j its two-dimensional spin block has eigenvalues
`2 pi j/L + 2 pi integer`. Consequently

```
tr_spin(q_j) = (2 pi/L) (2j + L n_j),     n_j in Z.
```

If the spin trace of any nonzero-separation position kernel vanished,
then the integer polynomial

```
A(z) = sum_(j=0)^(L-1) (2j+L n_j) z^j
```

would vanish at a nontrivial L-th root. For **odd prime L**, that root is
primitive and its rational minimal polynomial is
`Phi_L(z)=1+z+...+z^(L-1)`. For completeness, irreducibility follows by
Eisenstein at L after the substitution z=w+1: all nonleading coefficients
of `((w+1)^L-1)/w` are divisible by L, and its constant L is not divisible
by L squared. Since degree A is at most L-1, divisibility forces all its
coefficients to be equal. But their residues `2j mod L` are distinct.
This is a contradiction. The executable certificate retains arbitrary
symbolic integer n_j, performs the polynomial remainder, and checks odd
prime periods 3,5,7,11,17,31,97. It does not replace all branches by a finite
branch search.

Now define a range-R density sum precisely as `A_R=sum_X a_X`, where every
a_X is supported in a set of sites of torus graph diameter at most R.
Arbitrary many-body terms, number-changing terms and coefficients are
allowed. If sites x,y have distance greater than R, no summand contains
both. For each summand at least one outside-site number operator commutes
with it and distinguishes `|x,s>` from `|y,t>`, so its one-particle matrix
element vanishes. Thus every such sum has a banded one-particle kernel.
It cannot equal any of the number-preserving logarithms just classified
when the prime ring contains a separation beyond R. The number-preserving
assumption belongs to Q, not to each proposed local term.

This excludes a uniformly finite-range exact density for these translation
logarithms. It does **not** exclude a different conserved observable,
number-changing logarithms, algebraically decaying tails, or a continuum
identification after a controlled restriction of states and errors.

## Quantitative locality boundary for the principal branch

With `q_xy=<x|p|y>`, r=x-y modulo L, direct differentiation of the finite
geometric sum gives, for r nonzero,

```
L odd:  q_xy = (i pi/L) (-1)^r csc(pi r/L),
L even: q_xy = (pi/L) (-1)^r [1+i cot(pi r/L)].
```

The diagonal is zero for odd L and pi/L for even L. The same kernel acts
on each spin. In either case its magnitude is
`pi/(L sin(pi r/L))` for `1<=r<=floor(L/2)`. Therefore, for any integer
R>=0 with L>2(R+1), every range-R density sum satisfies the full-Fock bound

```
||P-A_R|| >= |<R+1,s|(P-A_R)|0,s>|
           = pi/[L sin(pi(R+1)/L)] >= 1/(R+1).
```

This bound uses only one-particle matrix elements, so arbitrary higher-body
local corrections do not evade it. It does not require A_R to preserve
number or be translation invariant. The quantitative bound is for the
**principal** P; the all-branch theorem above claims nonzero tails, not the
same uniform magnitude for arbitrary integer branches. A bound in full
operator norm is not automatically a bound on selected low-energy states.

## Controlled initially smooth sector and its interaction leakage

There is a positive restricted approximation. On any finite periodic
box, let Omega project onto N<=Nmax with every occupied momentum satisfying
`|k_j|<=kappa` along the axis considered, with **Nmax kappa<pi**. This
guarantees no wrap of the total initial momentum. The nearest-neighbor
quadratic S_j=dGamma(sin(p_j)) then satisfies

```
||(P_j-S_j) Omega|| <= Nmax kappa^3/6 = delta,
```

using `|k-sin k|<=|k|^3/6`. Its sign follows the translation convention;
it is the negative of the D_j lift in the previous momentum package.
H0 preserves Omega, but the contact generally does not. On N<=Nmax,
`||V||<=U floor(Nmax/2)` for U>=0. Duhamel comparison with the free
evolution yields, for every finite model time t,

```
||(I-Omega) exp(-itH) Omega|| <= ell
     := min(1, |t| U floor(Nmax/2)),
||(P_j-S_j) exp(-itH) Omega|| <= delta + (pi+Nmax) ell.
```

The second line splits the evolved state with Omega and uses
`||P_j-S_j||<=pi+Nmax` on the retained number sector. These are operator
bounds, hence also cover reference-entangled inputs in Omega. They are
uniform in volume at **fixed Nmax**, not at fixed particle density.
The rational evaluator uses the sufficient domain Nmax kappa<3 and
pi<22/7; at Nmax=2,kappa=1/4,t=1/100,U=1 it reports an evolved error below
0.057. No asymptotic scaling of U or physical model time is inferred.

The exact zero-momentum pair on the four-site ring has contact matrix
element U/4 to the pair `(1,-1)`, outside any cutoff below pi/2. This
negative control explicitly retains the leakage term. Small initial
momenta alone do not justify projecting the interacting action.

## Scope and replay

The complete ring also rejects deleting Fock permutation signs or using
the periodic charge on the old open action: their squared commutator
Hilbert--Schmidt norms are exactly 4 and 2 at c=1/8,U=1. Twelve ideal
adjacent FSWAPs implement its T on all 256 occupation columns. They use
the existing FSWAP primitive; this exact word identity is **not** a new
certificate of literal pulse error or a compiled measurement of P.
No additional reference preparation, source readout, repair covariance,
hardware precision or physical time is assumed to follow from it.

The finite global conservation premise is closed for the supplied
periodic action. Local momentum Ward columns, a symmetric Lorentz stress,
a physical chart/time and scaling law, geometric modular weights,
entropy/area stationarity, Newton normalization and an Einstein-register
inhabitant remain separate. The supplied action and boundary are explicit
model data, not unexplained assumptions being claimed as discoveries.

```sh
python -m model.periodic_momentum --check model/periodic_momentum/results.json
python -m pytest -q -W error model/periodic_momentum
cd Lean && lake env lean QFT/PeriodicMomentum.lean
```
