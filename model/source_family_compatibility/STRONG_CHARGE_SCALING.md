# Volume-uniform charged spectral localization

The actual centered ground and its original fermion operators have a
volume-uniform strong-contact spectral bound. For the ballistic spacing
law `c=v/a`, a family with `U/c -> infinity` sends the original charged
response out of every fixed physical-energy window. This conclusion uses
no projection onto a singly occupied source space and no bound containing
the total number of edges. It is a statement about original response
probabilities, not a hard spectral gap or a particle identification.

## Same-source response vectors

Use the invariant full centered ground, signed retarded clock, and exact
source operators of `COMPOSITE_CAR.md`. Write `alpha=U/2>0`. The isometries
V0 and V1 are the retarded vectors of the original annihilator a and the
actual composite `b=2(n_opposite-1/2)a`. Their identities are

    V0*V0=V1*V1=I,  V0*V1=0,
    L V0=V0 h+alpha V1,
    Y=L V1-alpha V0,  V0*Y=0.                         (1)

L is the selfadjoint signed retarded generator. Its two components use
the same positive physical excitation clock with opposite time signs.
They retain the separate addition and removal measures. No negative
physical excitation energy or additional independent fermion is supplied.
For a maximum-degree bound d and hopping convention c/2 per spin,

    ||h|| <= c d/2,
    ||Y|| <= 5 c d sqrt(1+d^2).                       (2)

The second estimate is the exact odd local-current anticommutator bound
in the composite proof. In its decomposition `Y=V1 C+R`, the two summands
are orthogonal. Thus the theorem retains the actual residual R; setting
it to zero is unnecessary.

The normalized full ground projectors satisfy (1)-(2). The source
translation averages and their stated thermodynamic ground subsequence
preserve them. Their constants do not depend on volume, U, or the chosen
unit lattice smearing f. Finite-support smearing proves the domain
identities first; the bounded maps in (1)-(2) extend them to all square
summable f by closure of L. Transporting the source chart transports the
state and the original/composite probes together. This argument makes no
identification of repairs in different charts.

## Spectral bands with original probabilities

Put `W_+=(V0+V1)/sqrt(2)` and `W_-=(V0-V1)/sqrt(2)`. These are orthogonal
isometries. Their exact residuals are

    (L-alpha)W_+=(V0 h+Y)/sqrt(2),
    (L+alpha)W_-=(V0 h-Y)/sqrt(2),
    [(L-sign alpha)W_sign]*[(L-sign alpha)W_sign]
       =(h^2+Y*Y)/2.                                 (3)

In particular let

    beta^2=c^2 d^2 [101/4+25d^2].                    (4)

Then `||h^2+Y*Y||<=beta^2`. A measured smaller upper bound on this actual
positive Gram can replace (4). For any Borel set F whose distance from
both alpha and -alpha is at least w>0, the spectral theorem gives

    ||1_F(L)W_sign f|| <= beta ||f||/(sqrt(2)w).

Since `V0=(W_++W_-)/sqrt(2)`, the triangle inequality gives the original
signed response bound, including its factor:

    mu_f(F)=||1_F(L)V0 f||^2
             <= min(1,beta^2/w^2) ||f||^2.            (5)

This applies to the complement of the two intervals centered at
`+/-U/2` with radius w. For a fixed physical window `[0,E]`, if E<alpha,

    mu_add,f([0,E])+mu_remove,f([0,E])
       <= min(1,beta^2/(alpha-E)^2) ||f||^2.          (6)

Zero-energy contributions from both components are included. Each
original branch obeys the same upper bound; no division by preparation
success is used. A local normalized a_(x,s) has exact original addition
and removal weights 1/2, so its bound can also be capped at 1/2 per branch.
For general smeared f only the sum of the two masses is exactly ||f||^2;
onsite density I/2 does not imply that each smeared mass equals one half.

Equation (5) supplies certified broad energy windows for an existing
complete source energy instrument, with its rejected bins and preparation
failures retained. It does not construct that instrument's reference
service or its finite resolution. A measured interval with a boundary
uncertainty must enlarge the admitted band or keep that boundary weight.

## Scaling consequences and their scope

The exact sufficient condition for loss of fixed-window original weight is

    alpha(a)>E eventually,  c(a)/(alpha(a)-E) -> 0.   (7)

For `c(a)=v/a`, v>0, and `U(a)/c(a)->infinity`, condition (7) holds for
every fixed E. It holds uniformly in the volume and the unit smearing,
including packets selected around any kinetic valley and momentum-dependent
smearings. The entire fixed-energy original canonical charge response
therefore vanishes on this family. This result does not exclude source-even
collective observables, nonlocal charged sectors, or another capacity law.
It does not establish a positive gap for the full many-body spectrum.

The ratio condition alone is insufficient if c tends to zero. At U=1 and
c->0 the same source approaches the atomic response at physical energy
1/2. Its weight stays inside any fixed window with upper endpoint greater
than 1/2. The complete dimer family below verifies this countercontrol.

There is also a controlled description in the interaction-normalized
clock. Set `W=[V0,V1]`, `X=[[0,I],[I,0]]`. The two columns of

    L W-alpha W X=[V0 h,Y]

are orthogonal by (1), so this map has norm at most
`gamma=5 c d sqrt(1+d^2)` for d>=1. Duhamel yields

    ||exp(-it L/alpha)W-W exp(-itX)||
         <= min(2,|t| gamma/alpha).                  (8)

In particular the original retarded correlation differs from
`cos(t) I` by at most this bound. When c/U->0 its signed response
converges to equal atoms at scaled frequencies -1 and +1. The two
original physical charge measures converge to mass 1/2 at scaled energy
1. This follows from (5), (8), and exclusion of weight near zero; it does
not require an exact half-mass assertion for finite smeared fields.
The statement concerns response vectors and correlations. It is not an
arbitrary-input quantum-channel approximation or a construction of a
regular interacting field hierarchy in that scaled clock.

An energy rescaling q(a)>0 rescales the same source clock. The physical
locality velocity bound is a fixed constant times `q(a)c(a)a`. If
`U/c->infinity`, a finite nonzero interaction-normalized band energy
`qU/2` and a->0 imply `qca->0`. Thus that energy convention makes the
proved source locality cone collapse at the same spatial scaling.
Keeping the ballistic clock instead gives (7). This is a constraint on
this particular common spatial/clock scaling, not on all source families.

A charge-clock shift has an additional exact limitation. Let
`K_mu=K-mu N` and keep the same source reference. Gauge covariance gives
addition energies `E_add-mu` and removal energies `E_remove+mu`. If
mu>=0, the shifted removal weight in [0,E] is contained in the original
removal weight in [0,E]; if mu<=0, the same statement holds for addition.
More sharply the branch whose energy increases uses the original window
[0,E-|mu|], which is empty when |mu|>E. Hence at least one conjugate local
branch obeys (6) for every choice mu(a). No scalar charge shift keeps both
original half-weight branches in bounded energy windows under (7).
Arbitrary mu need not preserve the ground property of the reference;
negative shifted excitations signal that failure. This argument uses the
exact charge-frequency shift and makes no ground-state assertion for K_mu.

The complementary weak-contact result `0<=U=o(c^(2/7))` in
`../ballistic_free_boundary/README.md` identifies the full canonical
hierarchy as free. Neither theorem selects a unique behavior in the
intermediate class between those sufficient bounds. The moment bounds
`C,R=O(c)` alone do not make the retarded two-block approximation converge
at a fixed physical time or fixed spectral resolution when c diverges.
In particular an uncontrolled two-block eigenvalue is not a mass readout.

## Cell-normalized contact coefficients

On a finite lattice region define the piecewise constant cell field
`Psi_a(x)=a^(-3/2) a_n` on the cell of volume a^3 at n. This is an exact
normalization of the lattice CAR algebra; it is not a claim that an
unsmeared continuum field operator exists. Direct cell integration gives

    integral (Psi_up* Psi_up)(Psi_down* Psi_down) dx
       =a^(-3) sum_n N_nup N_ndown,
    N=integral sum_s Psi_s* Psi_s dx.

Consequently the exact quartic coefficient in this cell notation is
`g(a)=a^3 U(a)`. The same centered clock includes the term `-U(a)N/2`;
it cannot be discarded while comparing its ground response. Equivalently
one may use centered cell densities and retain the scalar constant.

If `c=v/a` and `g(a)->g0>0`, then `U/c=g(a)/(v a^2)->infinity`. Equation
(6) gives the explicit original fixed-window bound

    mu_add,f([0,E])+mu_remove,f([0,E])
      <= min(1, v^2 d^2[101/4+25d^2] a^4
                    /[g(a)/2-Ea^3]^2) ||f||^2,       (10)

whenever the denominator before squaring is positive. Thus the direct
finite nonzero cell-contact scaling has O(a^4) original odd response in
every fixed energy window. Holding that coefficient nonzero does not
produce a finite-energy interacting original canonical field on this
selected centered ground and physical clock. This excludes that specific
trajectory; it does not exclude a critical finite U/c trajectory, a
renormalized different observable, or a source-even collective sector.

## Complete finite-source checks

For integers n>=3 and any rational scale r>0, take the actual full dimer
with

    c=2nr, U=(n^2-4)r, sqrt(U^2+4c^2)=(n^2+4)r.

Its full ground energy is -n^2 r. Each local charge preparation has the
two exact positive-energy atoms

    E_-=r(n^2/2-n+2), E_+=r(n^2/2+n+2),
    w_-=(1+4n/(n^2+4))/4,
    w_+=(1-4n/(n^2+4))/4.                            (9)

The implementation verifies the full positive shifted matrix, unique
ground, complete charge-sector spectral projectors, the actual original
weights, and both residual identities (3) by rational arithmetic. It
retains the preparation complement of weight 1/2 and all zero spectral
outcomes. These finite source values test the theorem's normalization;
they do not substitute a dimer for the cubic continuum.

The choice r=1/2 has c=n, U=(n^2-4)/2 and realizes ballistic spacing
`a=1/n` with diverging original charge energies. The choice
`r=1/(n^2-4)` has U=1 and realizes the atomic countercontrol. Neither is
fitted to experimental masses.
