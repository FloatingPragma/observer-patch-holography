# Dilution of bounded local ground observables

Let the connected fair binary source have V ports, fixed total occupation q,
and exchange generator G=sum_e(I-S_e). Its unique sector ground is
Omega_(V,q)=binom(V,q)^(-1/2) sum_(|x|=q)|x>. The positive transfer P and
counted Hamiltonian H=-2M log(P) have the same ground. This is the ground
of the specified source, not the faithful reference P/Tr(P). The construction
below applies to controlled H or G evolution; it does not identify unread
repair with that unitary evolution.

For a region S of s ports, let Q_S be its empty-occupation projector. Exact
counting gives

    p_s = 1-<Q_S> = 1-binom(V-s,q)/binom(V,q) <= s q/V.       (1)

When q>V-s the numerator is zero. An equivalent formula for <Q_S> is
product_(j=0..s-1)(V-q-j)/(V-j), with a zero factor interpreted exactly.
It needs only s factors, irrespective of the capacity.

## Bounded observables and spectral weight

Let A be supported in S and preserve its total occupation: [A,N_S]=0.
The empty local sector is one-dimensional, so A Q_S=c_A Q_S and
Q_S A=c_A Q_S, with |c_A|<=||A||. Therefore

    |<A>-c_A| <= 2||A|| p_s,
    <A* A>-|<A>|^2 <= 4||A||^2 p_s.                       (2)

The first inequality follows by restricting A-c_A to the occupied block.
For the second, the expectation minimizes ||(A-c)Omega|| over scalars c;
choose c_A and use ||A-c_A||<=2||A|| on that block.

For another such B supported on r ports, and any real controlled times,
Cauchy-Schwarz and invariance of Omega under H and G give

    |<A* alpha_t(B)>-<A*> <B>|
        <= 4||A|| ||B|| sqrt(p_s p_r).                    (3)

There is no restriction on t or on a capacity-dependent clock rescaling.
The total connected spectral measure of A is its variance, so it obeys
(2) as well. Thus, if q/V tends to zero, every uniformly bounded observable
on a uniformly bounded number of ports in this number-preserving dictionary
is asymptotically indistinguishable from its scalar empty-sector statistics;
any convergent limit is scalar, and its connected spectral weight is zero.
This includes the fixed-q=2 family. A putative nonzero particle pole in this
particular dictionary has residue tending to zero, whatever rescaling is
chosen for its energy location.

A necessary quantitative visibility condition follows without fitting any
mass. If ||A||<=K with K>0 and the connected spectral weight is at least
w>0, then

    s q/V >= w/(4K^2).                                   (4a)

For q>0 this requires at least ceil(w V/(4K^2 q)) source ports. At q=0,
positive connected weight is impossible within this dictionary. This is a
necessary condition only; meeting it does not construct a pole or identify
a particle.

This does not exclude fields with increasing support, norms growing with
capacity, a positive-density sector q/V->p>0, or charged operators on an
admitted sum of occupation sectors. Those are different observable/state
attachments. It also does not assert a vanishing spectral gap: poles and
their observable residues are distinct quantities.

## Complete local instruments, records and source state

Take a fixed finite protocol of complete local instruments. At step j and
for every recorded prefix h, require fine Kraus operators K_(j,h,a) supported
on at most s_j source ports, with

    sum_a K_(j,h,a)* K_(j,h,a)=I,
    [K_(j,h,a),N_(j,h)]=0.                                (4)

Every zero, failure and hidden fine Kraus slot is retained. Coarse outcomes
may be obtained only afterward. The condition in (4) is branchwise local
number preservation, not merely conservation of mean number or conservation
of a larger source-plus-charged-ancilla system. A supplied nonempty charged
ancilla can invalidate the condition and is not a free resource here.

On the empty local block write K_(j,h,a) Q=c_(j,h,a) Q. Completeness gives
sum_a |c_(j,h,a)|^2=1. Define the comparison instrument by these scalar
amplitudes acting on the unchanged source. Its record probabilities may
depend on the prefix; it is an ordinary complete classical protocol with
no source response beyond the local empty-sector amplitudes.

Let V_j and V_j^0 be complete Stinespring isometries for the actual and
scalar instruments. They agree on the empty block, and each has norm one.
For the common ground,

    ||(V_j-V_j^0)Omega|| <= 2 sqrt(p_(s_j)).               (5)

This bound also holds for a classically controlled prefix: its orthogonal
record blocks average the squared error, bounded uniformly by the right
side. Insert any record-controlled H or G unitary intervals between steps.
They preserve Omega exactly and do not affect the scalar-prefix source.
Replace instruments in a telescoping order for which the preceding prefix
is scalar. Every later complete isometry has norm one. The full purified
terminal vectors then differ by at most

    epsilon = 2 sum_j sqrt(p_(s_j)).                     (6)

Consequently their complete terminal source-plus-record states have half
trace distance at most min(1,epsilon), and every complete classical history
law has total-variation distance at most that value. An untouched reference
can be retained: the pure initial source ground factors from it. No
postselection, multiplication by the number of outcomes, or removal of
rare branches occurs. Adaptive stops can be padded by an identity instrument.
The bound is for a fixed declared finite protocol; it does not stay small
for arbitrarily many interventions or source regions growing with V.

The proof covers the ideal number-preserving B and diagonal source services
and literal endpoint C/E instruments on their stated local supports. Stored
non-TP numerical maps require their separate complete-error/gain budgets.
A service specified only on an occupied encoded subspace is not automatically
a complete instrument on the empty block required by (4).

## Exact native endpoint witness

One literal fair E attempt has two raw-zero branches P_even/sqrt(2).
In the source ground, unequal endpoint occupation has probability

    r = 2q(V-q)/(V(V-1)).                                (7)

The scalar empty-sector comparison retains both fair coins and raw zero.
The actual classical raw-record law differs from it by exactly r in total
variation. Retaining the quantum source gives the exact distance

    d = r/2 + sqrt(r-3r^2/4).                             (8)

For each accepted coin, the difference is one half of
(1-r)|e><e|-|Omega><Omega|, where
Omega=sqrt(1-r)|e>+sqrt(r)|o>. Its two-dimensional trace norm is
sqrt(r-3r^2/4). The two accepted records and all odd records are orthogonal;
the odd records contribute total trace r, proving (8). This also shows why
small classical event probability does not justify dropping source/reference
coherences from a complete comparison.

The finite checker executes the public native interpreter on the literal
H66 ground (V=12,q=2), checks all 25 raw slots and the full accepted
operators, and recovers r=10/33 exactly. Formulae (1)-(8) concern the same
binary endpoint law on each federation level. A capacity-independent Bell
experiment from a separately supplied known occupied A0 B1 input does not
contradict dilution of the delocalized Dicke ground.

## Density-mode frequency corollary

The exact density closure in `model/federated_observer/SCALING.md` states,
on every source occupation sector,

    G F_f = F_(L f),     F_f(x)=sum_i f_i x_i,

where L is the actual port-graph Laplacian. Because Omega has constant
configuration amplitudes, the same identity is

    G(F_f Omega)=F_(L f) Omega.

For a nonconstant mean-zero mode Lf=lambda f with sum_i |f_i|^2=1,

    H(F_f Omega)=h(lambda) F_f Omega,
    ||F_f Omega||^2=q(V-q)/(V(V-1)),
    h(lambda)=-2M log(1-lambda/(2M)).                     (9)

The residue follows directly from the hypergeometric two-point law. This
is an exact ground-sector spectral branch of the interacting exclusion
source, not an assumption of independent particles. On the source transfer
floor alpha, the ground comparison gives

    lambda <= h(lambda) <= lambda+lambda^2/(4M alpha).    (10)

Thus whenever the inherited spatial scaling supplies lambda/(2M)->0,
the controlled density frequency follows the same Laplacian scaling.
A common clock multiplier cannot replace it by a square-root dispersion or
change ratios between these modes. The explicit periodic-cycle control in
`federated_observer/SCALING.md` has its first two gapless density energies
in ratio tending to four, for every nontrivial q. That control graph is not
the selected federation. The selected-source conditional diffusive and
two-dimensional spatial limit is the separately specified dictionary in
`native_spatial_scaling/PROOF.md` and its inherited field theorems.
Equation (9) identifies the controlled ground frequency on that same
Laplacian dictionary; it supplies no different spatial continuum.

This prevents identifying these gapless density modes with a linearly
dispersing relativistic particle merely by taking a transfer logarithm or
changing clock units. It does not exclude hydrodynamic or nonrelativistic
subbranches of an interacting relativistic theory. A massive relativistic
particle also has a quadratic low-momentum excess energy; a nonzero mass,
other dispersion regimes and their source observables would require a
separate identification. These statements do not classify every protected
source sector as a density mode.

The code's dyadic square-root upper bounds have a declared finite precision,
default 128 bits and at most 192 bits. They are safe finite estimates with
an additive rounding floor of one grid unit per root; they do not compute
arbitrarily small asymptotic rates. The analytic bounds (1)-(10) do not have
that numerical floor.
