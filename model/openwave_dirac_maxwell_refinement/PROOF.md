# Pauli-root compact-field records

## Source and native carrier

The actual48-mode native lepton Bogoliubov frame supplies two spin coordinates
for each charged family, each antiparticle, and both internal cells. Its full
96-row particle/hole probe is retained. The charged energy ratios and electron
anchor are inputs from the frozen common-mass construction. The kinetic
prescription is a declared extension on the complete excitation Fock space;
the neutral and other native sectors retain their spectral and identity
complements. No strong-source or signed-Dirac unitary equivalence is assumed.

Every link has L²(U(1)) space with all integer electric fluxes. On a plaquette,
Gauss reduction with charge q=3Q and compensating quantum reference−q gives
E_e=n−q1_(e<x). Native hops alter exactly one link flux and preserve every
local constraint. The complete one-excitation carrier has x=0,…,3, both spin
states and every n∈Z. At each link configuration the covariant Pauli operator
is D=σxΠx+σyΠy. Its square contains iσz[Πx,Πy]; this term is retained exactly.
The matter energy is m+sqrt(m²+D²)−m+(r/2)L_cov, with r=1/16. The positive
field energy is (1/32)ΣE²+(1/8)(1−cos loop). Both components are positive,
so their closed form sum defines the autonomous source. The electric operator
is unbounded and no flux cutoff defines the source. Second quantization of
the positive matter operator preserves excitation number. An A-dependent
signed-Dirac diagonalization, pair currents and transformed field conjugates
are not asserted.

At zero threaded probe phase, the two spin blocks are complex conjugates and
the cross-spin blocks vanish exactly. Both spectra are retained by this
intertwiner. The executed initial spin is up; every opposite-spin measurement
outcome is retained with zero probability. The external threaded phase is
exp(iqθ) on one matter link. Its first and second derivatives include the
noncommuting derivatives of every power of D². Separately, conjugation by the
relative dipole supplies its exact current and contact term. Native-third
charge units and electron-charge units differ by3 for currents and9 for
contacts. A threaded flux is not a pure-gradient Ward control.

## Infinite flux and numerical errors

Here ||D²||≤c²=1/128. The degree-six binomial polynomial approximates
sqrt(m²+D²)−m with error at most mρ⁷/(1−ρ), ρ=c²/m², for m≥1. Before
compression,2degree|q|+2 padded flux rows retain every polynomial return path.
The prescription never squares a clipped Dirac operator.

With W_s=exp(s|n|), s=1/2, each charged hop gains at most exp(s|q|).
The weighted root series converges when ρ_s=ρ exp(2s|q|)<1. After placing the
unbounded electric diagonal and scalar Wilson/magnetic terms in the exact
interaction picture, the bounded generator has norm at most
B_s=mρ_s/(1−ρ_s)+(1/8)exp(s)+r exp(s|q|).
Thus a bare-flux input obeys ||Q_Cψ(t)||≤exp(−s(C+1)+|t|B_s).
Duhamel gives the full-to-Galerkin error at most
(1+|t|B_0)exp(−s(C+1)+|t|B_s). Rational Taylor sums and outward dyadic rounding
bound every exponential. The polynomial error is added separately.

The finite eigensystem is complete. Real and imaginary eigenvector coordinates
are rounded to denominator2²⁶, the Hermitian matrix and eigenvalues to2²⁴.
Integer-overflow guards precede exact complex Gram and residual products.
The exact polynomial-to-dyadic matrix error includes the displayed matrix
rounding and binary64 coefficient/summation bounds. At degree six the
Pauli-square powers are exact Gaussian-integer dyadics, with numerator
magnitude at most4⁶. Fixed outward Gram/residual/source budgets are accepted
only after their complete live checks. The bound
||exp(−itH)−Q exp(−itΛ)Q*||≤g+|t|r sqrt(1+g)
is applied with the additional source-matrix error. Rational phases, including
the parent's independently checked pi interval, are compared with every actual
complex output amplitude. The preparation-weighted total-variation budget
also includes probability squaring. Native probe weights are displayed
floating diagnostics; the certified table is position/flux plus failure.

The positive outside electric floor (C+1−|q|)²/8 and bounded crossing operator
supply a Schur-complement lower bound on each full infinite one-excitation
energy. The true full rotor vacuum is subtracted. These reference-dependent
sector energies are neither isolated physical poles nor new empirical fits.

## Complete apparatus and work

The input is matter vacuum with bare cycle flux0, field energy1/8.
Preparation has success1/2 and an explicit failure identity. Success creates
the Gauss-dressed native excitation; every native/Fock complement remains in
the instrument. The terminal detector absorbs the charge at its recorded
position and retains the charge reference and all link fluxes. Separate
apparatus states store every native96-row and position/flux outcome. Mean
preparation work is (m+initial kinetic+initial electric)/2; the field constant
cancels. Mean terminal work is (−m−kinetic+1/8−magnetic)/2. Their sum equals
half the final electric energy, the change in total source energy including
the unmeasured failure branch. The actual records verify this identity.

## Free spatial refinement

On the three-dimensional free lattice the same prescription has
E_a=√(m²+c²Σsin²(ap_i)/a²)+(r/a)Σ(1−cos(ap_i)),
while the quadratic Maxwell symbol is
ω_a=(2c/a)√Σsin²(ap_i/2). The common speed satisfies c²=2EB=1/128.
The absolute errors against √(m²+c²|p|²) and c|p| are bounded by
ca²|p|³/6+ra|p|²/2 and ca²|p|³/24 respectively. These follow directly from
|sin x−x|≤|x|³/6 and1−cos x≤x²/2. All eight corners are retained;
seven have positive Wilson lift2r times corner weight/a. At r=0 all eight
zero-kinetic corners are explicitly retained. This free-symbol convergence
proves no interacting boost symmetry, pair threshold or absence of a preferred
regulator frame.
