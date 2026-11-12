# Stress and dressed clocks of the autonomous wave source

## Same action, projected force

This package uses exactly the Hamiltonian, source ground, reference domain
and time unit of `openwave_dynamic_wave`. It does not replace that source
by its static mean geometry. Wave and dressed-clock experiments below have
separately specified preparations on the same full12800-state carrier.

In a simultaneous geometry-coordinate sector, write
Kq=exp(2q)Lexp(2q) and h=sqrt(Kq). Varying source coordinate Q_a changes
q by F_a. With D_a=diag(F_a), the unique derivative G_a satisfies

    h G_a + G_a h = 2(D_a Kq + Kq D_a).                 (1)

This follows by differentiating h²=Kq; h>0 makes the Sylvester map
invertible. The same spectral geometry conjugation as in the source
constructs the complete field force on all geometry and wave states.
Equation(1) is a derivative of the actual scalar action, not an imported
continuum stress formula. The harmonic force is the projection
G_a=sum_x F_xa(dh/dq_x). The uniform ambient-profile identities are

    sum_x dh/dq_x = 2h,    sum_a dH_native/dQ_a = H_native.       (2)

For the first identity, K(q+s)=exp(4s)K(q). The factor2 includes both the
declared lapse and spatial metric dependence; it is not a factor1
rest-clock response. Uniform ambient variation is an action derivative,
not a dynamical zero mode of the Dirichlet geometry. No equation of state
or Einstein stress tensor is inferred from(2).

## Exact finite Ehrenfest equation

Define the complete geometry momentum by
P_a=κ sum_b M_ab i[Hg,Q_b]. The identities [N,a]=−a and [N,a*]=a* hold
at a product cutoff, including its boundary completion. Therefore

    dQ/dt=(κM)⁻¹P,       i[Hg,P]=−κS Q.

Finite Q,P do not have exact canonical commutators. Let V=H_native+H_field,
F_a=dV/dQ_a and R_a=i[V,P_a]+F_a. Then the complete source obeys

    dP_a/dt = −κ(SQ)_a − F_a + R_a.                  (3)

R is evaluated and retained. Equivalently, its field part solves the
Sylvester equation whose right side is i[Kq,P_a]+dKq/dQ_a. This follows
by commuting h²=Kq with P, so no canonical commutator is silently used.
The boundary terms are appreciable at cap1. A cap2 numerical comparison
reports them separately and is not a proof of a cap-independent limit.

The geometry, native-coupling and field energy rates are their exact
commutators with H. Their sum vanishes identically. Executed wave and
vacuum controls have the same native preparation and Hamiltonian. At t=4,
the field changes mean geometry coordinates by approximately
(−0.00686739,−0.00778769); geometry energies are about0.0329070 with the
wave and0.00442191 with vacuum. No external geometry update supplies that
response. With couplings removed, the source factorizes and its geometry
vacuum is stationary.

## Complete stress and momentum observations

For each field-force component, its full spectral projectors give196 wave
outcomes, plus all4 vacuum/geometry outcomes at eigenvalue0. A separate
momentum instrument measures all4 geometry-momentum values and every
wave site or vacuum. Every complementary native sector is retained.
The instrument sums to the identity and extends to arbitrary quantum
references. Momentum and position readouts are different experiments.
Each complete apparatus can exchange at most30 source-energy units.

The force and instrument numbers are floating evaluations with explicit
Sylvester and commutator residuals. They are not asserted to be interval
enclosures of every force eigenprojector. The underlying unitary state
errors are bounded by the wave package's dyadic certificates. The exact
identities(1)–(3) do not depend on those floating residuals.

## Complete dressed native clocks

In the field-vacuum sector, let E_ij be the geometry ground energy of
Hg+E_i exp(Q0)+E_j exp(Q1), where the E_i are the native excitation gaps.
The producer certifies all four eigenvalues of all64 native sectors with
exact dyadic residuals and interval source entries. The clock at the
second cell has dressed gap Δ_j(i)=E_ij−E_i0. For j=1,...,7, the ratio
Δ_j(1)/Δ_j(0) is bounded by interval division of positive gaps.

The seven frequency shifts range from approximately−1069.0878ppm to
−981.1887ppm. Their spread is strictly in
[87.8991138,87.8991463]ppm. Thus one exact redshift factor does not describe
these finite-strength, backreacting native clocks. This is a finite-source
result, not evidence against an observed equivalence principle.

Each clock is prepared in the coherent sum of its two stationary dressed
branches. All native-X and geometry-number outcomes are retained; the
wave and other native complements have probability0 for this preparation.
The overlap of the two geometry ground vectors determines visibility.
These source-consistent preparations are distinct from the initially
unentangled geometry/wave packet. Spectral preparation is a declared
service. Complete record values are numerical evaluations; the frequency
gap and nonuniversality intervals have independent residual certificates.

## Quantitative test-clock domain

Fix background E_i with i=0 or1 and consider
A(s)=Hg+E_i exp(Q0)+s exp(Q1), 0≤s≤E_7. This auxiliary family describes a
continuously tunable probe coefficient; the native clocks use s=E_j.
The coordinate bound |Q_a|<0.067 and ω_min>349/500 imply a uniform gap
greater than9/50: the total coefficient is below77/20 and
exp(0.067)−exp(−0.067)<671/5000. Scalar shifts of the positive potential
give the Weyl gap bound directly. The radius of exp(Q1) about its interval
midpoint is below17/250.

For its unique ground energy f(s), finite-dimensional variational
concavity, Hellmann–Feynman and the reduced resolvent imply

    f'(s)=<exp(Q1)>_s,
    |f''(s)| ≤ 2(17/250)²/(9/50),
    |(f(s)−f(0))/s − <exp(Q1)>_0|
       ≤ [(17/250)²/(9/50)] s.                       (4)

The scalar midpoint cancels from the reduced resolvent. Integrating the
second-derivative bound proves(4). For the auxiliary coefficient s=1/1000,
the effective-lapse error is at most289/11250000, about2.56889×10⁻⁵.
That coefficient is a controlled probe intervention, not an observed
native mass or an additional empirical fit. Equation(4) states exactly
how a universal test-clock response is recovered as backreaction vanishes,
while retaining the finite-strength differences of the executed clocks.

## Scope

The package derives the dynamic force, its finite boundary correction and
a quantitative clock-response domain from one declared source. It does
not supply a relativistic conserved stress tensor, continuum geometry,
physical G, laboratory clock identification, a physical universality test
or nonlinear Einstein dynamics. The scalar metric choice and harmonic
geometry architecture retain their input status.
