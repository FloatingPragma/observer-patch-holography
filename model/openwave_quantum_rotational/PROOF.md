# Quantum angular momentum and calibrated charged records

The construction consists of fermionic matter in complete spatial oscillator
shells and a quantized transverse Maxwell mode space. Its finite Hamiltonian
and complete instruments define simulation-generated records. The native
Clifford representation, positive energy law, rotational regulator, trap,
charge assignment, coupling and preparation services specify this member.
No identity with the periodic cubic source or its three-lepton spectrum is
assumed.

## Carrier, spatial action and quantum field

Let V be the span of the three-dimensional oscillator functions with total
number at most two, oscillator length one. It has dimension ten. The actual
projected differential operators are p=P(-i grad)P and
L_i=-i epsilon_ijk a_j* a_k. The native sixteen-dimensional Clifford module
has the exact unitary decomposition into four complete Dirac blocks given
by the source intertwiner. All sixteen internal coordinates are retained.
There are160 fermion modes, with full antisymmetric Fock space of dimension
2^160. Their common declared charge is-1; Q=-N_ferm. This global charge
assignment does not supply a local longitudinal Gauss-law construction.

The transverse spatial modes are

    A_b(x)=exp(-|x|²/5) (b cross x)/2.

Their radius satisfies R²=5/2. Direct Maxwell integration gives
c0=pi^(3/2)R^5/4 and c1=5pi^(3/2)R³/8=c0. Set kappa=c0/32. The field
frequency sqrt(c1/c0) is one. With eta=1/8, use
b_i=eta(a_i+a_i*) and pi_i=(a_i-a_i*)/(2i eta) for the three vector
oscillators. The field carrier retains every state of total number at most
C, with C in{1,2,3}. Its dimensions are4,10,20. This complete total-number
cutoff commutes with rotations.

The field angular momentum is the actual mode Noether generator
Jf_i=-i epsilon_ijk a_j* a_k. It obeys the exact rotation algebra and makes
b a vector even at the cutoff. The projected free Maxwell energy, with its
vacuum constant subtracted, is Hf=Nf. Squaring already-truncated canonical
coordinates instead omits a boundary term. Precisely,

    eta² sum pi_i² + sum b_i²/(4 eta²) -3I/2
      +(C+3)P_top/2 = Nf.

The complete boundary term is retained. Finite oscillators do not satisfy
an infinite-dimensional canonical commutator on their highest shell.

## Positive kernel and full fermion extension

In one native Dirac block define

    B = sigma dot p + sum_i M_i b_i,
    M_i = sigma dot P[exp(-|x|²/5)(e_i cross x)/2]P,
    D = tau_3 m + tau_1 B,                 m=16.

The sign in M uses the declared charge-1. Since the two tau matrices
anticommute, the full absolute Dirac operator is

    |D| = I_tau sqrt(m²+B²).

Both signed Dirac branches and every native block are retained. Computing
one20-dimensional spin/spatial block tensored with the field is an exact
multiplicity reduction: the full one-particle carrier consists of eight
identical such blocks. Their multiplicity is not a physical species count.
Let K denote the complete160-by160 matrix of field operators supplied by
|D|. Define the full source Hamiltonian by its ordered CAR lift,

    H = sum_ab c_a* K_ab c_b + dGamma(Nm) + Nf.

The field Hamiltonian is included once. This is a linear self-adjoint
Hamiltonian on the full finite matter/field carrier. Its unitary evolution
and complete measurement instruments extend to every external quantum
reference.

Positivity of this operator-valued lift requires its block order. On the
one-particle/field space K-mI is positive. Write K-mI=B0*B0, with operator
blocks (B0)_ra, and put d_r=sum_a (B0)_ra c_a. Photon operators commute with
CAR operators, so the exact identity

    sum_ab c_a* (K-mI)_ab c_b = sum_r d_r* d_r >=0

holds on the full Fock space. Thus

    H >= m N_ferm + dGamma(Nm) + Nf.

The source vacuum is the unique zero-energy source state. Every source
sector with at least two fermions has energy at least32. A separate idle
compensator has its declared zero energy; this does not assert uniqueness
of a vacuum after adjoining an arbitrary idle reference. A two-mode full
Fock control verifies positivity, while transposing only the matter block
indices gives a negative eigenvalue-1/2.

## Exact dressed half-spin

Let Jm=dGamma(L+Sigma/2), and let J=Jm+Jf. The projected differential
operators and Gaussian multiplication coefficients obey exact vector
covariance. Hence [B,J]=0, and the functional calculus and ordered CAR lift
give [H,J_i]=0 for every component. The generator rotates actual matter
wavefunctions and the Maxwell mode amplitudes. It is distinct from a
canonical-spin label or a logarithm of a cubic permutation. The quantum
matter and Maxwell torques cancel as operators:

    i[H,Jm_i]+i[H,Jf_i]=0.

The exact free kinetic square on the20-row block satisfies
T(T-5I/2)=0. Also ||M_i||<=R/2<4/5, and
||b_i||<=2 eta sqrt(C). Therefore

    ||B|| <= a_C = 8/5 +(3/5)sqrt(C),
    0 <= sqrt(m²+B²)-m <= epsilon_C I,
    epsilon_C=a_C²/(2m)<1.

The receipts use rational upper bounds for sqrt(C). On a fixed native
block and tau sign, H0=m+Nm+Nf has a two-dimensional lowest subspace P0:
spatial vacuum, photon vacuum and the two native spin coordinates. All
other eigenvalues are at least m+1. For
H_s=H0+s(sqrt(m²+B²)-m), s in[0,1], min-max bounds isolate exactly two
levels below m+epsilon_C. Every higher level is at least m+1.

For a low eigenvector u, the complementary resolvent gives
||Q0 u|| <= epsilon_C/(1-epsilon_C) ||P0 u||. Thus P0 maps the isolated
subspace bijectively onto the initial spin doublet. It intertwines rotations,
since both projections commute with J. The dressed subspace carries exactly
the same irreducible spin-half representation. Consequently

    J²=3/4 I

on the entire dressed doublet, and rotational invariance makes its energy
one scalar E_C. Its Jz=+1/2 line is a full interacting matter/field state.
The full sixteen-dimensional lowest one-particle band consists of identical
copies. E_C<17<32 also makes it the lowest source charged energy among all
fermion-number sectors. The Maxwell part alone has integer angular spectra;
the half-integer Casimir belongs to the coupled matter and Maxwell fields.

The finite executions measure a nonzero Maxwell angular contribution,
approximately1.2056062e-5 in the cap-two up state, with the matter contribution
completing total Jz=1/2. The dressed Casimir does not rely on either separate
expectation being1/2. Arbitrary entangled inputs exhibit opposite nonzero
matter and Maxwell torques. Omitting the field torque fails the control.

## Preparation and complete localized records

A retained charge+1 reference fermion compensates the source charge. A
specified even pair-creation pulse prepares the source oscillator-vacuum,
spin-up seed and the opposite reference from their vacua. The reference
Hamiltonian is zero. The pulse and its work are declared apparatus inputs.
The full source energy projector onto[m,m+epsilon_C], together with its
complement, is a complete binary filter on the entire Fock/field carrier.
Because C<=3<m, no empty-source photon state lies in this window; sectors
with two or more fermions lie above it. Conditional success prepares the
exact dressed up line. The resolvent bound gives original success at least
1-[epsilon_C/(1-epsilon_C)]². Numerical evaluations give its actual success.
Every filter failure and the retained compensator are included.

For a spatial ball, the position effect is P chi_ball P on the finite
oscillator space, tensored with the other factors. Its complement completes
the observation. The code integrates the actual Hermite products using
incomplete gamma functions. For even monomial powers2p_i, the integral over
a radius-r ball of exp(-|x|²) product x_i^(2p_i) is
product Gamma(p_i+1/2) times the regularized incomplete gamma function with
shape sum p_i+3/2 and argument r². Odd integrals vanish.

The sum of squared oscillator basis functions through shell two equals
pi^(-3/2)exp(-r²)(2r^4+5/2). Cauchy-Schwarz therefore bounds every normalized
state's outside-ball probability by

    10 erfc(r) + (4/sqrt(pi)) exp(-r²)(r^5+5r³/2+5r).

At r=5, erfc(r)<=exp(-r²)/(r sqrt(pi)), 1/sqrt(pi)<4/7 and the positive
Taylor lower sum for exp(25) give the rational upper bound in the receipt,
less than1.1e-7. This is a localized complete-shell construction with a
supplied trap; it is not a free translating soliton. Radial position,
photon number and total angular observables commute with one another.
The position/photon success records and filter failure form2(C+1)+1 original
records. A separate complete angular instrument retains the up outcome, its
full angular complement and preparation failure. The photon number
measurement can disturb the stationary
state; its work is reported separately from premeasurement conserved energy.

## Certified finite energies and magnetic current

The source probe is the actual uniform-field derivative E=partial_B B,
with ||E||<4/5. For K=sqrt(m²+B²), X=partial_B K solves

    KX+XK=BE+EB.

Thus ||X||<=a_C(4/5)/m. This is the current derivative of the same joint
quantum Hamiltonian; no Pauli moment is inserted. In the dressed up line,
the static magnetic normalization is g_C=4 E_C |<X>|, using total Jz=1/2
and the same vacuum-relative rest energy.

All interval calculations use exact algebraic matrix coefficients and an
explicit trial vector of Gaussian integers divided by the square root of
their exact squared norm. Floating eigensolvers choose this trial vector;
the ensuing interval certificate does not assume that it is an eigenvector.
With S=B², the degree-six Taylor polynomial p6(S) of sqrt(m²+S) is a lower
operator bound. Its positive remainder is at most

    delta_C=(33/2048) a_C^14/m^13.

This follows from the positive seventh scalar derivative and spectral
calculus. The lower-polynomial Rayleigh value and variance are evaluated
with outward interval arithmetic. If their intervals give Rayleigh R and
residual norm r, then the actual Rayleigh value is in[R,R+delta_C], and its
residual is at most r+delta_C. The exact excited-band floor m+1 gives the
Temple lower bound

    E_C >= R_lower -(r_upper+delta_C)^2/(m+1-R_upper-delta_C),
    E_C <= R_upper+delta_C.

Degeneracy of the spin doublet does not change this argument, since its
lowest energy is scalar and the next distinct energy is at least m+1.
The energy residual and exact Jz residual bound the trial's distance from
the unique up line. For x=a_C²/m²<1, the derivative tail is bounded by

    (2 a_C (4/5)/m) (33/2048) x^6(7-6x)/(1-x)^2.

Indeed the absolute binomial coefficients decrease from the seventh term,
and differentiating each noncommuting power S^k yields k ordered products.
The convergent geometric derivative sum gives the displayed bound. The
code evaluates all ordered products, then adds this tail and the bounded
observable error due to the certified state distance. The resulting energy,
current and g intervals concern this finite source, with no physical
uncertainty interpretation.

## One calibration and finite cap comparison

The cap-two ground energy is used once with the empirical electron energy
0.00051099895069 GeV. Its certified energy interval gives an interval for
the GeV-per-native-energy unit. The cap-one and cap-three records use that
same unit; none is separately retuned. Energy and current intervals are
subtracted to bound finite changes between the complete carriers. Every
field occupation, boundary completion, native block and fermion sector
remains specified in each construction.

These comparisons do not prove an infinite-photon, spatial-continuum or
Lorentzian limit. They supply neither a local longitudinal gauge completion,
observed electron identity, source-selected fine-structure constant, nor
muon and tau branches. The declared global charge and empirical unit remain
inputs to the simulation-generated particle interpretation.
