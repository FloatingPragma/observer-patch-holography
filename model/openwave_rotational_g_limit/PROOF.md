# A controlled native magnetic coefficient and precession limit

The source uses the complete N=2 oscillator-shell carrier, native Clifford
matrices, Gaussian transverse field and actual angular generator J of
[the angular construction](../openwave_angular_torque/PROOF.md). Set the
oscillator energy coefficient to one, charge q=-1, field radius R=6/5,
mass m>=64, and self-field coefficient kappa=m^-6. The same energy functional
defines a stationary magnetic packet and a quasistatic precession experiment.
Quasistatic elimination of the field minimizes its instantaneous energy;
it is distinct from retaining propagating Maxwell mode momenta. The coupling
route controls the field error, and is not a physical fine-structure input.

Write D=beta m+A, where A=alpha dot (p+A_self+A_ext). Since {beta,A}=0,
the complete absolute kernel is sqrt(m²+A²). This identity includes both
Dirac signs. The code evaluates sqrt(m²+t)-m as t/(sqrt(m²+t)+m), avoiding
subtraction of large close numbers without changing the Hamiltonian.
All four native Dirac blocks remain present. The prepared beta-positive
block is an invariant sector of the positive kernel, not a discarded
negative-energy branch of D.

## Exact finite algebra

On the full forty-row block at zero field, T=A² satisfies T(T-5I/2)=0.
For E=partial_B A with a uniform z probe, Q=E² satisfies
Q(Q-I/4)(Q-I/2)=0. Both are positive squares. The executable verifies these
identities with exact radicals on every row, hence

    ||A||<=sqrt(5/2)<8/5,  ||E||<=1/sqrt(2)<4/5,
    ||partial_B T||<=sqrt(5)<3.

The oscillator vacuum with native spin up has
<T>=3/2, <partial_B T>=1 and <E²>=1/4. The spin-down derivative is -1.
These values are calculated from the minimally coupled native Dirac matrix.
They supply the L+2S magnetic coefficient; a Pauli term or a g value is not
inserted. Finite-shell boundary terms are retained, and no canonical
commutator identity is imposed on the highest shell.

Let H_0=sqrt(m²+T)+N. Its difference from m+N is a positive operator of norm
at most epsilon_0=5/(4m). In each beta-positive J_z=+/-1/2 sector, the
oscillator vacuum is the unique state below the number gap one. The full
lowest eigenvector lies within
eta_0=2 epsilon_0/(1-epsilon_0) of that vacuum, after a phase choice.
The complementary resolvent proves this bound and the low-level isolation;
it does not truncate the evolved source to the vacuum subspace.

The exact derivative X=partial_B sqrt(m²+T) solves KX+XK=partial_B T.
Comparing it with (partial_B T)/(2m), the positive Sylvester inverse gives

    ||X-(partial_B T)/(2m)|| <= ||T|| ||partial_B T||/(4m³).

The full stationary derivative therefore has, in the supplied-mass
normalization, |g_m-2|<=12 eta_0+15/(2m²). Here g_m=4m|mu_z| because the
actual axial angular momentum is exactly 1/2. Orbital and native-spin
contributions are measured with this same generator. Dividing instead by
the separate canonical-spin expectation would define a different quantity.

## Self-consistent static field

The Gaussian field derivative in a unit amplitude direction has operator
norm at most r=3/5. This follows by compressing the multiplication operator
alpha dot A_b and using sup_x |A_b(x)|<=R|b|/2. Its magnetic metric obeys
5<c_1<7. With a_B=8/5+(4/5)|B|, every normalized state's field-energy
minimizer lies in the ball

    |b| <= b_* = kappa r a_B/(5m-kappa r²).

Indeed its current norm is at most r(a_B+r|b|)/m. The field energy is
coercive; every minimizer satisfies the displayed inequality. In that ball
the field Hessian is bounded by
C=r²/m+r²(a_B+rb_*)²/m³. Thus kappa C/5<1 proves a unique minimizing field
for every normalized matter state, not merely convergence of one iteration.

The self-field changes the matter kernel by at most
delta=((a_B+rb_*)²-a_B²)/(2m). The axial lowest eigenline persists with gap
at least d=1-epsilon_0-2delta. Its stationary current map has derivative
at most C+2[r(a_B+rb_*)/m]²/d. Multiplication by kappa/5 is less than one.
This proves a unique self-consistent axial branch by contraction on its
invariant amplitude interval. Axial symmetry makes its transverse current
zero, so the scalar fixed point also solves the full three-component field
equation. The field momentum is zero, and total J_z is exactly 1/2.
At zero self-field the isolated low representation has J²=3/4: the source
commutes with J², and its nonzero overlap with the oscillator vacuum fixes
the representation. A fixed axial self-field can mix different j values.
Exact J_z=1/2 therefore does not assert that this finite-self-field matter
state has J²=3/4, or that a semiclassical field carries a quantum Casimir.
The complete static J² probabilities and the zero-field Casimir control
are reported separately.

The state change from the zero-field eigenline is at most 2delta/d.
The probe derivative norm is at most (4/5)(a_B+rb_*)/m, and its field
derivative at most (4/5)r/m+(4/5)r(a_B+rb_*)²/m³. These bound the static
moment shift. The rest energy is the actual empty-vacuum-relative matter
energy plus c_1|b|²/(2kappa), with additional shift at most
delta+7b_*²/(2kappa). The receipt evaluates all resulting static g bounds
as rational numbers. Numerical stationary eigenvectors are evaluations of
this proved branch, rather than an independent existence argument.

## Finite-time precession on the same member

Fix B=1/10 and duration T=64m. Prepare the normalized equal superposition of
the exact lowest beta-positive J_z=+/-1/2 eigenlines of the full zero-self-field
probe Hamiltonian. These are states of the full kernel; no two-state spin
Hamiltonian is inserted. Their relative phase makes the transverse angular
expectation positive. Angular selection rules make <J_x+iJ_y> a single
rotating term with frequency Delta E. All other source modes are retained.

Let epsilon_B=a_B²/(2m), rho_B=a_B^4/(8m³), and
s_B=epsilon_B²/(1-epsilon_B). A scalar square-root remainder plus the full
complementary Schur bound places each of the two exact low energies within
rho_B+s_B of its vacuum value. The vacuum values have difference B/m:
the B²/4 contribution is identical for both spins. Therefore

    |Delta E-B/m| <= 2(rho_B+s_B).

This is an analytic bound on exact eigenvalues; numerical root tolerances
do not enter it. The two vacuum-distance bounds give transverse radius at
least R_J=1/2-5[2epsilon_B/(1-epsilon_B)]>0. For the quasistatic source the
uniform kernel perturbation gives the full-state Duhamel bound T delta.
Using ||J_+||<=5, the angular expectation error is at most
5(2+Tdelta)Tdelta, and its phase error is bounded by twice that value/R_J.
This includes all unprepared modes and state-dependent field reaction.
The receipt combines these rational bounds with the same stationary rest
energy in g_T=2 E_rest |Delta arg<J_+>|/(BT).

At m=16384 the certified uniform result is |g_T-2|<0.004863. The static
bound also tends to zero. Numerical evaluations give static g about
1.99999999884 and frozen-field precession g about 2.00000000061; their
printed digits are not certified physical precision. Both static and
precession conventions use the actual total angular generator and the
same source member. Complete angular records and the supplied-mass
convention are reported alongside them.

This is a finite-regulator, low-momentum parameter theorem. It does not
derive the physical electron mass, alpha, anomalous magnetic moment, a
continuum limit of the cubic source, or equivalence to the interacting
Peter-Weyl Hamiltonian. Its quasistatic dynamics is semiclassical and does
not define a linear completely positive channel on arbitrary references.
