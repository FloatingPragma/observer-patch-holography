# Autonomous signed Dirac flow on the retained field carrier

This is a declared amendment, `native-local12-autonomous-signed-dirac-noncompact-v5-proposed`.
It keeps the nativefield25 one-particle carrier, Gauss reference, routed
integer currents, four spin coordinates and all 193 physical field modes.
It replaces the old ordered full-shift macro by a continuous joint generator.
The old source files and their frozen results are unchanged. This contract
does not claim that either source is a reducing sector of the other.

## Carrier, routes and common clock

The Gauss-reduced carrier is C^8_position tensor C^4_spin tensor
L2(R^193). Its electric field in position sector x is
E=P_x+F p, where P_x=B^T G(G_reference-n_x) and F is the frozen coexact
plus harmonic orthonormal frame. The 190 coexact coordinates have
lambda_j>0 and omega_j=sqrt(a b lambda_j), a=1/2, b=3. The three harmonic
coordinates retain their full noncompact free-particle Hilbert space.

For each actual positive axis route x to y, write J_xa for its integer edge
current and j_xa=F^T J_xa. Define the spin-independent unitary

    T_a |x,s> W(q,p)phi_u
      = |y,s> exp(i j_xa.q/2) W(q,p+j_xa)phi_u.

The positive route is a permutation of the eight positions. The owning
negative route at y is exactly the reversed path with current -J_xa.
Thus the negative operator is T_a^dagger, including its Weyl phase.
All 48 signed paths satisfy the literal integer inverse and Gauss equations;
the implementation also checks nonzero displaced packets at shear time0.73.
In unreduced notation the simultaneous matter/field shift preserves the
operator Gauss law, because B J_xa=delta_x-delta_y.

Keep the existing physical clock t=tau0 h u, with tau0=1, h=1/8 and mass m=1.
The autonomous generator in that clock is

    H_phys=H_field/h + m beta + sum_a alpha_a Pi_a,
    Pi_a=(T_a^dagger-T_a)/(2 i h),
    H_field=a(p_co^2+p_h^2)/2+b sum_j lambda_j q_j^2/2+H_Coulomb.

The equivalent cell-time generator is
H_cell=h H_phys=H_field+h m beta+sum_a alpha_a(T_a^dagger-T_a)/(2i).
A physical time interval dt therefore uses field interval dt/h and mass
angle m dt. The Coulomb term and coexact vacuum scalar are scaled together
with the field. This retains the optical gap sqrt(3)/h from ROW19 and the
earlier common clock; an unscaled H_field would change their relationship.
The retained `UNSCALED_PROTOTYPE_*` files are a rejected development control,
not the final convention or a source claim.

On amplitudes psi_x=exp(i p x), Tplus has eigenvalue exp(-i p h). The
adjoint-minus-forward definition gives Pi=sin(p h)/h. The reverse order
would give the negative symbol. The five-site scalar check fixes this sign
convention; it is an algebraic control, not an executed enlarged field model.
On the actual two-cell-per-axis torus every free corner has sin(k)=0.
All eight corners and both twofold-degenerate energies +/-m are retained.
No larger-volume relativistic dispersion is inferred from that aliased case.

## Joint generator and energy control

H0=H_field/h is self-adjoint on the domain of the 190 oscillator plus
three free kinetic operators. Its sector Coulomb potential is a bounded
nonnegative eight-position diagonal. The Schwartz core is sufficient for
the displayed route and commutator identities. Each T_a is unitary and
commutes with its spin alpha_a. Hence Pi_a is bounded self-adjoint with
norm at most1/h. The complete bounded matter term V obeys

    ||V|| <= M=|m|+3/h=25.

The bounded-perturbation theorem defines H_phys=H0+V self-adjoint on D(H0).
Its lower bound is E0/h-M, where E0=sum_j omega_j/2. All form-domain inputs
evolve under one strongly continuous unitary group with conserved expectation
of this same H_phys. If Erel=<H_phys>-E0/h initially, then at every time

    0 <= <H_field>-E0 <= h(Erel+M),
    <sum_j N_j> <= h(Erel+M)/min_j omega_j,
    <p_h^2> <= 2 h(Erel+M)/a.

These bounds concern this new autonomous source. The old macro changes the
expectation of this candidate joint H on the explicitly checked states, so
it is not silently identified with this flow. That comparison is separate
from the earlier observation that its field energy changes. There is no
inference that every possible Floquet logarithm equals the physical sum.
The bounds do not control the noncompact harmonic coordinate q_h, imply
spatial binding, or construct a stationary charged state.

Lie products

    (exp(-it H0/n) exp(-it V/n))^n

converge strongly to exp(-it H_phys). Splitting V into its mass and three
bounded axis terms has a separate operator-norm error at most
t^2(3|m|/h+3/h^2)/n, using the pair commutator bound and unitary telescoping.
The outer H0/V convergence is only strong here: its unbounded commutator
does not acquire a uniform operator-norm bound by declaring it small.
No numerical cutoff, native pulse iterator, all-input reference error,
preparation service, or finite-time native Trotter rate is established by
that strong convergence statement. Those require the separate implementation
and cutoff contract being produced by the charge/spin lane.

## Current and stationary means of this source

Let T_xa denote the positive route restricted to initial position x and
R_a=sum_x J_xa T_xa. The physical current operator is

    Jdot = sum_a alpha_a(R_a+R_a^dagger)/(2h).

The unreduced physical electric field satisfies [E,T_xa]=J_xa T_xa;
the same full current follows after adding the position-dependent
longitudinal field to the transverse canonical momenta. Its expectation
uses Re<alpha_a T_xa>/h, which can be negative and is not a branch
probability. The implementation evaluates these coherent matrix elements
and independently checks B<Jdot>+d<n>/dt=0 using 2 Im<psi|P_x V|psi>.
It keeps every route, spin and field cross term.

The coexact Ehrenfest equations are

    dq/dt=a p/h,  dp/dt=-b Lambda q/h+F_co^T<Jdot>,

with dq_h/dt=a p_h/h and dp_h/dt=F_h^T<Jdot>. Necessary stationary
means satisfy p=0, b Lambda q=h F_co^T<Jdot>, p_h=0 and
F_h^T<Jdot>=0. This continuous relation is distinct from the kicked
macro's p=-j/2 and cot(omega/2) kernel. It is a self-consistency condition,
not a construction of a solution or an intrinsic magnetic moment.

The three literal positive two-step winding cycles at root0 have harmonic
projection norm squared1/3 and rank3. A cell-dependent displacement cannot
remove all harmonic route factors: around a closed cycle such dressing
cancels, leaving its nontrivial Weyl holonomy. Consequently no common
normalized free harmonic reference has been proved to factor from this
charged Hamiltonian. A compact quotient or Bloch-fiber construction needs
its own precise sector and state transfer. An exact selected fiber of a
noncompact direct integral is not automatically a normalizable input.

## Signed bands, reference and existing results

All four spin components remain. No positive-band projection is applied.
For D=sum_a alpha_a Pi_a, its square retains the magnetic commutators
from noncommuting covariant momenta; no spin term is installed separately.
This one-particle signed operator is bounded below after adding the positive
field, but that is not a filled Dirac-sea vacuum, a many-body Fock stability
theorem, a positive-particle dictionary, or a pair-production construction.
An infinite-volume or continuum limit with its lower band would require
new domain, reference and energy analysis. The fixed neutralizing Gauss
reference is unchanged and is not a dynamically supplied antiparticle.

The existing RER `openwave_dirac_maxwell_refinement` proof supplies an
autonomous compact positive-Pauli-root architecture with full-flux bounds.
Its source and reference differ, so neither its ground bounds nor its
native preparation are transferred here. The existing
`openwave_dirac_field_connection` proof establishes U*EU=E+A and
alpha(E+A)^2 under field-dependent signed-band diagonalization. This
amendment avoids that hidden omission by retaining the original signed
spin basis. Any later band projection or change of basis must transform
the full field and current, including the connection and its square.

The ROW19 compact optical direction is orthogonal to every actual route.
Therefore the same exact dark oscillator still factors from this new
Hamiltonian, now with its retained physical factor1/h. This is a reuse of
the already proved route orthogonality, not a stationary charged excitation.
