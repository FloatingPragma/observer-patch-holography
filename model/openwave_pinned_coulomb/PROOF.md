# Exact pinned charge-pair minima of the quadratic-field source

## Hamiltonian and complete pinning instrument

Use the finite declared Hamiltonian of
[the native-current field construction](../openwave_current_coulomb/PROOF.md),

    H' = sum_x h_x + lambda T + 4 N_ref + (kappa/2) rho^Y G rho^Y,
    lambda=2^-34, kappa>0, G=(D*D)^+.

The local Hamiltonian h includes the actual quantum scalar and Yukawa
operators. Every cell retains all 30 CAR modes and five scalar states;
every link retains all 119 Peter-Weyl states. The original nonblank-link
penalty is absent. The complete reference Fock space is retained.
Physical states satisfy all twelve native local Gauss constraints, and
therefore total source-plus-reference hypercharge is zero. The periodic
Green operator removes the constant mode. It does not introduce a
neutralizing external charge or change the native charge normalization.

The accepted reference sector contains one native e^c reference particle,
with hypercharge +1 and trivial weak/colour representation. Let M_0 be
the complete local matter-plus-scalar gauge-singlet space. Let M_chi be
the complete local subspace with hypercharge -1 and trivial weak/colour
representation. It includes every allowed occupation, scalar state and
multiplicity, rather than one visible charged-energy block. Both spaces
are invariant under h. The actual local ground energy in M_0 is E_0.
Write e_chi for the minimum of h on M_chi. No numerical value or
identification of e_chi with a selected visible pole is assumed.

For each cell y, Pi_y imposes M_chi at y and M_0 at every other cell,
with the accepted reference and full joint physical link space. The
Pi_y are mutually orthogonal. Together with I-sum_y Pi_y they give a
complete instrument. The complement includes reference failures,
screening clouds, other local charge/weak/colour sectors and all states
with a different pinned-density pattern. These outcomes are not discarded
from the source or identified as failures of a Coulomb theorem outside
its domain. This instrument is more specific than an unscreened charge
projector that permits additional nonabelian arrival irreps.

## Exact hopping compression

Every native fermion hypercharge is nonzero. Each transport block has
one constant hypercharge, and every hopping word has distinct endpoints.
Thus a hopping word changes the total local matter-plus-scalar density
by opposite nonzero charges at its two endpoints. In particular,

    P_rho T P_rho = 0

for every fixed local hypercharge pattern, including every Pi_y. This
argument retains the scalar, reference and all quantum link factors; it
does not assume a one-particle truncation or a particular flux path.
The executable checks the complete native hypercharge list and transport
blocks, the actual finite-volume hopping paths, and a full physical
CAR/scalar/link column of H'.

## Tensor lower bound and its saturation

All gauge generators vanish on M_0. On M_chi the only nonzero generator
is the scalar Y=-1. Hence the Gauss constraints on Pi_y act on the flux
and reference factors without resolving the local gauge multiplicity
spaces. The physical pinned carrier has the form

    M_chi,y tensor (tensor_{x != y} M_0,x) tensor F_y,

with the reference spin multiplicity included in F_y. The finite native
flux space F_y is nonempty: a simple e^c string joins the reference to y,
and the actual endpoint current matrices satisfy all twelve Gauss
constraints. At y=0 the blank flux is allowed. The same string works
for every vector of M_chi, because its gauge representation is the same
one-dimensional character. The actual local spectral embedding and an
explicit two-L-occupation singlet supply independent nonemptiness
witnesses for M_chi.

On this complete pinned carrier the Green term is the scalar
kappa C_y, with C_y=G_00-G_0y. The compressed operator is

    Pi_y H' Pi_y = h_chi,y + sum_{x != y} h_0,x + 4 + kappa C_y,

tensored with the identity on F_y. Hopping contributes zero. The local
operator lower bounds therefore give

    inf spectrum(Pi_y H' Pi_y) >= (V-1) E_0 + e_chi + 4 + kappa C_y.

A local ground vector in M_chi, local ground vectors in M_0 and the
physical e^c string attain that lower bound. The exact constrained
minimum is consequently

    E_min(y) = (V-1) E_0 + e_chi + 4 + kappa(G_00-G_0y),
    E_min(y)-E_min(0) = kappa(G_00-G_0y).

The unknown complete charged minimum cancels exactly. Computing or
selecting its value is unnecessary for this separation-energy observable.
Flux loops and their multiplicities are retained. They contribute no
energy in this declared Hamiltonian and cannot lower the sum of the
complete local minima. The degeneracy does not obstruct the minimum
energy difference. It would matter for a vacuum-selection or autonomous
transport problem, which is not identified with this static calculation.

The equality holds at the declared nonzero hopping because it concerns
expectation minima on compressed subspaces. The Pi_y need not be
invariant under full-H' evolution. Pinning is an explicit apparatus
constraint, not a claim that free particles remain pinned without work.
The original linear-link Hamiltonian and its ground, clock and
preparation probabilities do not transfer to H'.

## Computed potential and static force

The Green operator is computed from the incidence Laplacian, with no
sampled 1/r potential inserted. Exact side-two and side-four checks give
nearest-cell minimum differences 7/48 and 21/128 at kappa=1. The finite
static force is the negative centered difference of the exact minimum
energy. Side lengths 64, 96 and 128 test its large-box form. It is
attractive and approaches the comparison -kappa/(4 pi r^2); potential
differences approach -kappa/(4 pi r). The native unit current fixes the
charge factors. The counted distance, coupling kappa and conversion to
physical length/energy are separate declared inputs. Keeping all native
current and Gauss identities while changing kappa from one to two
doubles the energy and force. These identities therefore derive no
physical fine-structure constant.

The floating Green profiles are numerical finite-volume comparisons,
without outward-rounded continuum error enclosures. The exact theorem
is the constrained-minimum identity on every specified finite source.
It is neither an unrestricted screened-ground force nor a dynamical
scattering calculation, and it supplies no transverse photon Hamiltonian.
