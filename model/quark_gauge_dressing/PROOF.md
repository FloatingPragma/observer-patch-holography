# A native quark charge island with positive electric energy

## Full-source symmetry and its reference domain

Let S be the signed permutation on the twelve native quark coordinates
with ordered pairs `(0,12),(1,13),(4,14),(2,3),(5,6),(7,8)`. The first
coordinate maps to the second with sign plus; the second maps to the first
with sign minus. On each spatial site apply the canonical particle-hole
map

`c_(x,s,i) -> D_x epsilon_st S_ij c_(x,t,j)^dagger`

on quarks, leaving the six lepton modes and every scalar coordinate fixed.
This is a full finite CAR automorphism with a unitary Fock implementation.
Exact normal ordering exchanges the entire native up and down quantum
Yukawa polynomials and preserves the lepton polynomial and scalar energy.
Every local pair contains two factors D_x, so its fixed site coefficient
is preserved. The exact complete spatial identity

`h C+C conjugate(h)=0`, `C=diag(D_x epsilon)`

proves invariance of the quark hopping and texture on all occupations.
The full joined Hamiltonian therefore has this symmetry; it is not
restricted to the texture-free cell.

The transformation sends baryon number B to -B and hypercharge Y to Y-B.
It preserves the span of the eight color and three weak generators by
checked signs. Thus the native singlet domain with B=0 is invariant.
The full one-cell physical ground has B=0 by its owning ground theorem,
so its up-removal and down-addition measures are related by this map.
The spatial ground's baryon sector is not computed, and its invariance
under this transformation is not assumed. The operator symmetry alone
is not a proof that every selected spatial reference has identical
up and down spectra.

## The declared electric-field extension

Place the entire native source at vertex x of the committed thirteen-vertex
cone and a compensating quantum reference at r=1, with x distinct from r.
The 343 internal sites, when that source member is used, are internal
coordinates of this charge island. This construction does not impose a
separate Gauss law at those internal sites or supply covariant spatial
hopping between them. The added field is longitudinal and electrostatic;
no transverse photon oscillator is included.

Use the inherited edge counting metric. With edge gradient D, curl C and
`L=D^T D`, the exact rational Green matrix is

`G=(L+11^T)^(-1)-11^T/13^2`.

The code checks the literal incidence matrices, `CD=0`, `rank D=12`,
`LG=I-11^T/13` and `G1=0`. Define

`rho_x=delta_x-delta_r`, `f_x=D G rho_x`,
`C_x=||f_x||^2/2=rho_x^T G rho_x/2`.

Every vertex and every edge coefficient is retained. Gauss' equation is
`D^T f_x=rho_x`, and the field has zero curl. For any other field with the
same divergence, write it as `f_x+e_T` with `D^T e_T=0`. The two terms are
orthogonal, so

`||f_x+e_T||^2/2=C_x+||e_T||^2/2`.

This proves that C_x is the strictly positive minimum electric energy of
one compensated unit charge. It is derived from the graph and positive
field energy, rather than assigned separately to each quark species.

For a V-site native member, retain the finite flux register
`ell=k/6`, `|k|<=36V+3`. The bound contains every native matter charge,
whose absolute maximum is 6V, plus the global scalar cap's charge one half.
The compensator retains direct and dual quantum representation sectors
covering all source charges. Its intrinsic Hamiltonian is zero. All
nonabelian reference coordinates required by charged probes are retained.
The electric operators and Hamiltonian are

`E_e=f_x,e ell`, `H_E=(g/2) sum_e E_e^2`, with one common primitive `g=1`.

On the full kinematic carrier impose all thirteen operator constraints

`G_v=(D^T E)_v-delta_vx Y_source-delta_vr Y_reference=0`.

Since x differs from r, these constraints are exactly
`ell=Y_source=-Y_reference`. The isometry J that appends the corresponding
flux eigenstate obeys

`G_v J=0`, `J^dagger H_E J=g C_x Y_source^2 >=0`.

This is a finite gauge-fixed longitudinal representation. Independent
compact link variables, local gauging of the internal native geometry and
nonabelian gauge-field dynamics are not claimed.

## Same-clock charged records

The source law commutes with total hypercharge. The local full-Fock
commutator and all thirty charged-generator identities are checked
exactly; the spatial kinetic operator is the identity on native internal
coordinates and the scalar frequency is the identity on doublet indices.
The physical clock is therefore

`exp(-it (H_native tensor I+H_E))`, with hbar one.

Prepare the native singlet ground reference with its complete retained
purification, the compensator blank and flux zero. The ground-vector
selection is the same source-only singlet preparation; the electric term
vanishes on every branch of that preparation. Its original spectral
preparation probabilities, including non-ground branches, are unchanged.
The selected neutral reference is stationary and its energy is unchanged.
No assertion identifies it as the ground of the enlarged physical carrier:
compensated charged sectors require their own comparison.

A charged native probe carries the opposite quantum reference and shifts
the flux register by its exact source hypercharge q. Its range is physical
by the operator Gauss constraints. Since `[H_native,Y_source]=0`, every
state in that complete target has the additional energy `g C_x q^2`.
Consequently the full spectral measure translates as

`R_dressed(z)=R_native(z-g C_x q^2)`.

Every residue, dark target and outside-support outcome is unchanged. This
is an exact intertwining on the entire target, not a selected energy fit.
The one-cell example consumes the complete native quark measures with an
explicit transfer proof. It adds no source occupation cutoff and does not
copy their distinct four-unit reference energy: this extension declares
zero intrinsic compensator energy and derives its field energy explicitly.
The original success weight one half and failure weight one half are
retained; the failure source and all quantum references remain outputs.

For up-singlet removal, down-singlet removal and quark-doublet removal,
the charges are `2/3,-1/3,-1/6`. Thus all corresponding up and down lines
split by `g C_x/3>0`, while every original spectral weight remains attached.
The field term changes under the quark symmetry from `Y^2` to `(Y-B)^2`,
which explains the splitting on the charged target sectors. No per-flavor
Yukawa coefficient or measured mass enters.

This energy depends on the declared source-reference geometry. It is a
finite electrostatic dressing energy. It does not establish an intrinsic
rest mass, an observed generation count, a GeV scale or a relativistic
propagating particle limit. It is a concrete source extension that breaks
the finite up/down energy equality through positive gauge-field energy.
