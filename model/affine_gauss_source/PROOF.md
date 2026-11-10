# One finite affine-Gauss source for static charge and transverse records

## Source, relation and finite carrier

This is a declared abelian construction on the repository's contractible
cone. It is distinct from the 119-state native link Hamiltonian. No native
vacuum, charge normalization, fermion representation or clock is transferred.
The charge unit, matter energy, reference energy and hbar are all one.
The physical meaning of those units is not calibrated.

`cone_whitney_bridge.load_geometry` and `cone_maps` supply the same oriented
13 vertices, 42 edges, 50 faces and 20 tetrahedra used by the Whitney cone.
Write D for the edge-by-vertex gradient and C for the face-by-edge curl.
They satisfy `CD=0`, `rank D=12`, `rank C=30` and `ker C=range D`.
The metric choices here are identity matrices on edges and faces. They
are counting metrics, not the geometric Whitney mass matrices. The
positive-normal-frame argument of `WhitneyConeModes.lean` applies with
these declared positive metrics. The executable independently checks the
full incidence ranks and polynomial spectral identities.

The existing discrete Green/Thomson result supplies the longitudinal
split. Here its dimension-independent rational inverse routine is applied
to the cone Laplacian `L=D^T D`; its twelve-port wrapper is not applied to
the thirteen-vertex cone. Define

```
G=(L+11^T)^(-1)-11^T/13^2,
LG=I-11^T/13,     G1=0,
P_T=I-D G D^T,   K=C^T C.
```

The code verifies these identities exactly, as well as `P_T^2=P_T` and
`KP_T=K`. The constant scalar mode is removed by the neutral-charge
constraint. There are no transverse harmonic zero modes on this cone.
The eight distinct positive eigenvalues of K and their multiplicities are

```
4-sqrt(5): 3,  3: 5,  6-sqrt(5): 3,  4: 4,
6: 4,  4+sqrt(5): 3,  7: 5,  6+sqrt(5): 3.
```

The remaining twelve eigenvalues are zero gauge directions. Exact
polynomial projectors in K give every displayed spatial residue; their
symmetric idempotence proves positivity. Write `omega_j=sqrt(lambda_j)`
on the 30-dimensional transverse space. In particular `1<omega_j<3`.

The transverse quantum carrier is finite:

```
F_3 = direct sum_(n=0)^3 Sym^n(C^30),
dim F_3 = binomial(33,3)=5456.
```

It uses total occupation at most three, rather than thirty independent
cutoffs. It is invariant under any change of frame within a degenerate
frequency space. Thus choosing an algebraic orthonormal normal frame
changes coordinates, not the finite source. Such a frame exists by the
finite real spectral theorem and can be constructed from the checked
projectors by Gram-Schmidt.

Matter consists of thirteen declared spinless fermionic modes with the
full `2^13` CAR carrier. The static compensator is at vertex r=1 and has
fourteen states `R=0,...,13`, charge `-R`, energy `R` and parity `(-1)^R`.
There are twelve longitudinal bits ell_y for y unequal to r. Define

```
L_charge=sum_(y!=r) ell_y (delta_y-delta_r),
E_L=D G L_charge,
G_v=(D^T E)_v-n_v+R delta_(v,r).
```

With the transverse coordinates below, `D^T E=L_charge`. The physical
space is exactly the common zero kernel of all thirteen commuting G_v:
`ell_y=n_y` for y unequal to r, and `R=N=sum_y n_y`. It has dimension
`2^13*5456=44695552`. The full declared kinematic carrier has dimension
`2^13*14*2^12*5456`. All source occupations, reference charges and cutoff
states are retained. The gauge-fixed longitudinal representation is a
finite design choice; no compact independent local-link realization is
claimed. These constraints are exact operator identities on the physical
space, not constraints imposed only on expectation values.

## The cutoff obstruction and its exact boundary completion

Let P project total boson occupation to at most three and set
`a_j=P b_j P`. These are finite matrices. Define

```
q_j=(a_j+a_j^dag)/sqrt(2 omega_j),
p_j=i sqrt(omega_j/2)(a_j^dag-a_j),
A=sum_j v_j q_j,
E=E_L-sum_j v_j p_j,       B=C A.
```

The minus sign gives the convention `dot B=-C E`. The finite a operators
do not satisfy the unrestricted canonical commutator. For occupation n,
`a_j^dag a_j=n_j`, while `a_j a_j^dag=n_j+1` below total cutoff and zero
at total occupation three. Consequently

```
(||E||^2+||B||^2)/2
 = ||E_L||^2/2 + Tr(omega)/2 + dGamma(omega)
   -(1/2) P_top [dGamma(omega)+Tr(omega)].
```

Longitudinal-transverse orthogonality removes all cross terms. Dropping
the boundary term changes the ground. In the bare quadratic expression,
a three-quantum state in a lowest-frequency mode has energy
`3 omega_min/2`, whereas the bare vacuum energy is `Tr(omega)/2`, at
least `15 omega_min`. Thus this uncorrected finite cutoff does not supply
the intended vacuum. The one-mode two-level example also gives a constant
bare quadratic Hamiltonian.

The declared finite Hamiltonian includes the positive boundary completion:

```
H=N+R+(||E||^2+||B||^2)/2
    +(1/2)P_top[dGamma(omega)+Tr(omega)]
 =N+R+||E_L||^2/2+Tr(omega)/2+dGamma(omega).
```

This is one finite electric, magnetic and boundary Hamiltonian. Its clock
is `exp(-itH)`. The boundary term is a disclosed spectral service, not a
local plaquette interaction. The executable checks the coefficient
identity on all 5456 photon occupation states and all thirty modes.

Every term commutes with all G_v. The ground on both the full kinematic
carrier and the physical carrier is unique: all matter, reference and
longitudinal registers are zero and every photon mode is in its vacuum.
All other terms in the final displayed expression are nonnegative, and
G is positive on nonzero neutral loads. Preparing these declared blank
registers has success exactly one. This assumes initialized finite blanks;
no autonomous cooling or native implementation is asserted.

Exhausting all 4096 longitudinal bit patterns gives the least nonzero
longitudinal energy 26/217. A lone unmatched matter or reference occupation
costs one, and every transverse quantum costs more than one. Thus the
full kinematic gap is 26/217. On the physical space every matter excitation
includes its reference and costs at least two; the actual physical gap is
exactly `sqrt(4-sqrt(5))`. These distinct gaps refer to the stated
kinematic and Gauss-physical domains.

The actual finite equations `[H,q_j]=-i p_j` and
`[H,p_j]=i omega_j^2 q_j` hold on the entire cutoff carrier. Hence
`dot B=-C E` and `dot E=C^T B` hold for the static-charge source, with
`D^T E=L_charge` constant. They require the completed H. A finite canonical
commutator is not inferred from these equations. Every vacuum field word
with at most six linear field insertions agrees exactly with the
unrestricted Gaussian comparison: a nonzero return path with six ladder
operations never exceeds total occupation three. No error claim for
unrestricted Weyl kicks or arbitrary highly occupied states follows.
For a frequency-one mode the finite vacuum eighth moment is 81/16,
whereas the unrestricted Gaussian eighth moment is 105/16. This exact
counterexample excludes an all-moment or unrestricted Weyl transfer.

## Static compensated charge and its energy

The physical embedding V maps a matter occupation n and photon vector to
the matched registers `R=N`, `ell_y=n_y`. A compensated creation operator is
`D_y=V c_y^dag V^dag` on the physical space. Its explicit full extension
raises the reference charge and, for y unequal to r, the y longitudinal
bit, together with the actual fermion creation. Its action has the CAR
sign `(-1)^(sum_(i<y)n_i)`; invalid boundary moves are zero. All Gauss
commutators vanish on the complete carrier, since its changes satisfy

```
delta L_charge=delta_y-delta_r,
delta n=delta_y,      delta R=1.
```

For y=r the longitudinal change is zero. The reference parity changes
with the odd matter operation, making the joint operator even. The
source operator is an odd CAR creation; physical spin and spin-statistics
identifications are absent.

On the actual ground, D_y creates the exact charged-sector ground at
fixed occupation y with the same transverse vacuum. Its electric mean is
`D G (delta_y-delta_r)` and its excitation energy is

```
E_y-E_vac=1+1+(1/2)(delta_y-delta_r)^T G (delta_y-delta_r).
```

The first one is matter energy and the second is retained reference
energy. The finite Green electric costs are

| Source position relative to boundary reference | Electric cost |
| --- | --- |
| Same vertex | 0 |
| Apex | 26/217 |
| Adjacent boundary vertex | 33/217 |
| Boundary surface distance two | 40/217 |
| Opposite boundary vertex | 6/31 |

The boundary surface distances refer to the inherited twelve-port
icosahedral graph; the cone has an additional path through the apex.
These are exact finite Green energies. They do not prove an inverse-distance
continuum law, a thermodynamic Coulomb phase or an itinerant charged
threshold. The reference-sector global ground places the matter at r,
where its electric cost is zero. Other rows are constrained static charge
preparations, not degenerate global grounds of the whole one-charge sector.

Use the contraction `D_y/2` as a charge pump. Its complete two-flag
instrument has success probability 1/4 and failure probability 3/4 on the
neutral ground. The success state has the displayed exact energy; the
failure state stays the ground. Both branches satisfy every Gauss
constraint. The pump is an explicitly supplied finite intervention and
can supply work. Its reference energy is included in the measured
transition frequency.

## Complete transverse and charged-field records

Let P_lambda be a checked projector of K. On the actual neutral ground,
the full finite spectral residues for a transverse electric edge e and
magnetic face f are respectively

```
r_E(lambda,e)=sqrt(lambda) (P_lambda)_(e,e)/2,
r_B(lambda,f)=(C P_lambda C^T)_(f,f)/(2 sqrt(lambda)).
```

They are nonnegative and supported at the actual excitation frequency
`sqrt(lambda)`. Their time transforms are exactly
`sum_lambda r(lambda) exp(-i sqrt(lambda)t)` for all real t on this H.
All thirty modes, including degeneracies, enter the polynomial projectors.
No floating-point eigenvector establishes a residue.

For the chosen total cutoff, `||a_j||<=sqrt(3)`. Cauchy-Schwarz over thirty
modes, `1<omega_j<3`, `||P_T e||<=1` and `||C_f||^2=3` give both field
operator norms strictly below `sqrt(540)<24`. Dividing the field by 24
therefore defines a contraction on the whole finite carrier. Its pump
success followed by energy measurement has probabilities `r/24^2`, with
pump failure `1-sum r/24^2` and the energy complement retained at zero
on this prepared input. The finite Kraus definitions retain that
complement on every other input.

A direct combined experiment first creates the compensated charge, then
pumps the *full* electric edge field E_e and measures H. Every longitudinal
bit contributes a dipole field of norm less than one, so
`||E_L,e||<12`, and `||E_e||<36<48` on the entire kinematic carrier.
Use `E_e/48` for this second pump. Let l_e be the static electric mean
at the selected y. Its elastic field residue is `l_e^2`; its transverse
residues are the same r_E as above, at total energies
`E_y-E_vac+sqrt(lambda)`. The twelve retained records are

```
charge-pump failure:                         3/4,
charge success, field-pump failure:          (1/4)[1-(l_e^2+sum r_E)/48^2],
elastic charged energy:                     l_e^2/(4*48^2),
each of eight inelastic charged energies:   r_E(lambda)/(4*48^2),
other energy:                               0.
```

They sum exactly to one. Every positive transverse residue gives a
strictly positive unconditional combined charged-field record. The
correlation in the same charged preparation is
`l_e^2+sum r_E exp(-i sqrt(lambda)t)`. Thus static charge, its compensator,
its electric field and the transverse spectral support use one declared
finite state and Hamiltonian clock.

For arbitrary inputs use the failure operators `sqrt(I-B^dag B)`, each
pump B, and the full spectral projectors of H for the stated energy bins
and their complement. Their products give the complete sequential Kraus
instrument, with summed squared effects I. The bins may contain states
from other matter sectors; this does not change the displayed probabilities
on the specified input. All Kraus operators commute with every G_v. A
retained twelve-state pointer suffices for the complete joint record; its
clock is the identity. Other control storage and work are separate declared
services. The full kinematic source, charge reference and pointer dimension
is `12*2^13*14*2^12*5456`.

This finite result supplies no motion of matter, scattering, compact-link
realization, continuum photon theorem or observed charge/spin/mass
identification. Adding matter hopping changes H and requires its own
actual-ground and response analysis.
