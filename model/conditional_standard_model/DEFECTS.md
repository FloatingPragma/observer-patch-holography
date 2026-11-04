# Defect classes and massive particle families

The conditional Standard Model action specifies its fields and mass matrices.
A family-to-defect correspondence is a separate premise: it requires a map
from quantized defect states to those fields that preserves gauge currents,
spin, statistics, Hamiltonian, observables and the chosen clock. The catalogue
in [defects.py](defects.py) records proven invariants and the domains in which
they are preserved. Every Standard Model defect assignment is null.

## The supplied onsite mass family

The [native mass classification](../native_mass_channels/README.md) exhausts
Hermitian onsite matrices anticommuting with the fixed two-spin Weyl hopping:

```
B_n=m0*(-1)^(n1+n2+n3)*I + sum_j m_j*(-1)^n_j*sigma_j.
```

For uniform coefficients put r=|m|. The gap is `abs(abs(m0)-r)`.
There are three connected gapped components. The cones `m0>r` and `m0<-r`
contract to `(1,0)` and `(-1,0)`. The middle component `|m0|<r` retracts to
`m0=0, |m|_1=1`, an octahedral sphere homeomorphic to S2. Explicitly,

```
scalar cone: (m0,m) -> ((1-t)*m0+t*sign(m0), (1-t)*m),
middle:      (m0,m) -> ((1-t)*m0, (1-t+t/|m|_1)*m),  0<=t<=1.
```

The cone inequality is preserved by convexity. For the middle component,
`|m(t)|=(1-t+t/|m|_1)*r >= (1-t)*r > |m0(t)|` for t<1,
and the endpoint has positive radius. The endpoints are fixed by these maps.
The code recomputes the complete 96-equation, 32-variable edge nullspace of
dimension four and exact rational retraction witnesses.

Within this family, walls separating distinct components require leaving
the gapped coefficient space. Each component has trivial fundamental group,
so this classification supplies no nontrivial line-defect winding. The middle
component has `pi2(S2)=Z`: one point-defect type with infinitely many integer
charges. The three components are not three generations. Opposite spin-mass
directions lie in the same component; the executable rational rotation joins
them with a constant uniform gap. This path makes no claim about the full
spectrum of a spatially varying interface.

The degree witness uses nonvanishing continuous piecewise-linear interpolation
on an octahedral sphere. The hedgehog has degree +1 and its negative has degree
-1. Coordinate perturbations below delta<1/3 keep every interpolated face
separated from zero by at least `1-3*delta`, giving an explicit protected
homotopy class. Bare finite-lattice vertex labels need an interpolation and
an admissible configuration space to define this invariant.

## Protected indices and ordered-spin textures

The [continuum radial operator](../topological_defect_spectrum/PROOF.md) has
internal grading `Gamma=i*gamma1*gamma2*gamma3*beta1*beta2*beta3`. Its explicit
bounded hedgehog and antihedgehog each have two normalizable zero modes, with
chiral indices +2 and -2. A continuous bounded chiral perturbation path
preserves this Fredholm index while zero stays outside the essential spectrum.
The calculation proves no index formula for arbitrary integer-degree profiles.

The [original radial lattice operator](../bounded_defect_index/PROOF.md) uses
the nonlocal spatial-inversion grading `R psi(n)=psi(-n)`. Its index is +2
for both orientations. Each nonzero inversion pair contributes an invertible
off-diagonal block, and the origin contributes two even zero rows. The
essential gap permits the homotopy from the atomic limit to the full hopping.
Compact self-adjoint inversion-odd perturbations preserve this index, as do
bounded inversion-odd perturbations of norm below m. An inversion-breaking
scalar shift moves the zero levels despite an essential gap. This index is
distinct from the sphere degree and the continuum internal chiral index.

The [strong-mass lattice theorem](../radial_defect_gap/PROOF.md),
`m>sqrt(3)*v/a`, isolates exactly two zero modes and gives a uniform nonzero
gap and quantitative localization. For a=v=1, m=64 and rational hopping bound
k=2, the certified nonzero gap is at least 62. Fixed positive physical m and v
with a tending to zero leave this strong-mass domain. None of these operator
indices identifies electric charge, Lorentz spin, defect formation energy or
native-repair preservation. A zero eigenvalue of a pinned operator is not a
massless propagating particle identification.

The [ordered-spin ground orbit](../spin_texture_energy/PROOF.md) is
`SO(3)/SO(2)=S2`. Point-boundary maps have degree `pi2(S2)=Z`; compact textures
with fixed direction at infinity have Hopf invariant `pi3(S2)=Z`. The explicit
compact SU(2) lift has absolute degree
`-(2/pi)*integral_pi^0 sin(f)^2 df=1`, and hence absolute Hopf invariant one.
The executable catalogue recomputes this identity and the exact rotation-energy
identities on the actual selected quantum-ground orbit.

For this compact family, `DeltaE(R(x/L))=L*E_cont+O(1)` with `E_cont>0` in the
fixed exchange clock. An isolated radial hedgehog in the same orbit obeys
`liminf DeltaE/L>=4*pi*J*kappa>0`. These macroscopic families do not provide a
finite-mass particle. The statements concern their energy expectations;
they neither exclude every finite-energy spectral component nor classify all
quantum defects. An unconstrained finite-lattice rotation assignment deforms
to the identity, so a microscopic conserved topological operator does not
follow from the continuum orbit alone.

## The canonical electroweak Higgs orbit

For the assumed single canonical complex doublet and positive symmetry-breaking
parameters, the minimum constraint is `H^dagger H=vEW^2/2`. Its homogeneous
orbit is S3. Write H=(z1,z2), w=vEW/sqrt(2). The matrix

```
U=(1/w)*[[conj(z2), z1], [-conj(z1), z2]]
```

has determinant one, is unitary on this constraint, and maps `(0,w)` to H.
This gives the explicit transitive SU(2) orbit used by the code. In this
canonical setting, `pi1(S3)=pi2(S3)=0`; the homogeneous Higgs vacuum orbit
does not require stable strings or point monopoles. A central quotient acting
trivially on the doublet does not change that orbit. Other global gauge groups
and stabilizer embeddings require their own analysis. The distinction is
discussed by [Gripaios and Randal-Williams, Topology of the Electroweak Vacua](https://arxiv.org/abs/1610.05623).
This statement does not exclude embedded strings, singular monopole sectors,
gauge saddles or composite defects. The integer `pi3(S3)` alone proves no
stable massive particle.

## The conditional family map

The candidate target channels are up-type quarks, down-type quarks and charged
leptons, together with W+, W-, Z and the Higgs field. A neutrino mass extension
is a separate conditional target requiring its own action. Photon, gluon and
minimal massless-neutrino fields are excluded. Any additional zero singular
values of a supplied charged Yukawa matrix are also excluded individually.
Degenerate or ordered Yukawa eigenvalues carry no assigned topological charge.

A proposed assignment must supply the invariant's physical configuration
space, its dynamical preservation or quantified lifetime, gauge representation,
spin and statistics, finite formation energy, motion, and the field/state
intertwiners. It must also satisfy the interaction constraints in
[defect_constraints.py](defect_constraints.py). For an exactly conserved additive
charge, a decaying parent transfers its total charge to its decay products.
Exact topology cannot independently make every unstable massive species a
stable isolated charged sector. Composite, nonadditive and approximate
topological mechanisms have different hypotheses.

This catalogue classifies the stated onsite-mass family and records specific
ordered-spin and operator invariants. It does not classify all possible OPH
defects or establish a Standard Model particle identification.
