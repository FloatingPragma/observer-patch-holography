# Lattice degree, source energy and finite deformations

## The actual ordered source and its conditional state

Use the fixed-spacing spin Hamiltonian and the same selected ordered
quantum ground as [spin_texture_energy](../spin_texture_energy/PROOF.md):

```
H=J sum_<xy>(S_x.S_y-1/4),                  J>0,
omega(S_x)=(-1)^x m_order e1,               m_order>1/1000,
C_ab=omega(S_x^a S_y^b)=diag(a,b,b),
kappa=-(a+b)>=m_order^2/2>1/2000000,         b<=0.
```

The ground-selection rule is an input inherited from the source theorem.
It is not derived from native repair. The spin Hamiltonian is the stated
strong-coupling limit of the complete centered contact source. The source
clock, exchange J and lattice spacing are retained, without identifying
spin directions with physical Lorentz spin or a supplied Dirac mass.

Choose a finite assignment of onsite spin rotations R_x in SO(3), equal
to identity elsewhere, and apply the corresponding local unitaries to
omega. Put n_x=R_x e1. The exact expectation above the uniform reference is

```
Delta E=J sum_<xy> Tr[C(R_x^T R_y-I)]
       =J sum_<xy>{kappa(1-n_x.n_y)
          +(-b)(1+R_11-R_22-R_33)},                  (1)
```

where R=R_x^T R_y in the summand. The second factor equals
2(1-cos theta)u_1^2 for a relative rotation about u. Every bond's excess is
nonnegative, and

```
Delta E >=(J kappa/2) sum_<xy>|n_x-n_y|^2.             (2)
```

The uniform reference has bond energy J(a+2b-1/4). The actual spin-1/2
bond has spectrum {-J,0}; consequently each excess in(1) is at most J.
This upper bound uses the source bond operator, not an assumed classical
product state. All statements below concern this explicit local-unitary
orbit and its nonnegative bond-energy sum. They do not minimize over all
quantum states with the same one-point directions.

## Cubic shells and an exactly computable integer

For an integer center c and positive integer radius r, triangulate the
boundary of c+[-r,r]^3. Divide each unit square with corners
p00,p10,p01,p11 into triangles (p00,p10,p01) and (p11,p01,p10), with outward
orientation. Both edges incident to the first vertex of each triangle
are actual nearest-neighbor source bonds. The diagonal is used only for
interpolation. Every surface source edge occurs in exactly two square
boundaries, including edges along the twelve cube ridges.

For each triangle with directions u,v,w impose the finite admissibility
condition that all three pairwise dot products are nonnegative. Normalized
barycentric interpolation

```
(lambda0*u+lambda1*v+lambda2*w)/|lambda0*u+lambda1*v+lambda2*w|
```

is then defined for every lambda_i>=0 with sum1. The triangles glue into
a continuous map from the shell sphere to S2. Its oriented area divided
by4pi is an integer degree q. This is a derived integer for the specified
interpolation, not a microscopic superselection operator.

The code computes degree without eigensolvers, floating angles or
rounding to an integer. A nonzero rational direction vector w_x may be
used before normalization. Choose a rational regular target y lying on
no plane spanned by a triangle edge. For each nondegenerate triangle,
solve [w_u,w_v,w_w]lambda=y over the rationals. If all lambda_i>0, count
the sign of its determinant. The sum is the degree: each such preimage
has that orientation sign, and all other triangles miss the target.
Degenerate triangles map to arcs and have no preimage of this regular
target. The implementation rejects a target on any edge plane, so no
boundary convention silently changes the count. Positivity of the
normalization factors leaves cone membership and determinant sign
unchanged.

The executable also requires a connected closed oriented triangulation,
one circular link at each vertex and Euler characteristic2. Omitting a
triangle or changing one orientation cannot produce an accepted degree.

## A shape-independent energy cost of shell charge

For an admissible spherical triangle its unsigned area is

```
A=2 atan(|det(u,v,w)|/[1+u.v+v.w+w.u]).
```

The denominator is at least1. Also
|det(u,v,w)|=|u.((v-u) cross(w-u))|<=|v-u||w-u|.
Using atan z<=z and2ab<=a^2+b^2 proves

```
A <= |u-v|^2+|u-w|^2.                                 (3)
```

The two terms are precisely the two source bonds at that triangle's
right-angle square corner. Summing the two triangles of a square counts
its four sides once. Summing all squares therefore counts each shell
source edge twice. The signed area and(3) give

```
4pi|q| <= sum_triangles A <=2 sum_shell_edges |n_x-n_y|^2.
```

Together with(2), the actual source energy on the shell obeys

```
Delta E_shell >=pi J kappa |q|.                        (4)
```

No radial profile, continuum Taylor approximation or field equation is
used in this inequality. The finite admissibility condition is required.
It cannot be discarded for rough lattice assignments.

Distinct concentric integer-radius shells have disjoint source edges:
every endpoint of a surface edge has the same infinity norm r. More
generally any declared collection of shells with disjoint edge sets can
be summed. Hence

```
Delta E >=pi J kappa sum_r |q_r|.                      (5)
```

For an infinite assignment, onsite rotations define an automorphism on
each local observable and hence a state; its energy excess is the extended
sum of the same nonnegative bond terms. An assignment with admissible
nonzero degree on infinitely many such shells has infinite source energy.
A finite-energy assignment can have nonzero charge on only finitely many
admissible disjoint shells, with their number at most
Delta E/(pi J kappa).

Admissibility at large radius also follows from finite source energy.
If a triangle is nonacute, one pair has dot product less than0 and hence
chord square greater than2. That chord is either a source edge or the
square diagonal. In the latter case its square is at most twice the sum
of the two actual apex-edge chord squares. Thus a nonacute shell has
source chord-square sum greater than1 and energy greater than J kappa/2.
Combining bad shells and admissible charged shells, whose edges are
disjoint, gives the finite bound

```
Delta E >=(J kappa/2) N_bad+pi J kappa sum_good |q_r|.   (5a)
```

Therefore finite energy permits only finitely many nonacute shells and
only finitely many admissible charged shells. Every sufficiently large
shell is admissible and has degree0. This conclusion requires no separate
asymptotic smoothness or admissibility hypothesis. The finite degree
calculation and bound(4) retain their local admissibility condition.

Within the finite rotation orbit, zero excess energy implies every
nearest-neighbor direction agrees by(2), hence n_x=e1 throughout the
connected lattice with its fixed identity exterior. This selects the
uniform direction among these energy minimizers without classifying
other quantum states or the local frame when b=0.

For separated opposite charges, take two shell families of outer radius R
whose enclosing cubes are disjoint. If every shell in each family has
degrees +q and -q respectively, their combined cost is at least
2pi J kappa R|q|. For center separation D in infinity norm, R can be any
positive integer with2R<D, provided these charge hypotheses hold. This
is a lower bound on explicitly assigned directions within the orbit.
It is not a derived particle potential or a proof of a metastable pair.

The radial direction n_x=(x-c)/|x-c| on any such shell is admissible and
has degree+1. Pairwise coordinate differences within a triangle are at
most1, so products of matching integer coordinates are nonnegative;
the fixed face coordinate contributes r^2>0. Radial projection of a
convex cube boundary is an orientation-preserving homeomorphism to S2.
Target reversal gives degree-1. These are exact finite examples without
an assumed continuum shape between lattice sites.

## What a change of degree costs

Consider any continuous onsite-rotation path from a configuration with
nonzero interpolated shell degree to the uniform assignment. The path
need not remain acute. Whenever every triangle's convex hull avoids0,
normalized barycentric interpolation is continuous in space and time,
so its degree cannot change. Some triangle must therefore have
0 in conv(u,v,w) at a crossing.

At the chosen square-corner apex u this implies u.(v+w)<=0. To prove it,
write alpha*u+beta*v+gamma*w=0 with nonnegative coefficients. If alpha>0,
dotting beta*v+gamma*w=-alpha*u with v+w gives
-alpha*u.(v+w)=(beta+gamma)(1+v.w)>=0. If alpha=0, either a coefficient
vanishes, which cannot yield0 from a unit vector, or v=-w. Thus

```
|u-v|^2+|u-w|^2>=4,
Delta E_at_crossing >=2J kappa >J/1000000.             (6)
```

Both chords are actual source edges of that triangle. This is an absolute
energy floor at a topology-changing configuration. It is not a barrier
above the initial defect energy, an activation time or a lifetime bound.
Indeed the charged-shell lower bound(4) already exceeds(6) for |q|>=1.

## A finite deformation with a source-energy bound

Every finite rotation assignment can be continuously erased by onsite
rotations. Place its support in an N by N by N cube and reset sites to
identity in lexicographic order, one site at a time. Choose a unit
quaternion lift with nonnegative scalar part at each site. Normalizing
(1-t)q+t(1,0,0,0), 0<=t<=1, gives a continuous rotation path: its squared
norm is at least (1-t)^2+t^2>=1/2.

For a completed prefix, all reset-reset bonds have zero excess and all
unreset-unreset bonds retain their initial nonnegative contributions.
Only cut bonds can add energy. A lexicographic cut through the cube has
at most N^2+N+1 edges: at most N^2 across planes, N across a partial row,
and1 within that row. The one site rotating during a step touches at
most six additional bonds. Since every bond excess is at most J,

```
max_path Delta E <=Delta E_initial+J(N^2+N+7).          (7)
```

The bound includes bonds to the identity exterior. The path is an
explicit controlled deformation, not an autonomous solution of the spin
Hamiltonian or the microscopic field Hamiltonian. It proves that this
finite lattice charge is not an exact obstruction to local controls.
It supplies neither a low-energy relaxation trajectory nor a sharp bound
on the least possible energy barrier.

## The native repair law does not select these spin directions

Use the same full-domain C or E repair and retain both fair coin values
and all sixteen raw slots. On the bare protected carrier with exactly one
fermion at every site, the native branch law gives current0, unchanged
occupation and spin, sign+1 and raw0. Thus for every edge

```
K_(coin,raw) P =delta_(raw,0) P/sqrt2.                  (8)
```

For arbitrary spin states entangled with an untouched reference, the
complete dilation is the identity source channel tensored with the fresh
coin |+> and raw0. Reading the coin gives a fair bit independent of the
spin state. Induction proves the same statement for fixed finite
repair-only schedules, or schedules depending only on their retained
repair records, with their entire records retained. Onsite spin
rotations preserve this carrier, so records from repair-only insertions
between specified rotations also contain no spin-direction information.
The controlled rotations themselves remain supplied operations.

The package checks(8) on every full-Fock input of the two-site and
three-site native maps, including a nonadjacent seam. It also retains an
explicit control in which the actual source hopping takes a protected
spin state out of the carrier. Therefore(8) is not a claim about
interleaved full Hamiltonian evolution, finite-c dressed ground states,
stationary states of the complete process or a selected winding sector.

The source energy gives finite constraints on the admitted ordered orbit;
repair-only records do not supply a texture preference there. The supplied
radial Dirac mass texture is a distinct operator coefficient. No mapping
from these rotated quantum grounds to that coefficient, no autonomous
formation law, and no observed-particle assignment follows. Responsive
Gram or curvature identities supply no independent winding or absolute
mass selection for this package.
