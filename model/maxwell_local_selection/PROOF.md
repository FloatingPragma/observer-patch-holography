# Local gauge action and conservative field evolution

The committed carrier has twelve ports, thirty oriented seams and twenty
oriented triangular faces. Write `D` for its seam-by-port gradient and `C`
for its face-by-seam curvature. The existing finite theorem stack proves
`CD=0`, `rank D=11`, `rank C=19` and covariance under the sixty committed
orientation-preserving port permutations. This package reads those same
tables; it introduces no geometric mesh or particle mass.

Two conditional selection results reduce the independent choices in the
finite Maxwell construction. They do not derive their constitutive
hypotheses from the native record process.

## Exhaustive classification of the local quadratic action

Consider a real quadratic action on seam potentials,

```
S(A)=constant+ell^T A+(1/2) A^T H A,    H=H^T.
```

Require the following properties.

1. Gauge invariance: `S(A+D chi)=S(A)` for every `A,chi`.
2. Strict face locality of the Hessian: an off-diagonal entry `H_ij` may
   be nonzero only if seams `i,j` occur in the same committed triangular
   face. This is stronger than arbitrary finite-range or vertex-sharing
   locality. The diagonal is unrestricted.
3. Covariance under all sixty committed rotations, acting on seam
   coordinates with their induced orientation signs.

Then, and only then,

```
S(A)=constant+(kappa/2)||CA||^2.
```

The coefficient `kappa` is an arbitrary real number. Nonnegative action
relative to its constant is equivalent to `kappa>=0`; a nonzero such
quadratic form is equivalent to `kappa>0`. The constant and normalization
are not fixed by these hypotheses.

There is a short proof directly from the triangles. Gauge invariance of
the quadratic part is equivalent to `HD=0`; its linear part requires
`D^T ell=0`. Fix an edge `e=(u,v)` of a face `f`, whose third vertex is
`w`. In the `(e,w)` entry of `HD`, strict locality permits exactly the
other two edges `a,b` of that face. Their signs in `CD=0` give

```
H_ea/(C_fe C_fa)=H_eb/(C_fe C_fb).
```

Symmetry of `H` propagates one common weight `w_f` to the three
off-diagonal pairs of this face. Each allowed off-diagonal pair belongs
to exactly one triangle. Subtracting the weighted face forms leaves a
diagonal matrix which still annihilates `D`. Each edge has a nonzero
endpoint incidence, so each remaining diagonal entry vanishes. Thus

```
H=sum_f w_f C_f^T C_f=C^T diag(w_0,...,w_19) C.
```

The weights are unique because each is read off any of its face's
off-diagonal pairs. The producer checks all sixty edge/third-vertex
identities and the unique-face incidence of all sixty pairs on the actual
carrier. No positivity of the individual face weights was assumed.

An independent exact elimination verifies the classification. Strict
face locality leaves ninety symmetric coefficients: thirty diagonal
entries and sixty unordered off-diagonal pairs. The coefficient matrix
of `HD=0` has rank seventy. The twenty face forms have embedding rank
twenty and their constraint product is zero. These are rational ranks,
not numerical singular-value cutoffs or a search through action choices.

For every committed rotation its signed seam permutation `U` and face
permutation `V` satisfy `CU=VC`. Uniqueness of the weight representation
therefore turns invariance of `H` into invariance of `w`. The rotation
action is transitive on faces; equivalently the exact weight constraint
matrix has rank nineteen. Thus all weights equal `kappa`. Finally the
sum of the sixty signed seam action matrices is exactly zero. Group
averaging consequently annihilates every invariant seam vector, so
`ell=0`. This also removes a possible invariant linear term.

Every hypothesis matters. Unequal positive face weights retain gauge
invariance and face locality but lose rotation covariance. The operator
`(C^T C)^2` is gauge invariant, positive and rotation covariant but is not
strictly face local. Both are checked as exact countercontrols. General
nonquadratic local functions of curvature are outside this classification.

The result strengthens `LocalFaceMaxwellAction.lean`: the sum of squared
face curvatures is the unique form under these explicit assumptions,
rather than an unrestricted choice of a gauge-invariant quadratic action.
The new classification is exact executable mathematics, not a new
Lean-kernel theorem.

## Deriving Ampere evolution from energy and reversal

Let electric seam data be `E` and closed magnetic face data be `B`, with
`1^T B=0`. This magnetic space equals `range C`, because `rank C=19`.
Take positive, supplied energy weights `epsilon,beta` and

```
W(E,B)=(epsilon ||E||^2+beta ||B||^2)/2.
```

Assume real linear autonomous evolution. Require the Faraday identity
`dot B=-CE`, conservation of `W` for every initial field, and time reversal
with `E` even and `B` odd. These requirements force

```
dot E=(beta/epsilon) C^T B,
dot B=-CE.
```

This uniquely determines evolution on the closed-flux carrier. An
extension to an unphysical constant-flux input is not selected by a
theorem restricted to closed flux; the producer uses the natural same
formula on the full ambient space, where `C^T 1=0`.

To prove uniqueness, write the unknown electric derivative as
`dot E=A E+F B`. Faraday fixes the entire other row. Energy conservation
is the polynomial identity

```
epsilon E^T A E+E^T(epsilon F-beta C^T)B=0.
```

The independent choices of `E` and closed `B` imply `A^T=-A` and
`F B=(beta/epsilon)C^T B` on that space. Time reversal means
`R G R=-G`, where `R(E,B)=(E,-B)` and `G` is the evolution generator.
Its electric-electric block gives `A=-A`, hence `A=0`. This yields
Ampere evolution without imposing it independently.

Without the specified reversal, an arbitrary skew electric term survives
energy conservation. The executable constructs a nonzero such term from
two actual curl-sector vectors. It preserves Faraday, energy and Gauss,
but violates time reversal. Thus energy conservation alone is insufficient.

The generated evolution preserves Gauss data because `D^T C^T=0`; it
preserves closed flux because `1^T C=0`. In the ambient fifty-coordinate
space the generator has rank thirty-eight and twelve zero modes: eleven
electric gradient modes and one constant magnetic mode. On source-free
Gauss data and closed flux there are nineteen oscillator pairs. This
finite mode count is not the continuum photon polarization count.

The temporal-gauge potential equation is

```
epsilon ddot A+beta C^T C A=0.
```

It is also the Euler--Lagrange equation of
`epsilon||dot A||^2/2-beta||CA||^2/2`. More generally the gauge-covariant
definitions `E=-dot A-D phi`, `B=CA` give Faraday directly. They remain
specified field and potential identifications. The electric energy
coefficient, magnetic coefficient and their ratio are not fixed by gauge
invariance or symmetry. Replacing `beta` by `2 beta` changes the exact
generator and every positive squared frequency by a factor of two while
retaining all structural hypotheses.

## Evolution parameter and discrete implementation

This result derives one continuous linear law from its stated energy,
field definition and reversal requirements. It does not select a native
repair word or a finite-step integrator. At a common nonzero step the
implicit midpoint and velocity-Verlet maps of a supplied oscillator are
distinct, reversible, symplectic maps with the same continuous generator.
The receipt supplies an exact rational example. In particular, the result
does not identify the leapfrog law of `TemporalMaxwellEvolution.lean` with
the integrated-prism law of `WhitneyMaxwellDynamics.lean`.

The Whitney mass/stiffness pair uses a separately supplied geometric
metric and interpolation. It is not constrained to the strict seam
locality class above. Existing Whitney action variations already derive
the recurrence, Noether identity and Gauss propagation from that action;
this package neither replaces its metric nor promotes it to source data.

To obtain a native physical Maxwell derivation still requires a source
field/readout map supporting the quadratic and strict locality hypotheses,
an identification of this positive field energy and its two weights,
selection of the reversible evolution and its clock, and the refinement,
current and laboratory attachment. The package does not identify native
repair dissipation or modular energy with this conservative field energy.
