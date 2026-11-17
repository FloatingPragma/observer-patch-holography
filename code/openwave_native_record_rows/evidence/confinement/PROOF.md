# Native centre records and the full-colour string bound

This is a new unpinned N+ candidate, `A-native-triality-string-v1`.
It changes memory, preparation and observation architecture, not the pinned
endpoint rule. The result below is an exact dimensionless static potential.
It does not establish the physical Cornell tension or identify the endpoint
records with physical quarks. No measured coupling enters this construction.

## 1. Finite theorem stack and source distinction

The row-12 backbone supplies the conditional fundamental SU(3) colour
representation and its centre on the native fifteen-coordinate representation.
The complete Fock-space identity in
`model/openwave_native_centre_string/PROOF.md` is
`exp(2 pi i Q) exp(2 pi i C)=I`: the fundamental endpoint carries nonzero
triality. This identity is an interface for the present construction; this
packet does not upgrade its physical-source identification.

`Screen/DiscreteGauss.lean` proves the complete neutral incidence fibre and
its cycle ambiguity. We retain the analogous complete mod-three fibre below;
we do not select a spanning-tree string as the entire state space.
The existing full-SU(3) Haar chain and positive-cone proofs correctly establish
linear energy for their supplied electric action. Their parameter kappa and
retrospectively selected 1109.457... do not enter here. Instead the positive
transfer, and hence the energy coefficient, are reconstructed from native
movement records.

## 2. Native finite local instrument

Use one persistent two-port source with one unit and the pinned
`native(loads,(0,1),coin)` with independent fair placement coins. Retain both
coin and signed current at every call. Let `M=abs(current)`. For either
current location of the unit, exactly one placement coin gives M=0 and one
gives M=1. Therefore successive M values are independent fair bits, also
conditional on the complete past. This assertion is checked on all initial
locations and all two-coin words; it uses the endpoint rule, not just coins.

An addressed edge has a classical counter `z_e in Z/3`. Two consecutive
movement bits write

    delta = 0                  if M1=0,
    delta = +1                 if M1=1, M2=0,
    delta = -1                 if M1=1, M2=1.

The observer updates `z_e <- z_e+delta mod3`. The counter is declared
apparatus. Its one-update transition is exactly

    P_e = (2 I + X_e + X_e^-1)/4.

It has eigenvalues 1 on the trivial character and 1/4 on each nontrivial
character. It is a positive symmetric stochastic matrix on the whole
counter space, rather than a selected signed sector alone. On a finite
graph with m edges, one layer addresses every edge once in a frozen order.
The conditional fair-bit result makes the counter transition the exact
product `T=product_e P_e`, despite the persistent bit-cell position.
The complete native tape and that position remain retained; positivity is
claimed for this closed counter-observable algebra, not for an unaudited
larger native-plus-memory transfer.

The count-time positive generator and the layer-time generator are different
declared clocks:

    H_count = sum_e (I-P_e) = (3/4) sum_e N_e,
    H_layer = -log T = log(4) sum_e N_e,

where N_e projects to nonzero edge character. The first uses one independent
Poissonized update per edge in one unit of reconstruction time; no Poisson
timing claim establishes spacetime symmetry. The second uses exactly 2m
primitive calls per layer. The spectrum and both clock normalizations are
reported separately. No physical-time identification is inferred.

If the native endpoint mutation holds its loads and reports zero current,
all M vanish and every counter update is the identity, giving H=0. Feeding
the unchanged coin directly into the counter would miss this control.

## 3. All charge and flux sectors, without an imposed path

Orient every edge of any finite graph G. The complete Fourier basis is
indexed by `t in (Z/3)^E`. Its vertex centre charge is `q=boundary(t)` and its
energy is `(3/4)|supp(t)|` in count units. Every t is retained, including
all q, closed cycles, and zero-energy t=0. The physical endpoint sector in
this theorem has q(x)=+1, q(y)=-1 and q=0 elsewhere.

For any connected component K of the nonzero-edge support, summing q over
its vertices cancels each internal oriented edge and gives zero mod3.
Consequently a component containing x must also contain y. Its support
contains a path from x to y, and therefore

    |supp(t)| >= d_G(x,y).

A shortest path, with flux +1 oriented from x to y and flux zero elsewhere,
is an admissible vector attaining equality. Thus minimizing over the FULL
Gauss fibre, rather than a chosen line, gives

    V_count(x,y) = (3/4) d_G(x,y),
    V_layer(x,y) = log(4) d_G(x,y).

The vacuum minimum is zero. Neutral cycles and branching three-triality
configurations remain in the spectrum. The finite audit enumerates every
flux vector of the declared chain, cycle, ladder and cube, retaining the
complete charge/energy multiplicity table. Disconnected endpoints have an
empty charged fibre; this theorem assumes connected x and y.

## 4. Full SU(3) extension and a uniform bound on all irreducibles

Let omega=exp(2 pi i/3), and describe an edge record by U_e in SU(3). Its
update multiplies U_e by omega^delta I. On the full Haar space L2(SU(3)^E),
the SAME source transition is

    P_e = (2 I + R_(omega I) + R_(omega^-1 I))/4.

It commutes exactly with both left and right SU(3) actions and hence with
every local gauge transformation. On any irreducible representation (p,r),
the centre acts by omega^(p-r). The eigenvalue is therefore 1 for triality
zero and 1/4 for triality nonzero, with NO bound on p or r. This is a uniform
statement about the complete Peter-Weyl space. It does not infer an
operator-norm gap from weak convergence of a discrete Haar sampler.

Add a fundamental colour space at x and its conjugate at y and take the
complete invariant Gauss space. The centre subgroup already forces the
support argument of section 3 in every nonzero invariant component. Thus
the same lower bound holds on the full SU(3) Gauss space. A fundamental
matrix-product path, contracted with the endpoint colour indices, is a
nonzero invariant spin network: its edge representations are fundamental
or conjugate fundamental along the path and trivial elsewhere. It attains
the bound. This proves the exact potential in section 3 on the full colour
space, uniformly across all higher representations.

In particular, the proof does not replace SU(3) by a truncated list of three
irreducibles. It does have an extensive zero-energy complement: ALL
centre-neutral irreducibles, including the adjoint, have zero local energy.
This is a triality projector energy, not the Casimir action of the older
colour chain. It does not establish a Yang-Mills continuum, finite thermal
partition trace, asymptotic freedom, or an adjoint mass gap. Those claims
cannot be borrowed from the older supplied action.

## 5. Native preparation limit for colour records

The counter process is finite and exact without any Haar sampler. To realize
the full colour reference as a limit of actual classical records, rather
than install a Haar random oracle, define this explicit optional preparation
family. For each of 18 uniforms, read B fresh native movement bits as an
integer K and use the midpoint `u=(K+1/2)/2^B`. Nine Box-Muller pairs produce
18 real numbers; group them into a complex 3x3 matrix. QR with positive real
diagonal in R gives Q. Multiply its last column by conjugate(det Q) to
obtain determinant one. A singular/ill-conditioned draw, with a Gram-Schmidt
residual at most 2^(-B/4), writes the identity and an explicit failure flag.
No draw is retried or discarded. Supply an independent uniform initial
centre counter by the three equally weighted, disclosed preparations.

As B increases, the midpoint uniform vector converges in distribution to
continuous independent uniforms. Box-Muller and the QR map are continuous
outside a measure-zero set. The complex Gaussian matrix is nonsingular
almost surely; the failure threshold tends to zero. The positive-diagonal
QR distribution is Haar on U(3), and determinant correction is SU(3)-left
equivariant, hence gives Haar on SU(3). This establishes convergence of the
full preparation law and all fixed bounded continuous record moments. The
QR fact and phase convention are described by Francesco Mezzadri,
"How to generate random matrices from the classical compact groups",
https://arxiv.org/abs/math-ph/0609050 .

At finite B the colour preparation is centre invariant after averaging its
three disclosed centre preparations, but generally not fully SU(3)
invariant. We do not assert exact finite-B full-gauge symmetry. Its centre
record transfer nevertheless retains the exact positive spectrum above;
the full SU(3) assertion is a source-limit theorem. The uniform triality
lower bound prevents a sequence of high irreducibles from reducing the
limiting charged energy below it. Every finite collection of colour-record
correlations converges, since the centre updates are continuous and the
native centre law is independent of the preparation history.

The numerical preparation audit is a finite diagnostic, not the proof of
Haar convergence or its rate. Raw bit tapes, failed draws, finite-B laws and
all sample moments are retained. The compiler's logarithm, trigonometric
functions and QR are specified classical decoder operations, with numerical
roundoff bounded in the audit by unitarity and determinant residual checks;
no claim of a physical quantum state-preparation device follows.

## 6. Observable, endpoint dictionary and physical scale

For any oriented path gamma the record score
`chi_gamma=product_e omega^(orientation_e*z_e)` has one-layer correlation
`4^(-length(gamma))`. Its logarithmic decay gives the layer energy. The
complete Gauss minimization in sections 3 and 4 identifies which endpoint
sector that path probes and establishes its ground energy, instead of
assuming a path score alone is a potential. The source control gives
correlation 1 and energy zero. Signed and complex record scores are retained
as scores, not treated as accepted probability weights.

To attach this to physical quarks one must join the source's endpoint
records to the row-12 colour/current representation on the same candidate,
including the full native charge instrument and the physical reference.
The supplied fundamental representation plus the existing quotient identity
is a precise algebraic interface; it is not a completed physical identity.

If one count-energy unit is E0 GeV and one graph edge is a fm, the tension
is `(3/4) E0/a` GeV/fm. With the layer clock the coefficient is log(4).
Neither E0 nor a is fixed by the present source. Choosing their ratio to
obtain 1 GeV/fm would calibrate the requested result. In a common relativistic
dictionary with dimensionless limiting speed c0, the dimensionless ratio
`sigma hbar c / m_e^2` would be `(3/4)c0/m_e,source^2` for the count clock;
the source electron energy and common speed must be supplied by an actual
joined measurement, not selected to match the tension. Graph distance also
requires the common physical ruler; a cubic nearest-neighbour graph alone
has directional lattice distances, not an isotropic continuum norm.

The packet therefore discharges a constructive native linear-potential
component with a full-colour uniform bound. It leaves the physical endpoint,
common clock/ruler, physical tension, one-candidate join, and external G6
review open. A Coulomb correction from the same charged source, where the
full Cornell form is required, is also owed. No full row award is claimed.
