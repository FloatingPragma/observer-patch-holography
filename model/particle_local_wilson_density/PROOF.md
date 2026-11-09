# Local Wilson spectral weight from plaquette energy density

The declared periodic gauge-coordinate field admits a volume-independent
bound on the **expected local weight** of its Wilson modes with absolute
eigenvalue below one half. The initial density and Hamiltonian are exactly
those of the [referenced local gauge state](../particle_local_gauge_state/PROOF.md).
The bound holds for every program time, on the coordinate comparison parent
and as a compressed effect on its actual finite word codes. It does not
assert that every configuration has a Wilson gap, that low modes are
spatially localized, or that an interacting relativistic limit exists.

## 1. The complete source kernel and a positive decomposition

Let V be the number of cells in a finite periodic three-dimensional cubic
lattice. Keep the actual fifteen-dimensional representation R of
G=S(U(3)×U(2)) and all four Dirac indices. On the one-body space of
dimension 60V, write the covariant color/site shifts as T_mu. Each is unitary.
Set D_mu=T_mu-I and Dbar_mu=I-T_mu^dagger=-D_mu^dagger. With the source's
Hermitian Clifford matrices gamma_mu, mu=1,2,3, and beta=gamma_tuple[3],

```
A = I + (1/2) sum_mu [(I-gamma_mu)D_mu-(I+gamma_mu)Dbar_mu],
H_W = -beta A,           H_W²=A^dagger A.
```

This is the Wilson parameter s=0, equivalently negative Wilson mass one,
of the [owning spatial field](../interacting_extension/chiral_lattice/PROOF.md).
It uses beta, not the source's distinct gamma5. Spatial shifts act on
color/sites and commute with the Clifford matrices. No fermion occupation,
internal representation row or Brillouin-zone corner is deleted.

The algebra in Appendix C of
[Hernández, Jansen and Lüscher](https://arxiv.org/html/hep-lat/9808010v2)
has a useful trace-norm version. Here is its complete three-dimensional
form, with every pair ordered. Put sigma_mu,nu=(i/2)[gamma_mu,gamma_nu] and

```
B_mu,nu = Dbar_mu D_mu Dbar_nu D_nu,
C_mu,nu = (i/2) sigma_mu,nu [Dbar_mu+D_mu,Dbar_nu+D_nu],
J_mu,nu = -gamma_mu [Dbar_mu+D_mu,Dbar_nu-D_nu].
```

Expanding the Clifford products gives

```
H_W² = I + (1/4) sum_(mu != nu) (B_mu,nu+C_mu,nu+J_mu,nu).       (1)
```

Because `Dbar_nu D_nu=-(Dbar_nu-D_nu)`, the first term splits exactly as

```
B_mu,nu = (D_nu D_mu)^dagger(D_nu D_mu)+R_mu,nu,
R_mu,nu = -Dbar_mu [D_mu,Dbar_nu-D_nu].                         (2)
```

Consequently `H_W²=I+P+E`, where P is the sum of the positive squares in
(2), divided by four. E is the remaining sum. The total E is Hermitian,
although its individual ordered-pair terms need not be. These identities
use only unitarity and the Clifford relations; they require no spectral
gap or pointwise plaquette admissibility. The code checks them on the
complete native representation, including a configuration with zero modes.

## 2. A trace estimate that counts bad modes

For a fixed positive plaquette orientation mu<nu, let

```
L_mu,nu = sum_x ||I-R(U_mu,nu(x))||_1.
```

Here the norm is the fifteen-dimensional matrix trace norm. Each of the
four signed-shift commutators `[T_mu^sign,T_nu^sign]` is a unitary shift
times a block-diagonal difference of the two transports around one
plaquette. Unitary invariance of singular values gives its trace norm
exactly L_mu,nu, with base points translated or orientation reversed as
necessary. Both transformations preserve the complete periodic sum.
The corresponding norm after retaining all Dirac indices is 4L_mu,nu.
These facts also cover short periods and repeated links.

Use `||Dbar_mu||<=2`, `||gamma_mu||=||sigma_mu,nu||=1` and the ideal
property of the trace norm. Relative to the bound 4L_mu,nu on each
spin-lifted elementary commutator, the coefficients for R, C and J are
respectively 4, 2, 4. There are two ordered pairs for each positive
orientation and the prefactor in (1) is 1/4. Thus

```
||E||_1 <= 20 sum_p ||I-R(U_p)||_1.                            (3)
```

This sum is extensive. It is not replaced by a volume-divergent operator
norm advertised as a local estimate.

Let `Pi=1_{|H_W|<1/2}`, the spectral projector of the full coordinate
symbol, including every degenerate eigenvector. Since P>=0,

```
(3/4) Tr Pi <= Tr[Pi(I-H_W²)]
             = -Tr(Pi P)-Tr(Pi E) <= ||E||_1.                 (4)
```

The first inequality can be strict for nonzero Pi; the weak form is
sufficient. Define the same positive magnetic cost as the field action,

```
s_p = 1-Re Tr R(U_p)/15 = ||I-R(U_p)||_F²/30,
S = sum_p s_p.
```

The trace/Frobenius inequality in dimension 15 gives
`||I-R(U_p)||_1<=15 sqrt(2s_p)`. There are 3V positive plaquettes.
Cauchy–Schwarz in the plaquette sum and (3)–(4) therefore yield the
pointwise mode-count inequality

```
Tr Pi <= min(60V, 400 sqrt(6 V S)).                            (5)
```

It is compatible with arbitrarily small nonzero S producing no low mode,
and with a sparse set of rough configurations carrying low modes. It
does not assert a lower bound on the smallest eigenvalue.

## 3. The same quantum state supplies a local bound

Let rho be any state of the complete coordinate/Fock/reference system
whose coordinate marginal is translation invariant. It may have gauge
coherences, matter entanglement and an arbitrary retained reference.
Expectations of coordinate multiplication operators are integrals against
its positive diagonal coordinate measure. Translation covariance of H_W
and Pi makes the expected trace of the site-diagonal block of Pi the same
at every site. For a set X of cells, with the full Dirac/M15 projector P_X,

```
<Tr(P_X Pi P_X)> = (|X|/V)<Tr Pi>
 <= min(60|X|,400|X| sqrt(6 <S>/V)).                           (6)
```

The last step uses concavity of the square root. Agreement of the gauge
diagonal measure under translations suffices; no particular fermion number
or single-particle preparation is assumed.

The [local gauge-state planner](../particle_local_gauge_state/state.py)
uses one per-link degree L and one action coefficient kappa for all
volumes, retains the complete Fock carrier and proves

```
<S>_t/V <= s_0+(e_E+2b)/kappa <= eta/48000
```

at every model time, where `eta=2^(-target_bits)` and
`b=60(2+epsilon_m)/a`. Substitution into (6) gives

```
<Tr(P_X Pi P_X)>_t <= min(60|X|,2 sqrt(5 eta)|X|).              (7)
```

All time and volume quantifiers refer to that same prepared polynomial
Hamiltonian and unchanged clock. The density is not a ground state or a
field-only Gauss singlet; its conjugate reference supplies the compensating
charges. This result does not transfer the bound to another preparation.

## 4. A local sign-approximation consequence

The source Legendre polynomial q_m obeys the bounds
`|q_m(x)|<=1+epsilon_m` on the full Wilson spectral interval and
`|sign(x)-q_m(x)|<=epsilon_m=16(15/17)^(m+1)` on `|x|>=1/2`.
Choose sign(0)=0; another bounded convention at zero obeys the same
estimates. Spectral calculus, followed by the finite trace, gives

```
||[sign(H_W)-q_m(H_W)]P_X||_HS²
 <= 60|X| epsilon_m²+(2+epsilon_m)² Tr(P_X Pi P_X).             (8)
```

The coefficient2+epsilon_m is conservative. Taking expectations and
using (6) or (7) gives a volume-independent local Hilbert–Schmidt error.
The executable stores exact rational coefficients and the squared
radicand, so no rounded square root becomes a proof premise.

Equation(8) controls a spatial block of the one-body sign function in
the actual quantum coordinate state. It does not control the full-Fock
norm of an extensive Hamiltonian difference, or its commutator with a
time-evolved local detector. Inferring a complete finite-duration
interacting record comparison from (8) requires the separate local
Duhamel estimate stated below.

## 5. Complete finite effects and exact controls

On the actual finite word code with projection P_K, K>=3VL, the same
translation-invariant density, compressed Hamiltonian and energy argument
apply. Compression of the nonnegative coordinate inequalities preserves
them. Define a complete binary effect by

```
E_low = P_K [Tr(P_X Pi P_X)/(60|X|)] P_K,
E_other = I_code-E_low.                                      (9)
```

Both effects are positive and sum to identity on every code input and
every Fock occupation, including the unused directions of the preparation.
The square-root instrument with both pointer labels is a permitted
finite apparatus. Its coefficients are finite integrals of a bounded
spectral function. Its exact operation may be declared a primitive; a
compiled approximation needs a complete-channel error. The retained
reference is untouched. Neither outcome is discarded or renormalized.
The low-label probability has the explicit bound from (7) divided by
60|X|. This normalized probe is a diagnostic of Wilson spectral weight;
it is not a detector of an observed particle or a particle mass.

Equation(9) is the compression of a coordinate-symbol effect. It is
not the spectral projector of a separately truncated Wilson matrix.
The pointwise eigenvalue control and the finite observer effect have
different domains, joined by the stated compression.

The exact native control uses weak links
`(I+i sqrt3 sigma1)/2`, `(I+i sqrt3 sigma2)/2`, I, lifted into the actual
fifteen exterior basis vectors. All six ordered identities (1)–(2)
hold. The full sixty-dimensional squared Wilson kernel has ranks 16, 28, 16
at eigenvalues 0, 1, 3. Its complete uniform auxiliary spectral law is
4/15, 7/15, 4/15, and S=3/5. The exact projector polynomials sum to identity.
The flat native symbol also checks (1), with H_W²=I. These coordinate
controls verify algebra and multiplicities; a point coordinate is not
called a normalizable gauge-field preparation.

## 6. The remaining physical interface

The result supplies local low-mode **weight**, not decay of each low-mode
projector between distant regions. A volume-uniform finite-duration
comparison of the full interacting exact-sign and polynomial programs
must bound the same state's local response to their Hamiltonian difference,
including the dynamical electric field and the retained reference. This
is the [local Duhamel target](../particle_field_preparation/REFINEMENT.md#the-missing-local-operator-estimate).
Equation(8) supplies its local spectral input but is not that theorem.
Intermediate interventions need their own energy and symmetry accounting.

An interacting chiral relativistic limit additionally needs the stated
refinement, state, current and connected-response estimates. The Wilson
kernel's auxiliary gap is not a particle mass gap; the free overlap field
has a zero-energy mode even when its Wilson kernel is gapped. No measured
mass, source-to-species assignment or Lorentz-invariant scattering amplitude
is inferred from this construction.
