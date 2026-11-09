# An isolated particle band of the declared high-field spin source

The declared cubic spin Hamiltonian has an isolated lattice particle band
for sufficiently small positive `J/B`. Its first momentum-dependent term
has order `J²/B`, and its two sublattice charges have the same dispersion.
These statements concern the ground representation of this Hamiltonian.
They do not identify a particle observed in nature or derive this action
from an unobserved native repair schedule. A refinement with bounded
positive rest energy and bounded `J/B` has vanishing physical propagation
speed on the same ground clock.

## 1. Source, field normalization and reference

Let `epsilon_x=(-1)^(x1+x2+x3)`. On a finite bipartite graph the physical
source is

\[
 H_B=J\sum_{xy}(S_x\cdot S_y-I/4)
       -B\sum_x\epsilon_x S_x^1,\qquad J,B>0.             \tag{1}
\]

The actual integer matrix in `native_ordered_ground.source.matrix` uses
the Marshall frame. A product of Hadamard matrices then gives exactly

\[
 \widetilde H_B={J\over4}\sum_{xy}
       (X_xX_y-Y_xY_y-Z_xZ_y-I)-{B\over2}\sum_xZ_x.       \tag{2}
\]

Both frame changes are onsite unitaries. They preserve spatial support,
spectra, complete outcome probabilities and the evolution parameter. Put
`n_x=(I-Z_x)/2`, `a_x^dagger=(X_x-iY_x)/2=|1><0|` and `lambda=J/B`.
Adding the scalar `B|Lambda|/2` and dividing by `B` gives

\[
 N+\lambda V,\qquad N=\sum_xn_x,
\]
\[
 V=\sum_{xy}\left\{ {a_x^\dagger a_y^\dagger+a_xa_y\over2}
       -{1\over2}+{n_x+n_y\over2}-n_xn_y\right\}.       \tag{3}
\]

The particle generator is `K_B=H_B-E_0(B,J)`, or its positive
infinite-volume ground-representation generator. The scalar in (3) is a
calculation convention; the excitation energies always subtract the
interacting ground energy. Flip number is not conserved. Its parity and
the staggered charge described in Section 4 are conserved.

For the infinite simple cubic lattice group the three positive-axis bonds
at each site into one interaction supported on
`Lambda_0={0,e1,e2,e3}`. Each shifted bond in (2) has norm `J`, so the
dimensionless interaction satisfies

\[
                      \|\phi_x\|\le3\lambda.            \tag{4}
\]

The onsite Hamiltonian `n_x` has a unique product ground and the simple
excitation `mu=1`. All sums of two or more positive onsite energies are
integers at least two. Thus the nonresonance condition is satisfied.

## 2. The infinite-volume band and the imported smallness condition

[Yarotsky, Theorems 1--4](https://arxiv.org/abs/math-ph/0411042)
apply to (3)--(4): for sufficiently small interaction norm the product
ground has a unique gapped continuation, and the spectral component near
one is a translation-covariant particle subspace with real analytic
momentum dispersion. The constants depend on the finite interaction
range. The source substitution supplies `3J/B` as the norm bound; the
paper supplies no evaluated numerical threshold for this source.

In particular there is a number `lambda_*>0` such that for
`0<lambda<lambda_*` the positive generator has a nonzero band projection
`P_1`, and

\[
 P_1\mathcal H\simeq\ell^2(\mathbb Z^3),\qquad
 K_B|_{P_1\mathcal H}\simeq E_{B,J}(k),\quad
               k\in(\mathbb R/2\pi\mathbb Z)^3.         \tag{5}
\]

The torus in (5) uses radian momentum. No finite choice of `J/B` is
certified to satisfy the imported threshold by a small-box diagonalization.

The boundary convention does not change the selected infinite state.
Yarotsky's complete-star restriction omits every interaction whose entire
support exits the box. The source restriction retains every internal
bond. These are positive-weight bipartite graphs with identical couplings
on each fixed interior region. The cut-bond estimate in
[the fixed-field locality proof](../spin_fixed_field_locality/PROOF.md)
applies to both exhaustions, so their fixed-local-support ground marginals
converge to the same state. This identifies the state in (5) with the
declared source ground, without imposing a boundary field or an additional
state average.

Write `delta=3c_2 lambda`, where `c_2` is the qualitative spectral constant
in Yarotsky's Theorem 1. Shrink `lambda_*` so that `delta<1/3` and all band
construction bounds hold. The normalized spectrum lies in neighborhoods
of the nonnegative integers of radius `delta` times the integer. Hence the
band near one is separated from zero and from the spectral sets starting
at `2(1-delta)`. Its width is at most `2delta`, smaller than the positive
gap `1-delta`. These are inequalities involving an unspecified theorem
constant, not a numerical admission test.

## 3. Ground subtraction and the first dispersive coefficient

The finite calculation admits positive integer edge multiplicities
`w_xy`. Let `A_xy=w_xy`, `d_x=sum_y w_xy`, and `s_x=sum_y w_xy^2`.
All sums over distinct edges below use these weights once per undirected
edge. Direct perturbation theory for (3) gives

\[
 {E_0\over B}=-{|\Lambda|\over2}
       -{\lambda\over2}\sum_{xy}w_{xy}
       -{\lambda^2\over8}\sum_{xy}w_{xy}^2
       +O_\Lambda(\lambda^3).                           \tag{6}
\]

The first-order one-flip block is diagonal. Its absolute coefficient is
`-sum_edges w/2+d_x/2`. A one-flip state can leave its subspace only by
creating a pair on an edge disjoint from the flipped site. Every such
intermediate state has three flips and denominator `1-3=-2`.
The pair-creation amplitude is `w_xy/2`; it is added before squaring
when several listed bonds have the same endpoints. Before ground
subtraction the second-order diagonal is
`-(sum_edges w^2-s_x)/8`. For `x!=y` the off-diagonal entry is
`-sum_z w_xz w_zy/8`. Subtracting (6) therefore gives

\[
 {K_{\rm eff}\over B}
   =I+{\lambda\over2}\operatorname{diag}(d)
       +{\lambda^2\over8}\{2\operatorname{diag}(s)-A^2\}
       +O_\Lambda(\lambda^3).                           \tag{7}
\]

Equation (7) is an exact statement about finite Taylor coefficients.
`algebra.expansion` evaluates the original full source matrix and all
complementary states with rational arithmetic. It also checks the full
order-two intertwining equation; the second wave correction includes
five-flip states where they occur. It does not estimate a numerical
finite-coupling remainder.

For the simple infinite cubic graph `d_x=s_x=6`. The second coefficient
has finite support: its diagonal is `3/4`, its six axial distance-two
entries are `-1/8`, and its twelve mixed-axis entries are `-1/4`.
These follow by counting all 36 two-edge walks on `Z^3`, independently of
any finite torus. With `gamma(k)=sum_(j=1)^3 cos(k_j)` its symbol is

\[
 {E_{B,J}(k)\over B}
   =1+3\lambda+{\lambda^2\over2}\{3-\gamma(k)^2\}
       +O(\lambda^3).                                  \tag{8}
\]

The uniform meaning of the remainder in (8) requires the following
analytic-kernel argument, in addition to the finite matrix calculation.
Fix an exponential weight `alpha>0` and then choose a fixed cluster
parameter `epsilon<exp(-2alpha)`. In the proof of Yarotsky's Lemma 1 the
ground coefficients are the fixed point of an analytic contraction in the
corresponding weighted cluster norm. Explicitly write
`Omega(v)=exp(S(v)) Omega_0`, `S(v)=sum_(I nonempty) vhat_I`, where
`vhat_I` creates the excited cluster `v_I` from the onsite vacuum. Let
`P_I` extract its cluster coefficient and let `H_(I,0)` be the positive
restriction of `N` to that excited cluster space. The creation operators
commute, so the nonempty coefficients of the eigenvector equation are

\[
 v_I=F_\lambda(v)_I
   =-\lambda H_{I,0}^{-1}P_I
           e^{-S(v)}V e^{S(v)}\Omega_0,\qquad I\ne\varnothing.
                                                               \tag{8a}
\]

The scalar eigenvalue contributes only to the empty coefficient.
The weighted norm is
`sup_x sum_(I containing x) ||H_(I,0)v_I|| epsilon^(-d_I-1)`,
where `d_I` is the minimum connecting-graph length. Expanding (8a) in
commutators gives the absolute cluster majorants in Lemma 1. They bound
the map and its Lipschitz constant with `|lambda|`, for complex as well
as real coupling. Choose a smaller complex disk on which this analytic
map preserves a fixed ball and has contraction factor strictly below
one, uniformly in volume. Iteration from zero converges uniformly there;
its limit is analytic.
This step uses the contraction construction, rather than inferring
coupling analyticity from momentum analyticity.

Choose the fixed Riesz contour `|z-1|=1/4`, shrinking that disk if needed.
The relative resolvent series in the proof, Eqs. (14)--(15), converges uniformly
on the contour in the weighted excitation coordinates. Contour integration
therefore supplies analytic projected local excitation coordinates. Their
Gram and Hamiltonian kernels are analytic in the weighted convolution
algebra

\[
               \|t\|_\alpha=\sum_x e^{\alpha|x|_1}|t_x|. \tag{9}
\]

For this real source the real-coupling inner products have an analytic
continuation using the analytic left and right ground coordinates; complex
conjugation of a variable is not used as a holomorphic operation. The
cluster bounds and projected-vector estimates, Eqs. (19)--(24), control
these local pairings uniformly. Their complex continuation is an analytic
calculation device, not a positive state or orthogonal physical projection
at nonreal coupling. The Gram kernel
approaches the identity in (9), so its inverse square root is given by
its norm-convergent power series. This produces an analytic band kernel
in (9). All these operations commute with the local thermodynamic limit
by the same uniform estimates.

After the ground subtraction, the order-one and order-two coefficients
in (7) depend only on paths within distance two of their entries. They
therefore stabilize in the bulk. A near-identity change from the finite
intermediate-normalization basis to the orthonormal band basis cannot
change the second bulk coefficient: the preceding coefficients are
`I` and `3I`, which commute with that change. Cauchy estimates in (9)
give a bound `C|lambda|^3` on the remaining kernel on a smaller disk.
Since `|x|^2` is controlled by an exponential weight, the remainder in
(8) is uniform in `C^2(k)`. This argument proves existence of a radius
and a remainder constant. Neither is evaluated by the finite certificate.

For sufficiently small positive `lambda`, the band is nonconstant:
the order-two difference between `k=(pi/2,pi/2,pi/2)` and `k=0` is
`9lambda^2/2`. Cubic symmetry and inversion make zero stationary; the
positive leading Hessian and the uniform remainder give

\[
 E_{\min}=B+3J-{3J^2\over B}+O(J^3/B^2),\qquad
 D^2 E(0)={3J^2\over B}I_3+O(J^3/B^2).                 \tag{10}
\]

The minima are zero and `pi(1,1,1)` in the full torus for sufficiently
small coupling. Away from neighborhoods of these two points, the
order-two minimum is separated by a fixed positive amount; in those
neighborhoods its Hessian is positive. The uniform `C^2` estimate and
the exact symmetries preserve both statements. The second-order quartic
term about zero is

\[
       -{J^2\over8B}\left\{|k|^4+\sum_jk_j^4\right\}.  \tag{11}
\]

The cubic anisotropy in (11) and the lattice period are retained. A
positive curvature does not establish a relativistic dispersion law.

The eight-site open source cube has degree three. The declared periodic
`2x2x2` control counts both positive-axis bonds on each length-two axis:
it has `d_x=6`, `s_x=12`, and coefficient `(24I-A^2)/8` in (7).
Its parallel-edge amplitudes cannot be identified with the bulk
unit-edge coefficient `(12I-A^2)/8`.

## 4. Internal charges and an original-probability local preparation

In the transformed frame the physical uniform spin generator is

\[
 C={1\over2}\sum_x\epsilon_xZ_x,\qquad
 C-C_0=-\sum_x\epsilon_x n_x,\qquad
 C_0={1\over2}\sum_x\epsilon_x.                         \tag{12}
\]

Every created pair has opposite sublattice signs, so (2) commutes with
`C`. The finite ground has charge `C_0`: it is the unique continuous
continuation of the product ground and cannot change its discrete charge
without a degeneracy. In infinite volume the charge denotes the onsite
`U(1)` automorphism and its ground-fixed representation, rather than an
unconverged sum of local generators. One-flip states have relative charge
`-epsilon_x`. The band projection preserves this charge exactly.

Ordinary one-site translations are symmetries of (2), and exchange the
two relative charges. In the physical frame they include the onsite
rotation which compensates the staggered field. The even-parity lattice
translations preserve each charge. Thus each charge has one band over
the even-translation momentum torus, with equal energies. In full-torus
notation this is the exact identity

\[
                 E(k+\pi(1,1,1))=E(k).                 \tag{13}
\]

A full-lattice momentum eigenvector and a definite charge are different
bases of this folded description. The labels are internal spin charges;
no electromagnetic charge or natural species assignment is made.

At a specified site use the two Kraus operators

\[
 K_s={X_x-iY_x\over2},\qquad K_f=n_x,\qquad
 K_s^*K_s+K_f^*K_f=(1-n_x)+n_x=I.                       \tag{14}
\]

The success operator has relative charge `-epsilon_x`, and the other
operator is neutral. On the same ground, the unnormalized band weight is

\[
 q_x(\lambda)=\|P_1 K_s\Omega_\lambda\|^2.              \tag{15}
\]

The projected-local-vector Gram estimate in Yarotsky's construction,
Eq. (22), gives `q_x(lambda)->1` as `lambda->0`. The same estimate bounds
the off-diagonal Gram kernel exponentially. Consequently a sufficiently
small coupling has a positive original band weight and a stable spatial
Gram matrix. There is no numerical floor asserted at a specified value
of `J/B` without an evaluated smallness bound.

The complete probabilities are `q_x`,
`||K_s Omega||^2-q_x`, and `||K_f Omega||^2`; they sum to one. Dividing
by the success probability is unnecessary. The spectral split in this
statement is a mathematical projection; it does not declare an exact
native implementation of an infinite-volume spectral projector.
Finite source instruments, clock resolution, local marginal error and
native service-failure outcomes require their corresponding complete
readout bounds. The positive-field preparation and locality proofs
provide finite ground and local-state approximations with those scopes.

In the Marshall source frame (14) is the complete `X` pointer followed by
`Z` on its positive outcome: `K_s=Z(I+X)/2`, `K_f=(I-X)/2`. Embed these
spin blocks into the full onsite Fock space with ordered states
`empty,up,down,double`. Set the success block to zero on empty and double,
and set the failure block to the identity on those two states. The two
full four-state effects sum to the identity. Both Kraus operators commute
with the original fermion grading `Gamma=(-1)^N`, and successful operators
on disjoint sites commute. Flip-parity oddness in (3) is therefore distinct
from original fermion-parity oddness. Under a grading-preserving local
dictionary this source operator cannot create a single odd physical
fermion from a reference of fixed grading, by the argument in
[the exact grading exclusion](../particle_identification/README.md#exact-grading-exclusion).
Even grading alone determines neither Lorentz spin nor a physical species.
`readout.creation_certificate` checks the complete four-state isometry
and these grading statements. This abstract completion is distinct from
the literal pointer circuit outside the protected spin code: the circuit
does not perform the added occupancy-failure test. Their equality on the
protected code is sufficient for the finite source calculation below.

### A finite product residual and an executed native instrument

Let `n=|Lambda|`, `e=sum_edges w_xy` and `S=sum_edges w_xy^2`.
The product vacuum in the frame (2), or its staggered `X` image in the
physical frame, has the exact Rayleigh value and squared residual

\[
 r_0=-Bn/2-Je/2,\qquad
       \|(H_B-r_0)\Omega_0\|^2=J^2S/4.                 \tag{15a}
\]

Each distinct bond creates a different two-flip basis vector with
amplitude `Jw_xy/2`, which proves the residual identity. The exchange
operator satisfies `-Je I<=H_J<=0`. The onsite field alone has its
second eigenvalue `-Bn/2+B`. By the min-max principle the second
eigenvalue of the full source, counting multiplicity, obeys

\[
 \lambda_1(H_B)\ge-Bn/2+B-Je=r_0+D,\qquad
                         D=B-Je/2.                    \tag{15b}
\]

If `D>0`, the variational bound `E_0<=r_0<lambda_1` proves uniqueness
of the full finite ground without a smallness theorem. Every excited
component in the spectral expansion of `Omega_0` has
`|E-r_0|>=D`. Its total weight is therefore at most
`J^2S/(4D^2)`. For the normalized ground `Omega_g`, the pure-state
half trace distance satisfies

\[
 {1\over2}\|\,|\Omega_0\rangle\langle\Omega_0|
              -|\Omega_g\rangle\langle\Omega_g|\,\|_1
       \le \min\{1,J\sqrt S/(2D)\}.                    \tag{15c}
\]

For the actual open eight-site cube with `J=1,B=400`,
`e=S=12`, `r_0=-1606`, `D=394`, and the squared residual is three.
Using `sqrt(12)<=4` gives the rational upper bound `1/197` in (15c).
This is a direct finite-source certificate; it uses no fitted eigenvector
and no long adiabatic ramp.

`preparation.execution_certificate` executes the literal native program
from a fresh blank input. The zero-ramp staggered-product preparation uses
136 pulses. The complete `X` pointer at site zero uses 25 pulses, followed
by six pulses implementing logical `Z` conditional on pointer outcome
zero. There are 167 pulses, 17 work bits and 131072 work-space rows. The
execution retains both pointer outcomes and all 608 native first-failure
slots. The pointer's negative outcome has zero probability on the ideal
prepared product; its slot is retained. It can have positive weight for
the exact interacting ground.

Let `epsilon_native` be the sum of the exact angle and port bounds
returned by the native compiler. This is a complete-isometry bound,
including the first-failure branches. It is independent of the numerical
execution residual. Channel contraction of (15c) through the ideal
complete instrument and the compiler bound give

\[
 {1\over2}\|\rho_{\rm executed}
          -\rho_{\rm ideal\ ground\ instrument}\|_1
       \le {1\over197}+\epsilon_{\rm native}<{1\over100}. \tag{15d}
\]

No success branch is renormalized. The host execution is a floating-point
calculation and its norm and normalization residuals are diagnostic;
the rational residual certificate and analytic angle/port bounds justify
(15d). The source and pointer ideal maps preserve the protected code at
their completed-word boundaries. The stated bound includes native service
failure and does not assert that the literal native circuit equals the
abstract four-state occupancy completion on arbitrary out-of-code inputs.

The finite coupling `J/B=1/400` is not certified to lie below the
unspecified infinite-volume Yarotsky threshold. Equations (15a)--(15d)
certify a complete finite ground-instrument approximation and its actual
native execution; they give no infinite-volume band-residue error or
physical particle assignment. The separate adiabatic resource plan has
its stated finite-ground scope and is not the executed tape in this example.

## 5. Lattice scattering and the speed bound

For sufficiently small positive coupling the ground is unique, the
interaction is translation invariant and finite range, and (8) is a
nonconstant analytic isolated band. Its nonzero Hessian at zero also
shows that its analytic Hessian determinant is not identically zero.
Its zero set therefore has measure zero. Independently, the band-width
estimate in Section 2 gives
`(h-h) intersect Spec(U)={0}`, where `h={(E(k),k)}`: every positive
energy in the difference set is smaller than the ground gap. The only
zero-energy spectral point is the ground at zero momentum.

These verify the regular and the pseudo-relativistic alternatives in
[Bachmann, Dybalski and Naaijkens, Definition 4.1 and Section 4]
(https://arxiv.org/abs/1412.2970). Their wave-operator theorem supplies
lattice scattering states on this sector. The word pseudo-relativistic
denotes the difference-set condition above; it does not assert Lorentz
symmetry. Nontrivial scattering amplitudes and asymptotic completeness
do not follow from the band construction.

For disjoint supports with `d(X,Y)>=1`, the source locality estimate is
independent of the onsite field:

\[
 \|[\tau_t(A_X),B_Y]\|
 \le \|A_X\|\|B_Y\||X|\,2^{-d(X,Y)}
                     (e^{48J|t|}-1).                  \tag{16}
\]

Passing to the product onsite-field interaction picture leaves every
bond support and norm unchanged. The derivation of (16) therefore has
no factor growing with `B`. Since `log(2)>=2/3`, the Euclidean
Lieb--Robinson cone may be taken to have speed `72J`; graph distance
dominates Euclidean distance. Proposition 5.2 of the same paper bounds
the velocity of an isolated smooth observable particle shell by any
such cone speed. Its proof uses the cone commutator decay, which (16)
supplies. Its Section 4 premises have been checked above, so

\[
                        |\nabla_k E(k)|\le72J.          \tag{17}
\]

The exact finite source gap is at least `B/8` by the checked coefficient
substitution in [the finite-field gap proof]
(../spin_twist_closure/FINITE_FIELD_GAP.md), which imports
[Rayudu and Takahashi, Theorem 11](https://arxiv.org/abs/2607.10765v1).
It also holds for the positive-field infinite-volume ground generator:
the local ground limit exists, and the spectral-gap passage is the one
in Bachmann--Dybalski--Naaijkens, Proposition 5.4. In particular

\[
                           E_{\min}\ge B/8.             \tag{18}
\]

## 6. Bounded rest energy under spatial refinement

Consider a family of these sources with lattice spacing `a->0`,
couplings `J_a,B_a>0`, and a fixed ceiling
`J_a/B_a<=lambda_bar<lambda_*`. Use one positive clock
factor `s_a` for the ground generator and all sectors. The physical
dispersion near any carrier momentum `k_a` is

\[
 e_a(p)=s_a E_{B_a,J_a}(k_a+a p),\qquad
 M_a=s_a\min_k E_{B_a,J_a}(k).                          \tag{19}
\]

Suppose `M_a<=Mbar<infinity`; a finite positive limiting rest energy is
a special case. Equations (17)--(18) give

\[
 s_aB_a\le8\overline M,\qquad
 |e_a(p)-e_a(q)|\le72a s_aJ_a|p-q|
            \le576a\overline\lambda\,\overline M|p-q|. \tag{20}
\]

Thus every finite pointwise limit on a connected physical-momentum
window is constant. In particular the family cannot converge there to
`sqrt(M^2+v^2|p|^2)` with `M>0,v>0`. This conclusion requires no
convergence of derivatives and no uniform perturbative Taylor estimate.
If a band is flat, the same conclusion is immediate without the
nonconstant-shell hypothesis used to derive (17).

The locality premise is an observed local or quasilocal field map with
the spatial embedding `x->a x` and the same positive source generator.
The onsite frame changes preserve it. Passing to the even translation
lattice changes a fixed basis convention and does not remove the factor
`a`. A global vacuum-energy shift leaves every excitation energy and
(18) unchanged. A sector-dependent subtraction changes the common
generator and is outside (19).

The exact premise escapes are consequently explicit. A fixed-cutoff
lattice interpretation does not send `a` to zero. A trajectory approaching
another regime can violate the bound on `J/B`. A different dynamical
generator or a spatial identification with no uniform local or quasilocal
control changes the clock or locality premise. Subtracting the rest energy
and magnifying a quadratic band can describe a nonrelativistic sector,
with an altered energy convention. None of these alternatives is supplied
as a relativistic or empirical identification by this construction.

The high-field result is distinct from the declared joint family
`p=n^-2, B=J n^-5`, for which `J/B` diverges. The high-field theorem
does not apply to that family. Its failure to produce a nonzero-speed
massive refinement neither transfers a sharp band to the joint family
nor excludes a sector built under different source premises.
