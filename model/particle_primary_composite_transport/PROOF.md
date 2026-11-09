# Exact transport of physical three-defect composites

Two pairs of odd physical occupation sectors have positive full-space
energy thresholds on the primary eight-cell source. Their localized
singlets also propagate under its unchanged interacting Hamiltonian.
The first nonzero nearest-neighbor transport moment is exactly
\(-3/4\), at order three. It gives a controlled positive short-time
transition probability and distinct complete energy laws for coherent
spatial preparations. An energy band, a relativistic dispersion relation
and a particle rest mass require additional spectral information.

## 1. Source and physical preparation

Fix member zero of the [primary contract](../particle_common_contract/PROOF.md):
all sixty fermionic modes per cell, the actual fifteen-dimensional native
gauge representation, degree-one product link cutoff, all Gauss
constraints, the selected polynomial density \(\rho\), its idle conjugate
reference and source clock \(e^{-itH}\). The Wilson and magnetic
interaction coefficients retain their specified value one. Use the
[physical sector witnesses](../particle_primary_sector_threshold/PROOF.md)
with native occupation increments

\[
 \Delta=\pm(0,2,0,0,1),\qquad
 \Delta=\pm(1,0,0,2,0),                       \tag{1}
\]

in charge order \((-4,-3,1,2,6)\). The first is a weak epsilon
contraction times the charge-six singlet. The second is the color epsilon
contraction of one charge-minus-four and two charge-two entries. The
actual twelve gauge generators annihilate these local vectors. The
signed vectors contain three positive-beta particles or three
negative-beta holes above the fully filled negative-beta anchor.
Their integer-mask exterior chart is related to the frozen source's
lexicographic subset chart by the exact row permutation
`(0,1,4,2,5,7,10,3,6,8,11,9,12,13,14)`. The native matrix representation
is conjugated by this permutation before acting on the witness masks.
The code checks both the source central charges and every matrix entry
against the exterior-minor construction in that chart. This is a basis
isomorphism of the same source, not a change of gauge representation.

Let \(\phi_x\) be a normalized local witness at cube vertex x, the anchor
at every other vertex, and constant functions on every link. All eight
vectors are physical, orthonormal and in the same conserved occupation
sector. Their source mean is

\[
 b=\langle\phi_x,H\phi_x\rangle=-720+9=-711.  \tag{2}
\]

The preparation is an explicit intervention. The paired-memory
construction in the sector proof applies to any normalized physical
target in this sector, including \(\phi_x\) or a coherent combination of
its eight translates. An opposite-sector target is retained in the
second memory; the original \(\rho\) and every reference correlation move
to the first memory by the parity-corrected graded exchange. Neither the
source action nor its time unit changes. No probability calculated for
this prepared target is assigned to the original selected density.

## 2. Exact third moment on the actual cutoff

Every hopping term moves one fermionic occupation across one interior
edge. Electric, onsite and magnetic terms move no matter occupation.
Transferring (1) from distinct vertices x to y requires three net
occupation transfers. Therefore

\[
 \langle\phi_y,H^j\phi_x\rangle=0\quad(0\le j<3).
                                                               \tag{3}
\]

At order three, nonzero transfer requires x and y to be adjacent. All
three factors must be directed hops on that one edge. A nonadjacent
pair has graph distance at least two and requires at least six hops.
Pair creation or a spin-changing hop cannot contribute at minimum
order: the initial and final occupation masks differ in exactly three
annihilated and three created modes, each of which must be used once.

Compression causes no unproved intermediate-state assumption. Let P be
the actual product link cutoff and \(H_{\rm par}\) the complete parent
operator. Each \(H_{\rm par}\phi_x\) is retained: the electric action on
constant links vanishes, onsite action leaves the vector fixed, a single
hop has one native matrix entry, and a magnetic plaquette has one entry
on each of four links. The Gauss support is preserved. Thus

\[
 \begin{aligned}
 \langle\phi_y,(PH_{\rm par}P)^3\phi_x\rangle
 &=\langle H_{\rm par}\phi_y,
                  PH_{\rm par}P H_{\rm par}\phi_x\rangle\\
 &=\langle H_{\rm par}\phi_y,
                       H_{\rm par}^{\vphantom3}H_{\rm par}\phi_x\rangle.
 \end{aligned}                                                  \tag{4}
\]

This identity uses admission of the first actions on both sides. It
does not assert admission of the second action, an invariant small
subspace, or equality of higher moments under compression.

For the surviving one-link term, local gauge invariance removes the
link coordinate pointwise. Conjugating the fermions at one endpoint by
the actual native R(g) changes the bond link g to the identity and leaves
both endpoint singlets fixed. The fully filled negative-beta modes
have determinant one. The complete group is connected, and the exact
local gauge kernels therefore give invariance under every group element.
The integrand is constant, so its normalized Haar integral equals its
identity-link value.

In the actual beta chart the three Wilson blocks obey

\[
 P_+ S_i P_+=-\tfrac12 I_2,\qquad
 P_- S_i P_-=+\tfrac12 I_2.                    \tag{5}
\]

The three extra-particle hops can occur in six orders. Exact fermionic
signs give \(3!(-1/2)^3=-3/4\) between normalized translated singlets.
For holes, moving the occupied negative-beta modes in the opposite
direction supplies the fermionic hole sign and gives the same answer.
The code executes the actual 120-mode two-site signs, all three Wilson
axes, both signed sector pairs and a nonreal nondiagonal native group
element. The pointwise gauge argument proves the full Haar statement;
the group sample is an independent exact control of the chart.

Consequently, with K=H-bI and cube adjacency matrix A,

\[
 \langle\phi_y,K^3\phi_x\rangle=-\tfrac34 A_{yx}
 \quad(x\ne y).                              \tag{6}
\]

More generally, replacing either local singlet by any vector of the
same all-particle or all-hole gauge-singlet spin space gives
\(-3/4\) times their local inner product. Equation (5) transports both
beta spin rows identically. The
[complete local gauge kernel](../particle_primary_composite_spin/PROOF.md)
has dimension six, with two source spin-Casimir values: \(3/4\) with
multiplicity two and \(15/4\) with multiplicity four. The code calculates
the entire six-by-six transport matrix on all three axes and both signed
sector pairs. It is exactly \(-3I_6/4\). Thus both binary-cube components,
G1 and H, are transported, without assigning Lorentz spin to them.

## 3. Controlled full-source transport

The first-action packets in the sector proof give, for all four signed
choices in (1),

\[
 \|K\phi_x\|^2=v=\frac{132301}{735},\qquad
 \langle\phi_x,K^3\phi_x\rangle
       =d_3=\frac{369869271}{42875}.            \tag{7}
\]

The second identity evaluates the exact first-action packet compression:
\(\sum n_{r,k}(105C_r/2+3k)+g\mathbf1^TA\mathbf1\), where
\(g=6/22050\) and A in this formula is the five-loop packet matrix from
the sector proof. The use of this compression computes an inner product
of admitted vectors, without claiming its invariance.

Write \(H=H_0+V\), with the actual dimensionless action at spacing one.
The full local operator bounds give
\(\|V\|\le J=360+6/7\). Exact orthogonal packet norms give

\[
 \|(H_0-b)V\phi_x\|^2=
 \begin{cases}2219349/5&\text{weak pair and singlet},\\
 4433973/10&\text{color triple}.
 \end{cases}                                                   \tag{8}
\]

If R denotes an outward rational upper bound on
\(\sqrt{(8)}+J\sqrt v\), then \(\|K^2\phi_x\|\le R\). The same
bound holds at every cube corner. Taylor's integral formula and
self-adjointness give

\[
 \left|\langle\phi_y,e^{-itK}\phi_x\rangle+
                       \frac{i t^3}{8}\right|
 \le \frac{R^2|t|^4}{24}                       \tag{9}
\]

for adjacent vertices. Indeed the fourth derivative is
\(\langle K^2\phi_y,e^{-isK}K^2\phi_x\rangle\), whose magnitude
is at most \(R^2\). No global spectral truncation enters this estimate.
Since the scalar phase from b changes no probability,

\[
 0<|t|\le\frac{3}{2R^2}
 \quad\Longrightarrow\quad
 |\langle\phi_y,e^{-itH}\phi_x\rangle|^2
                         \ge\frac{|t|^6}{256}>0.        \tag{10}
\]

The certified time bounds are approximately \(4.94489\,10^{-8}\) and
\(4.94553\,10^{-8}\) in the original source time unit. This deliberately
conservative estimate proves actual interacting propagation. It supplies
no useful macroscopic transport rate or extrapolation to late times.

The detector is the physical rank-one projector
\(|\phi_y\rangle\langle\phi_y|\) and its full complementary outcome.
It tests the entire declared bare vector, including constant links and
the other cells; it is not an unproved local position detector. Both
outcomes and both physical memories and the reference are retained.
Equation (10) is its unconditional prepared-state success probability.

## 4. Spatially coherent preparations have different energy laws

For cube character \(k\in\{0,1\}^3\), prepare

\[
 \Phi_k=8^{-1/2}\sum_x(-1)^{k\cdot x}\phi_x.
\]

All eight preparations have the same mean b and centered variance v.
The cube adjacency eigenvalue is \(3-2|k|\), so (6) and (7) give

\[
 \langle\Phi_k,(H-b)^3\Phi_k\rangle
          =d_3-\tfrac34(3-2|k|).              \tag{11}
\]

The four values for \(|k|=0,1,2,3\) are distinct. Their largest
difference is \(9/2\). Hence their complete source-energy probability
laws are different, even though their first two moments agree. These
are moments of the full interacting H on normalized physical
preparations, not eigenvalues of an eight-state Hamiltonian. The bare
compression has \(PHP=bP\), yet (7) gives positive variance and (6)
nonzero transport. Replacing the full dynamics by that compression
would remove precisely the propagation proved here.

The cube character is a preparation label. This open finite cube has no
proved conserved spatial momentum or relativistic mass shell from
(11). The positive sector thresholds, spatial transport and distinct
spectral laws supply three compatible properties of the same physical
sectors. Isolated branches, physical statistics under exchange,
scattering, controlled refinement and a rest-energy calibration have
separate mathematical requirements. The GeV mass table is empty.

## 5. Exact weak-interaction transport coefficient

The interpolation \(H(\lambda)=H_0+\lambda V\) has the primary action
only at \(\lambda=1\). Within either conserved sector in (1),
\(H_0\) has lowest energy b. Its eigenvectors at b have exactly three
defects and constant links. Every other eigenvalue is at least b+6:
the defect number has the same odd parity and is at least three, while
a nonconstant link costs at least \(63/4>6\). A neutral particle-hole
singlet in an unused cell attains the additional defect cost six.
Let \(P_0\) be this entire ground eigenspace and
\(R_0=(H_0-b)^{-1}\) on its orthogonal complement, zero on \(P_0\).
Haar means give \(P_0VP_0=0\).

The third-order effective coefficient between adjacent localized
composites is therefore the exact finite matrix element

\[
 t^{(3)}=\langle\phi_y,V R_0 V R_0 V\phi_x\rangle.     \tag{12}
\]

The full eigenspace is retained in this formula. At minimum transfer
order the one-hop intermediate has electric energy \(e_r=105C_r/2\)
in the first transferred native representation. The two-hop intermediate
component which can return to constant links under the last hop is in
the dual of that last native representation. Its energy is the same
\(e_s\). All other link components have zero overlap after the final
hop. These duals belong to the actual degree-one alphabet. Thus each
of the six constituent orders in the exact CAR calculation receives
the factor \(1/(e_r e_s)\).

For the weak and color triples respectively, (12) is

\[
 \begin{aligned}
 t^{(3)}_{\rm weak}
 &=-\frac{1}{4(189/8)^2}
   -\frac{1}{2(189/8)(63/4)},\\
 t^{(3)}_{\rm color}
 &=-\frac{1}{4(147/4)^2}
   -\frac{1}{2(147/4)(42)}.                  \tag{13}
 \end{aligned}
\]

Both signs of each sector give these same exact coefficients on all
three source axes. The code independently computes the actual CAR
paths with their first and last native denominators, including a
nontrivial native link control. A complete link-coordinate proof uses
the invariant contraction and its uniquely matching last-hop dual,
rather than evaluation of a general resolvent at a single group point.

The finite norm bound separates this decoupled spectral cluster from
its complement for \(|\lambda|<3/J\), with \(J=360+6/7\).
The primary value one is outside that certified interval. Equations
(12) and (13) provide exact coefficients about the decoupled action;
their values are not extrapolated to a primary dispersion or energy.
The tensor-product, isolated-level and smallness requirements of
lattice quasiparticle theorems are not supplied by this coefficient.
