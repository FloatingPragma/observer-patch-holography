# Event cost, family multiplicity and canonical quantum-link dressing

The [executable certificate](families_chirality.py) establishes two finite
connections. The squared displacement of the existing native near-rotation
service selects the rank-three slow band without a supplied dimension
window. A unitary on the complete fermion/link carrier supplies canonical,
gauge-covariant dressed annihilation operators, avoiding the failure of a
linear operator-valued projection. Neither result identifies three physical
generations or an interacting chiral continuum.

## Reused constructions

The [native family selector](../common_source/native_family_selection.py)
constructs a 36-state protected pair from the actual twelve-port incidence.
Its complete spectrum, Reynolds intertwiner and recorded cooling instrument
are reused directly. Its neutral family pair is separate from the pair
carrying gauge currents. A rank-three spatial frame is not identified with
three generations by this construction. The [finite family theorem](../../Lean/Screen/A5FamilyBand.lean)
proves the exact coefficient-space projectors and their strict Laplacian
cost order. The [family attachment certificate](../../code/a5_closure/family_band_attachment_certificate.py)
binds the finite cost and pole-residue interfaces and distinguishes a
rank-45 tensor attachment from physical generation identification.

The [native quantum-link attachment](../standard_model_source_map/gauge_attachment.py)
supplies all twelve currents, the 60 link levels and their endpoint action.
The [interacting chiral certificate](../particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md)
supplies the exact raising-link counterexample to naive projected CAR. The
[overlap regulator](../common_source/overlap_regulator.py) supplies the
gauge-covariant Euclidean construction for admissible classical links,
including a free single-zero symbol and a volume-independent ideal locality
tail. Its classical-link hypotheses do not hold for the raising links.
The [chiral electric construction](../interacting_extension/chiral_electric/PROOF.md)
supplies a connection and the positive geometric correction on a
fixed-rank abelian gauge-coordinate patch. These results are compatible
within their separate domains; none licenses replacing quantum links by
classical unitary matrices or suppressing mirror records.

The [spatial-overlap and quantum-flux construction](../interacting_extension/chiral_lattice/PROOF.md)
also supplies a conditional coordinate completion. Its isometry sends the
native vacuum to the constant function and each native matrix-coefficient
level to \(\sqrt{d_r}D^r_{ab}(g)\) in
\(L^2(S(U(3)\times U(2)),\mathrm{Haar})\). Multiplication on that space
has additional representation transitions; it is a different completion
from the finite raising operator. At fixed finite spatial volume, its
electric energy plus bounded spatial-overlap and plaquette terms gives a
self-adjoint gauge-covariant vectorlike Hamiltonian. It supplies neither
uniform locality near Wilson-kernel gap closures nor a selected chiral
state. This route can supply dynamical unitary coordinate links if the
source selects the completion. The exact finite charged-shift trace
obstruction does not prohibit that infinite coordinate realization.

## 1. Per-event cost selects the finite band without a dimension window

Let \(R_g\) be the 36-dimensional permutation representation of a native
proper rotation. The incidence rule selects the twelve order-five rotations
that carry ten ports to adjacent ports. This set is closed under inverses.
Write \(C=\sum_gR_g\), and sample its twelve members uniformly. For every
coefficient vector \(v\), including an arbitrary untouched tensor factor,

\[
 \mathbb E_g\|(R_g-I)v\|^2
 =\langle v,\mathcal C v\rangle,
 \qquad \mathcal C=2I-\frac C6.                         \tag{1}
\]

This follows by expanding each squared norm, using unitarity and
\(C=C^\dagger\). The code also recomputes the left side from all twelve
complete matrices. The existing coefficient isometry \(W:\mathbb C^{12}
\to\mathbb C^{36}\) satisfies

\[
 W^\dagger W=I,\quad CW=W(2A+2I),\quad
 W^\dagger\mathcal CW=\frac{5I-A}{3}.                 \tag{2}
\]

Thus the finite theorem's Laplacian comparison is the same quadratic
comparison, up to a positive factor. The full 36-state cost spectrum is

| Event cost | Multiplicity |
| --- | ---: |
| \(0\) | 2 |
| \((5-\sqrt5)/3\) | 3 |
| \(2\) | 20 |
| \((5+\sqrt5)/3\) | 3 |
| \(5/2\) | 8 |

The strictly smallest positive eigenspace has dimension three and is
exactly the retained native cooling projector. This statement uses neither
the interval \([3,5]\) nor a preselected three-dimensional coordinate
frame. The two zero-cost invariant states are retained. Excluding them is
an explicit non-invariant-response condition; unconstrained cost
minimization selects the invariant space.

The incidence-defined near-rotation rule is also a specified choice. If all
60 rotations are sampled uniformly, the same calculation gives
\(\mathcal C_{\rm all}=2(I-P_{\rm inv})\), with one positive cost of
multiplicity 34. Symmetry by itself therefore does not select the
near-rotation service or its rank-three minimum.

Equation (1) is the cost of the specified rotation/probe service. It is
not an identification of a squared Hilbert-space displacement with heat,
elapsed time, a particle mass, or the number of passive C/E repair events.
If this service is triggered with a state-independent probability
\(q>0\) per original attempt, retaining the no-trigger result gives
\(q\mathcal C\). The ranking is unchanged, but the physical trigger and
its clock require their own source attachment.

The existing complete cooling law is

\[
 \mathcal E_P(\rho)=P\rho P+
 \operatorname{Tr}[(I-P)\rho]\frac P{\operatorname{rank}P}. \tag{3}
\]

It preserves every coherent state in the selected band, with the reset
labels and complementary input information retained in the record/bath.
This is an added controller. Corrected passive C/E branches on the code
have the form \(\sqrt{w_h}I\), so their complete unread map is the identity.
In particular it fixes \(I_{36}/36\), whereas (3) sends that state to
\(P/3\). The two channels cannot be identified.

The selection requirement is concrete: derive the choice of the
non-invariant minimum-cost controller from the full observer/repair/memory
program, then identify its protected output algebra with the multiplicity
factor of the propagating matter currents. Covariance and complete
probability alone allow (3) for other invariant projectors, including the
rank-20 eigenspace. The external-copy ambiguity also survives:

\[
 \operatorname{Tr}[(\rho\otimes I_n/n)(A\otimes I_n)]
 =\operatorname{Tr}(\rho A).                           \tag{4}
\]

The same identity holds for every source operator word. Source observables
that act trivially on an added copy factor cannot determine its dimension.
An actual matter attachment must break that ambiguity; a rank-three
carrier calculation alone does not do so.

## 2. Canonical CAR through a complete gauge-invariant unitary

Select the charged singlet component of the existing fifteen-state matter
representation at two endpoints, with annihilators \(c_x,c_y\). Retain
all four occupation states and all 60 link levels. The active link is the
existing raising matrix \(U=|e,e\rangle\langle0|\), not a classical
unitary. Define the Hermitian dressed exchange

\[
 S=c_x^\dagger Uc_y+c_y^\dagger U^\dagger c_x,
 \quad S^3=S,\quad
 W=e^{i\pi S/4}
   =I+(1/\sqrt2-1)S^2+iS/\sqrt2.                      \tag{5}
\]

The 240-dimensional operator is identity on the vacuum, double occupation
and inactive link states. Every component of the twelve native endpoint
Gauss generators commutes with \(S\), hence with \(W\). The code checks
all 24 identities using the complete link generators, including their
action on the other 58 levels.

Set \(b_i=Wc_iW^\dagger\). Unitary conjugation gives

\[
 \{b_i,b_j\}=0,\qquad
 \{b_i,b_j^\dagger\}=\delta_{ij}I.                   \tag{6}
\]

These are identities on all occupations and link levels. Since
\([G,W]=0\), the dressed operators retain the same endpoint charge
transformation law. Their coefficients include occupation/link dependence;
they are not replaced by a linear one-particle expression
\(\sum_jP_{ij}(U)c_j\). That linear expression has the additional
commutator terms established in the retained counterexample.

This is a canonical finite construction using the same link and current
algebra. The angle \(\pi/4\) is a declared demonstration control. It is
unrelated to the charged-lepton Koide phase and contains no measured mass
input. No native pulse tape is asserted executed for (5), and no physical
time calibration is inferred from a dimensionless pulse area.

## 3. The dark-complement rank cannot be hidden

On the one-fermion space, let \(P_0=|x\rangle\langle x|\otimes I_{60}\).
The canonical projector \(P=W_1P_0W_1^\dagger\) has rank 60. For the
retained dressed-edge matrix

\[
 A=\begin{pmatrix}0&iU\\-iU^\dagger&0\end{pmatrix},
 \qquad P_{\mathrm{old}}=(A^2-A)/2,
\]

the exact rank is one. Direct calculation gives

\[
 P=P_{\mathrm{old}}+P_{\mathrm{dark},x},\qquad
 \operatorname{rank}P_{\mathrm{dark},x}=59.            \tag{7}
\]

A unitary cannot turn rank 60 into rank one. The rank-one sign projector
therefore cannot be identified with one canonical fermion mode tensored
with the unchanged 60-level link in this construction. Equation (7)
exhibits the required completion rather than deleting it. This rank
argument concerns this tensor-factor route; it is not a theorem against
every gauge-dependent fermion bundle.

For \(K=I-2P\) and \(\gamma_5=2P_0-I\), the finite operator
\(D=I+\gamma_5K\) satisfies the Ginsparg-Wilson identity exactly. Its 118
dark one-particle columns have \(D=0\). It is consequently a sign and CAR
control, not a relativistic kinetic regulator or a species count. The
displayed \(\gamma_5\) is an endpoint grading; it is not identified with
spacetime chirality. A
unitary change of variables preserves a parent spectrum and cannot remove
the naive lattice action's unwanted low modes.

## 4. Keeping and controlling the complementary occupations

The canonical heavy occupancy \(N_H=b_y^\dagger b_y\) is a projector
that commutes with every Gauss generator. A declared positive parent
\(H_0=I-S\) and penalty \(\Lambda>1\) give

\[
 H_\Lambda=H_0+\Lambda N_H.                            \tag{8}
\]

All 240 states are retained. In the \(N_H=0\) subspace, the compressed
parent is the identity and the complementary mixing has norm one. The
complete spectrum is

\[
 1\ (119),\quad 1+\Lambda\ (119),\quad
 1+\Lambda/2\pm\sqrt{\Lambda^2/4+1}\ (1\text{ each}).
\]

The shift of the lowest branch from the compressed value obeys

\[
 0\le1-E_-\le\frac1{\Lambda-1}.                      \tag{9}
\]

This is a finite gauge-covariant, canonical complement-control example,
with exact checks at \(\Lambda=2,4,8\). The penalty changes the
Hamiltonian and is supplied. It neither derives a mirror gap of the
spatial source action nor selects physical masses.

## 5. Minimal operator data for the interacting chiral connection

For a gauge-covariant, fixed-rank projector path \(P(t)\), Kato transport
has anti-Hermitian one-particle generator \([\dot P,P]\). A path from a
constant spin projector must preserve its ranks in each gauge-charge
sector, not only the total rank. On a common quantum-link carrier the
many-body generator can be formed as
\(\widehat K(t)=\sum_{ij}c_i^\dagger K_{ij}(t)c_j\).
Its unitary transport is canonical on full Fock space. Its restriction to
one fermion realizes the specified projector transport. Gauge covariance
of \(K\) makes \(\widehat K\) commute with the total Gauss action.
Higher-occupation terms generated by conjugation must be retained.

The executable source must select such a local projector path and
controller, its compatible electric connection and positive geometric
term, and a kinetic/complement Hamiltonian with bounds uniform in volume
and refinement. It must also supply the anomaly-compatible measure and
positive physical state on the same operational clock. The admissible
classical overlap construction and abelian electric patch are useful
starting points for these operators; their hypotheses do not prove that
the complete nonabelian quantum-link source supplies them. The alternative
coordinate completion requires its own source-selection proof and
uniform control of the gauge configurations that approach kernel-gap
closures. No target
particle masses or mixing angles enter the finite results above.

Run the exact certificate and tests with:

```sh
python -m model.particle_calibration_selection.families_chirality
python -m pytest -q model/particle_calibration_selection/test_families_chirality.py
```
