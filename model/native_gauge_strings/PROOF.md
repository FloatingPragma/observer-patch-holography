# Gauge-invariant propagation of one dressed fermion

The configured Hamiltonian of
[native quantum Gauss constraints](../native_quantum_gauss/PROOF.md) has
an exact spatial reduction on a rooted tree. Its complete one-fermion
Gauss sector consists of one contracted matter-and-flux state at every
vertex. This includes the root. The reduction uses the literal sixty-state
link coordinates and the owning Hamiltonian's electric projector, with
all five matter blocks retained before imposing the constraints.

The graph, link register, coherent coupling, static dual reference and
clock are specified inputs. The result describes a dressed excitation
bound to that reference when the string energy grows with distance.
It makes no assignment to an observed particle or asymptotic scattering
state.

## The full configured carrier

Let \(\mathcal T\) be a finite rooted tree. Orient each edge from child
\(x\) to parent \(y\), and retain the actual native matter module
\(M=\bigoplus_s V_s\), of dimension fifteen. The owning twelve
anti-Hermitian generators give five blocks. In their actual central-current
normalization their eigenvalue/dimension pairs are

\[
 (-3,1),\quad(-1,3),\quad(-1/2,6),\quad(3/2,2),\quad(2,3).
\]

These numbers label the source representations, without an empirical charge
calibration. Every nonvacuum block has nonzero central charge. The exact
commutator equation on \(\operatorname{End}(V_s)\) has rank
\(d_s^2-1\) for all twelve actual generators. Its kernel therefore has
dimension one. Because the representations are unitary, a proper invariant
subspace would give a nonscalar commuting orthogonal projector. Each
\(V_s\) is consequently irreducible. Distinct central eigenvalues make
\(\operatorname{Hom}_G(V_s,V_t)=0\) for \(s\ne t\), including the two
three-dimensional blocks. This verifies representation identity rather
than inferring it from dimensions.

Every edge retains

\[
 \mathcal L_{xy}=\mathbb C|0\rangle\oplus
          \bigoplus_s(V_{s,x}^*\otimes V_{s,y}),\qquad
 \dim\mathcal L_{xy}=60,
 \qquad U^s_{ab}=|\bar a,b;s\rangle\langle0|.
\]

Fix one block \(r\) and one static dual reference \(V_r^*\) at the root.
The entire one-fermion Hilbert space is

\[
 \mathcal H^{(1)}_r=
 \left(\bigoplus_{v\in\mathcal T}M_v\right)\otimes V_r^*
               \otimes\bigotimes_{e\in E(\mathcal T)}\mathcal L_e.
\]

Its dimension is \(|\mathcal T|\,15d_r60^{|E|}\). No spatial spin or
family register is part of this owning model. The sector fixes total
fermion number to one, a conserved quantity of its Hamiltonian; it makes
no assertion that other Fock sectors are absent from the full theory.

At each vertex, \(G_v(T)\) is the matter current plus the incident link
currents, with \(-T^T\) at a child endpoint and \(T\) at a parent
endpoint, and the dual-reference current at the root. These are the same
local operators whose full-Fock Ward identity is proved in the owning
package. The Hamiltonian is

\[
 H=J\sum_{(x,y),s,a,b}
       \big(c^\dagger_{x,a}U^{s,xy}_{ab}c_{y,b}+\text{adjoint}\big)
       +\epsilon\sum_e(I-|0\rangle\langle0|)_e,
 \qquad J\in\mathbb R,\quad\epsilon\ge0.
 \tag{1}
\]

The electric term in (1) is the native quantum-Gauss occupation-projector
law. The separate gauge-Higgs source has a Casimir-weighted electric law;
on the string in block \(r\), that law would give its block energy times
the path length. The coefficient \(\epsilon\) in (1) is not silently
replaced by that different configured energy. Every one of the 1,800
matter/link input columns on one edge is checked against the owning
Hamiltonian implementation. Both operators act as the identity on the
reference, so this comparison extends to its full factor and untouched
ancillas.

## Exhaustion of the Gauss kernel

Fix the mobile fermion's position \(v\). Local Gauss operators preserve
that position and each link's block. The invariant space is thus an
orthogonal sum over complete matter/link block assignments.

Consider a leaf other than \(v\) or the root. It has no matter or
reference. A nonvacuum incident link would contribute a nonzero scalar
central charge, so the leaf's Gauss equation forces that link to be in
its vacuum. Remove the leaf and repeat. Every branch off the unique
root-to-\(v\) path is removed. This argument also handles \(v\) equal to
the root, when every link is forced to vacuum. A tree has no residual
cycle on which a closed flux could survive.

Along the path, an intermediate vertex has only two nonvacuum endpoint
factors. Its invariant tensors are exactly
\(\operatorname{Hom}_G(V_s,V_t)\). The actual-current calculation above
forces the two blocks to coincide and leaves a single contraction. The
root reference selects \(r\), and the final matter factor must be
\(V_r\). In a fixed block assignment, the representation factors at
different vertices can be reordered into a tensor product of vertex
spaces. Their invariant subspaces therefore factor. The one-dimensional
contraction at each vertex exhausts the global kernel; it leaves no
unexamined intertwiner or flux multiplicity.

If the path is \(v_0=0,v_1,\ldots,v_k=v\), define

\[
 A_v=\sum_{a_0,\ldots,a_k=1}^{d_r}
 |v,a_k\rangle|\bar a_0\rangle_{\rm ref}
 \bigotimes_{j=1}^{k}|\bar a_j,a_{j-1};r\rangle_{(v_j,v_{j-1})}
 \bigotimes_{e\notin\text{path}}|0\rangle_e.
\]

Then \(\|A_v\|^2=d_r^{k+1}\), and different positions are orthogonal.
The complete Gauss kernel of \(\mathcal H^{(1)}_r\) is

\[
 \mathcal K_r=\operatorname{span}\{A_v:v\in\mathcal T\},
 \qquad \dim\mathcal K_r=|\mathcal T|.
 \tag{2}
\]

The executable example uses a four-vertex branched tree. It exhausts
every central matter/link block assignment, checks every actual local
Gauss generator on every contracted state, and compares the full ambient
Hamiltonian output with the claimed reduction. The general exhaustion
in (2) follows from leaf elimination and the exact irreducibility proof,
without a numerical rank extrapolation from that example.

## The same Hamiltonian propagates the dressed object

For an edge from child \(x\) to parent \(y\), raising the vacant final
link sends \(A_y\) to \(A_x\). Its adjoint sends \(A_x\) to
\(d_r A_y\), because the child index is summed. Edges with no matter at
either endpoint act as zero. Moving back toward the root annihilates the
last occupied link; moving into a child appends a link. Thus no term in
(1) leaves the complete kernel (2).

In the normalized basis
\(|v;r\rangle=d_r^{-(\operatorname{dist}(0,v)+1)/2}A_v\),

\[
 H\big|_{\mathcal K_r}
  =J\sqrt{d_r}\,\operatorname{Adj}(\mathcal T)
     +\epsilon\operatorname{diag}(\operatorname{dist}(0,v)).
 \tag{3}
\]

The factor \(\sqrt{d_r}\) follows from the unnormalized literal raising
operator. The normalized partial-isometry transport in
`model/gauge_defect_transport` uses a different normalization and is not
substituted for (1). The code retains exact rational coordinates in
\(A_v/\sqrt{d_r}\), whose Gram matrix is
\(\operatorname{diag}(d_r^{\operatorname{dist}(0,v)})\), and checks
\(GH=H^\dagger G\). Dropping this Gram factor changes transition
probabilities.

The position projectors \(P_v\) sum to the identity on the entire
one-fermion carrier, commute with every local Gauss generator, and retain
all positions including the root. They restrict to the complete rank-one
menu on (2). The declared Cayley step

\[
 U_\tau=(I+i\tau H/2)^{-1}(I-i\tau H/2)
\]

is unitary and preserves the energy of that same \(H\). The certificate
retains all position histories, including repetitions of the previous
position and zero-probability outcomes. Readout generally changes the
hopping energy. At every round its work is the mean energy after readout
minus the mean energy before readout, and these terms telescope to the
total energy change. A physical detector reservoir or reset-cost model
is not supplied.

## Volume control and interpretation

For an additional locally finite infinite rooted tree of bounded degree
\(D\), the diagonal operator \(H_0=\epsilon\operatorname{dist}(0,v)\)
has compact resolvent when \(\epsilon>0\). Every finite-radius ball is
finite, and the potential tends to infinity outside those balls. The
hopping perturbation has norm at most \(D|J|\sqrt{d_r}\). It follows
that (3) is self-adjoint on the diagonal domain and has compact resolvent.
This produces discrete reference-bound levels in this specified graph
model; it does not construct freely propagating asymptotic particles.

For continuous Hamiltonian evolution from support inside radius \(R_0\),
compress to the ball of radius \(R\ge R_0\). Leaving it requires at
least \(m=R-R_0+1\) hopping insertions. The interaction-picture Dyson
series gives, for total evolution time \(t\),

\[
 \|\psi(t)-\psi_R(t)\|
 \le 2\sum_{n\ge m}\frac{(D|J|\sqrt{d_r}\,t)^n}{n!}.
 \tag{4}
\]

For several position readouts, expand all evolution intervals jointly in
the Dyson series. The diagonal record isometries preserve radial support
and norm, including with untouched references. Terms with fewer than
\(m\) total hops coincide in the full and compressed experiments. The
ordered insertion volumes, summed over all distributions of \(n\) hops
among the intervals, equal \(t^n/n!\), where \(t\) is the total evolution
time. Thus (4), capped at one, bounds the total variation of the entire
position-record law, without multiplying a one-time estimate by the number
of readouts. It concerns the same initial state with support in the stated
ball; a ground-state preparation requires its own approximation bound.

The Cayley experiment is a specified discrete evolution. For a finite
ball its comparison with the continuous evolution has the separate bound
\(\|U_\tau-e^{-i\tau H}\|\le|\tau|^3\|H\|^3/12\), obtained from
\(|2\arctan z-2z|\le2|z|^3/3\). The per-step bounds add for complete
histories. A growing ball with a fixed step therefore receives no
automatic continuous-time limit from (4).

The derivation selects the complete invariant sector of a configured
source. It supplies no choice of a physical spatial graph or refinement,
no relativistic dispersion, chiral regulator measure, source-selected
coupling or laboratory energy scale. Extra matter can create additional
intertwiners and screen the string; cycles admit closed flux. Those
carriers require separate kernel analyses. Equation (3) is not extended
to them by omitting their states.
