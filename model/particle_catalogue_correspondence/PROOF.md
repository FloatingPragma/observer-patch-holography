# What a centre-class catalogue determines

Let the declared internal space be the exterior algebra of
`C^3 ⊕ C^2`, with occupation numbers `n3,n2`. Its diagonal current is

\[
6Y=3n_2-2n_3,\qquad T_3=(n_{w+}-n_{w-})/2,
\qquad Q=T_3+Y.
\]

The code enumerates all 32 basis states and independently recomputes the
fifteen selected states from `native_matter_attachment.representation()`.
The latter occupy exterior degrees two and four. Their auxiliary even
degree is independent of the odd canonical anticommutation operators used
for matter in the interacting Hamiltonian.

For `k=6Y mod 6`, reduction modulo three gives `k=n3 mod 3`; reduction
modulo two gives `k=n2 mod 2`. Thus `k` determines the two centre characters.
The multiplicities are `(4,6,6,4,6,6)`. Class zero contains hypercharges
`-1,0,+1`; class three contains `-1/2,+1/2` and electric charges
`-1,0,+1`. Neither integer hypercharge nor electric charge can be recovered
from the class alone. Complementing all five bits conjugates every charge,
while class three is its own inverse. A class-three label therefore does
not imply a Majorana particle or a self-conjugate representation.

## Full Gauss constraints

On `C^3 ⊗ (C^3)*`, the eight normalized Hermitian generators act by
`T_a X-X T_a`. The program constructs all nine-dimensional matrices and
verifies exactly

\[
\sum_a \operatorname{ad}(T_a)^2
 =3(I-P_1),\qquad
P_1=|I\rangle\langle I|/3.
\]

The singlet has dimension one and Casimir zero; the octet has dimension
eight and Casimir three. Both have trivial centre character. Centre
neutrality therefore cannot replace any of the eight non-Abelian Gauss
constraints. The same issue occurs for weak singlet and triplet sectors.

## A gauge match with different statistics

Two left-quark representations contain

\[
(\mathbf3,\mathbf2)_{1/6}\otimes
(\mathbf3,\mathbf2)_{1/6}
\supset(\overline{\mathbf3},\mathbf1)_{1/3}.
\]

The antisymmetric colour weights and antisymmetric weak singlet are
verified explicitly. This component has the same gauge representation as
the internal `d^c` component. A two-fermion creator has even physical
fermion parity; a single `d^c` creator has odd parity. Their baryon-number
labels are respectively `2/3` and `-1/3` under the usual quark convention.
Consequently, a class-conserving register rule `1+1→2` cannot be identified
with two quarks turning into one antiquark through a parity-preserving
evolution with no further record or excitation. A retained fermionic
reference can carry the missing parity, but belongs to the complete event.
This is compatible with an even diquark channel. Spatial and spin factors
must also satisfy the antisymmetry of its complete two-fermion state.

## Transfer to observations

An internal representation is a useful candidate assignment. A transfer
from a classical centre register to an interacting quantum field must also
supply an operator map for every full gauge current, the physical grading,
the preparation and boundary reference, and their common clock. Complete
spectral and detector records then determine whether a candidate is an
elementary-field channel, a composite, a boundary excitation or a resonance.
The finite gauge identities neither determine its pole energy nor replace
its propagation calculation. The exact finite interacting spectra in
`coherent_defect_wave/chiral_joint` retain their own source domain.
