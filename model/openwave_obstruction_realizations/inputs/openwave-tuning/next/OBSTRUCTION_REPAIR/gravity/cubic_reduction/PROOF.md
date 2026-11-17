# Exact G1 reduction of this neutral matter change

The digital ball and the flat improved source are invariant under the24
proper cubic rotations. Their spin lifts act on the first Pauli factor;
the second Dirac factor is retained. Let S_g be a consistently used lift
in the fundamental two-dimensional projective representation G1 and U_g
the spatial permutation times S_g tensor I. The choice of the sign of a
lift cancels between S_g conjugate and U_g.

The matrix units

    P_ij=(2/24) sum_g conjugate((S_g)_ij) U_g

satisfy P_ij P_kl=delta_jk P_il. On each spatial orbit, the script finds an
orthonormal basis B1 for the image of P11 and sets B2=P21 B1. Taking every
orbit retains every multiplicity copy of G1. No parity projection is made,
so both parity sectors are retained. On their combined image the actual
vacuum covariance has the form I_2 tensor C_G1 by Schur's lemma.

The neutral preparation uses ALL four normalized projected radial spinor
seeds. Its momentum envelope is cubic invariant and its normalization
matrix commutes with the spin representation. Consequently its mode matrix
transforms as X(gx)=S_g X(x) S_g^dagger. Its four columns furnish two G1
copies. Both particle and hole subspaces, and hence their covariance
change, lie entirely in the complete G1 isotypic space.

All other irreps are exactly unchanged, not numerically discarded. Therefore

    delta S=2 Tr[h(C_G1+delta C_G1)−h(C_G1)],
    delta S_first=2 Tr[log((1−C_G1)/C_G1) delta C_G1].

The factor2 is the irrep dimension. The script never substitutes a smaller
number of multiplicity copies. `verify_small.py` compares these formulas to
the unreduced492-mode source: entropy-change error2.39e-12, first-derivative
error2.42e-11, matter leakage outside G1 below1.22e-17, and all covariance
intertwining checks pass.

For a large ball the reduced matrix is constructed orbit by orbit directly
from the actual Fourier covariance. The unreduced matrix is never allocated.
Radius9 has12,284 full modes but a1,050-dimensional multiplicity matrix,
17.64MB. The complete run uses35.75s CPU. The saved C_G1.npy is a checkpoint,
not an approximation.

This proof applies to the named symmetric matter preparation. A geometric
coframe perturbation generally changes EVERY irrep and cannot use the
G1-only entropy difference without a separate support theorem.
