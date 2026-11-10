# A prepared stationary background of the joined source

Let `G` be the ten retained link words, `a` the sixty oriented address values,
and `j` any of the three internal triplet components. In the exact same
joint source and tick `U` as the parent construction, define

\[
 \Omega_j={1\over\sqrt{7200}}
 \sum_{g\in G}\sum_{p=0}^{11}\sum_{a=0}^{59}|g,p,j,0,a\rangle.
\]

## Stationarity and constraint

The address reflection fixes the uniform address. All central links act as
the identity on the triplet. A seam swap fixes the uniform port amplitude
with blank parity memory. Thus every `W_g` fixes its wave factor. For each
wave port, the conditional interior move permutes `G`, so all five controls
fix the uniform link factor. Consequently `U Omega_j = Omega_j`, including
the actual reverse tick. The interior Gauss involution only permutes links
on the triplet, hence `Gamma Omega_j = Omega_j`. The three states are
orthonormal; every normalized superposition of them is also fixed.

This is a degenerate pure stationary background. It is not a claim of a
unique vacuum, a minimum of an independently defined Hamiltonian, or absence
of boundary flux. Its odd regional flux and compensating exterior are the
specified boundary sector.

## Complete finite preparation

Start from the definite initial link word, uniform ports and addresses,
blank parity memory and a chosen triplet component. These wave preparations
are declared apparatus services. The complete interior Gauss measurement
has probabilities `1/2,1/2`; retain the negative result. The normalized
positive link state is `(ket(g0)+ket(g0 xor star))/sqrt(2)`.

The connected ten-word move graph supplies a finite word `V_k` of admitted
interior moves taking `g0` to `g_k`. Every move commutes with the star flip,
so it also takes `g0 xor star` to `g_k xor star`. Order the words with five
representatives first and their star partners second. Prepare a uniform
ten-level ancilla, apply `V_k` controlled by its label, inverse Fourier
transform that ancilla and read **all ten** outcomes. For outcome `ell`,
the unnormalized link state is

\[
 {1\over10\sqrt2}\sum_{k=0}^9 e^{-2\pi i\ell k/10}
 (|g_k\rangle+|g_k\mathbin{\mathrm{xor}}\star\rangle).
\]

Odd outcomes have probability zero. Each even outcome has probability `1/5`.
Outcome zero is exactly the uniform ten-word link state after normalization,
so it prepares `Omega_j`. Its success probability is `1/5` conditional on
positive Gauss preparation, or `1/10` per initial fresh trial. The complete
instrument retains all failed outcomes. Fresh independent retries have
failure probability `(9/10)^k`; a retry requires a fresh input, not the
unjustified reuse of a disturbed one. Exact coefficients are checked in the
cyclotomic field with polynomial `z^4-z^3+z^2-z+1`.

The controlled words are compiled from the parent endpoint-controlled
repair, rather than assumed to be unconditional idle operations. Copy the
wave port into a blank twelve-level scratch register by modular addition;
conditionally exchange that port with port 5; apply the specified `C_e`;
undo the exchange and the copy. Every allowed interior edge contains port
5, so this acts as `T_e` on the link, returns the original wave port and
cleans the scratch. On arbitrary scratch inputs each substep remains a
permutation. All `5*10*12*12=7200` basis cases are checked. These substeps
preserve the interior constraint on the triplet preparation domain. The
claim does not extend that substep property to the doublet.

The coherent ancilla, conditional routing, port scratch and Fourier readout
are explicit optional apparatus primitives. Idle U alone is not claimed to
prepare its own background. Pure-state initialization is also a supplied
preparation resource; no entropy is discarded in the retained instrument.

## Stationarity is not positive energy

For any finite Hamiltonian `H` and a faithful stationary density `rho` on
the full matrix algebra, the two-sided transition frequencies occur with
positive weight whenever `H` has distinct eigenvalues. To see this, choose
joint eigenvectors with `rho_i>0` and energies `E_i != E_j`. The matrix units
`|i><j|` and `|j><i|` have opposite frequencies under `A -> exp(itH) A exp(-itH)`;
their squared norms in the state are respectively `rho_j` and `rho_i`, both
positive. Thus the faithful finite reference cannot itself be a positive-energy
vacuum on the full algebra unless the dynamics is trivial there. This does
not rule out a pure reference, a restricted algebra or a controlled limit.

Even a pure fixed vector of U does not select the logarithm of U. The
[energy-branch control](../spectroscopy/PROOF.md) gives positive generators
with the same ground vector and identical ticks but different gaps. The
present preparation solves a stationary-state input; the physical energy
and relativistic vacuum obligations require their own source specification.

## Connected preparation and excitation readout

The [tensor-wave spatial construction](../spatial/PROOF.md) fixes identical
products of these backgrounds. It supplies exact local records independent
of remote boundaries outside the finite causal cone. This property belongs
to that explicitly declared circuit, not automatically to a CAR occupation
extension.

The [same-U spectroscopy](../spectroscopy/PROOF.md) measures all five faces
on `Omega_0`, applies both branches `(I+/-U^2)/2`, then retains all four
Fourier-clock bins. The active branches are orthogonal to the stationary
background and carry total probability `283/7200` conditional on its
preparation. Their finite frequencies are not rest masses.
