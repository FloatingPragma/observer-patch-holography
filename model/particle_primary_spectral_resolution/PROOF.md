# Prepared scalar-odd energy windows

The source is member zero of the
[primary contract](../particle_common_contract/PROOF.md), with its
unchanged physical Gauss carrier, selected polynomial preparation,
conjugate reference, action and Hamiltonian clock. A complete intervention
prepares a specified physical scalar-odd vector and retains the displaced
source state in additional idle memory. Its measured energy window has
original probability greater than one quarter. This is a different
experiment from the selected mixed-state scalar detector.

## 1. The physical packets and their residual form

Use the twelve orthogonal physical packets from the
[scalar-odd calculation](../particle_primary_odd_residual/PROOF.md).
The normalized anchor is \(\phi=T\Omega/4\), with onsite energy
\(-714\). It has cubic type \(A_1\), odd spatial inversion and even CAR
parity. The six hopping packets have squared norms \(n_i\), diagonal
onsite/electric relative energies \(d_i\), and include both charge-six
occupation sectors. The five loop packets each have squared norm
\(g=6/22050\). Their relative compressed matrix is \(A\), with
magnetic correction \(C=-N/420\). All these data are regenerated from
that source's actual producers.

For real rational weights \(t_i\) and \(u_j\), define
\[
 v=\sum_{i=1}^6t_i h_i+\sum_{j=1}^5u_j\ell_j,\qquad
 s=\|v\|^2=\sum_i n_it_i^2+g\sum_j u_j^2,\qquad
 q=\sum_i n_it_i+g\sum_j u_j.
\]
The normalized physical vector is
\(\psi=(\phi-v)/\sqrt{1+s}\). Every packet has the same cubic and
inversion type and the same five global native occupations as the
beta-filled constant-link blank \(\Omega\). The blank is orthogonal to
every packet: their charge-six beta-defect count is at least two.
The Rayleigh shift relative to \(-714\) is exactly
\[
 e=\frac{-2q+\sum_i n_id_it_i^2+g\,u^TAu}{1+s}.       \tag{1}
\]

Let \(V\) be the full hopping and magnetic interaction, and let
\(M(t,u)=\|Vv\|_{\rm parent}^2\). The owning exact Haar/CAR calculation
supplies this quadratic form, including the charge-six second-action
Gram \(\left(\begin{smallmatrix}93/2&123/2\\123/2&348\end{smallmatrix}\right)\).
All first-action packets are in the actual degree-one product cutoff.
The projection of \(Vv\) into their span has squared norm
\(q^2+g\|Cu\|^2\). The exact retained part of the residual therefore
has squared norm
\[
 (q+e)^2+\sum_i n_i[1-(d_i-e)t_i]^2
             +g\|\mathbf1-(A-eI)u\|^2 .                \tag{2}
\]
Its orthogonal complement is bounded by
\(M(t,u)-q^2-g\|Cu\|^2\). Parent-to-source compression contracts that
norm. Combining (1)--(2) gives a full-space residual bound for arbitrary
weights; no invariant-subspace assumption or Ritz eigenvector equation
is used.

The exact rational coefficients in the executable are
\[
 t=10^{-6}(17264,24063,14217,18842,31291,28137),\qquad
 u=10^{-6}(5697,9729,4390,6466,13859).
\]
They specify the probe preparation and were chosen from the source
residual form, without measured masses. They do not alter the action.
Global optimality is not asserted. Exact rational arithmetic certifies
center \(c=-714+e\) and variance cap \(\sigma^2\) satisfying
\[
 c\simeq-717.75357,\qquad \sigma^2<49/2,\qquad
 \sigma<4.947619.                                       \tag{3}
\]
The source has some scalar-odd eigenvalue in \([c-\sigma,c+\sigma]\).
Using the same source's certified ground bracket \([L_0,U_0]\) gives
an excited energy in approximately \([1.0666,14.6923]\); the receipt
contains the exact rational outward bounds. This encloses an eigenvalue,
without isolating a unique line or bounding every excitation away from
zero.

## 2. A physical preparation with retained memory

Introduce an identical physical source carrier A with idle Hamiltonian
zero, in the specified ready state \(\Omega_A\). Its occupancy is an
explicit apparatus resource. It is not taken from, or silently substituted
for, the original selected source density \(\rho_S\). The source's
conjugate reference R is retained. The additional memory has the same
physical-dimension upper bound as the source; at this member the declared
allocation uses at most 696 qubits. The exact Gauss dimension is not
computed.

The admitted exact algebraic matrix service implements on A
\[
 W_A=I-|\Omega-\psi\rangle\langle\Omega-\psi|,\qquad
 W_A\Omega_A=\psi_A .                                  \tag{4}
\]
Orthogonality and normalization give
\(\|\Omega-\psi\|^2=2\), proving unitarity. Both vectors are physical,
even, and have identical global native occupations. The operator is
identity outside their span, hence preserves Gauss support, total CAR
parity and those occupation operators. It need not preserve spatial
inversion: the apparatus prepares an odd spatial state.

Apply the graded exchange of the complete source and memory carriers.
On homogeneous vectors it acts as
\[
 S_g(|x\rangle_S|y\rangle_A)=(-1)^{p(x)p(y)}
                                |y\rangle_S|x\rangle_A.
\]
The exchange preserves each native occupation summed over S and A;
individual source occupations can change. This is a complete even
unitary, including every odd sector. Because
\(\psi\) is even, its sign is one on the prepared input. For every
possibly entangled original source-reference density,
\[
 S_g\bigl(\rho_{SR}\otimes|\psi\rangle\langle\psi|_A\bigr)S_g^\dagger
   =|\psi\rangle\langle\psi|_S\otimes\rho_{AR}.          \tag{5}
\]
Tensor factors in (5) are ordered with S and A adjacent. The displaced
state and all original reference correlations survive in A and R. No
trace-out, successful-outcome normalization or measurement of the
reference is performed. Reversing the two unitaries recovers the input.
The exact small-carrier control contains both even and odd sectors and
an entangled reference; it also distinguishes graded from ordinary
exchange on odd-odd inputs.

The source Hamiltonian, Gauss constraints, energy zero and clock are
unchanged. The source marginal is deliberately changed by the stated
intervention. It is not claimed that the selected polynomial density
was this trial, that the original scalar detector created it, or that
this preparation occurs by autonomous source evolution. The source work
mean is controlled directly:
\[
 c-U_0-\epsilon\le
 \operatorname{Tr}H_S(\psi_S-\rho_S)\le c-L_0,
 \qquad\epsilon=1/16.                                  \tag{6}
\]
Here \(\psi_S\) in a trace denotes its rank-one density. A closed,
autonomous energy-conserving controller is not derived. The specified
Householder and exchange are apparatus services, with all input memory
and output registers retained. Any physical compiler supplies its own
complete error and resource bounds; no dense primary matrix, preparation
tape or primary histogram is executed by this receipt.

## 3. Original probabilities and finite energy records

After (5), the source spectral law is the probability measure
\(\mu_\psi(B)=\langle\psi,1_B(H)\psi\rangle\). Its total mass is one,
with mean \(c\) and variance at most (3). Set \(r=23/4\). Then
\[
 \mu_\psi([c-r,c+r])\ge1-\sigma^2/r^2
                         >137/529>1/4.                 \tag{7}
\]
The lower endpoint exceeds \(U_0\), so all energy in this window is
strictly excited. The original excitation window is
\([c-r-U_0,c+r-L_0]\), with the same weight. Every exact spectral vector
is in the reducing physical \(A_{1u}\) subspace. This attaches a complete
prepared experiment to an excited scalar-odd window; it does not attach
the selected mixed-state \(K_{00}\) response to one eigenvalue.

The [complete source energy clock](../particle_primary_joint_spectroscopy/PROOF.md)
is valid for this pure input as well: its tail bound is spectral and does
not require a stationary density. Use energy resolution \(\eta=1/64\),
tail allowance \(1/1024\) and complete implementation allowance
\(1/1024\). Every J-bin output, the source, memory and reference is
retained. The original probability of a reported energy in
\([c-r-\eta,c+r+\eta]\) is at least
\[
               137/529-1/512>1/4.                      \tag{8}
\]
This reported interval also lies above the ground upper bound.
The same joint bound controls the event that the true spectral energy
is in (7) and its recorded value is within \(\eta\). No tail is
conditioned away. The source clock uses its original signed energy chart,
without shifting the zero to c. Its exact and rational Cayley operations
are functions of H and preserve the scalar-odd support.

If the actual preparation/exchange output has half trace distance at
most \(d\) from (5), subtract \(d\) from (8). The receipt supplies the
exact declared matrix service, not a zero-cost or efficient compiled
implementation. This member-zero experiment proves no refinement-uniform
particle, Lorentz dispersion, Standard Model identity or rest mass in GeV.
