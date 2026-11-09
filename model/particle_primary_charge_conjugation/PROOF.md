# Equal spectra and complete records in conjugate source sectors

The construction uses the actual [primary Wilson family](../particle_common_contract/PROOF.md).
The gauge group, all sixty CAR modes per cell, product link cutoff,
polynomial preparation, idle reference and source Hamiltonian clock are
unchanged. It supplies a finite particle/hole symmetry of the full source,
not an identification with a measured particle or a relativistic CPT theorem.

## 1. The native coordinate and spin involutions

Let \(\kappa(g)=\overline g\) on \(G=S(U(3)\times U(2))\).
This is a Haar-preserving group automorphism. In the actual native basis,
\(R=\Lambda^2(\mathbb C^5)\oplus\Lambda^4(\mathbb C^5)\);
every entry is a minor with real coefficients. Thus
\[
 R(\kappa(g))=\overline{R(g)},\qquad \det R(g)=1.       \tag{1}
\]
The determinant follows either from the exterior-power determinant formula
and \(\det g=1\), or from the traceless native generators and connectedness.
In particular the five central charge/dimension pairs obey
\(-4(3)-3(2)+1(6)+2(3)+6(1)=0\). Dropping any charged block can remove
this cancellation. The proof concerns the full native representation.

Pullback by \(\kappa\) is a **linear unitary** on the link-function
space; it does not conjugate wavefunction coefficients. It permutes the
native matrix-entry/conjugate-entry alphabet and preserves word degree,
so it commutes with the actual finite per-link projection. It also preserves
the native trace metric, hence the electric Casimir. Parent products are
formed before this projection throughout.

Use the source Clifford chart, with \(\beta=\gamma_4\), and set
\(M=\gamma_2\), meaning `gamma_data()[1]` in code. It is the real matrix
\[
 M=\begin{pmatrix}0&0&0&-1\\0&0&1&0\\0&1&0&0\\-1&0&0&0\end{pmatrix},
 \quad M=M^T,\quad M^2=I,\quad \det M=1.
\]
For each actual forward Wilson coefficient
\(S_i=(-\beta+\beta\gamma_i)/2\), exact multiplication gives
\[
 -M^T\beta^TM^*=\beta,\qquad
 -M^TS_i^TM^*=S_i^\dagger.                              \tag{2}
\]
These identities involve the complete Dirac matrices, including both beta
signs. Link conjugation without the CAR transformation fails (2)'s purpose.

## 2. Full CAR, Gauss support and the action

On the finite CAR space define the unitary particle/hole transformation
\(C_Fc_{x,\alpha r}C_F^\dagger=\sum_bM_{\alpha b}c^\dagger_{x,br}\).
It preserves the CAR and is unitarily implementable in finite dimension.
Choose \(C_F|0\rangle=|\mathrm{full}\rangle\) in the declared mode order.
There are \(F=60V\) modes. The square has phase
\((-1)^{F(F-1)/2}\det(I_{15V}\otimes M)^*=1\), so \(C_F^2=I\).
It preserves CAR parity, since \(F\) is even. This transformation acts
on every occupation sector, not just the filled-beta anchor.

For a same-site one-body matrix h, CAR gives
\[
 C_F(c^\dagger h c)C_F^\dagger
   =\operatorname{tr}h-c^\dagger M^Th^TM^*c.             \tag{3}
\]
The onsite trace is zero. On different sites, reordering supplies the
minus sign without a trace term. Combining (1)--(3) with link pullback
therefore maps every forward hopping term to its Hermitian adjoint:
\[
 -S_i\text{ transformed by transpose and }M
     =S_i^\dagger,\qquad
 R(\overline U)^T=R(U)^\dagger.
\]
Thus the sum of both orientations is unchanged. Electric terms are
unchanged by their metric isometry, and each real plaquette trace is
unchanged by complex conjugation. Let \(C\) be the combined CAR/link
unitary. It preserves the full parent Hamiltonian and product cutoff.

At every gauge vertex the full matter representation has determinant one.
Particle/hole conjugation of its Fock lift consequently has no extra
determinant character. Together with the link pullback it sends the gauge
action at g to the action at \(\kappa(g)\). Haar averaging gives
\[
 CP_GC^\dagger=P_G,\qquad CHC^\dagger=H.                 \tag{4}
\]
Ghost endpoints and their constant-link physical factors are retained.
No additive energy shift or counterterm occurs in (4).

## 3. Occupations and finite spatial transformations

For each native block of dimension d and every region with V cells,
\[
 CN_qC^\dagger=4dV-N_q,\qquad
 C(N_q-2dV)C^\dagger=-(N_q-2dV).                         \tag{5}
\]
The same calculation at one cell and the charge-trace cancellation give
\(CQ_{6,x}C^\dagger=-Q_{6,x}\). With \(P_\Delta\) the physical
projector for the five deviations from filled-beta occupations,
\(CP_\Delta C^\dagger=P_{-\Delta}\). These reducing sectors have exactly
equal energy spectra, including multiplicities and the same energy zero.
This does not assert that either sector contains a single particle.

For every actual binary cube spin lift S, matrix multiplication gives
\(M^TS^*M^*=S\). Site permutations and link inversions commute with
the corresponding combined transformation, hence C commutes with proper
cube rotations. Under spatial inversion \(\mathcal I\), the spin lift
is beta and \(M^T\beta^*M^*=-\beta\). The full Fock lift therefore obeys
\[
 C\mathcal I C^\dagger=(-1)^{N_F}\mathcal I.              \tag{6}
\]
Odd conjugate states have opposite inversion eigenvalues when those are
defined; even conjugate states have the same eigenvalue. This is a finite
cube statement, without continuous spin or a Lorentz particle assignment.

## 4. One state, clock and complete record law

The source preparation is a real polynomial in H on its physical support,
so (4) gives \(C\rho C^\dagger=\rho\) and \([C,e^{-itH}]=0\).
This unitary symmetry does not reverse time or replace the polynomial
state by a Gibbs state. In particular the original occupation-sector
weights of \(\Delta\) and \(-\Delta\) are equal.

For any complete physical instrument \(\{K_r\}\), its conjugate menu
\(K'_r=CK_rC^\dagger\) is complete. All source waiting times and finite
energy-clock operators, being functions of H, intertwine with C. Each
full history operator therefore obeys \(L'_r=CL_rC^\dagger\) and
\[
 \operatorname{Tr}(L'_r\rho L_r'^\dagger)
       =\operatorname{Tr}(L_r\rho L_r^\dagger).             \tag{7}
\]
Failures and all energy bins are included. Opposite charge-transfer
labels can be attached by the [coherent meter](../particle_primary_charged_readout/PROOF.md),
with every original probability retained. For an explicitly prepared
vector psi, C psi has exactly the same complete energy distribution.

The idle reference is not discarded or transformed during either readout.
For the canonical purification, \((C\otimes\overline C)\Psi=\Psi\),
so the two complete outputs are related by this unitary on their source
and reference descriptions. This is a comparison identity, not an extra
reference operation. For an arbitrary correlated input, the corresponding
identity requires conjugating its source input too; invariance is not
asserted for arbitrary source-reference states. Source-only conjugation
need not fix the canonical vector.

Equations (4)--(7) supply actual symmetry and opposite-sector energy
matching at every finite capacity. Nonzero detector access, a propagating
branch, relativistic dispersion, continuum spin and observed antiparticle
identity require their own constructions. No numerical primary spatial
matrix or GeV rest-mass spectrum is computed here.
