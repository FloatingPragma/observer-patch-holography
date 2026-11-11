# Translation and energy readouts of periodic defects

The geometry is an even cubic supercell of side L≥6, repeated periodically.
Each cell contains the two compact opposite-degree cores at (0,0,±1).
The on-site term is M n(x)·σ, with M≥8 and the normalized direction from
`model/compact_defect_pair/topology.py`. Every forward nearest-neighbour
bond is iσ_j/2. A bond crossing a cell face carries exp(iθ_j); its reverse
is the adjoint. The texture, bond strength and clock are common to every
cell. This periodic member differs from the open finite compact source.
The scalar and interacting native fermion sectors are absent from this
quadratic propagation control.

## Complete band count

Write h(θ)=h_0+T(θ). The on-site operator has four zero eigenvalues, one
positive and one negative eigenvalue of magnitude M on every other site.
Each directional hopping has norm at most one, since it is a Pauli matrix
times the imaginary part of a unitary shift. Hence ||T||≤3 for every θ.
The ordered-eigenvalue perturbation inequality gives exactly four roots
in [-3,3], with every other root of absolute value at least M−3>3.
This counts the entire isolated band at every Bloch phase. It neither
requires a zero root nor supplies a count of physical generations.

Discrete Fourier transformation of the supercell coordinate gives the
blocks h(θ) for θ_j=2πr_j/N_j in a ring of N_j cells. This is a unitary
change of basis on the complete source. The executable separately constructs
every edge of a two-cell ring and checks the intertwiner on all rows.
At general phase it diagonalizes every row of the selected block, retaining
its complete complement as well as the four low roots. The selected phases
0 and ±π/2 belong to four-cell rings. The calculation does not enumerate
every momentum of a single three-dimensional ring.

For a homogeneous texture n=(0,0,1), Fourier transformation on individual
sites gives the independent full-spectrum control

\[
 E_\pm(q)=\pm\sqrt{\sin^2q_x+\sin^2q_y+(M-\sin q_z)^2},\qquad
 q_j=(2\pi m_j+\theta_j)/L.
\]

The shift in the third momentum component is characteristic of this Pauli
texture. An on-site two-component Pauli term is not an anticommuting Dirac
mass, as established by the complete Clifford test in
`model/particle_candidate_atlas/DEFECT_FIELD_JOIN.md`.

## Complete observer instrument

Fill all negative-energy eigenmodes of the quadratic Hamiltonian. For each
normalized eigenmode a use both Kraus branches c_a†/sqrt(d) and c_a/sqrt(d),
where d=2L³. The canonical anticommutation identity makes the sum of their
effects I on the full Fock space. In the filled sea exactly d branches have
weight 1/d; the other d have zero weight. Each nonzero branch has positive
energy transfer |λ_a|. Addition carries the mode's crystal momentum and
removal its negative, relative to the filled sea. A dual probe reference
carries the opposite conserved labels. The orthogonal classical record
contains the branch and energy, including zero-probability branches in the
instrument definition. No conditioning on the four low modes is performed.
A direct sum over the selected blocks with equal prior weight supplies a
complete finite comparison experiment. This does not replace the spectrum
of the interacting defect, scalar and thirty-mode native source.

The exact joint spectral readout is an instance of the commuting energy and
supercell-translation algebra used in
`model/particle_joint_translation_readout/PROOF.md`. Position and momentum
need not possess simultaneous sharp records for this instrument to exist.
No source or ground-state transfer from that separate gauge program is used.

## Rest-energy test

For each positive band calculate the secant

\[
 c_a^2(\theta,j)=
 \frac{E_a(\theta e_j)^2-E_a(0)^2}{(\theta/L)^2}.
\]

A positive isotropic relativistic branch with its rest point at zero
momentum requires positive secants approaching a common limiting speed
squared. The M=8 controls contain negative secants and strong directional
variation. The M=64 bands are much flatter; values comparable to numerical
residuals do not support a sign claim. The receipt labels the secants as
floating diagnostics, not interval-certified derivatives. The calculation
supplies neither a relativistic refinement nor a physical rest-mass
assignment to the zero-momentum energies. It excludes treating those
energies as sufficient evidence for that assignment. A band minimum at a
different momentum, a moving texture, or an interacting propagator requires
its own calculation.

Changing the source energy unit multiplies every energy by the same positive
number. It cannot remove negative secants, their anisotropy, or alter mass
ratios. No per-band GeV factor is introduced.
