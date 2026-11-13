# Certified full-source spatial-refinement spectrum and clock moments

The coarse source is the full six-face cube at spacing1 and g²=64; the refined
source is the full36-face cube at spacing1/2 and g²=16. The coefficients are
(alpha,beta)=(32,1/32) and(16,1/4). They obey the exact electric matching imposed
by the existing Haar edge subdivision: a coarse link becomes two fine links.
This is a declared comparison, not a derived running coupling or a fitted mass.

## Complete source, complete frame, and infinite complement

The coarse proof frame contains every state through electric energy6, and the
existing complete trivalent census proves the next energy is13/2. The refined
201-state frame contains every state below5 by the arbitrary-valence proof in
`openwave_refined_haar_source/PROOF.md`. Thus their full Hamiltonian complements
are bounded below by208 and80, respectively. Nonnegative magnetic potentials
are retained in these bounds; no representation cutoff defines either source.

For every frame column, the exact CG multiplication rule enumerates every
Wilson output, including repeated intertwiner multiplicities. Rows outside the
frame form the complete B=QHP. The maximum absolute row sum times maximum
absolute column sum bounds ||B||². Algebraic coefficient intervals are outward
rational enclosures; omitted rows cannot be hidden in a numerical tolerance.

The fixed integer full eigenframe is checked against the actual exact matrix
using the existing complete-frame Gram/residual theorem. In particular it
certifies every finite eigenvalue, not selected Ritz values alone. If d is the
positive full-complement floor and lambda1^P<d, the Schur/min-max theorem gives
lambda1(H)>=lambda1^P_lower-||B||²/(d-lambda1^P_upper), with the finite Ritz value
as upper bound. The ground trial is a fixed integer vector. Its Rayleigh value
and complete residual include every output channel. The Temple inequality,
using the independent excited lower bound, encloses the true ground energy;
the residual divided by the separated excited threshold bounds the sine of
the ground-state angle. Hence the resulting positive gap belongs to the full
Haar source. The trial frame and its fixed rational coordinates are proof
witnesses, never a source definition or a fresh platform-dependent checksum.

## Local fields, variance, and physical-clock information

The fine calculation reads an interior fundamental square and an embedded
coarse boundary square traversing eight fine links. The latter generally lies
outside the201-state frame even when applied to its Haar vacuum; its entire
Wilson action is therefore evaluated and retained. Projecting it away would
spuriously eliminate the coarse observable. The coarse calculation reads the
corresponding four-link Wilson square. Ground-state angle bounds and operator
norms enclose each mean and second moment, hence its connected variance.

A useful exact finite theorem removes a large extra calculation. Write a
fundamental SU(2) matrix as x0 I+i x.sigma, with sum x²=1, and use generators
sigma_a/2. A link derivative of normalized trace is -x_a/2, so the carré du champ
of a length-L Wilson loop is Gamma(W)=L(1-W²)/4. Since the magnetic multiplication
term commutes with W, the full double commutator gives

m1=< (W-<W>)Omega,(H-E0)(W-<W>)Omega >
   =alpha <Gamma(W)>=alpha L(1-<W²>)/4.

The domain follows from bounded smooth Wilson multiplication preserving the
Casimir form domain. This identity applies to the interacting ground state;
Omega here means that ground state, not the Haar constant. The normalized
first spectral moment is m1/Var(W). It is a clock-sensitive average excitation
energy, not an asserted pole mass or GeV assignment.

The centered Euclidean correlation is the Laplace transform of a positive
spectral measure of total mass v=Var(W), first moment m1, and support above the
certified gap gamma. Jensen and the gap give, for t>=0,

v exp(-t m1/v) <= C(t) <= v exp(-t gamma).

Outward rational exponential bounds produce the stated intervals at source
times0,1/128,1/64 and1/32. These are actual full-source correlation enclosures;
they do not reconstruct a continuum correlator by assuming refinement
agreement. The output publishes any coarse/refined difference and whether
intervals overlap. Overlap is compatibility of these bounds, not proof of
equality, convergence, Lorentz symmetry or a common renormalized continuum.
Fine interior low modes remain present even when boundary electric moments
agree. Uniform repeated refinement, volume control and Yang–Mills UV
identification remain open, as do physical GeV calibration and Nature validation.
