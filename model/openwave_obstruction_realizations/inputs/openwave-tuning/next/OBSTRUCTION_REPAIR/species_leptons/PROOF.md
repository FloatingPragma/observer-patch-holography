# Neutral-family and beta-emission repairs

This package changes no earlier producer, acceptance run or source. `audit.py`
pins every consumed numerical object and writes its own complete result.
Architectures may be declared: their uniqueness is not a prerequisite. A
target-informed coefficient remains a calibration, however. The construction
rule in `reverse-engineering-reality/AGENTS.md` makes both points explicitly.

## Row 10: exact obstruction of the single-phase construction

The existing uncalibrated neutral coefficient is `M = V* Y V†`, where `Y` is
real symmetric positive and `V = I + (exp(i phi)-1) P` for a real rank-one
projector P. If F is the real charged eigenframe, the current matrix is

    C = Fᵀ V F = I + (exp(i phi)-1) u uᵀ.

It is complex symmetric for every phase, real projector and charged spectrum.
Consequently its modulus-squared matrix B is symmetric. Rephasings do not
change this condition. This is stronger than the failure of the one chosen
angle calculation, but it applies only to this family of constructions.

Write x = sin² theta12, y = sin² theta23, z = sin² theta13, and
K = 2 sqrt(x(1-x)y(1-y)z). In the standard three-angle convention the
off-diagonal equality B12=B21 requires

    cos(delta) = [xy - xz - (1-x)yz] / K.

The repository comparison triple (0.308, 0.470, 0.02215) gives 1.90614,
outside [-1,1]. These are consumed comparison inputs, not a new measurement
lookup. The script also checks all six relative neutral assignments. Four
require cos(delta) equal to 1.90614, -2.16218, -4.47188 or 5.07128. The
remaining two have a nonzero antisymmetric part independent of delta.
Thus relabelling the neutral columns does not repair the comparison at
these values. No experimental uncertainty interval is asserted here.

### Executed architectural repair

The repository already declares the native pulse alphabet (2,3,10), phase
-i and ascending primary word (2,3,10) before its CP evaluation. We reuse
that word, without searching for a fit. Products of noncommuting real
projector phases need not be symmetric. The resulting current has
max |B-Bᵀ| = 0.0123213 and removes the proved one-pulse obstruction.

For any unitary V, M=V* Y V† remains symmetric and

    (VF)ᵀ M (VF) = Fᵀ Y F = diag(lambda).

The script checks V's unitarity and M's symmetry symbolically, reconstructs
the Takagi identity numerically, and computes all 36 complete five-outcome
current instruments (three momenta, four times, three initial families).
Failures are retained. The same internal Higgs/Ward factorization applies
to every symmetric family coefficient; no new internal charge is introduced.

This architecture predicts (sin² theta12, sin² theta13, sin² theta23) =
(0.0443755, 0.123449, 0.178037). It also fails the comparison. Removing the
single-pulse obstruction is therefore a solved structural subproblem, not a
solution of the physical PMNS row. Another declared dynamics must produce
different output. A search selecting words by these measured angles would
be calibration and must be recorded as such.

## Row 11: the asymptotic residual cannot be removed by more rounds

The independent calculation uses the declared circulant eigenvalues

    lambda_k = 1 + sqrt(2) cos(2/9 + 2 pi k/3),   m_k = lambda_k².

It gives the ordered ratios 206.7703159727 and 3477.4728371046. The first
differs from the consumed muon comparison by 442.0158 of that comparison's
quoted standard deviations. Increasing the finite round count converges to
this discrepant value. The script retains the finite-round ladder and the
local phase sensitivity. The phase change that would fit the muon is
reported solely as a diagnostic, labelled target-fitted and not adopted.

An interacting mass-shell correction on a frozen source could in principle
change this conclusion. It must be calculated from that source, including
its electron reference, rather than chosen to cancel the residual. No such
correction is derived in this package. The original retrospective ancestry
and the missing common-source transfer remain in force. Architectural
non-uniqueness alone is not a reason to reject the construction.

## Row 25: a three-dimensional physical-direction repair

Keep the existing complete 120-state hadron source, actual local weak
currents, all three neutral masses and their weights, and both contact cells
at (0,0,0) and (1,0,0). Change only the declared free lepton regulator to

    p_j = 2 sin(theta_j),
    m(theta) = m0 + sum_j (1-cos(theta_j)),
    E = sqrt(|p|² + m²),
    v_j = sin(theta_j) [4 cos(theta_j) + m] / E.

This is Wilson coefficient w=1 for both charged and neutral leptons. It is
an architectural amendment, not a claim about the unchanged literal repair
law. The weak operator and its full-Fock norm bound are unchanged.

The script independently projects the stored full local currents onto the
actual neutron/proton quartet and verifies their vector/axial matrix
elements. It then evaluates the existing exact proton-channel spin trace

    [1+3g²+(1-g²)a.b+2g(1-g)a_z+2g(1+g)b_z]/2,

where a=p_e/E_e and b=-p_hole/E_hole. The two-cell form factor
cos²((theta_e,x-theta_hole,x)/2) and the finite-time sinc² mismatch remain.
All neutral weights remain; no continuum angular distribution is substituted.

### Why physical direction now has a strict sign

For w=0, theta_z -> pi-theta_z preserves p and E and reverses velocity.
The source is invariant under the corresponding staggered reflection.
Therefore every velocity-odd emission statistic vanishes, even though its
Dirac-momentum-label statistic need not vanish. This is independently
reproduced on cubes of side 4 and 8. Side 8 has nonzero free velocities;
its cancellation is not the zero-velocity accident of the side-4 grid.

For w=1, pair momenta with the same sin(theta_z) and opposite cos(theta_z).
Write c=|cos(theta_z)| and a0=m0+1+sum_(x,y)(1-cos(theta_j)) > 0. The two
masses are a0-c and a0+c. The positive-c member has lower energy, and its
velocity has the sign of sin(theta_z). If the negative-c member has the
opposite velocity, its contribution is smaller in magnitude; otherwise
both contributions have the same sign. Thus the pair's contribution to

    S = sum_theta [p_z/E] sign(v_z)

is nonnegative. For L divisible by 4, the two theta_z=+/-pi/2 planes are
strictly positive. On those planes E<7 for both species in the consumed
mass range, giving S/L³ >= 4/(7L).

The x-contact form factor does not depend on either z coordinate. The
linear electron term therefore gives the instantaneous directional bound

    -O_e >= (Gt)² g(g-1) 2/(7L),
     O_antinu >= (Gt)² g(g+1) 2/(7L).

All other terms vanish under the appropriate z reflection. The neutral
weights sum to one. These bounds require only 1.2752<g<1.2754, the inherited
source interval, and neutral masses below 10^-6 in electron units. They do
not claim that this interval was predicted independently.

### Complete Fock-space error, not a first-pair truncation claim

Use t=1/100000 and the repository's unchanged G=1/64. The declared Wilson
kinetic term changes the interaction picture by a unitary; it does not
change the inherited full weak norm bound 54. The Dyson state error is

    delta <= (54 G t)²/2.

The selected first-order amplitude is at most 4Gt, so any selected
signed probability differs by at most (8Gt+delta)delta. All other hadrons,
flavours, helicities and multipair sectors are in this bound.

The mismatch satisfies |Delta|<14. The exact inequality
|sinc²(Delta t/2)-1| <= Delta² t²/12 bounds the finite-time weighting change
by 16(Gt)²(14t)²/12. `audit.py` subtracts both errors using rational
arithmetic and obtains strictly positive margins at L=4 and L=8. Divide
by 120 for the original preparation, retaining failure probability 119/120.
The strict directional signs therefore do not rely on numerical quadrature
or unbounded floating roundoff. The printed numerical probability intervals
are diagnostics; only the displayed rational sign margins are certificates.

The controls are w=0, opposite chirality and an equal classical mixture
of the two mirror currents (with its fair coin retained). This last control
is an incoherent mixture, not the coherent sum of the currents. Opposite
chirality reverses the signal; the mixture and w=0 give zero. The script
retains all twelve records.

### On-shell doubler exclusion in all three dimensions

If any electron momentum component lies in a wrong half-zone
|theta_j|>=pi/2, its energy is at least sqrt(8). A neutral in such a zone
costs at least 2, and an accompanying electron costs at least 1. The
consumed energy release is Q=2.53099<sqrt(8). Thus on the free energy shell
both outgoing particles lie in the central half-zone in every coordinate,
where each velocity has its momentum's sign. Off-shell states at finite
time have not been discarded; they remain in the complete calculation.

### Scope still open

This repairs the physical-direction parity observable in three dimensions
on the declared Wilson-amended source, with a complete short-time Fock-space
error account. It does not derive the neutron-proton gap, the selected weak
chirality, the axial magnitude or a neutron lifetime. Those empirical and
architectural inputs remain explicit. Transfer to a common uncalibrated
colour/weak/electromagnetic source remains a separate obligation.
