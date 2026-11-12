# One autonomous source for waves and native clocks

## Geometry on the same finite graph

Let L be the positive Dirichlet Laplacian on the7×7 planar grid. Fix source
sites (2,3),(4,3), and let B insert their two coordinates. Define

    F=L⁻¹B(B*L⁻¹B)⁻¹, S=F*LF, M=F*F.

The producer constructs these matrices over exact rationals. F takes the
specified values at the two sites, is harmonic at every other site, and
uniquely minimizes q*Lq for those values. If q=FQ, the quadratic geometry
action has stiffness κS and inertia κM, with κ=100. Thus spatial profiles,
elastic coupling and inertia follow from the same grid action. Restricting
geometry to this harmonic two-coordinate architecture is declared; this
is not an invariant-subspace or adiabatic theorem for an unrestricted49-mode
geometry. Reflection symmetry gives S=[[s0,s1],[s1,s0]] and
M=[[m0,m1],[m1,m0]]. Their normal frequencies are
ω±=sqrt((s0±s1)/(m0±m1)).

With product occupation cap C, let a± be the two complete truncated
lowering operators and c±=[4κω±(m0±m1)]^(-1/2). Put
Q0=c+(a++a+*)+c−(a−+a−*), Q1=c+(a++a+*)−c−(a−+a−*), and
Hg=ω+N++ω−N−. The two Q operators commute exactly. The number Hamiltonian
includes each finite oscillator's upper-boundary completion; it is not
obtained by discarding that term from a naively truncated quadratic action.
At C=1, |Q_a|<0.067 and |q_x|<0.1.

## Matter, field and a common clock

Both native cells retain all eight eigenlevels of the positive excitation
matrix K from `openwave_energy_lapse.source.native`. The field carrier is
the vacuum plus all49 planar one-quantum modes. For each commuting geometry
coordinate value q, define h(q)=sqrt(exp(2q)Lexp(2q)). On the full carrier,

    H = Hg + exp(Q0)K0 + exp(Q1)K1 + (0 direct-sum h(FQ)).

The primary dimension is64×4×50=12800. Every term is positive; the unique
source ground is native00, geometry vacuum, field vacuum, with energy0.
External references may have their own zero-energy degeneracy. Both native
level operators and the field vacuum/one-quantum decomposition commute
with H. These exact reducing blocks replace a full12800-state
diagonalization while retaining every complementary block and outcome.

H generates one autonomous linear unitary group. The field acts on geometry
and geometry acts on the field, with no geometry refresh or classical
profile inserted during evolution. Total energy is conserved; separate
geometry, native coupling and field energies need not be conserved.
Finite recurrence is retained. The scalar metric choice exp(2q), the
restricted geometry architecture, κ, cap, field-number carrier and
external preparation/readout schedule are declared. The field is not a
proved three-dimensional Maxwell or relativistic photon theory.
The instantaneous positive-square-root interaction conserves field number
by construction. Its equivalence to canonical quantization of a scalar on
a time-dependent metric, including possible pair creation, is not asserted.

## Exact finite numerical certificate

For C=1 the common geometry-coordinate eigenframe is an exact Walsh matrix
with entries ±1/2. Each of its four q profiles has interval-enclosed
algebraic coordinates. The positive wave matrix Kq=exp(2q)Lexp(2q) has
λ_min≥8 exp(−0.4)/64>1/16. For a proposed dyadic square root A, exact
integer products and outward interval source entries enclose
R=A²−Kq. Positive Gershgorin bounds establish A>0. The Sylvester equation
for A−sqrt(Kq) gives

    ||A−sqrt(Kq)|| ≤ ||R||/(λ_min(A)+1/4).

Every root matrix is rounded to denominator2^46. Exact Walsh conjugation
assembles their direct sum at denominator2^48. Exact interval native
energies, exponentials and geometry frequencies bound the remaining
source-to-dyadic error. These bounds use the spectral minimum and
Frobenius bounds from exact integer square-root ceilings; no eigenvalue
solver accuracy is presumed.

For the196-dimensional wave block and each4-dimensional vacuum block,
the eigensolver proposes a dyadic O,D. With
η≥||O*O−I||<1/2 and R≥||HO−OD||, its polar factor W is unitary,
||O−W||≤η and ||O||≤2. Every source block used here has norm below30.
Consequently ε=2R+120η bounds ||H−WDW*||, and

    ||exp(−itH)−O exp(−itD) O*|| ≤ |t|ε+3η.

The exact source errors are included before applying this bound. Interval
phases are rounded with checked error below2^(-48) in each component.
The arithmetic contract is IEEE754 binary64 round-to-nearest, no overflow,
ordinary complex products and dot sums of dimension at most196. Writing
u=2^(-53), γ_k=ku/(1−ku), the conservative bound
16dγ_(32d)+2^(-40)<2^(-28), d≤196, covers forming and applying a propagator,
including phase and normalization arithmetic on the stated unit vectors.
The producer adds that guard to every step. This contract is explicit;
it is not a hardware-independent proof of a numerical library.

The initial packet is defined by the emitted30-bit integer complex
coefficients normalized by their exact integer squared norm. A Gaussian
formula proposes those integers. The normalization evaluation error is
bounded by2^(-40). If δ bounds a state-vector error, a record probability
differs by at most2δ+δ². These bounds also hold with arbitrary quantum
references. For sixteen retained-history steps, the norm error is at most
16δ+256δ² plus the initial normalization error, covering the slight norm
growth of approximate propagators. The small additional decimal guard in
the output covers final bound evaluation.

For the initially pure geometry/wave branch, partial trace gives a geometry
density error at most epsilon=2delta+delta squared. Its purity error is at
most epsilon[1+(1+delta)²], plus an explicit2^(-35) trace-evaluation guard.
The certified purity at time4 is below one, establishing geometry/wave
entanglement on the same autonomous source.

## Complete wave, clock and arrival instruments

The input geometry is its vacuum, with native levels either(0,0) or(1,0),
and an equal mixture of field vacuum and the full-grid packet. Thus the
control and source are different preparations of one H. The zero-native
control includes the wave's own backreaction; it is not a flat geometry.
Each single-time instrument records every native level, every geometry
number and every field site or vacuum. The63 complementary native sectors
have zero probability for these preparations and are retained explicitly.
Both wave and vacuum geometry records are included.

A separate coherent input uses native(|0>+|1>)/sqrt(2) at the first cell.
The exact block decomposition supplies complete native-X/geometry/field
records on this single input. Its clock visibility depends on the field
and geometry overlaps. This connects coherent native clock observations
to the dynamical wave source without identifying records from different
initial states. These preparations retain their distinct labels.

The arrival instrument applies the same U(1/2), then records every receiver
site x≥5 and geometry number on a first hit. The no-hit branch continues;
all final wave-site/geometry and vacuum/geometry timeout records are kept.
At sixteen steps, their sum is one for the exact instrument. This is a
scheduled first-arrival measurement, with its quantum disturbance included,
not an unobserved classical trajectory. Each source evolution conserves H;
the readout can exchange up to30 source-energy units per step.

The operational mean arrival time is conditioned on a hit, and both hit
and timeout probabilities are reported. Bounds on its numerator and
denominator give a native-source delay interval
[0.010034899,0.010075386] source-time units. Hit probabilities are about
0.1680870 and0.1682521 on the initial wave branch; timeouts are not removed.
At t=4 the transverse mean relative to the zero-native control lies in
[−0.001846753,−0.001846538] sites. Reflection about y=3 commutes with the
complete H and reverses this deflection. The C=2 numerical regulator
control retains both signs, with mean delay about0.01053894. It is not a
cap-independent limit certificate.

## Scope

The certificate supplies actual wave bending and an operational delay on
an autonomous quantum geometry, sharing one source, retained reference
domain and time unit. It does not import the31³ static lapse profile or
the earlier two-cell geometry state. Physical G, laboratory time, observed
lensing, Einstein dynamics, a full photon spectrum and irreversibility are
not inferred from the finite construction.
