# The OPH Consistency Stack

How the strange-loop hypothesis generates a universe, and why it generates exactly
one if its equations close. This file enumerates every consistency
requirement the program invokes, the freedom each one removes, the mathematical
carrier of each step, and the current status. It is the load-bearing spine of
the argument: every paper cites the stack instead of re-deriving it, and every
open item in it maps to a closure ledger row or a named gate in the gap register of
[PROOF_SPINE.md](PROOF_SPINE.md). Companion to
[STRANGE_LOOP_PRINCIPLES.md](STRANGE_LOOP_PRINCIPLES.md), [STRANGE_LOOP.md](STRANGE_LOOP.md),
[CLOSURE_LEDGER.md](CLOSURE_LEDGER.md), and
[COMPRESSION_SCORECARD.md](COMPRESSION_SCORECARD.md).

## The selection chain

The argument has one shape from end to end: begin with the space of all
self-consistent worlds, and let each consistency requirement quotient that
space. The claim structure is a chain of forcing steps.

**C1. Self-reading (SL-0).** A world that is the fixed point of its own
description must contain observers, records of those observers' readings, and
a mechanism that keeps records consistent under continued reading. This removes
all worlds without internal record-keeping. Carrier: the observer-patch normal
form (consensus paper); the record/repair architecture (screen microphysics
paper). Status: principle plus theorem package on the declared branch.

**C2. Overlap consensus.** No observer reads the whole world; descriptions
must agree where patches overlap. The public world is the quotient normal form
that survives agreement. This forces the algebraic-quantum language (patch
algebras, states, trace/Born readings) and, on the S² substrate chart (SL-1),
forces the kinematics: Conf⁺(S²) ≅ SO⁺(3,1), spatial chart
H³ = SO⁺(3,1)/SO(3), dimension 6 − 3 = 3. Removed freedom, conditional on the SL-1 chart: given S², signature, dimension, and
local symmetry group follow; the chart itself is a declared principle and is counted
where principles are counted. Carrier: compact
paper, Lorentz branch; consensus paper, quotient normal forms. Status: theorem
on the declared branch; realized-branch nonemptiness open.

**C3. Modular self-consistency.** A closed system with no external clock must
generate its own dynamics from its state–algebra pair. Tomita–Takesaki modular
theory supplies exactly this: the restricted state carries the flow, thermality
is the KMS consistency condition, and time is internal readout. Removed
freedom: no external time parameter exists to tune. Carrier: compact paper,
modular-flow sections. Status: established mathematics installed on the
declared branch.

**C4. Transportable-charge consistency.** Internal labels must transport
coherently across overlaps. The category of transportable charges, under the
compact-gauge receipt, reconstructs a compact gauge group (Tannaka–Krein). The
explicit one-Higgs matter package and the Minimal Admissible Realization rule
select SU(3)×SU(2)×U(1)/Z₆ with the hypercharge lattice, N_c = 3, N_g = 3.
Removed freedom: the gauge menu. Carrier: compact paper, gauge reconstruction;
MAR clauses declared as inputs in the scorecard. Status: theorem conditional on
the declared receipt and clause inputs.

**C5. Entropic-gravitational consistency.** Local horizon thermodynamics of
the record surfaces, made consistent patch by patch, yields the Einstein branch
(Jacobson's argument installed in the patch setting). The recovery leaves
exactly one global ambiguity, +Λg_ab. Removed freedom: the form of the gravity
law; remaining freedom: one global number. Carrier: compact paper, five-axiom
recovered core. Status: theorem on the declared branch.

**C6. Record-existence (SL-2).** A screen at exact golden-ratio balance
carries no events; records require detuning, and the declared detuning is the
Gaussian-normalized width of one electromagnetic observation:
P = φ + √π·α. Removed freedom: the pixel ratio is pinned to one equation;
remaining freedom: one local number. Carrier: fine-structure paper; SL-2.
Status: principle (declared law), with its closure tested in C9.

**C7. Inner/outer identity (SL-3, SL-4).** Under self-simulation the measured
fine-structure constant and the substrate pixel readout are one quantity, and
the record capacity is one number shared by substrate and world. A simulated
physicist's measurement of α *is* a substrate readout. The identity gives the
closure equations their test: the values completed C8 and C9 maps force must
equal the corresponding physical observables the world measures. Removed freedom: the theory keeps no authority
over which universe it describes; the measured world is the arbiter. Carrier:
STRANGE_LOOP_PRINCIPLES.md, three-layer statement. Status: principle; this is
the step that turns a completed common-observable residual into a measurement
of the hypothesis. Incomplete maps and mismatched observables do not qualify.

**C8. Capacity self-reading.** The universe must read back its own boundary
capacity without deficit or slack: N = F(N), with F the readback map
Cap∘Obs∘nf. Carrier: compact paper D6; flagship synthesis sections. Status:
open, reduced. The coupling theorem G2-GAP-1 conditionally constructs F: modulo
CP-1 (balance), CP-2 (inversion form), CP-3 (averaging carrier), the readback
fixed point is the bridge capacity, certified at 3.5321315434×10¹²². The
nominal 6.6% difference between that conditional bridge value and the
Λ-located value is not an unconditional capacity test. F and CP-1 to CP-3
remain open, and the joint Planck posterior has not been propagated through
the Λ-to-N map. CL-3 becomes evaluable only after those obligations close.

**C9. Loop closure.** The detuning law, run through the full forward chain
A_T(P), must return the pixel it started from: P = φ + √π/A_T(P), exactly.
Carrier: ../code/P_derivation; fine-structure paper. Status: open. The source
chain contracts to the certified unique fixed point α⁻¹ = 136.994835177413…
and the gauge-width chain to 137.035660136946577…. These are incomplete
declared-map outputs because the Ward-projected hadronic transport is absent
(CL-1/CL-2). The externally timestamped historical v2 target in
`../falsification/frozen_targets/` is algebraically invalid and specifies no
verdict. Corrective v3 is a permanently inactive post-target-access erratum
scaffold. A valid future test needs a detached successor whose complete method
is frozen before withheld data or produced by an audited clean-room operator.
The 2026-07-16 V1 grid cannot execute such a test: its source embeds canon
target values, the directing session had target access, and the grid envelope
is exploratory.

**C10. Unit self-consistency (SL-5).** One clock anchor (cesium) connects
substrate units to laboratory units; every SI display is bookkeeping downstream
of it. Removed freedom: unit conventions cannot carry physics. Carrier:
STRANGE_LOOP_PRINCIPLES.md; scale-bridge sections. Status: declared bridge; G rides it as an
input.

**Endpoint.** C1–C6 remove every structural freedom except a two-real-parameter
family (P, N) together with the declared discrete selections the scorecard
counts. C7 identifies the two parameters with measured readouts. C8–C9 are two
equations on the two-parameter family. If both equations hold and their fixed
points are unique, the stack admits exactly one universe. Whether that universe
is ours is decided by one reading: the closure ledger at zero. "Consistency alone generates exactly our universe" is a
theorem-shaped claim whose proof obligation is an experiment. Its target is
frozen, but no qualifying target-blind payload has been executed.

## Uniqueness lemmas

The selection chain needs "the fixed points are unique" stated as mathematics.
Three lemmas carry it.

**L1 (one P).** On any interval I where the interval-arithmetic
evaluation of the closure map g(P) = φ + √π/A_T(P) certifies g(I) ⊆ I and a
derivative bound |g′| ≤ L < 1, the Banach fixed-point theorem gives existence
and uniqueness of P in I. Status: discharged on an explicit interval, and
globally as at-most-one on the declared physical domain.
`../code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json`
proves g(I) ⊆ interior(I) and L ≤ 0.0724 for both readout maps (mpmath.iv
outward rounding, derivative enclosure by forward-mode interval AD, SU(2)/SU(3)
edge-sum tails bounded by geometric majorants so the enclosure covers the
infinite-cutoff sums). Certified unique fixed points: source map
α⁻¹ = 136.994835177413… (enclosure width 7.2×10⁻²⁴); gauge-width map
α⁻¹ = 137.035660136946577…. Stage 2 is reached for these two incomplete
declared maps. Stage 3 is not evaluable until Ward-projected transport
completes the map (CL-1/CL-2). The global at-most-one statement for each
incomplete declared map is discharged on the declared physical domain by
`../code/P_derivation/runtime/p_global_uniqueness_certificate_2026-07-14.json`:
sup |g′| < 1 is certified on every piece of an interval subdivision of
α ∈ [0.005, 0.01] (α⁻¹ ∈ [100, 200], the declared solver scan window of
`paper_math.solve_closure`), 256 pieces per readout map at the declared
cutoffs su2=120/su3=90 with the tail majorants folded in, worst piece
L ≤ 0.3041, exceptional set empty. Each declared map therefore has at most
one fixed point on the declared domain, and with the stage-2 existence
certificate exactly one. Monotonicity of A_T in P is no longer required for
the domain-global statement; it remains open only as an analytic refinement
beyond the declared domain.

**L2 (one N).** Under CP-1 to CP-3 the coupled readback map is affine in the load
coordinate with contraction 1/2, and its unique fixed point is the bridge
capacity (G2-GAP-1, certified enclosure at relative width 1.6×10⁻²⁵). The lemma
is conditional on the three premises; discharging CP-1, the counting theorem, is
the highest-value open mathematics in the program.

**L3 (no-landscape corollary).** Under L1 for a completed registered P map and
L2 after CP-1 to CP-3 are discharged, the SLH principle set admits at most one
(P, N). A target-blind completion landing outside its registered basin then
falsifies that formulation (STRANGE_LOOP_PRINCIPLES.md rules 2 and 7). The
present certificates do not discharge this corollary: L1 applies to two
incomplete P maps, and L2 is conditional. The no-landscape conclusion is a
proof obligation, not an evidential result.

## Dependency structure of the closure ledger

The ledger rows are not independent. They reduce to four generating objects:

| Generator | Object to construct or compute | Ledger rows it closes or moves |
|---|---|---|
| G1 | Complete the Ward-projected hadronic transport and run a target-blind payload under a detached future successor whose full method-selection chronology is blind or clean-room audited. The 2026-07-16 V1 grid is exploratory only: its source embeds canon target values, its directing session had target access, and its sampled min/max envelope is not an interval certificate. Its `S_hadronic` range contains the zero-EW CL-1 point diagnostic 0.8954 but is 1.17×10⁸ times wider than v2's historical invalid tolerance; v3 registers no active threshold. The V1 report also conflates the closure-residual name `Delta_source_residual_inverse_alpha_units = 0.0414658…` with the payload's total `Delta_source`; the correct total endpoint target is `Delta_required_P = 8.7280337…`. No promotion follows (../code/particles/hadron/ward_projected_payload/PAYLOAD_STATUS.md) | CL-1, CL-2 |
| G2 | The capacity readback map F and its contraction certificate. Status 2026-07-16: G2-GAP-1 is a conditional theorem (../code/capacity_readback/G2_GAP_1_COUPLING_THEOREM.md). Its fixed point equals the EW bridge capacity N = π·exp(6π/(P·α_U)) only modulo CP-1, CP-2, and CP-3. The certified conditional fixed point at P_fwd encloses 3.5321315434e122 with relative width 1.6e-25 (../code/capacity_readback/runtime/F_candidate_coupled_certificates.json). CL-7 remains open. After F and CP-1 to CP-3 close, CL-3 also requires propagation of the joint Planck posterior through the Λ readout (../code/capacity_readback/F_CONSTRUCTION_2026-07-14.md) | CL-3, CL-4, CL-7 |
| G3 | Define the physical electroweak readout and complete its scheme map before comparing it with data. The emitted pair is a running/tree chart coordinate; the quoted references are stale PDG 2025 mass-dependent-width Breit-Wigner parameters, not pole masses, so CL-5 is not physically evaluable. The 96-entry one-loop sweep tests only the declared chart menu. The 2026-07-16 two-loop result applies SM two-loop increments to an MSSM one-loop baseline and excludes that hybrid prescription, not two-loop repair in general. The pole packet is a scale-dependent partial PRTS/Feynman-gauge prescription with open tadpole/FJ and vev schemes; its JKV cross-check audits a slope term rather than the complete finite packet. No raw W/Z pull from these packets is evidence | CL-5 |
| G4 | Solver hygiene: matched printed pair at certified precision | CL-6, closed 2026-07-14 (converged precision-100 reruns; identity to 35+ digits; CI test) |

Beyond the ledger, the mass sector waits on one further object: a
source-derived flavor-orbit selector (the proven two-modulus
non-identifiability of the light-quark equations names exactly what is
missing). The program is therefore finite and enumerated: four generators and
one selector stand between the current standing and full closure. Every one of
them is a computation or a construction; reinterpretation closes nothing.

## Reading rules

- A C-row's "removed freedom" claim holds on its declared branch and carrier;
  branch nonemptiness and receipt gates stay visible at the carrier.
- No C-row is evidence for any other. Evidence enters only through the closure
  ledger and the compression scorecard.
- The stack may grow rows (a newly identified consistency requirement is a new
  C-row with a new carrier), and rows never silently strengthen: changes here
  are dated, like ledger verdicts. Adopted 2026-07-14.
