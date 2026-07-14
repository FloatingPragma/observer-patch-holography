# The OPH Consistency Stack

How the strange-loop hypothesis generates a universe, and why it generates exactly
one if its equations close. This file enumerates every consistency
requirement the program invokes, the freedom each one removes, the mathematical
carrier of each step, and the current status. It is the load-bearing spine of
the argument: every paper cites the stack instead of re-deriving it, and every
open item in it maps to a closure ledger row or a named gate. Companion to
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
closure equations their test: the values C8 and C9 force must equal the values
the world measures, exactly. Removed freedom: the theory keeps no authority
over which universe it describes; the measured world is the arbiter. Carrier:
STRANGE_LOOP_PRINCIPLES.md, three-layer statement. Status: principle; this is
the step that makes every residual a measurement of the hypothesis.

**C8. Capacity self-reading.** The universe must read back its own boundary
capacity without deficit or slack: N = F(N), with F the readback map
Cap∘Obs∘nf. Carrier: compact paper D6; flagship synthesis sections. Status:
open. F is not yet constructed (CL-7); the one-capacity requirement already
yields a live 6.6% contradiction between the electroweak-bridge value and the
Λ-estimated value (CL-3).

**C9. Loop closure.** The detuning law, run through the full forward chain
A_T(P), must return the pixel it started from: P = φ + √π/A_T(P), exactly.
Carrier: code/P_derivation; fine-structure paper. Status: open. The source
chain contracts to the certified unique fixed point α⁻¹ = 136.994835177413…
and the gauge-width chain to 137.035660136946577…, both outside the SL-3 basin
by the missing hadronic transport term (CL-1/CL-2); the blind completion
against the externally timestamped frozen target in
`falsification/frozen_targets/` is the decisive armed experiment.

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
theorem-shaped claim whose proof obligation is an experiment, and the
experiment is armed.

## Uniqueness lemmas

The selection chain needs "the fixed points are unique" stated as mathematics.
Three lemmas carry it.

**L1 (one P).** On any interval I where the interval-arithmetic
evaluation of the closure map g(P) = φ + √π/A_T(P) certifies g(I) ⊆ I and a
derivative bound |g′| ≤ L < 1, the Banach fixed-point theorem gives existence
and uniqueness of P in I. Status: discharged on an explicit interval, and
globally as at-most-one on the declared physical domain.
`code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json`
proves g(I) ⊆ interior(I) and L ≤ 0.0724 for both readout maps (mpmath.iv
outward rounding, derivative enclosure by forward-mode interval AD, SU(2)/SU(3)
edge-sum tails bounded by geometric majorants so the enclosure covers the
infinite-cutoff sums). Certified unique fixed points: source map
α⁻¹ = 136.994835177413… (enclosure width 7.2×10⁻²⁴); gauge-width map
α⁻¹ = 137.035660136946577…. Stage 2 of the basin-then-contract protocol is
reached for P; the stage-3 landing verdict is unchanged (outside the SL-3
basin, CL-1/CL-2). The global at-most-one statement is discharged on the
declared physical domain by
`code/P_derivation/runtime/p_global_uniqueness_certificate_2026-07-14.json`:
sup |g′| < 1 is certified on every piece of an interval subdivision of
α ∈ [0.005, 0.01] (α⁻¹ ∈ [100, 200], the declared solver scan window of
`paper_math.solve_closure`), 256 pieces per readout map at the declared
cutoffs su2=120/su3=90 with the tail majorants folded in, worst piece
L ≤ 0.3041, exceptional set empty. Each declared map therefore has at most
one fixed point on the declared domain, and with the stage-2 existence
certificate exactly one. Monotonicity of A_T in P is no longer required for
the domain-global statement; it remains open only as an analytic refinement
beyond the declared domain.

**L2 (one N).** Once F is constructed (CL-7), the same schema applies: a
certified contraction interval for F yields existence and uniqueness of N_CRC.
The monotone-and-bounded shape of Cap∘Obs∘nf is the design target; the lemma is
conditional until F exists.

**L3 (no-landscape corollary).** Under L1 + L2, the SLH principle set admits at
most one (P, N). There is no vacuum landscape to relocate into: a blind closure
computation landing outside its basin falsifies the formulation, permanently
(STRANGE_LOOP_PRINCIPLES.md rules 2 and 7). Maximal falsifiability is here a corollary of
uniqueness. This is the sharpest structural contrast
with landscape frameworks, and it is the reason the program can afford to
freeze targets before computing payloads.

## Dependency structure of the closure ledger

The ledger rows are not independent. They reduce to four generating objects:

| Generator | Object to construct or compute | Ledger rows it closes or moves |
|---|---|---|
| G1 | Ward-projected hadronic transport, blind, vs the frozen target. Status 2026-07-14: payload harness built; the first-principles internal-mass bracket S_eff in [0.558, 1.054] contains the required 0.8954; closing the window needs 4×10⁻⁹ relative precision on the hadronic moment, beyond every current method with or without measured data (code/particles/hadron/ward_projected_payload/PAYLOAD_STATUS.md) | CL-1, CL-2 |
| G2 | The capacity readback map F and its contraction certificate. Status 2026-07-14: every declared-structure candidate certified and excluded (fixed points at 1 to 10³ nats); the unique 10¹²²-scale reading couples the port load to the D10 observation step and reproduces the EW bridge expression, so CL-7 is the coupling theorem G2-GAP-1 (code/capacity_readback/F_CONSTRUCTION_2026-07-14.md) | CL-3, CL-4, CL-7 |
| G3 | Repaired forward electroweak chain from the strange-loop principles. Status 2026-07-14: preregistered 96-entry sweep executed; the certified-pixel menu is exhausted at one loop; the SM triple admits no fixed point (Landau ceiling disjoint from the fixed-point window) and β_EW in {3, 5} land three orders of magnitude off, so those selections are consistency-forced; the single landing at the rounded legacy carrier is a counted two-fold selection awaiting its own frozen test (falsification/preregistered/ew_repair_results_2026-07-14.json) | CL-5 |
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
