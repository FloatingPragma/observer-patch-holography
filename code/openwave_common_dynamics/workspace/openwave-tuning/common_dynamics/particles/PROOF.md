# Particle, hole and coupled readout contracts

All free correlators in this package use the frozen common W or the explicit
local mass-coin proposal in `MASS_PROOF.md`. All coupled records call the
actual `vector_field/autonomous.py` API. No independent Hamiltonian,
particle trap, observed mass, or force array is supplied to a consumer.

## Correlations and persistence

For the supplied CAR mode extension Γ(U), the empty reference gives
⟨empty|c_f Γ(U)^n c_g†|empty⟩=f*U^n g. This is an exact one-particle
correlator of the same source. It does not establish that empty is an
energy ground state. Finite-dimensional unitary correlators have poles on
the unit circle and recurrence; these facts alone are insufficient for an
infinite-volume particle pole, decay width or measured lifetime.

An oriented graph cycle defines a compact state: place opposite amplitudes
on the two directed moving slots of each edge. Its divergence vanishes,
so A*ψ=0 and Rψ=−ψ, while Sψ=−ψ. Hence Wψ=ψ exactly. The executable uses
the native triangle(0,1,5); all six nonzero slots are retained. Under the
mass proposal its phase becomes e^{iμ}. This is an exact localized charged
mode with zero base quasienergy, not a positive-mass mobile species.

The held-event state (|port0,event0⟩−|port0,event2⟩)/√2 has R=−1 and S=+1.
Its base phase is−1 and its massive phase is−e^{iμ}, with no motion. It is
a control against identifying a finite quasienergy or apparent localization
with inertial rest mass. It remains localized even under the autonomous
field word because it never reaches a moving slot; the field still evolves.

Point/event preparation, the compact triangle and this held state are each
propagated for16 complete source steps. Source and event padding, norms,
correlations and all96 position marginals are recorded. Base native3-bit
QPE distinguishes the two compact modes at phase bins0 and4 using the
frozen energy apparatus and its full unconditional error certificate.

Finite-box acoustic spacings at development sizes2,3,4 decrease with box
size; the zero-momentum base acoustic phase is zero. These spacings are
not a rest-mass determination. The supplied local mass term has a different,
explicit source origin, detailed separately.

## Hole reference and global phase

The fully filled1152-physical-mode state is a separately supplied reference.
On its one-hole subspace Γ(U) acts as det(U) conjugate(U), in the standard
signed missing-mode basis. Base W is real and det(W)=1, so its hole and
particle spatial propagators coincide while the transported charge/current
has opposite sign. The corresponding Gauss backgrounds are different.
This does not construct simultaneous particle and antiparticle species,
annihilation, pair creation, or a stable interacting vacuum.

For the mass proposal, Tr(R) on physical modes is−960, so
det(Uμ)=exp(i960μ). Relative to the same filled reference, the hole
propagator is conjugate(Uμ)=U_{−μ}. The executable checks that conjugation
and the physical determinant. A controlled full-source experiment must
retain the determinant phase; it cannot inherit the base-W probe by
relabelling the particle as a hole.

## Coarse native detector without event dephasing

A full occupation read measures too much to implement a coarse intermediate
position measurement. The supplied instrument instead creates seven fresh
blank flags. For address bit b it reversibly XORs occupations of every
encoded mode whose position address has bit b. There are7168 CNOT calls.
It measures the seven flags with fresh native C/E pointers and reverses
the7168 CNOTs. On the one-particle2048-mode sector the effects are exactly
P_xρP_x for all128 position addresses. Coherence among the16 event labels
at a fixed position and all unmeasured field/reference coherence survives.

On arbitrary Fock sectors the same circuit measures seven occupation
parities; it is not interpreted as one particle's position. The filled
physical reference has zero address parities, so the same coarse circuit
also reads one-hole position without resolving its event label.

The optional joint terminal instrument additionally reads every1152 field
bit. That field dephasing is an explicit intervention. For the field source,
flags3203..3209 are outside2048 matter wires,1152 field wires and the three
producer scratch wires3200..3202. There are7 fresh blank flag resets and
1159 fresh ready-pointer resets. All reset timeouts and native fine pointer,
padding and export outcomes remain in the unconditional law.

The apparatus-only bound is the sum of14336 native CNOT operator errors,
one pointer error per measured bit, and a reset-failure union bound using
exact rational worst-initial-state success lower bounds. These are used as
uniform conditional bounds: reset independence is not assumed and a product
of marginal success probabilities is not asserted to be the exact law.
The actual first-failure law uses each history-dependent conditional factor.
In the actual reset implementation every conclusive fine Kraus word is a
scalar multiple of |target⟩⟨x|. Its success effect is therefore diagonal;
the minimum listed success probability is its minimum eigenvalue, including
inputs entangled with a reference. Private reset tapes remain in the actual
instrument. Comparing a successful blank carrier with the original reduced
reference marginalizes those private tapes; it does not assert their
independence from that reference. The terminal particle/field TV claim uses
that observable marginal with failures retained. The finer tape law is
specified by the actual reset kernels, not a product-reference ansatz.
The consumed `vector_field/RESET_DOMAIN.md` supplies the stronger fine-tape
comparison: replace each successful weight |a_w|² by |a_w|²/s_x. That
ideal fine instrument keeps its resolved initial label and reference
correlations. Its half-diamond distance from the actual reset, with failure
flags retained, is at most1−min_x s_x. This proves the conditional union
and composition bounds without assuming an independent carrier ensemble.
It excludes source preparation and evolution until explicitly composed in
`massive_coupled.error_contract`. The complete three-step combined bound
adds3203 initial source blank resets, actual one-particle preparation, and
three full matter/current/field/mass words. The architecture's coherent
services and source response retain their declared ancestry.

`massive_phase_contract` separately compiles the controlled encoded
one-particle Wμ word, including every scalar phase. It adds the native
source, control-wrapper, preparation, pointer and24 blank/one ready reset
errors. This probe is not transferred to the full Fock/field source without
controlling that entire augmented word and allocating its actual wires.

## Exact source continuation and coupled measurements

The finite development member has electric phase α=.17 and plaquette phase
β=.41 on face0; the other207 couplings are zero. It is a one-face member,
not a homogeneous Maxwell limit. `coupled_response.py` executes particle
and separately prepared hole states for three steps with the actual API.
Every288-edge current is computed from the same reflected moving amplitudes.
The unconditional law obeys charge-weighted continuity

\[
q\,[p_{n+1}-p_n]+B\langle J_n\rangle=0.
\]

This statement is checked against all96 position changes and all288 current
components. The complete electric state remains in every branch; field
energy is a moment of that state, not a supplied force.

The actual coarse position instrument at step1 branches the source without
normalizing any branch. Each branch receives two more identical source
steps, and the full two-time law sums to one. Its final position distribution
differs from the uninterrupted experiment, so measurement backaction is
visible rather than implicitly erased by a success normalization.

`COUPLED_JOINT_LAWS.json.gz` and `MASSIVE_JOINT_LAWS.json.gz` emit every
nonzero position/field block and all12 event probabilities for each terminal
experiment. Omitted ideal blocks are zero. The first file also records every
two-time coarse-position branch with its original probability. Native
complement outcomes are controlled by the factorized complete circuit and
reset certificate, not assigned zero by a conditional success law. SHA256
pins cover compressed and uncompressed complete laws.

The saved probability tables are binary64 evaluations of the complete finite
linear program. Normalization, continuity and conservation residuals are
numerical controls. The native operator/error certificate compares the
emitted native program with the mathematically exact ideal circuit law;
it is not a rigorous rounding enclosure for every saved decimal entry.
An interval evaluation of those tables would add a separate numerical
certificate. No physical-error claim is inferred from their small residuals.

The mass intervention uses the actual autonomous API after the source-local
kick, with all current and field records on the same source. Its full native
Γ(kick) word and error are emitted. The observed finite position responses
establish interaction and mass-coin effects in this development member.
They do not themselves establish a two-charge force law or interacting
renormalized mass. The logical step uses t=nτ₀h; bare-W pulse durations or
clock detector budgets are not transferred to this longer word.

No reserved acceptance case is run and no official row is changed here.
