# Blinded Generative Repair Kernel for IBM Quantum

## Status

This document defines the fixed-cutoff acceptance protocol for GitHub issue
`#509`. It is a declared, falsifiable hardware model and implementation
contract, not a hardware-run receipt. It is not a new OPH theorem and it does
not claim that repair dynamics select a physical vacuum.

The source boundary matters. The finite-quotient theorem surface says that
normal-form settlement does not select a probability law and that a repair
generator relaxes only after a reference law has been declared. The
edge-sector handoff likewise assumes a generator with stationary law

\[
\pi_\beta(\alpha)\propto d_\alpha e^{-\beta C_2(\alpha)}.
\]

The experiment below therefore tests an engineered observer-like system with
local state, a proposal boundary, readback, durable records, conditional
feedback, and a public evidence bundle. It does not use the IBM device to infer
that nature must choose the programmed controller.

## Primary Kernel

Let \(G\) be a finite group with a symmetric generator set
\(S=S^{-1}\). The observer-visible mismatch state is

\[
h=x_A^{-1}x_B\in G,
\qquad
\Phi(h)=\ell_S(h,e),
\]

where \(\ell_S\) is Cayley word distance and \(e\) is the consistent normal
form.

The conventional heat step is

\[
H_G(h,h')=
\frac{1}{|S|}
\#\{s\in S:h'=hs\}.
\]

To avoid periodicity and nonunique parity-class stationary laws, the benchmark
uses the preregistered lazy heat step

\[
L_G=\frac{I+H_G}{2}.
\]

For a proposed repair generator \(s\), the system writes the record

\[
a(h,s)=\mathbf 1[\Phi(hs)<\Phi(h)].
\]

It commits the update only when that public record is one:

\[
R_s(h)=
\begin{cases}
hs,&a(h,s)=1,\\
h,&a(h,s)=0.
\end{cases}
\]

Averaging uniformly over repair proposals gives \(R_G\). The primary ongoing
process and the depth-matched open-loop null are

\[
K_{\rm repair}=L_G R_G,
\qquad
K_{\rm heat}=L_G^2.
\]

Both begin from ordinary computational-basis states. No heat-kernel or OPH
target amplitudes are prepared. The held-out observable is the joint process
law of the heated-state record, repair-decision record, and final state, plus
the stationary law reconstructed from the complete empirical transition
matrix.

## Exact Frozen Predictions

For \(\mathbb Z_5\) with generators \(\{+1,-1\}\), the record-gated stationary
law is

\[
\pi_{\mathbb Z_5}^{\rm repair}
=\frac{1}{33}(21,5,1,1,5),
\qquad
\mathbb E[\Phi]=\frac{14}{33}.
\]

The open-loop lazy-heat null is uniform and has
\(\mathbb E[\Phi]=6/5\). The stationary one-step KL separation is
`0.3076363756104924` nats per transition.

For \(S_3\) with all three transpositions as generators, in the state order
\(e,(01),(12),(02),(012),(021)\),

\[
\pi_{S_3}^{\rm repair}
=\left(\frac12,\frac16,\frac16,\frac16,0,0\right),
\qquad
\mathbb E[\Phi]=\frac12.
\]

The open-loop lazy-heat null is uniform and has
\(\mathbb E[\Phi]=7/6\). The stationary one-step KL separation is
`0.3748900964125387` nats per transition.

These values are consequences of the frozen transition tables, not values
placed in an initial quantum state.

## Stationary-Law Diagnostic

The companion diagnostic makes the dimension convention explicit:

\[
\pi_{\kappa}(R)
\propto d_R^{\kappa}e^{-\beta C_R}.
\]

With a reciprocal proposal \(Q\), the quotient-lumpable Metropolis--Hastings
acceptance rule is

\[
A_\kappa(R,R')=
\min\!\left[
1,
\frac{d_{R'}^\kappa e^{-\beta C_{R'}}Q(R',R)}
{d_R^\kappa e^{-\beta C_R}Q(R,R')}
\right].
\]

The frozen hypotheses are:

- `kappa = 1`: one quotient-visible dimension factor, the OPH edge-law
  diagnostic;
- `kappa = 0`: an unweighted sector Boltzmann null;
- `kappa = 2`: a full-block or Plancherel multiplicity null.

This diagnostic is identifiable only when some proposal edge joins unequal
dimensions. `S3`, `A4`, and the seeded random spectrum have positive KL
separation. `Z3` and `Z5` have \(d_R=1\) in every sector and therefore give
exactly zero separation. They are mandatory negative controls: no analysis may
claim dimension-exponent evidence from an abelian run.

Every Metropolis--Hastings edge circuit receives `192` shots on each backend
slot. The `kappa = 1` comparison against `kappa = 0` and `kappa = 2` is a
secondary conditional-likelihood diagnostic with Holm correction. It does not
turn a sampler into a source-law derivation.

## Dynamic-Circuit Patch

The reference circuit uses:

- three qubits for the bounded group or sector state;
- one decision/readback qubit;
- a mid-circuit state record;
- a decision record;
- record-conditioned state feedback; and
- a final state record.

Random lazy-heat and repair proposals are expanded into balanced circuit
variants. Every random choice is therefore visible in the circuit catalog and
receipt. Invalid three-qubit codes are leakage states and are never discarded.

The frozen Cayley catalog contains exactly five arms:

- `record_gated`, the contemporaneous self-reading repair instrument;
- `open_loop_heat`, the matched second-lazy-heat null;
- `delayed_record`, which gates on the pre-disturbance state;
- `shuffled_record`, which assigns the record to the wrong proposal; and
- `inverted_record`, which reverses the contemporaneous decision.

Each balanced Cayley circuit variant receives `192` shots on each backend
slot. Separate four-qubit computational-basis diagnostic circuits receive
`512` shots per prepared code on each slot. Those shot counts are sealed in the
public manifest and cannot be changed after seeing outcomes.

The primary circuit performs:

1. prepare a computational-basis mismatch state;
2. apply the frozen lazy disturbance variant;
3. read the disturbed state into the public record;
4. compute the descent decision from that record;
5. record the decision;
6. apply the repair update only when the decision is one;
7. reset the decision qubit; and
8. read the final state.

The open-loop arm uses the same state, record, slot count, and final readback,
but its second lazy heat step ignores the decision. Transpiled circuits must be
padded to a common duration envelope before hardware submission.

This is the OPH technology differentiator in operational form: a bounded patch
reads its own local state, writes a record, uses that record to choose a repair
move, and publishes enough evidence to verify the resulting process.

## Blinding, Sealing, and Decoys

Before any QPU job is submitted, one private reveal file must bind:

- semantic group and element labels to physical bit strings;
- record-gated, open-loop, delayed-record, shuffled-record, and inverted-record
  arms to opaque identifiers;
- physical layout assignment;
- development versus held-out backend;
- circuit family, horizon, and seed; and
- cyclic, nonabelian, and seeded-random decoy identities.

The seal is a two-phase hash chain. Phase one creates a production-random
32-byte secret, the complete opaque circuit catalog, normalized-OpenQASM 3
logical-circuit hashes, backend roles and layouts, shot counts, and a catalog
precommitment. Phase two
derives the anonymous candidate probability tables and hardened analysis lock
from that already sealed catalog. It then binds the analysis-lock hash into the
manifest core, commits the secret and private payload, and computes the final
public-manifest hash. Rebuilding any circuit, mapping, role, layout, candidate
table, or analysis rule must break at least one link.

The private reveal is written outside the checkout with mode `0600`; the
analysis lock and public manifest are written without overwrite, with the
public manifest committed last. Only opaque identifiers, hashes, resource
envelopes, and anonymous candidate tables are available to blinded analysis.
Reveal occurs only after every submitted job, failed submission, layout,
exclusion, blind report, and event-journal hash has been recorded.

Required decoys and controls are:

- cyclic `Z3` and `Z5` controls;
- nonabelian `A4` and `S3` families;
- a seeded random spectrum;
- shuffled and inverted readback controllers;
- a label/layout-only model; and
- a calibration-convolved noise model.

Conjugate relabelings of the same kernel are label controls, not independent
physics decoys.

Hardware execution is role-scoped and sequential across the two independently
selected, jointly sealed backend/layout slots. Development is planned,
explicitly submitted, harvested, and verified first. A held-out submission is
rejected until every expected development group is registered and complete.
Each submission invocation addresses only one role and one job group, and
submission and harvest events form append-only SHA-256 chains. The held-out
chip uses the same sealed manifest and analysis lock; it is not a
post-development redesign.

## Likelihood Contract

The primary endpoint is the pooled `S3` joint process likelihood across the
two sealed backend/layout slots for all heated-state, decision, final-state,
and leakage outcomes. The comparison is a conditional multinomial likelihood
ratio given a preregistered positive contamination model. It is not a Bayes
factor or posterior-odds calculation: finite-sample calibration uncertainty is
not integrated. The point channel and its frozen lower- and higher-contamination
sensitivity channels are all reported, and the worst-case sensitivity bound
must clear the relevant likelihood-ratio threshold.

The preregistration target is:

- conditional likelihood ratio strictly greater than `10` against every frozen
  null on each independent backend/layout slot, including every frozen
  contamination-sensitivity channel;
- pooled conditional likelihood ratio strictly greater than `100`, including
  every frozen contamination-sensitivity channel;
- the correctly revealed semantic mapping uniquely preferred after label and
  layout multiplicity correction;
- all transition rows inside the simultaneous `99%` calibrated prediction
  envelope; and
- no unlisted exclusion or post-reveal layout choice.

`Z5`, stationary-law, alternate horizon, and individual-layout results are
secondary. Their family-wise p-values use Holm correction. The primary pooled
endpoint is declared once and is not selected from the secondary family.

No shot is dropped. Invalid codes, incorrect readbacks, failed feedback, and
leakage remain explicit outcome categories. The separate `512`-shot diagnostic
circuits test the complete four-bit basis readout, with exact one-sided binomial
tests and Bonferroni correction against a `15%` maximum joint basis-code error
fraction. They are diagnostic-only validity controls: they are not promoted to
a factorized causal channel for the primary dynamic-circuit likelihood.

The primary contamination model and its sensitivity alternatives are frozen
before held-out execution and are never refit on held-out outcomes. A failed
joint basis-calibration control or more than `10%` leakage in any locked row
makes the benchmark invalid rather than unfavorable. No mitigation or
postselection can rescue either gate.

## Failure Rule

After the calibration and leakage validity gates pass, the reduced
repair-kernel claim fails if either of these occurs:

- both independent backend/layout slots favor the same frozen null over the
  record-gated kernel by at least `100:1` throughout the frozen contamination
  sensitivity range; or
- calibration controls pass while the complete record-gated transition law is
  outside the simultaneous prediction envelope on both systems.

The run is invalid, rather than negative, if the reveal commitment fails,
submitted jobs or layouts are omitted, a target-dependent exclusion is made,
or a post-reveal layout is selected.

Random labels or decoys matching the semantic controller as often as the true
mapping falsify specificity. Dependence on a single layout blocks replication.

## Unavoidable Quantum-Mechanics Boundary

The dynamic instrument is a standard finite measurement-and-feedback quantum
channel. Unrestricted ordinary quantum mechanics predicts the programmed
instrument. Therefore the experiment can distinguish record-gated repair from
the frozen open-loop heat, state-preparation-only, noise, label, and decoy
models. It cannot produce a likelihood ratio for OPH against unrestricted
quantum mechanics because unrestricted quantum mechanics contains the
programmed controller.

A passing result supports only this statement:

> IBM hardware realized a blinded finite self-reading repair instrument, and
> its preregistered process law outperformed matched open-loop, calibrated
> noise, label, and decoy protocols.

It does not establish an OPH physical vacuum, the Einstein branch, the
Standard Model branch, or cosmology.

## Executable Surfaces

- `programs/generative_repair_kernel.py`: exact kernels, stationary laws,
  identifiability, likelihoods, and deterministic receipts.
- `programs/record_gated_cayley_circuits.py`: local dynamic-circuit preflight
  for the primary and open-loop kernels.
- `programs/generative_repair_circuits.py`: secondary dimension-exponent
  Metropolis--Hastings instrument.
- `programs/blind_preregister.py`: opaque catalog, commit/reveal, and digest
  verification.
- `programs/blind_preregister_orchestrate.py`: two-phase production seal and
  atomic artifact writer.
- `programs/cayley_blind_likelihood_analysis.py`: frozen conditional-likelihood,
  calibration-sensitivity, validity-gate, and reveal-independent analysis.
- `programs/record_gated_cayley_runtime.py`: explicit role-scoped plan, submit,
  and harvest path with append-only event journals.
- `programs/runtime_analysis_packet.py`: fail-closed adapter from complete
  submission and harvest journals to the blinded analysis packet.
- `tests/test_generative_repair_kernel.py`: pure finite-matrix tests.
- `tests/test_record_gated_cayley_circuits.py`: optional Qiskit/Aer circuit
  tests.
- `requirements-ibm.txt`: pinned optional IBM test environment.

Only the explicit `submit` action in `record_gated_cayley_runtime.py`, with the
sealed manifest and analysis hashes confirmed, can submit an IBM job. The
preregistration, preflight, packet-building, and analysis surfaces do not
submit hardware work. This protocol itself records no hardware outcome.
