# Identifiability and calibration checks

These executable controls establish what the supplied action, detector and
geometry benchmarks can identify. They retain exact counterexamples,
independent finite-matrix checks and a retrospective empirical comparator.
They do not select a physical action or establish OPH-specific empirical
agreement. The [source-selection guide](../SOURCE_SELECTION.md) states the
corresponding physical obligations.

## Repair and consensus in an empirical comparison

The scored oscillator programs compile a supplied preparation, Hamiltonian,
environment and detector. Their predictions do not take a native C/E law
argument. Executing the same program with local history repair can test that
solver's implementation, but cannot by itself identify a microscopic repair
instrument. Three controls make the boundary executable.

### What a current-controlled work observer can distinguish

[rule_observability.py](rule_observability.py) checks the actual native maps,
including every matrix unit of the full 78-state carrier. For a common
coin/current record, each C Kraus map is a partial injection of load states;
E splits that same map by an additional contrast label. Extracting the native
diagonal blocks therefore gives the same unnormalized work/reference state
for both candidates. This is an operator identity, not an assumption that
the input is classical or uncorrelated.

Work-only instruments commute with that block extraction. Induction extends
the equality through any finite protocol whose seam choices, work controls
and stopping decisions use only common records. Initial native/work/reference
entanglement is allowed. The complete common-record and work/reference output
laws agree if the terminal native coherence remains inaccessible. This
includes work feedback selecting later native seams; merely adding that loop
cannot make a work-only detector select C or E.

The contract excludes coherent native/work coupling, native coherence reads,
coherent native recovery and feedback using E's extra label. The existing
native interference read is a negative control: it gives `5/8` versus `3/8`
and falls outside the work-only interface. The
[retained audit](rule_observability_results.json) does not assign a likelihood
to laboratory records with an unidentified native apparatus map.

### Completing the supplied apparatus history locally

[consensus_observability.py](consensus_observability.py) applies the existing
[finite acyclic history-repair argument](../common_source/HISTORY_REPAIR.md)
to unnormalized density operators and the actual anharmonic program. The
compatibility equations are fixed before solving:

```text
rho_root = prepared_state
rho_child = declared_CP_branch(rho_parent).
```

Each repair reads one parent and writes one child. A fair covering schedule
recovers the unique compatible branching law; a topological pass suffices.
Uniqueness concerns that law, not a selected random outcome. Intermediate
solver assignments are candidate histories, not altered experimental records.
No measured response or solved terminal state enters the local update.

The [retained execution](consensus_observability_results.json) repairs all
130 existing primary endpoints, retaining all four spin/detector outcomes at
each. Each endpoint is a separate prepared experiment; shared unread computed
prefixes do not turn these into repeated measurements of one evolving state.
The maximum full output-matrix discrepancy from chronological execution is
below `1.8e-14`. Scoring the same 125 positive-duration observations gives
RMSE `0.2876427868` in the declared position-estimator units, with no refit.
This is the existing empirical comparison through another solver, not an
independent dataset or additional empirical confirmation.

At the existing case-a `0.4 ms` endpoint, forward, reverse and shuffled
covering orders give the same four-outcome law; an unfair schedule does not
close. The primary and harmonic-control programs both satisfy every
compatibility equation, but their output laws differ by total variation
`0.16139288`. Agreement therefore supplies no criterion for choosing the action.
This is local completion on the supplied event graph, not an implementation
of the full observer consensus tower or a physical repair-clock calibration.
Its edges use the source-executed apparatus maps; individual native C/E seam
events have not been shown to implement the global history solver.

### The laboratory dephasing is not a parity-record rate

The experiment's [Appendix D](https://arxiv.org/html/2603.04744v1) supplies
number dephasing with `L=sqrt(18 s^-1)*N`. Its decay rates on occupation
coherences `(0,1)`, `(1,2)` and `(0,2)` are `9`, `9` and `36` per second.
An extra occupation-parity record leaves the last coherence unchanged.
Identifying these two noise mechanisms would change the observation model.

[dephasing_observability.py](dephasing_observability.py) checks a stronger
restriction: nonnegative mixtures of fixed occupation-diagonal binary
pinching generators obey `Gamma_02 <= Gamma_01 + Gamma_12`. Number diffusion
violates this inequality. This excludes that direct generator identification,
including arbitrary positive rates. The smallest possible maximum error among
these three rates is exactly `6 s^-1`, attained by `(15,15,30)`. This is a
generator comparison, not a statistical error bar. It does not exclude weak measurements,
noncommuting operations, hidden-state reductions or an effective emergence
of number diffusion. The [exact control](dephasing_observability_results.json)
does not reject either complete native rule using the oscillator data.

## Action and geometry ambiguities

The compressed four-level calibration is exactly insensitive to

```text
B(N) = N(N-1)(N-2)(N-3)/24.
```

Adding this local positive term leaves every operation inside that calibrated
space unchanged but changes its eight-level extension. Greater precision on
the same calibration cannot remove the ambiguity. The
[identification experiment](../identification_design/README.md) addresses it
with a coherent preparation involving level four and a phase-sensitive read;
its declared two candidate laws differ by total variation `0.7071067812`.
The action family and coherent services are supplied.

Likewise, the unit-edge cube graph admits a planar embedding with the same
existing graph dynamics and observation records. Inferring rank three needs
the declared orthogonal-cube metric convention or additional metric-sensitive
observations. Connectivity rank and cycle count alone are not a physical
spatial dimension.

The scalar sampling grids retain the same seven modes. Increasing their
sampling density changes a gradient coefficient but does not enlarge the
field-mode family. The exact-gradient endpoint is appropriate for fixed-band
comparisons; omitted modes require a separate enlargement and shared
observable domain.

[architecture_counterexamples.py](architecture_counterexamples.py) executes
these controls; [architecture_results.json](architecture_results.json)
retains their exact statements and numerical contrasts.

## Independent empirical comparator

The [anharmonic transfer](../anharmonic_transfer/README.md) tests the supplied
spin/pulse model against 125 positive-duration position-estimator values from
[McGarry et al.](https://arxiv.org/html/2603.04744v1), using their
[public deposition](https://zenodo.org/records/18752287). An independent,
ancilla-free comparison uses

```text
H = 2*pi*500*N + 4000*cos[(pi/3)*(a+a_dagger)+potential_phase],
L = sqrt(18)*N.
```

It uses the same transferred preparation, potential phases, literal `h=.4`
detector and published simulated-duration coordinate, with no additional
fit. DOP853 integrates the finite Lindblad equation without normalizing or
projecting the state. This is a supplied effective-Hamiltonian comparator;
it does not execute native source words.

Its aggregate estimator RMSE is **0.317681**, compared with **0.287643** for
the detailed pulse model: an improvement of **9.46%** for the latter. The
pulse model improves four of the five cases; the simpler model is better
in case e. The zero-prediction control has RMSE **0.645911**. The comparison
is retrospective, with observed responses available when it was designed.
It is neither a reserved-data test nor a distinction between OPH and
ordinary quantum mechanics.

The 32-versus-48-level target calculation changes an estimator by at most
`1.62e-7`; trace residuals are below `2.6e-15`. These are numerical checks,
not unrestricted-occupation certificates or apparatus-error bounds. The
[occupation certificate](../occupation_control/README.md) supplies a separate
controlled comparison for the specified pulse model.

The empirical inputs and complete output are retained in
[empirical_probe.json](empirical_probe.json). The compressed
[transfer reference](data/transfer_reference.json.gz) fixes the preparation
and comparison payload, independently of a live receipt regeneration. Its
decompressed SHA-256 is
`71403f474c39705d20a6c30dc2f36b7556943dd4f1734a5a1a395184744d2398`.
The manifest binds the compressed bytes; the producer validates the
uncompressed digest before parsing. No empirical values are modified.

## Measurement provenance and interpretation

All 130 deposited position estimates, including the five zero-duration
values, satisfy

```text
|63.6*x - round(63.6*x)| < 3.13e-8.
```

This printed quantization is compatible with hundredth-sized signals divided
by `.636`, whereas the documented conversion uses `sqrt(2)*.4`. The pattern
does not authenticate a detector gain, displacement coordinate, accepted
trial count or extraction procedure. It provides no basis for rescaling the
observations after inspecting residuals. The published `h=.4` comparison is
retained with its uncertainty explicitly unresolved.

Preparation transfer is also an empirical assumption: the largest
zero-duration estimator residual is **0.340550**. The largest positive-duration
residual is about **1.0823**, exceeding the selected ideal finite-model radius
of **0.001517**. That radius does not include apparatus uncertainty or removal
of the occupation cutoff. Refining the time splitter alone therefore cannot
settle the observed discrepancy.

Selection dwell, detector-pulse durations, phase tracking and preparation
validity need independent laboratory metadata. The current comparison uses
the stated ideal selection and readout interfaces. Aggregate means without
authenticated accepted counts do not support a statistical rejection.

For normalized characteristic functions, Hermiticity requires
`chi(-beta)=conjugate(chi(beta))`. For a real-displaced thermal state,
`chi(u+iv)=exp[-(nbar+1/2)(u²+v²)+2i*alpha*v]`: conjugating beta conjugates
the characteristic value. Exact reflected copies in deposited grids are
not independent observations. The [readout controls](../readout_controls/README.md)
and [two-order Gram constraints](../characteristic_constraints/README.md)
apply these identities with their stated detector and preparation premises.

## Reused mathematical results

The [finite environment derivation](../anharmonic_transfer/DERIVATION_ENVIRONMENT.md)
uses positive Gaussian-phase quadrature and CPTP contraction. The
[splitting and conditioning bounds](../anharmonic_transfer/DERIVATION_BOUNDS.md)
retain the finite oscillator commutator
`[a,a_dagger]=I-d|d-1><d-1|`. Replacing it by the unrestricted identity
invalidates a uniform finite-matrix bound. Complete accepted/rejected outputs
give conditional TV at most `epsilon/q`, where `q` is the reference
acceptance; small acceptance amplifies error.

[test_math_probes.py](test_math_probes.py) checks these algebraic ingredients,
the exact cutoff-boundary residual, finite forced-oscillator counterexamples,
and complete-source number/parity word reductions. The implemented source
words are documented in the apparatus derivations. The
[Etienney–Robin–Rouchon occupation theorem](https://arxiv.org/html/2501.09607v6)
is attached with a controlled trajectory integral in the
[occupation derivation](../occupation_control/DERIVATION.md). Agreement
between two cutoffs alone does not substitute for that theorem.

Finite instrument realization reuses
[Stinespring dilation](https://community.ams.org/journals/proc/1955-006-02/S0002-9939-1955-0069403-4/S0002-9939-1955-0069403-4.pdf)
and [universal finite gates](https://arxiv.org/abs/quant-ph/9503016), conditional
on admitted coherent services. Compilation shows that the chosen law can be
implemented; it does not choose that law over other compilable laws.
The [quantum-network framework](https://arxiv.org/abs/0904.4483) and
[operational Markov criterion](https://arxiv.org/abs/1801.09811) require
comparisons under the same allowed interventions, including surviving memory
and reference outputs. The [joint clock construction](../joint_clock_process/README.md)
provides an explicit native example where a reduced instantaneous field state
is insufficient but its history-dependent description is exact.

Relativistic reconstruction separately requires a compatible correlation
hierarchy with the hypotheses of the corrected
[Osterwalder–Schrader theorem](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-42/issue-3/Axioms-for-Euclidean-Greens-functions-II-with-an-Appendix-by/cmp/1103899050.pdf).
A finite GNS construction or positive transfer matrix does not establish
those hypotheses. The known [four-dimensional Ising/phi4 scaling result](https://arxiv.org/abs/1912.07973)
also constrains the specified conventional scalar families; it is neither
a general exclusion of OPH nor a proof of a nontrivial continuum limit.
See [observer-relative recovery](../OBSERVER_RELATIVE_RECOVERY.md) for the
current finite construction and remaining compatibility obligations.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.identifiability_checks.architecture_counterexamples
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.identifiability_checks.empirical_probe
python -m model.identifiability_checks.rule_observability --check model/identifiability_checks/rule_observability_results.json
python -m model.identifiability_checks.dephasing_observability --check model/identifiability_checks/dephasing_observability_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.identifiability_checks.consensus_observability --check model/identifiability_checks/consensus_observability_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/identifiability_checks
```

The tests retain exact architecture comparisons, source/data bindings,
corruption controls and the complete empirical replay. Numerical diagnostics,
ideal mathematical bounds and physical uncertainty have separate meanings.
