# Finite particle response and field-action selection

The finite constructions below determine response operators and complete
observation laws on stated source domains. Their joint mathematical content
does not identify a physical particle spectrum or an early-Universe population.
Each construction retains its own action, preparation, instrument and clock
hypotheses. A physical composition must realize those hypotheses on the same
process and compare the same records.

## Born weights and the executable model

The [finite effect-valuation theorem](../../Lean/EventAlgebra/FiniteBuschGleason.lean)
derives the Born trace formula and its unique density matrix from normalized,
nonnegative weights additive on every coexistent pair of complex matrix
effects. It includes dimension two and assumes no continuity axiom. The
formula itself is a conclusion under these hypotheses.

The [operational bridge](../../Lean/EventAlgebra/OperationalBornBridge.lean)
starts with complete finite outcome weights and derives their disjoint-event
additivity by finite sums. If equal-effect implementations have equal weights and each coexistent
pair has a disjoint common implementation, those weights descend to the valuation
required by the theorem. Refining the pair `(E,0)` supplies an implementation
of every effect. These universal operational conditions are explicit
hypotheses, rather than consequences of finitely many calibration scores.

The [native record construction](../../model/native_born_bridge/PROOF.md)
computes record probabilities from endpoint repair and supplied fair coins
without calling a trace-probability runtime. Complete record coarse graining
selects a unique diagonal state on the record algebra. It does not determine
the input state or select a full complex state.
The [implemented tomography](../../model/native_interface/tomography.py)
does include real and imaginary analyzers spanning all 4356 Hermitian
coordinates of the 66-state sector, with aborts and timeouts retained.
Its coherent routing and quantum state/instrument semantics are supplied.
The distinction is between an implemented complex measurement system and
an independent derivation of its universal operational probability laws.
The complex representation and universal operational effect-additivity
obligations retain their declared status in the premise register.

There is also a [source-derived full matrix functional](../../model/native_quantum_reconstruction/PROOF.md).
The complete-attempt kernel `P` is obtained from all 3960 native branches and
satisfies `2I/3<=P<=I`. The twelve occupation fields and connected transfer
paths generate all 4356 matrix units. Scaled closed-history counts determine
their moments, whose linear extension has the unique positive density
`rho_L=P^L/Tr(P^L)`. Its canonical complex-linear extension is additive by
definition and positive because `P` is positive, so it meets the hypotheses of
the Born theorem on this reconstructed algebra. Its entries are classical
endpoint probabilities divided by the return probability, every connected
positive-definite symmetric stochastic matrix gives the same result, and native
C and E evolution converges to the uniform state without preparing it.

The [current-only return adapter](../../model/native_thermal_return/PROOF.md)
implements the required configuration comparison with one initial restored
census and one final census; signed currents supply intermediate configurations.
At the stated reconstruction settings, at most 11398 attempts per experiment
give complete-law error at most `7.41e-18`; charging the rare-path and return
normalizations gives systematic error at most `3.76596e-6` for every
reconstructed effect statistic. Sampling uncertainty is separate. Arbitrary
effects require signed or complex combinations of statistics across
experiments. They are not all single native acceptance events, whose direct
configuration effects are diagonal. Choosing the complexification and
identifying reconstructed unitary time with physical time are distinct
requirements.

The [modular-history construction](../../model/native_modular_dynamics/PROOF.md)
reconstructs ordered modular correlations from finite combinations of these
same native path counts. Exact binomial polynomials have explicit truncation,
coefficient-amplification and return-normalization bounds. The displayed
two- and three-time examples have systematic polynomial bias below
`2.868e-10` and `4.302e-10`; sampling error is separate. All references use
`H=-log(P)`, with modular parameter `s` corresponding to conventional
Heisenberg time `-Ls`. Proper coordinate corners do not intertwine the full
modular flow. Retaining a chosen sector label supplies compatible block
flows for an explicitly completed reference. This changes the reference
and does not implement a native unitary operation or identify physical time.

## Thermodynamics from native repair

The [finite four-law core](../../Lean/Thermodynamics/CofinalSpectralTailFamily.lean)
derives finite thermodynamic laws from a conditional-resampling repair law
and its energy calibration. Cofinal low-temperature control is required for
the stronger regulator-uniform extension, not for this finite core. The laws
are conclusions of these hypotheses. A physical interpretation must establish
the hypotheses on the same source process.

The [binary native attachment](../../model/native_thermodynamic_bridge/PROOF.md)
obtains a direct conditional-resampling instance. On a fixed binary-load
sector, averaging the native seam's two fair-coin outcomes gives
`P_e=(I+swap_e)/2`. This is precisely uniform resampling conditional on the
outside-seam loads. The [finite involution theorem](../../Lean/Thermodynamics/NativeBinaryThermodynamics.lean)
connects this form to the conditional-repair library. The connected seam
ensemble has one uniform stationary reference, which also uniquely maximizes
entropy in that sector. One seam alone leaves the masses of its separate
fibres undetermined.

For the unread population, `S(pP_e)-S(p)=D(p||pP_e)>=0`. Fixed seam words
telescope this identity. Independent random scan includes an additional
nonnegative mixing term; its kernel need not be idempotent. Retaining the
complete endpoint and coin/current tape gives
`H(X_n,J)=H(X_0)+n log(2)` for a fixed word, so the marginal entropy increase
does not imply loss of the initial information from the complete record.
The source, binary sector and fair-coin/exogenous scheduling assumptions are
explicit. Unbalanced configurations and record-dependent feedback require
their own analysis.

This closes the common-reference and local-resampling attachment on this
restricted population branch. Its modular energy is constant, so it cannot
identify a nontrivial temperature, heat reservoir, physical erasure cost or
cooling protocol. The physical energy/temperature and cofinal-refinement
conditions of the physical four-law interpretation retain their separate status.

A [native observer-memory construction](../../model/native_record_thermodynamics/PROOF.md)
supplies a nonconstant internal calibration. A stopped coin-prefix program
records `0`, `10` or `11` after respectively one, two or two native attempts.
Its weights are `w=(1/2,1/4,1/4)`; the label reads only coins, so they ignore
the hidden configuration by construction. Define the stored-prefix cost
`E=n log(2)`. At fair coins this equals label surprisal for every complete
prefix code, so `w=exp(-E)` holds by construction, and overwriting the memory with a fresh prefix
is exactly its global heat bath in at most two native attempts. The
[typed finite-theorem instance](../../Lean/Thermodynamics/NativeRecordFourLaws.lean)
constructs `RepairLawData` and `CalibrationData` from these counts and applies
`finiteFourLawCore`. That core holds for every full-support law with
`E=-log(ref)`, and it formalizes no contact relation or energy exchange with a
second system. Its gap is `log(2)` with a unique ground label.

Agreement of `b` sequentially generated independent prefixes gives weights
`w^b`, hence the Gibbs family of the same cost at integer multiplier `b`, as it
does for any law and its surprisal.
The coin law is independent of source state, so separate source preparations
are unnecessary. Keeping the old memory after disagreement gives a lazy
reversible kernel with exactly that stationary reference. Bounded retries
retain all rejection outcomes and give an explicit geometric error; almost
sure retry gives the exact heat bath with finite expected cost. This supplies
finite thermodynamic premises for the declared observer-memory protocol.
Stored-prefix cost is not total implementation work, laboratory energy or a
physical temperature calibration.

Independently, complete-attempt closed histories give
`rho_L=exp(-L H)/Tr(exp(-L H))` for `H=-log(P)`, which is reconstructed before
the loop length `L` is chosen. The return adapter also supports a distinct
incident-clock transfer and its permitted regional field dictionary.
Both comparisons include census failures, nonreturns and any waiting cap.
These Gibbs functionals describe reconstructed history ensembles; they are
not coherent working-state preparations or stationary Gibbs laws of the
unchanged open source process.

## Regional entropy and the Einstein premises

The [regional response construction](../../model/native_gravity_response/PROOF.md)
uses the same native length-two reference `sigma=C^2/171720`. Midpoint
filters give 65 independent normalized, faithful state directions. These
are postselected history ensembles with complete return and filter weights;
no causal quantum preparation is inferred. Each native seam has regional
sector dimensions `1,2,1` and exterior multiplicities `45,10,1`. Its canonical
tracial completion obeys the exact bulk/edge entropy split, with the central
Shannon term in the bulk. The modular first law has a quadratic remainder
bound for the same states, including noncommuting perturbations.

The full reference differs from that completion and has cross-sector
coherence, which rules out a reference-preserving retraction onto the chosen
regional algebra. A native pair-occupation filter preserving all twelve mean
occupations has a strictly positive bulk-entropy slope, while the fixed
geometric cut count has zero derivative. Thus these finite entropy premises
do not supply geometric stationarity. Physical stress, responsive area,
fixed proper volume, Newton normalization and a common small-ball refinement
require their own derivations.

## Carrier from directional readout

The [weighted carrier selector](../../model/carrier_selection_reduction/PROOF.md)
starts with a supplied real three-dimensional response and directional
quadratic measurements `tr(u u^T X)` of every symmetric response `X`.
Its positive resource weights sum to one. Informational completeness needs
at least six distinct unoriented channels. The smallest eigenvalue of the
weighted information operator is at most `2/15`; at the minimum channel
count, attaining that value forces equal weights and squared inter-axis
overlaps `1/5`. Exposing both signs gives twelve icosahedral ports. Positive
overlap adjacency has degree five, thirty edges and twenty triangular faces.
Those counts are derived without a spherical triangulation, Euler-defect
cost, preassigned port count or group action.

A bound within `10^-7` of the information optimum forces the same discrete
incidence through an integer sign-matrix gap. It does not force exact angles
or exact metric symmetry. The response dimension, quadratic measurement law,
minimum-resource and worst-direction information criteria, and inverse-port
convention are premises. The native rank-three construction depends on a
declared carrier and cannot independently discharge this alternative
selector's dimension premise.

## Source-counted localization and spectrum

The [native spectrum and localization search](../../model/native_quantum_defects/PROOF.md)
retains all 66 transfer modes and classifies every one of the 4096 port
cores. Six ten-port cores support three-dimensional confined amplitude
spaces with transfer eigenvalues `4/5`, `5/6` and `13/15`; no such space
fits on nine or fewer ports. These modes fail the reference-like exterior
test, and confinement is lost under admitted changes of selector weights.
Complete C/E record commitment also changes their transfer-sector weights.

A distinct construction starts from the faithful reconstructed
`rho=P/Tr(P)` and changes its one-port phase by `D=(-1)^n_a`.
The functional `D rho D` has exactly the same eleven-port exterior and
configuration probabilities, but positive `H=-log(P)` energy excess in
`[1/1008,1/672]`. Every complete native history preserves equality of the
exterior-plus-tape laws relative to the correspondingly evolved reference.
The contrast decays by `5/6` per C attempt or `7/10` per E attempt under
uniform scheduling. This is a count-derived localized functional; no native
phase gate, prepared working state, protected moving particle or physical
mass is established.

## Mass comparison and inverse coordinates

The [particle comparison](../../model/particle_spectrum_audit/PROOF.md)
independently evaluates the specified charged circulant and reproduces the
multiloop quark, W, Z, Higgs and top outputs through pinned SMDR. Its
[receipt](../../model/particle_spectrum_audit/results.json) and
[reference table](../../model/particle_spectrum_audit/pdg2026_reference.json)
state the distinct comparison observables. The charged coefficients are read
under an explicit pole-identification premise. Numerical resemblance does not
establish that premise, and the architecture has target-informed ancestry.

All twelve nonzero comparison outputs are outside their printed marginal
intervals. The charged masses differ by about minus 70.5 ppm; the five non-top quark
relative discrepancies range from about 23% to 209%. A common normalization
does not place the electron and muon in their printed intervals. Neither this
arithmetic comparison nor the printed experimental intervals supply a total
perturbative error bound or a covariance-aware joint significance.

For a coefficient `m=A*6^(-n)` and ratio `r=output/reference`, holding its
transport factor fixed gives `n_required=n+log_6(r)`. The u, d, s, c and b
coordinates are approximately 6.63023, 6.19375, 4.52837, 3.11653 and 2.25334.
Threshold dependence and self-scale mass definitions require a coupled rerun
for an actual changed boundary. These diagnostic coordinates are not a
selected replacement texture. Matched five-flavour running gives
`alpha_s(91.1879 GeV)=0.1194406595324884` for the stipulated branch; a
reference-ratio estimate is a different operation.

## Exact finite reductions

The [finite response proof](../../model/particle_finite_response/PROOF.md)
reduces the neutral transfer to a bright-sector quartic and its dark sector,
with detector backaction and transient-memory bounds retained. For its
specified unshifted amplitude-to-mass attachment, `0<=p<=1/12` and `0<h<=1/2`,
the normal-order squared-gap ratio is at most `0.004987197763...`, below the
comparison value near `0.03108`. This excludes that attachment. It does not
exclude arbitrary neutral actions or determine absolute neutrino masses.

The same proof reuses harmonic extension to obtain the passive seam metric
`Z>=I`. Its conserved family direction has zero action before and after
positive normalization. The passive construction therefore cannot generate
three nonzero charged masses. Scalar fixed-algebra conclusions retain their
particular occupation and autonomous-channel hypotheses.

The [source-selection proof](../../model/particle_source_selection/PROOF.md)
shows that complete classical record instruments do not select individual
coherent Kraus phases. Its explicit same-rate counterfamily has identical
classical instruments and different coherent spectra; it is separate from
the native 128-record law. Independently, a real symmetric three-dimensional
kernel averaged over a regular real three-cycle has a repeated eigenvalue.
Neither result excludes a separately specified coherent source and reference.

## Complete records, memory and preparation

The [native record proof](../../model/particle_record_equivalence/PROOF.md)
constructs a four-address encoding of the twelve-address family memory on
its recurrent support. Every native raw outcome and old-memory label is
retained. The branch maps intertwine for arbitrary admitted memory, work and
untouched-reference inputs, proving equality of every finite joint history.
Initial address zero uses a retained first-event startup rather than a
stationary substitution. This construction establishes one exact capacity
encoding, not independence of arbitrary physical capacities or preparations.

Ideal recorded recovery preserves the logical process up to known branch
weights on both matching codes, both repair laws and all thirty seams. It
cannot shift logical poles. Always-duty and C-duty kernels have different
normalized quadratic traces, so changing the duty is a change of action.
The nine-transition preparation bound `1/2579890176` and the two-clock
cancellation control a first-order kernel approximation, not a physical
mass uncertainty.

The [all-lag memory proof](../../model/particle_memory_response/PROOF.md)
uses the exact group inverse `R=(I-T+Pi)^(-1)-Pi`. Its complete one-sided
correlation sum is `K=(diag(w)-w w^T)/2+E R diag(w)`. The antisymmetric part
gives a coherent ordering potential; the symmetric part gives positive
diffusion. Exact three-ray reduction replaces a lag cutoff for stationary
product preparation and unread output. The corrected completely positive
semigroup has an explicit rational `O(h^2)` finite-time diamond-norm bound,
including an arbitrary untouched reference. It is not a full-record
semigroup equivalence.

The [prepared readout proof](../../model/prepared_quantum_readout/PROOF.md)
uses nineteen sweeps in a fixed seam order, 570 native attempts, to remove
independently supplied exact `I/66` preparation and initial binary support.
It admits arbitrary total-two inputs on the full 78-state carrier, including
entanglement with an untouched reference. Its preparation error is below
`3.72e-8`; every one of the 8,195 cap-64 suffix records, including failed
heralds and timeouts, participates in the comparison. The warm-up tape is
ignored without conditioning. The separate compiler allowance yields a
complete comparison bound below `3.03719e-6`. The total-two sector, quantum
instrument semantics, coherent services, work action and clock remain inputs.

The [finite calibration theorem](../../model/finite_instrument_calibration/PROOF.md)
replaces projective-effect equality and invariance of every certain state by
finite joint tests. For a complete Kraus instrument, prepare each vector of
a trusted basis and one uniform coherent vector per nontrivial outcome block.
If the assigned outcome and a trusted return-to-input test occur with unit
joint probability, every outcome map is its Lüders compression. A rational
error bound controls the complete instrument with arbitrary untouched
references, and hybrid replacement controls adaptive finite histories.
Wrong labels and failed returns enter the errors without postselection.
Exact examples with nonprojective effects exercise the approximate result.
The trusted preparations/readbacks and quantum probability semantics are
assumptions; finite zero-failure counts do not establish zero probability
error or authenticate the apparatus.
The [single-channel core](../../Lean/EventAlgebra/FiniteProbeInstrumentSelection.lean)
is Lean checked: unit returns on the coordinate and uniform pure states force
the identity channel. The complete multi-label and stability results have
analytic proofs and exact controls.

## Gauge bands and local field dynamics

The [gauge-band construction](../../model/particle_gauge_band/GAUGE_BAND_PROOF.md)
retains the full 343-site classical-link action and pairing. Its separated
band has 180 canonical annihilation modes. Exact bounds establish polar
transport between thirteen raw branches; numerical checks cover 91 unordered
projector comparisons and two representative polar transports. Localization
bounds control exterior/core probability.

A supplied sixty-state rotor completion verifies 1,416 covariance identities,
while leaving five link-energy differences unselected. It changes the native
raising-only link law. The [native quantum-link witness](../../model/particle_gauge_band/NATIVE_LINK_CAR_WITNESS.md)
has noninteger naive Fock occupation even on a Gauss singlet. Canonical
classical-band occupation therefore needs a separate construction for
operator-valued links.

The [local Maxwell selection theorem](../../model/maxwell_local_selection/PROOF.md)
classifies every real symmetric gauge-invariant seam Hessian with strict
face-sharing support on the committed complex. It has the form
`C^T diag(w) C`, with twenty independent face weights. Covariance under the
proper carrier rotations forces a common weight and removes an invariant
linear term. A nonzero, nonnegative quadratic part gives
`S=constant+kappa ||C A||^2/2`, with `kappa>0`. A direct triangle proof and
exact rational rank checks establish completeness. Strict face locality is
stronger than generic finite-range locality, and nonlinear actions are
outside the theorem.

Among real linear autonomous evolutions, positive specified weights
`epsilon,beta`, Faraday kinematics and
conservation of `(epsilon ||E||^2+beta ||B||^2)/2`, together with the stated
electric-even/magnetic-odd reversal, uniquely give
`dot E=(beta/epsilon) C^T B` on closed magnetic flux. Gauss data are preserved.
Removing reversal admits additional skew electric terms. The coefficients,
field interpretation, native implementation and physical clock are not
selected, and the theorem does not choose a discrete integrator.

The [smooth weak-field theorem](../../model/maxwell_weak_field/PROOF.md)
derives that quadratic form and linear generator as leading response without
requiring an exactly quadratic action or linear law. For a smooth autonomous
field law it assumes a smooth rotation-invariant energy with nonzero diagonal
electric susceptibility and strict face locality of the pulled-back magnetic
vacuum Hessian, Faraday kinematics, conservation and the stated reversal.
The positive Hessian also forces the vacuum to be an equilibrium. In its
energy norm, a supplied derivative bound and energy confinement give
`||x(t)-exp(tG)x0|| <= (L/2) J0 r^2 |t|`. A local quartic model improves
the amplitude dependence to cubic and preserves exact Gauss data. General
nonlinear conservative reversible flows need not preserve Gauss: the exact
counterexample distinguishes the leading Maxwell constraint from a full
nonlinear conservation law. The finite result supplies no native physical
clock, field identification or continuum limit.

The [full-sector interaction test](../../model/native_modular_dynamics/interaction.py)
uses the unchanged native occupation dictionary and all 66 states. Every
quadratic one-body operator has zero four-pair diagonal contrast, even with
arbitrary hoppings, on-site potentials and a scalar offset. The native
`J=I-P` has contrast `-1/15`, giving operator-norm distance at least `1/60`
from this entire quadratic class. Exact logarithm bounds give distance
greater than `0.01984562` for `H=-log(P)`. The unit-rate Poisson generator
`J` and discrete-transfer logarithm `H` use different clocks. The result
establishes finite interaction in this dictionary, without an interacting
chiral relativistic limit or an exclusion of arbitrary re-encodings.

## Physical identification

Exact equality of complete observations is a stronger requirement than
matching response matrices, isolated spectra or marginal outputs. Compatible
preparation and per-outcome intertwining give a sufficient finite criterion;
determinism alone does not construct the required comparison map.

Pressure likewise requires a common energy and operational strain law.
The existing [equation-of-state theorem](../../Lean/Thermodynamics/EquationOfState.lean)
exhibits fixed-volume energy histories compatible with different derivatives:
`H_alpha(V)=Q(V/V0)^(-alpha)` gives `p/rho=alpha`. Constant changes of volume
and energy units leave `p/rho` invariant. Consensus about identical
fixed-volume records cannot select their unobserved strain derivative.

The [operational strain theorem](../../model/operational_strain/PROOF.md)
identifies that derivative from three geometrically spaced energy readings
within the held-state family `H(V)=C+A(V/V0)^(-alpha)`, with `A>0` and
`alpha>=0`. Successive energy differences determine the contraction and
exponent; they also recover the amplitude and constant offset when the work
is nonzero. Three settings are generically minimal when all three parameters
are unknown. Exact logarithm intervals propagate meter error into pressure.
For general response, a supplied curvature bound controls a two-point
pressure stencil; without a response class or regularity bound, positive
smooth counterexamples defeat every finite collection of energy readings.

A complete two-attempt native E census reads the initial charge projectors
on arbitrary coherent/reference inputs. It thereby measures independently
calibrated charge-diagonal generalized forces. It does not supply the strain
interventions or their energy readings. An additive energy offset changes
`p/rho=alpha A/(C+A)` while leaving pressure unchanged, so identifying an
absolute energy density requires more than identifying work.

These finite results neither supply a source-selected physical mass spectrum
nor a common interacting chiral 3+1-dimensional field and Einstein
realization. They state the controlled constructions and the same-process
attachments needed to make those comparisons meaningful.
