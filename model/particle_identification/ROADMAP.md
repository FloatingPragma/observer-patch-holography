# From the executable source to a particle spectrum

The target is one source-derived effective field theory whose prepared
experiments reproduce particle production, propagation, scattering and decay.
An eigenvalue list, a fitted mass formula, or several compatible finite
examples does not establish that target. The current interacting extension
provides a complete finite gauge–lepton–Higgs control and its detector records.
It does not identify its response lines with nature's particles.

Reading the simulator, observer, scheduler, repair and memory specifications
is allowed. Those specifications define the mathematical model and can be
used in a derivation. A laboratory identification, a selected initial state,
or a continuum limit must still be proved from them; knowledge of the code
does not turn a separately chosen effective Hamiltonian into its output.

## Required interfaces and their order

The interfaces below concern the same source family, observer algebra,
state and clock. A result for a different prepared state or action cannot
silently supply a missing interface. Some work can proceed in parallel,
but a physical particle claim depends on the entire relevant chain.

| Interface | Established starting point | Required next result and acceptance condition |
| --- | --- | --- |
| 1. Source and capacity family | [Capacity definitions](../observer/CAPACITY.md), the executable repair law, and source-bound finite packets distinguish record capacity, global capacity and Hilbert dimension. | Give the actual family of source programs and states indexed by physical capacity. Derive every dependence of local couplings, preparation and clock on that parameter. A chosen lattice size or number of repeated trials is not such a derivation. |
| 2. Observer state and complete experiments | [Native quantum reconstruction](../native_quantum_reconstruction/PROOF.md) constructs a positive matrix functional. [The Born bridge](../native_born_bridge/PROOF.md) separates finite consistency from the all-effect theorem. The interacting extension supplies complete Kraus maps. | Show which effect/instrument family the observer implements, including failures and back-action, and prove its equivalences and composition laws on the same state. Supply a source preparation or a controlled approximation to its local restriction. |
| 3. Quantum observables and uncertainty | [Robertson](../../Lean/EventAlgebra/Robertson.lean) proves the finite-state commutator inequality. [Prepared-clock fluctuations](../predictions/README.md#shared-state-fluctuation-field) supply a conditional quantitative canonical limit. | Bind actual source observables to the finite theorem. For physical position and momentum, additionally identify translations, spatial readout, the limiting domain and the common action unit. Merely assigning two matrices those names is insufficient. |
| 4. Common reference, clock and energy | [Complete interacting records](../interacting_extension/PROOF.md) include the reference energy and cancel a common gauge offset. [Native records](../interacting_extension/source_selection/PROOF.md) give exact correlated attempt counts. | Derive the reference and clock used by all sectors, including preparation and comparison records. Show how dimensionless rates relate to one another before assigning units. A capacity count alone does not fix seconds or GeV. |
| 5. Spatial field and relativistic limit | [The continuum CAR net](../continuum_field_net/README.md), [canonical spin field](../spin_canonical_field/README.md), and [free mass shells](../source_mass_shell/README.md) provide different controlled ingredients. | Recover local observables and translations in three spatial dimensions with a common time, positive energy and the required relativistic covariance. Prove the relevant convergence of joint records, not just agreement of a low-momentum dispersion curve. |
| 6. One interacting matter and gauge theory | [The finite Higgs join](../interacting_extension/higgs_join/PROOF.md) retains the complete 438/1596-state lepton sectors. [The local electric completion](../interacting_extension/chiral_electric/PROOF.md) includes the changing-projector kinetic term. | Extend the common carrier to the required colored matter and scalar sectors; derive its action and coupling placement from the source. Establish consistent chiral gauge dynamics, locality, anomaly cancellation and a nontrivial interacting limit on that same carrier. Separate free, abelian and local-projector controls do not compose automatically. |
| 7. Controlled finite calculations | [Local simulation](../full_fock_field/LOCAL_SIMULATION.md) bounds complete local histories independently of ambient volume for its fixed supplied action and initial local restriction. | Control initial-state, local-coefficient, clock, volume, lattice-spacing and occupation errors together. Prove uniformity on the capacity family or an explicit finite-capacity correction. Increasing a cutoff without an error theorem is not an extrapolation certificate. |
| 8. Physical particle definitions | Full finite response measures retain signed frequencies, residues and all complementary states. Conditional representation results constrain possible sectors. | Identify observable species through charge, spin, propagation and interactions, using the appropriate spectral or scattering definition. Explain missing and extra sectors. A topological defect is a particle candidate only if the same dynamics and records establish that identification. |
| 9. Masses, mixing, widths and scattering | [Contact scattering](../contact_scattering/ONSHELL.md) has complete lattice wave operators and a unitary scattering matrix. [Moving contact scattering](../moving_contact_scattering/README.md) has a controlled nonrelativistic two-particle continuum. [Packet readout](../scattering_packet_readout/README.md) connects an incoming packet to a bounded finite detector contrast and executes a separate small source control. | Extend these results to the selected relativistic interacting theory and its common vacuum; the existing fine-spacing packet schedule has not been executed. Calculate gauge-consistent correlation and scattering observables with controlled errors, as specified in [the scattering interface](SCATTERING.md). Fix renormalization conventions where needed. Derive masses, mixing and decay laws from the selected action/state; do not choose couplings or labels to reproduce the target values. |
| 10. Experimental comparison | The repository distinguishes simulated records, conditional calibrations and measured inputs. | Compare the predicted joint detector distributions with specified experiments using the same observable definitions, units, resolution and uncertainties. Collider rates also require the incoming state and luminosity; hadron collisions require the relevant hadron structure. Agreement of selected masses alone does not establish agreement of the theory. |

Interfaces 1–4 constrain what counts as an admissible field construction.
Interfaces 5–7 constrain which finite calculations approximate it. Interfaces
8–10 determine whether its excitations are the measured particles. A
source-derived cosmological abundance additionally needs the early state,
expansion and subsequent dynamics; a spectrum by itself is not a population
distribution. Recovering responsive Einstein dynamics is a further common-state
geometry task, with the same stress and reference data, rather than a
consequence of finding particle masses.

## The central common-model test

For each admitted finite experimental menu `M`, specify the source
preparation/encoding `E`, the complete source process `S_M`, the proposed
field process `T_M`, and the observation map `R`. Both outputs include every
retained record and the remaining accessible state. A sufficient strong
comparison on the declared input domain is

```math
\frac12\|R\circ S_M\circ E-T_M\|_\diamond\leq\epsilon_M.
```

An energy-constrained or explicitly restricted input domain may be used
instead, but must be stated and preserved by the experiment. A theorem for
one prepared state is weaker than the diamond-norm statement and should be
reported as such. Adaptive choices require a bound for the complete history,
not just its individual one-time marginals. The maps, common clock and error
budget must come from the same source family. For a continuum claim the
admitted menus must cover the claimed observable algebra and the relevant
errors must tend to zero in the specified limit. A finite-energy experimental
claim can instead retain an explicit nonzero error budget.

This is the missing bridge between the current native process and the chosen
interacting gauge–Higgs Hamiltonian. The
[repair-action calculation](../repair_action_audit/README.md) already proves
that the reduced native repair channel is dissipative, so simply naming it
a Hamiltonian does not establish this bridge. Retained environments and
observer memory can matter; their operations must be included explicitly.
Similarly, the [phase-selection proof](../particle_source_selection/PROOF.md)
shows that complete classical instrument maps alone do not determine a
coherent lift's phases. Full simulator specifications may supply additional
coherent structure, as allowed here. The required argument must point to
that structure and prove the comparison, rather than infer it from classical
event rates.

The [passive-clock bound](PASSIVE_CLOCK.md) makes one failed shortcut
quantitative. After 200 unread repair attempts, any static encoding and
decoding of a logical qubit through one native binary carrier remains more
than `0.74` away in half-diamond distance from every unitary channel. This
uses the full quantum contraction, not just population mixing. It excludes
neither access to the retained records nor active recovery or coherent
memory. Such services must be part of the actual common-model construction.

The fixed program's data flow makes this distinction concrete.
[`record_clock_member`](../native_program/members.py#L39) constructs the
66-state native clock experiment and calls `learned_flow.source_program`
with record-inferred relative weights. In the
[`interpreter`](../native_program/interpreter.py#L328), `_choices` selects
seams from the declared weights and prior current, and `_apply_native`
executes the endpoint law. These entry points do not emit the interacting
Higgs carrier, its Hamiltonian or a preparation of its ground state.

The [Higgs response entry point](../interacting_extension/higgs_join/run_joined_response.py#L183)
instead passes electric, magnetic, hopping and scalar-energy coefficients
equal to one, with fermion and reference energies equal to zero. It invokes
`block_eigh` and uses the lowest spectral columns as its prescribed initial
state. In the [carrier builder](../interacting_extension/higgs_join/interacting_carrier.py#L421),
`primary_family_matrices` reads `Y,D,J` and sets the ordered-control angle to
`log(5/4)`; `build_join` accepts the family matrices and scalar controls and
reports `source_selected_scalar_parameters=False`. The
[family-reference construction](../interacting_extension/dependencies/source_reference.py)
binds its native geometry and memory data while declaring its auxiliary
preparation and ordered feedback services. These readable specifications
determine the configured benchmark. A bound on finite detector information
does not prevent an observer from inspecting and calculating from them.

Three source maps are therefore required:

1. From the continuing observer/scheduler/memory program to this carrier,
   coupling placement and action, including the scalar cap and reference.
2. From the actual continuing source state to the common interacting initial
   state, with its preparation records and a controlled error.
3. From an emitted complete source pulse/readout history to the same
   Hamiltonian and probes, with a common clock and an error bound for every
   retained outcome.

The [active native code](../common_source/ENCODED_MEMORY.md) supplies exact
ideal record-assisted recovery and explicit coherent source services on
its stated code. It demonstrates why the unread-channel obstruction does
not exclude active or record-retaining implementations. Storage capacity
and the ability to program a finite simulation do not establish that the
fixed observer program emits this particular interacting action or state.

## State and clock selected by the native member

There is a derived state and a compiled dimensionless clock for the named
native observer. For its connected binary total-load-two carrier, positive
port weights `theta_i`, fixed positive seam selection and declared biased
endpoint law, the [source mixing proof](../local_repair_family/README.md)
gives the unique quantum stationary state

```math
\rho_* = Z^{-1}\sum_{x:\,\sum_i x_i=2}
             \left(\prod_i\theta_i^{x_i}\right)|x\rangle\langle x|,
\qquad
-\log\rho_* = (\log Z)I-\sum_i(\log\theta_i)n_i.
```

The sum is over binary occupation vectors. The proof controls arbitrary
reference-entangled inputs, including decay of off-diagonal blocks. It
therefore supplies more than a stationary population fit. For the declared
`25/8` member, the fixed-law preparation bound after 1600 unread attempts is
below `10^-6`. This concerns the source marginal and an untouched reference,
not replacement of the full growing archive by an independent state.
Stationary finite-FIFO posteriors retain their actual correlations and have
the restricted modular action on their support corners.

The [record-inferred implementation](../local_repair_family/record_clock.py)
prepares the continuing terminal training carrier and compiles the inferred
potential `-sum_i log(theta_hat_i)n_i`, with all future outcomes retained.
Its retained complete-output bound against the true stationary-reference
clock is below `0.093762`; the comparison with the actual finite-prepared
marginal's clock is a distinct statewise bound below `0.104425`. These
dimensionless error statements include the stated inference and apparatus
budgets. They do not equate modular parameter, repair-attempt count and
laboratory time.

The [joint count clock](../native_clock/clock_bridge.py) derives the
stationary ratio of `5/33` local movement ticks per incident hit and a
correlated finite-count error bound. This relates specified record counters;
it does not identify an interacting Hamiltonian with the modular generator.
The [energy-identification criterion](../energy_identifiability/README.md)
states which coherent comparisons determine relative generator offsets.
Its full-record invariances concern the declared experimental menu. Reading
a supplied Hamiltonian in the simulator specification remains allowed and
fixes its configured offsets without providing a physical field identification.

Thus the missing selection is precise: this source-derived state and
occupation-potential clock must be related to a common interacting field
state and action by an explicit observation/encoding theorem. The Higgs
benchmark's scalar, hopping and gauge terms do not appear merely by taking
the logarithm above. A chosen compiler program that simulates those terms
establishes an implementation conditional on its program input; it does
not establish that the continuing native observer selects that input.

The [direct-map tests](../source_state_bridge/PROOF.md) make two necessary
conditions explicit for the declared member. Its stationary eigenvalues
are `625/7505`, `200/7505` and `64/7505`, with multiplicities 1, 20 and 45.
The full isometric image of this state has trace distance at least
`1376/1501` from every pure full-state target. A pure field subsystem can
instead be prepared while its entropy is carried by another quantum
system; that system and the preparation operations require an explicit
source realization.

There is also a separate all-input carrier bound. If all quantum information
crosses a single cut through a 66-dimensional carrier, with arbitrary
classical records but no additional quantum path or entanglement assistance,
the decoded channel has half-diamond distance at least `62/73` from the
438-dimensional identity and at least `255/266` from the 1596-dimensional
identity. This is a test for faithful full-carrier encodings, not a lower
bound for a selected ground-state experiment or restricted measurement menu.
The named observer's declared quantum work services are outside this
66-dimensional bottleneck when they carry additional quantum information.
These results therefore reject specific direct identifications while
leaving active preparation and larger encoded constructions to be tested.

## Spatial constructions and excluded shortcuts

The [native density-field analysis](../physical_space_time/OCCUPATION_FIELD.md#a-complete-linear-readout-preserves-the-decay-spectrum)
restricts one proposed bulk reconstruction. For the complete centered
Gaussian density/Green observation space on the spherical refinement, a
fixed linear readout whose composite with the covariance square root is
bounded, injective and has dense range preserves all covariance-normalized
lag eigenvalues. A continuous nonconstant common clock change preserving
the stationary semigroup must be a positive linear rescaling. In this
domain the native decay rates grow linearly with their ordered index,
whereas the complete Euclidean free scalar on a smooth compact
three-dimensional space has cube-root growth. Such a readout and clock
cannot identify these two complete processes.

This result does not concern the separately prepared finite Higgs carrier.
It does not exclude a restricted boundary algebra, a nonlinear dictionary,
an interacting limit, another scaling family or a holographic reconstruction
with different hypotheses. Each route needs its own map of complete
observations and errors. Adding a radial coordinate or matching only an
equal-time covariance does not prove that map.

The [active product-field implementation](../observer_time_fields/REFINEMENT.md#one-product-field-state-and-active-clock)
executes an actual screen readout and a declared angular/radial field clock
with retained residuals and complete outcomes. Its radial grid, spectral
encoder, prepared state and protected services are declared inputs. The
retained cell-visible encoder is a finite packet with a nonzero spectral
residual; no growing nonconstant-mode encoder is executed. The
Higgs benchmark does not consume this spatial readout or clock. The
[twelve-null-direction criterion](../native_lorentz_transition/PROOF.md)
and [principal-cone calculation](../field_cone_compatibility/PROOF.md)
test candidate charts and field coefficients; neither selects the missing
common dictionary.

## Capacity does not equal the size of the calculation

No current theorem establishes that the three interacting response frequencies
are functions of physical global capacity, or that they become the observed
lepton masses at a particular capacity. The current calculation has declared
couplings, scalar cutoff, graph, preparation and model-time unit. Changing
those inputs is a change of model unless the source supplies the change.

Keep these parameters distinct:

| Quantity | Meaning |
| --- | --- |
| Global capacity | The source's global information-capacity parameter, when its defining antecedent is supplied. |
| Record capacity | The logarithm of the relevant number of admissible records. |
| Hilbert dimension | The dimension of a specified representation or retained carrier. |
| Spatial volume | The extent of a chosen spatial description. |
| Lattice spacing and occupation cutoff | Regulators of a particular effective calculation. |
| Number of copies or shots | Resources of a preparation or measurement protocol. |

It is therefore not necessary *in principle* to store an entire universe to
predict a finite local experiment. The existing locality theorem already
proves such a reduction for a specified action and state restriction. To use
it for the physical capacity family, one must additionally derive or bound
the local state and coefficients of that family. A common large-
capacity limit, or controlled corrections to it, would make a much smaller
calculation sufficient. If those local data retain an uncontrolled capacity
dependence, locality alone does not remove it.

The [capacity calculation](CAPACITY.md) gives a finite-history composition
bound for these distinct errors. It states sufficient conditions for a
controlled local prediction; it does not assume those conditions hold for
the physical OPH family.

The [native FIFO theorem](NATIVE_CAPACITY.md) derives one exact source-family
relation. With the interpreter's `export_oldest` policy and all exported
records retained, changing FIFO capacity among its admitted values 1 through
8 changes only the partition of the chronological tape between the archive
and memory. The complete native instrument, including its adaptive choices
and remaining quantum state, is identical after that reversible partition
change. Its error is zero for every input and reference. The `halt` policy
does not have this property. Neither policy identifies FIFO capacity with
physical global capacity or selects the coefficients of a field Hamiltonian.

## What an uncertainty statement establishes

For a density matrix and Hermitian observables the existing finite theorem
gives

\[
\Delta A\,\Delta B\geq\tfrac12|\operatorname{Tr}\rho[A,B]|.
\]

The [compensated-probe construction](UNCERTAINTY.md) binds this statement to
actual interacting observables, including their blocked outcome. Physical
position–momentum uncertainty is a further identification: it requires the
translation generator and position readout to have the appropriate canonical
relation and common action normalization. In finite dimension,
`Tr([X,P])=0`, so nonzero `i hbar I` cannot be an exact commutator on the whole
carrier. Canonical limits must specify their states and domains or a
controlled sector restriction. This is compatible with finite approximations
to QFT; it rules out one incorrect proof of that compatibility.

The [dense-effect argument](DENSE_EFFECTS.md) concerns a different interface:
why a countable operational effect domain can be sufficient for the existing
finite Born theorem. It does not infer such a domain from a handful of
tomographic settings.

## Particle energies need species-specific definitions

For a stable isolated massive sector, a suitable joint energy–momentum
spectrum can establish a rest energy and its dispersion, with nonzero
observable overlap. A finite transition difference relative to a declared
ground is only the starting quantity. Its clock and reference must be the
same as those used in the comparison.

An isolated one-particle pole is not a universal test for every Standard
Model label. Charged sectors with long-range massless fields require an
infrared-aware definition: under the hypotheses of
[Buchholz's Gauss-law result](https://www.sciencedirect.com/science/article/pii/037026938691110X),
charged states cannot be mass eigenstates. Unstable particles require
resonance/decay observables. Confined quark mass parameters require a specified
renormalization scheme and scale and are inferred through observable hadron
physics, as explained in the
[PDG quark-mass review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-quark-masses.pdf).
Massless propagation, bound states, thresholds and continuum backgrounds
must also survive the regulator and observation limits.

The final comparison is a prediction for complete experimental records,
including preparation and detector limitations. Simulator specifications
may help calculate that prediction. They do not replace the demonstration
that its quantities are those measured in the experiment.
