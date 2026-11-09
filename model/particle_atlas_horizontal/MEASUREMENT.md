# Preparation, measurement and shared records

The four-bin experiment is an exact illustration of a supplied quantum
instrument. Its bins, state and Fourier analyser are specified inputs. It
does not identify a native carrier coordinate with physical particle
position, or provide one position projector for every particle family.

## One preparation, two apparatus settings

A preparation is a reproducible procedure that establishes the state used
to predict the next records. In this illustration it prepares
`psi=(1,1,0,0)/sqrt(2)`. The two occupied bins have the same phase. This
describes the amplitudes for later records; it is no claim of two smaller
classical particles.

An apparatus couples a chosen observable to a pointer. The position setting
routes each spatial bin to its corresponding detector. The momentum setting
first applies the unitary four-bin Fourier transform, then uses the same
four detector ports. Momentum is represented by a spatial phase pattern,
so amplitudes from different bins can reinforce or cancel at an output.
Changing apparatus changes which question is asked of the same preparation.

A detector correlates one output with a distinguishable record. For a
complete instrument `{K_r}`, every outcome, failure or zero-weight label
belongs to the experiment. The corresponding effect is `E_r=K_r† K_r` and
`sum_r E_r=I`. The Born weights are `p(r)=Tr(rho E_r)`. The specified
projective instrument gives the conditional state
`K_r rho K_r†/p(r)` only when `p(r)>0`. Zero-weight records are retained
with no normalized conditional state. Probabilities alone do not determine
an instrument's state update.

Consensus compares, transmits and verifies the detector record. Distinct
classical record labels can be copied; their probabilities are unchanged by
sharing. Agreement does not freely choose an outcome or clone an arbitrary
unknown quantum state. A stable record fixes the outcome used in conditional
predictions. It does not create the particle's mass, spin or conserved charge.

## Exact graphical example

For the same prepared state, the complete distributions are

| Experiment | Four detector probabilities |
|---|---|
| Position apparatus | `1/2, 1/2, 0, 0` |
| Momentum apparatus | `1/2, 1/4, 0, 1/4` |
| After the position record `x0`, momentum apparatus | `1/4, 1/4, 1/4, 1/4` |
| After the momentum record `k0`, position apparatus | `1/4, 1/4, 1/4, 1/4` |

Each highlighted first record has probability `1/2`; it is one branch of a
complete experiment. The code retains all sixteen sequential histories for
each order. Reversing the phase of the second occupied bin leaves the
position distribution unchanged but changes momentum probabilities to
`0,1/4,1/2,1/4`. This is interference in a directly checked finite example.

A sharp position record broadens this state's momentum distribution.
Momentum has a well-defined probability law; the state has no single sharp
momentum value. A selected momentum mode similarly spreads the position
distribution. This concerns preparation spread and measurement backaction,
not a universal claim that every position measurement completely erases all
momentum information. Approximate joint measurements can be unsharp in both.

The code verifies the finite Robertson inequality for its actual matrices.
Finite matrices have `Tr[X,P]=0`, so they cannot satisfy
`[X,P]=i hbar I` with nonzero `hbar`. The continuum `hbar/2` bound requires
the appropriate common operators, domains and states.

## A double-path experiment with an actual readout

The same four-bin model gives a small, exact double-slit-style experiment.
The preparation supplies equal coherent alternatives at two apertures,
`psi=(|0>+exp(i phi)|1>)/sqrt(2)`. The supplied Fourier apparatus carries
those amplitudes to four screen detectors D0 through D3. Each run records
one screen port. Repeated preparations give a histogram of the four ports.
The diagram is an illustration with a specified apparatus, rather than a
derivation of a laboratory double slit from native repair.

Without a path marker, amplitudes from the two apertures arrive in the same
output channel and add before their squared magnitude is taken. At phase
zero, they reinforce at D0 and cancel at D2. The complete screen law is
`(1/2,1/4,0,1/4)`. D2 is a retained zero-weight detector outcome. Changing
the preparation phase to pi interchanges the bright and dark ports,
giving `(0,1/4,1/2,1/4)`.

A path marker is an additional physical coupling. It correlates aperture0
with one distinguishable marker state and aperture1 with another. The
joint state keeps the two path/marker alternatives, while the unread screen
law sums over the marker. With orthogonal markers the alternatives cannot
interfere in that screen marginal: its distribution is uniform. A person
does not have to inspect the marker. The retained path correlation is
sufficient. The complete joint law retains all eight records `(path,Dk)`,
each with probability `1/8`. Conditioning on either path, whose probability
is `1/2`, gives a uniform screen distribution. No path or detector record
is selected away.

For partial path marking, use normalized marker states d0 and d1 and define
`gamma=<d0|d1>`, with `|gamma|<=1`. The exact specified output state at port k
is proportional to `d0+exp(i phi)(-i)^k d1`. Thus

```math
p_k(\phi,\gamma)={1+\operatorname{Re}[\gamma e^{i\phi}(-i)^k]\over4}.
```

The full marker coupling is an isometry on all four input paths, including
the two paths with zero preparation weight. Taking its partial trace gives
the same four probabilities. Scanning phi at a fixed screen port gives
`p_max=(1+|gamma|)/4`, `p_min=(1-|gamma|)/4`, and visibility
`(p_max-p_min)/(p_max+p_min)=|gamma|`. This visibility is defined by a phase
scan; four sampled screen bins at one arbitrary phase need not attain both
extrema. At phase zero and gamma=1/2 the screen law is
`(3/8,1/4,1/8,1/4)`. Complex gamma shifts the fringe phase, with the
convention tested explicitly in the code.

The example separates a particle-like record from wave-like propagation:
each detector event has one classical label, while amplitudes for
indistinguishable alternatives interfere in the probabilities of those
labels. A path record changes the joint state and its reduced interference
law. Consensus shares the detector outcome after that physical process; it
does not choose the outcome, supply an amplitude, or require consciousness.

## What the final ledger means

The record starts with a known preparation procedure and apparatus setting.
The screen detector supplies the observed label. These fields define what
was measured and make repeated runs comparable:

| Field | Meaning |
| --- | --- |
| Run | A label distinguishing this preparation-and-detection attempt |
| Preparation phase | The supplied relative aperture phase phi |
| Apparatus setting | Unmarked paths, orthogonal path marker, or a specified partial marker |
| Path record | Aperture0 or aperture1 only when the apparatus supplies that record |
| Screen record | One of D0, D1, D2 or D3 |

For example, `(run42, phi=0, unmarked, no path record, D0)` is one possible
record, with probability `1/2` under that complete experiment. Its shared
copies have the same fields. It does not establish a sharp aperture value
that the apparatus did not record. In the marked experiment,
`(run43, phi=0, orthogonal marker, path1, D0)` is one possible joint record,
with probability `1/8`. These are illustrative ledger entries, not selected
simulation results. The certificate supplies every endpoint probability.

The screen label names an actual detector port within this specified model.
Calling it a physical momentum or position measurement requires the
apparatus-to-field dictionary. The Fourier-mode interpretation used here
does not supply a universal position detector for every particle family.

## The source theorem connections

The [native reconstructed functional](../native_quantum_reconstruction/PROOF.md#native-count-moments-and-their-positive-density-representation)
is derived from complete native transfer and closed-history statistics:
`rho_L=P^L/Tr(P^L)`, with `2I/3<=P<=I`. Its full complex matrix functional
has Born form. The counts used to reconstruct a general effect can include
scaled or signed combinations; they are not all individual native detector
acceptance events. Direct native definite-record effects are diagonal.

The finite effect theorem
[`finite_busch_gleason`](../../Lean/EventAlgebra/FiniteBuschGleason.lean#L1133)
derives a unique density representation from nonnegative normalized weights
additive on all coexistent matrix effects. The operational theorem
[`operational_record_probabilities_unique_born`](../../Lean/EventAlgebra/OperationalBornBridge.lean#L214)
connects record probabilities to this result when equal effects have equal
probabilities and every coexistent pair has disjoint implementations in one
experiment. These full complex operational premises are separate from
classical record consensus. The [flagship quantum section](../../flagship/from_observer_consensus_to_standard_physics.pdf#page=14)
states the finite public-record identities and this source boundary.

The [same-source position proof](../spin_covariant_localization/CONTINUUM_POSITION.md#exact-covariance-and-the-uncertainty-relation-on-this-subspace)
constructs `X_i=V i partial_i V*` and `P_j=V hbar k_j V*` on
`V Schwartz(R^3)` for the specified ordered-spin source and selected state.
Here `V` comes from that source's field Gram form, and the same spatial
translations generate momentum. Hence
`[X_i,P_j]=i hbar delta_ij I` and `Delta X_i Delta P_i>=hbar/2` on the
stated domain. The full detector retains a no-position complement; the
vacuum is in that complement. The time-invariant
[dynamical response-space extension](../spin_dynamical_localization/PROOF.md)
supplies a covariant position measure with explicit moment assumptions.
These spaces need not be isolated particle sectors.

The [finite source compiler](../spin_operational_position/PROOF.md#3-every-computable-finite-instrument-has-an-explicit-source-word)
uses admitted coherent load, pointer and refinable controls. Those controls
are inputs; the construction does not derive them from bare occupation
repair. The [carrier displacement identity](../carrier_position_translation/PROOF.md)
is a distinct coordinate readout. The rank-three carrier frame does not by
itself locate a field excitation. The
[finite Robertson theorem](../../Lean/EventAlgebra/Robertson.lean#L295)
applies to actual finite observables with their state-dependent commutator.
The particle paper's [probe and position discussion](../../paper/deriving_the_particle_zoo_from_observer_consistency.pdf#page=59)
keeps these scopes separate.

Photons require an appropriate detector position observable; a sharp massive
particle position projector is not supplied by the word "particle".
Quarks and gluons are compared through hadronic and jet observables, unstable
species through resonance and decay records, and neutrino flavor through
production, oscillation and detection. Their readout maps require the same
interacting state, physical clock and calibration as their field action.
