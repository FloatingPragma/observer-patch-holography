# OPH and the executable-model hypothesis

The executable model of OPH runs the OPH law, the simple local rule by which
neighbouring observers settle their disagreements. It moves each constituent
by reflected routing, binds two constituents that meet into a pair through a
motion-coupled contact law, forms that pair through a capture vertex that
sends away a marker, and opens the hop of a newborn pair over 128 updates.
The [specification](../simulator/SPEC.md) states the law, the
[simulator](../simulator/README.md) executes it, and the
[model page](../model/README.md) lists what is proven and measured on it.
The hypothesis is that the agreed records of its observers reproduce the
observable physics of our Universe.

The theory around the model is finished. On the OPH law the field of
one constituent has a proved 3+1D Klein–Gordon limit, moving clocks keep
proper time at long wavelengths, a bound pair forms in an encounter, and
formed matter obeys the long-wave relation omega^2 = m^2 + c^2 |k|^2. The
other observables of the theory on this model are
work in progress.

The outputs of the model meet some premises of finite theorems of the theory
as they stand. Its field responds in three independent spatial directions.
A presence detector at one address is the projective instrument of that
event. A laboratory's clock counts the updates of its own address, and
observers relate their counts through shared records. The gauge, gravity and
particle premises that other theorems consume are work in progress on this
law.

## Closure and a specified construction

The overall OPH hypothesis is closure: a complete history is both a timeless
self-referential structure and its internal computation. Observers recover
and realize their world's specification inside that history, including its
realizing hardware. Space, time and matter describe compatible observations.
The [basis](BASIS.md) and [interpretation contract](INTERPRETATION_AND_TERMS.md)
state the proposal and its adopted mathematical premises.

The creator may specify the carrier, the law, preparations, references,
instruments and finite resources. Simulated observers need not invent or
install the law. Neither uniqueness nor derivation of architectural choices
is a prerequisite for execution. The resulting records must satisfy their
physical comparisons with all calibration inputs visible.

## The model

A carrier is a torus with three spatial axes and a private cyclic register at
every address. One constituent lives on eight slots per address and private
position. Each update reflects the equal-slot part of the eight amplitudes
and sends every slot one address along its axis into the partner slot. Two
constituents that share an address in matched slots receive a contact coin
and step together as a bound pair; one constituent carries two rails whose
contact phases differ, and their relative phase is the pair's clock hand.
Where the two constituents meet head-on, a reversible vertex turns the
collision into a matched pair plus a marker that leaves on the address's own
wire. A newborn pair carries an age counter of 8192 states and opens its hop
by a reversible partial swap on the profile (1 - cos(pi u))/2 over 128
updates, then hops fully.

Exact statements about this construction, each checked by the simulator
against the literal update: the field of one constituent released in the
equal-slot state obeys the free Klein–Gordon equation with signal speed 1/2
under refinement, with error of order h^2 at fixed chart time; the matched
coincident states of the contact law form an invariant subspace, so a bound
pair never separates under the contact law; the mature pair's frequency obeys
cos(omega) = cos(m) (2 w_private + 2 w_spatial sum_i cos k_i), so
omega^2 = m^2 + 0.24919629 |k|^2 + O(|k|^4).

### Supplied

- the carrier and its sizes, and the private sector;
- the preparations and the detectors;
- the scalar reference standard of the unbound walk.

## How the model is tested

The state of the simulation does not contain the physical Universe or its
spacetime. Addresses, slots and wires are hardware of the computation and no
physical space; no test asks the carrier to be a three-dimensional lattice. Physical space, time, fields and particles are
effective descriptions that the observers inside the computation use for the
records they agree on. Every test is therefore a comparison between
distributions of admissible observer records and what observers in our world
record. A difference between microscopic state encodings is no failure of
the model, and a resemblance between them is no success.

The bridge is proved in the theory. The flagship's proposition on smooth
effective field descriptions shows that for every finite admissible menu of
observer readouts a smooth (1+3)-dimensional Lorentzian field theory
reproduces their complete predictions, and the spacetime paper gives the
rank-three spatial carrier and the causal and count limits. Observers who
inhabit the simulation can describe their experiences as relativistic 3+1D
continuum spacetimes. The executable model supplies part of the premises of
those theorems with its own local clocks and records: a field with a proved
3+1D continuum limit, complete instruments and the count clocks of its
observers.
The readouts in field form, the order of events with its placement and the
refinement family that the theorems consume in addition are work in
progress.

## Work in progress

| Observable | What the model has to supply |
| --- | --- |
| Clocks and energy | A portable clock of formed matter, whose rate one stationary component keeps at every incoming width, and a clock that starts with no supplied phase. |
| Signals and frames | One causal speed for every operation inhabitants can trigger, and a travelling laboratory that ranges by echoes and measures a persistent rod in its own frame. |
| Particles | Formation in an ordinary history, the class of an object as a record, charge, and the identification of formed matter with particles. |
| Interacting fields | A record of currents that later events read, and the gauge forces. |
| Gravity | A refinement of the carrier, and the regional response that the Einstein theorem consumes. |
| Independent comparison | Complete joint observations and interventions, input ancestry, uncertainties and rejection criteria compared with independent measurements. |

These requirements are conjunctive. An exact identity, a finite transition,
a fitted parameter or a conditional tensor implication cannot replace the
other hypotheses. A source-specific proof on another construction keeps its
own domain until an explicit map connects it.

## Gravity through the matched-regional theorem

The [matched-regional theorem](../model/matched_einstein_limit/README.md)
consumes ten entropy/area/stress observations in one orthonormal frame. Its
finite residual bound propagates three stated errors, with common positive
area coefficient and radius. The same regional data must enter all three
estimates. A continuum attachment requires errors negligible relative to the
fourth power of the shrinking proper radius, convergent tensor data and a
common nondegenerate Lorentz metric. Independent conservation and a common
vacuum reference give the corresponding constant-coupling and
cosmological-reference conclusions.

## Scope of physical quantities

Construction units, laboratory seconds and measured couplings are distinct
attachments. An input value named speed, field coefficient, reference energy
or vacuum strength does not derive physical `c`, the fine-structure constant,
Newton's constant or the cosmological constant. A finite matrix element is not
an on-shell decay rate. Frequencies of the model's own formed matter are
stated in radians per home tick of an observer's clock.

The intended photon-containing physical theory has no strictly positive
global vacuum gap. Sector gaps and pure Yang-Mills statements keep their
separate hypotheses and limits.
