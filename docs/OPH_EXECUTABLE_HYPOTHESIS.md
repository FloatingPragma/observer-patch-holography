# OPH and the executable-model hypothesis

The official executable model is
[polynomial-sv](../model/polynomial_sv/README.md). Its law, the contrast-memory
law, is the candidate native repair law of the theory: the simple local rule
by which neighbouring observers settle their disagreements. The hypothesis is
that the agreed records of its observers reproduce the observable physics of
our Universe. The [specification](../model/polynomial_sv/SPEC.md), the
[binding theorem](../model/polynomial_sv/PROOF.md) and the
[executable-model paper](../paper/executable_observer_construction.pdf)
state the law and what its observers read.

The theory around the model is finished. The model reproduces particle
stability: bound objects persist in free evolution, read on the observers'
agreed records. Its coupling is supplied. The other observables of the
theory on this model are work in progress.

## Closure and a specified construction

The overall OPH hypothesis is closure: a complete history is both a timeless
self-referential structure and its internal computation. Observers recover
and realize their world's specification inside that history, including its
realizing hardware. Space, time and matter describe compatible observations.
The [basis](BASIS.md) and [interpretation contract](INTERPRETATION_AND_TERMS.md)
state the proposal and its adopted mathematical premises.

The creator may specify the carrier, the law, preparations, references,
instruments and finite resources. Simulated observers need not invent or
install the laws. Neither uniqueness nor derivation of architectural choices
is a prerequisite for execution. The resulting records must satisfy their
physical comparisons with all calibration inputs visible.

## The model

The carrier is a federation of complete twelve-port icosahedral cells (the
theory's Lean library calls the cell the echosahedron) joined by external
seams that form a matching on ports. Every port holds a load of 0, 1 or 2.
Every seam keeps a private memory of eleven states, coherent working memory
of its two port observers and no public record. The memory writes the
contrast of the two loads, gives it a phase and unwrites it; on the loads this
echo acts exactly as

    H_kappa = (1/2) sum_e (I - SWAP_e) + kappa (6 sum_v n_v^2 - 2 sum_e n_a n_b),

an exchange of loads plus an agreement energy with coupling `kappa`. A local
emitter writes a pair of loads on one seam. A reader on every seam admits a
record only when the copies of its two port observers are valid and agree.

On every federation of complete cells, for `kappa > 3`, a pair of loads
written on one seam is read as a pair on some seam with probability at least
`1 - 40/(2 kappa - 6)^2` at every reading, independent of the number of
cells. The pair is read on seams other than the first, nothing holds it but
the law, and the same law without the agreement energy loses it. Run as local
events with no common clock, the pair holds when every port observer serves
each of its seams once per round, in any order.

### Supplied

- the coupling `kappa`;
- the coherent generator on the private seam memory;
- the reference registers of the degree term on ports with five seams;
- the emitter with its fuel, coupling and duration;
- the evolution parameter of the runs with a continuous parameter;
- the routing of the external seams of the federation.

## How the model is tested

The state of the simulation does not contain the physical Universe or its
spacetime. Ports, seams, cells and federations are hardware of the
computation and no physical space; no test asks the carrier to be a
three-dimensional lattice. Physical space, time, fields and particles are
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
continuum spacetimes. What the executable model supplies are the premises of
those theorems, with its own local clocks and records.

## Work in progress

| Observable | What the model has to supply |
| --- | --- |
| Clocks and energy | An observer clock normalized to proper time, from each port observer's own count of events. |
| Interacting fields | A record of currents that later events read, and the twelve port currents from the law's own services on the bound object. |
| Particles | A selection of `kappa`, a second sign of load, the class of an object as a record, and one limiting speed for all objects. |
| Gravity | A refinement map that preserves seams, and the regional response that the Einstein theorem consumes. |
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
an on-shell decay rate. Energies of the model's own bound objects are stated
in units of the exchange.

The intended photon-containing physical theory has no strictly positive
global vacuum gap. Sector gaps and pure Yang-Mills statements keep their
separate hypotheses and limits.
