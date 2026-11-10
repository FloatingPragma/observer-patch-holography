# Coupled ground and positive macrocell completion

The unchanged finite native Wilson–Maxwell source has its active physical
ground in the retained native complement, with multiplicity96 before
spectator references. Exact Coulomb/Peierls contractions and a weighted
Schur bound place its entire interacting code above that ground.

The separately declared positive completion
`H_new=H_old+m c(g) sum_x(I-P_x)` selects a unique active coupled code
ground. Each `P_x` is the complete code projector of one four-constituent
macrocell. The rational coefficient is computed from the same source's
vacuum contractions; no measured parameter is used. The old source and
receipt are unchanged. The completion preserves absolute code energies
and all original pattern-preserving complete count laws. It changes
complement energies and coherent phases between patterns.

The full physical gap is greater than `2.52078609e-10` source units,
distinct from the internal code gap greater than `1.7170018m`. A complete
finite energy-agreement preparation uses four source trials and three
retained pointers of dimension `157589952376326706`. Its all-pass
probability exceeds `.9995943286660` and the conditional main-source
ground trace-distance error is below `1.417923e-7`.

An actual stationary-ground pair/antiparticle/photon count has positive
weight, with a preparation-inclusive joint lower bound above `1.41692e-6`.
This is virtual occupation in the declared free-field frame. It is not
ground radiation, asymptotic scattering, an observed particle mass or
an interacting continuum theorem.

From the repository root:

```sh
python -B -m model.native_wilson_maxwell_ground --write
python -B -m model.native_wilson_maxwell_ground --check
python -B -Werror -m pytest -q model/native_wilson_maxwell_ground
```

`PROOF.md` contains the full physical complement comparison, exact
projector and Gauss intertwining, complete failure records, raw energy
and common-clock convention, controller work bounds and capacity.
The source coefficient and its positive completion are declared
architecture choices; primitive repair is not proved to select them.
The kinetic/current source retains its preexisting global code guards.

`research_reference.json` is an immutable local scientific comparison;
the producer recomputes every scientific field and checks exact equality.
The public package has no research-workspace or external absolute-path
dependency. `PROVENANCE.md` states source and import closure.
