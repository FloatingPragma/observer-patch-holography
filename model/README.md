# The executable OPH model

The official model is **CCF-LT-2** (`A-ccf2-local-capacity-v2`), the common finite matter, field and
response construction. A specified quantum source carries native matter,
charged and neutral transport, electromagnetic links and faces, weak
currents, and a response to total energy. Its instruments produce finite
observer records, with references and all outcomes retained.

This is an early calibrated test model. It is the main subject of current research and is under heavy
development.

The [specification](SPEC.md) and [design](DESIGN.md) fix the full law, source data, parameters,
preparations and execution scope. The [executable package](ccf_lt2/README.md)
contains its generator and calculations. The
[construction contract](CONSTRUCTION.md) explains what counts as a result
of the common model.

## Run the model

From the repository root:

```sh
python -m pip install -r model/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model
python -m pytest -q model/ccf_lt2
```

The default command regenerates and verifies the complete operator checks.
A shorter component replay is available with `python -m model --headline`.
The complete command checks both the two-vertex and three-vertex native
sources and their local energy response. It does not allocate the
exponentially large full native state vector.

## Matter and electromagnetic fields

Every vertex carries all 96 native fermion modes. Signed Dirac transport
retains particles, antiparticles, neutral modes and pair terms. The actual
calibrated matter coefficients are preserved. Electric flux, magnetic loops
and matter charge obey one oriented Gauss constraint; the current comes
from the same action. Complete finite checks establish these algebraic
relations and allowed channels. Physical particle spectra and a common
relativistic limit require their own calculations.

## Response to energy

Every vertex has independent response and capacity/memory registers,
coupled to its complete energy shard, including interaction and field terms. The full time-independent
Hamiltonian conserves total energy. Separate finite detector experiments
show a response to neutral energy; literal native basis-column calculations
verify that every source term participates in the combined operator.
An Einstein or Newtonian interpretation requires the source and observation
maps in the [proof dependency map](particle_joint_reconstruction/DEPENDENCIES.md).

## Proofs and physical comparisons

The [physical scope](PHYSICAL_AUDIT.md) separates complete operator checks,
finite diagnostic evolution and physical claims. The
[source contract](SOURCE_SELECTION.md) fixes input ancestry and the
[observation contract](OBSERVER_RELATIVE_RECOVERY.md) fixes complete records.
Other packages under `model/` contain mathematical components with their
own Hamiltonians, states and domains. Their theorems apply to this model
only through checked hypotheses; their numerical outputs are not silently
used as its predictions.
