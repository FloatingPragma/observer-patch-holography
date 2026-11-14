# Regional charge instrument on the joint record construction

The apparatus addresses three operational prerequisites for row 1 on the
actual finite `JointRecords` and selected `InteractingRecords` channels. It
provides no physical species assignment or electric-charge unit calibration.
It is a measurement apparatus on the bounded integer instrument, not a change to its dynamics.

A region is a declared subset of vertices of the subdivided port/seam graph.
Occupied transit seams are assigned as whole vertices to the region or its
complement. The resulting sharp regional load operator has eigenvalues 0 and
1. It is distinct from the half-endpoint effective charge in the record-energy
transit extension. Assigning a whole transit unit to a region is a supplied
observer geometry, without a continuum-region interpretation.

`RegionalCharge.instrument` implements the complete Lüders instrument
`rho -> P_q rho P_q` for q=0,1. It returns unnormalized and conditional branch
states, branch probabilities and the nonselective state. Optional refinement
also resolves port versus transit occupancy, retaining all four outcomes,
including impossible outcomes with zero probability. The complementary regional
operator is I-Q. A charge measurement can preserve coherence within a region;
the port/transit refinement can disturb additional coherence. Neither is
silently identified with the existing full-source native monitor.

Let C be the incidence of each half-seam oriented from its port to its seam
midpoint, n the one-hot occupancy of these subdivided vertices, and F the
actual integer half-seam record. The exact Gauss operator is `G=n+C F`.
For a region indicator chi, all terms internal to the region cancel:

```text
Q_R + boundary_R(F) = G_R;
Q_R = chi^T n, boundary_R = chi^T C F, G_R = chi^T G.
```

Charge, boundary and Gauss are jointly measurable diagonal observables.
`readout` returns their complete joint probability census. The identity holds
in every resolved outcome. For every retained native event, the code checks
`delta Q_R = -delta boundary_R`. The boundary change is computed independently
from the ledger's applied current. Rejected attempts have zero applied current
and unchanged charge/boundary even when their proposed current is nonzero.
Transit noops and zero-current events are retained. For coherent evolution,
`i Tr(rho[H,Q_R]) + i Tr(rho[H,boundary_R]) = 0` is checked on the same selected
Hamiltonian. No classical current is assigned to an unmeasured coherent history.

The fragmentation control prepares an equal joint-state superposition across
region and complement in the same one-load Gauss sector. It measures mean 1/2
and variance 1/4; its sharp regional outcomes are 0 and 1 with equal
probability. Charge monitoring removes cross-region coherence. This distinguishes
a fractional expectation from a sharp fractional load. Absence of a half-unit
eigenvalue is a property of the supplied one-load sector. It does not rule out
physical half-charge species, multi-load processes, other sectors or energetic
fragmentation channels. Those require a physical charge/species dictionary and
the corresponding dynamical construction.

The executed cases use the oriented triangle `((0,1),(1,2),(2,0))`, integer
cutoff 1, preparation port 0 with zero records, propagation mass 0.23, time
interval 0.37, commit probability 0.29, and region containing port 0 and
transit seam 0. The variants are the free channel and the selected first-order
and massive record-energy channels with coupling 0.29 and energy mass 0.23.
Tests also use the existing one-edge control. No cutoff convergence, clock
limit, many-load generator or reserved Fable configuration is studied.

Run from any directory:

```sh
oph-openwave regional_charge --output /tmp/regional-charge.json
python -m pytest -q candidate/tests/test_regional_charge.py
```

The apparatus and preparation are explicit supplied choices. A physical row-1
claim additionally needs the electric-charge dictionary and full admitted
sector/intervention contract. No official cell or full criterion follows.
