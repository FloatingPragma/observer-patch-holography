# Energy and work on the same spatial CAR law

This package supplies local energy operators, currents, drive power and
complete-record energy accounting on the **same** twelve-port,
1,080-mode ring and channel-correct Higgs action. It does not borrow an
energy law from the earlier unrelated full-Fock geometry. The 36 closed
internal components keep the largest Nambu calculation at dimension144;
all 2,160 Nambu directions and all CAR occupation sectors are retained.

```sh
python3 -m model.defect_field_energy.packet --check
python3 -m pytest -q -W error model/defect_field_energy
```

## Exact quadratic balance

For the centered CAR quadratic `F(B)=Psi† B Psi/2`, the CAR relations give
`[F(A),F(B)]=F([A,B])`. In normal order its identity term is
`c−Tr(h)/2`. That term is kept when converting commutators and work
operators; a scalar is not silently discarded because a Gaussian map
cannot observe it.

Local energy contains its complete onsite pairing and half each incident
kinetic bond. Thus `H=sum_x h_x` and

\[
J_{y\to x}=i[h_y,h_x]=-J_{x\to y},\qquad
\frac{d h_x}{dt_{\rm eff}}=
  \sum_yJ_{y\to x}+\dot\phi\,\partial_\phi h_x.
\]

The full ring has currents through range two under this allocation.
More distant local energies have disjoint CAR support and commute. All
internal currents cancel in the total energy derivative. This is an
operator identity on every occupation sector, including the collision
states and all massless modes. Independent three-mode commutator and
three-site full-Jordan–Wigner checks include the identity terms.

The default `component_energy` reference is `bare_propagation`, whose
scalar zero matches the literal propagation source factors. The optional
`common_static` convention subtracts eleven copies of the single-cell
vacuum energy from twelve replicated cells. The packet executes the
additional overall scalar source phase explicitly. Changing that fixed
reference changes neither currents nor energy differences; it is not
present in a Majorana matrix unless separately recorded.

## What the executed controls conserve

The source propagation law is an ordered product, not the exponential of
the simultaneous Hamiltonian sum. The package keeps that distinction.
For its selected four-mode factor, a declared two-segment interpolation
uses half the proper-time interval for each factor and twice its
Hamiltonian, reproducing exactly the existing ordered source word.
The common clock advances once for the entire interval. Both generators
use the same `1+j/4` multiplier through proper time.

Each constant segment conserves its own generator. Switching from the
reference Hamiltonian to the first control, to the second control, and
back to the reference supplies explicit work. Their sum accounts for the
reference-energy change, with the numerical segment residuals retained.
The following complete repair instrument can transfer additional energy.
All coin/raw energy contributions are retained before summing them; this
transfer is not labeled heat or controller work without a further
controller model.

On a source-prepared neutral test state, the coherent C policy transfers
`0.25 GeV`, while coherent E and both flagged controls transfer
`0.5 GeV`. These are consequences of the specified hopping coefficient
and preparation, not measured particle rest energies. Repair is therefore
not assumed to conserve a Hamiltonian merely because it preserves
particle labels or total charge.

## A neutral source readout

The two-site oracle used by the joint and modular packages has a
nontrivial Hamiltonian but **zero mutual energy current** under this
symmetric allocation. The packet instead reads the nonzero full-ring
operator `J_0to2` on its exact four-mode endpoint CAR support. This is an
operator-subalgebra identification, not a projection onto chosen field
states. The associated global CAR embedding uses the established
fermionic routing. The explicit source witness prepares the intermediate
and other omitted modes empty, with an opposite-charge spectator retained.

Separate Pauli-term measurements are inadequate as neutral energy
detectors: individual `XX` or `YY` pieces can change charge even when
their sum preserves it. A negative-control test keeps this issue visible.
The physical readout instead groups the exact symbolic quadratic into
number, same-charge hopping and opposite-charge pair observables. Each
group is diagonalized by an actual source Gaussian rotation, copied into
two neutral pointer bits, and uncomputed with fermionic routing intact.
Every one of the 64 input columns and all four outcomes is retained.
Every resulting field Kraus operator commutes with the actual electric
charge. No numerical visibility threshold selects terms.

Energy, local energy, drive power and ring-current **expectations** are
reconstructed from these groups on independent preparations. This is not
a projective measurement of total H or a measurement of its full energy
distribution when the groups do not commute. All spectator/reference and
native-record columns remain in the output factors. No detector shots or
hardware-accuracy claim is made.

## Kinematics and remaining scope

The same supplied lattice symbol has
`E(k)^2=sin(k)^2/a^2+E0^2`. For positive rest energy, its analytic curvature
is `E''(0)=1/(a^2 E0)` in dimensionless k, or `1/E0` in `p=k/a`, with
`c=1`. Lattice corrections and both kinetic zeros remain. This is
consistency of the chosen kinetic and mass terms, not an independent
mass derivation or a continuum inferred from twelve samples.

These finite energy/work identities and neutral readouts supply a
same-carrier input for further thermodynamic and gravitational tests.
They do not supply a3+1-dimensional stress tensor, a geometric null
modular kernel, a continuum limit, or Einstein dynamics.
