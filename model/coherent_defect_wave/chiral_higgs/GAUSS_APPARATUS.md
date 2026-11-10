# Compensated charged energy apparatus

The [complete quantum energy instrument](PROOF.md#4-complete-charged-energy-records-on-the-full-quantum-reference)
uses an e^c reservoir. That reservoir carries hypercharge one and is a weak
and color singlet. Its initially occupied setting has a nonzero gauge charge
unless its source boundary or apparatus contains the compensating charge.
A naked charged source creator does not preserve the source-only Gauss
kernel. This construction supplies the missing apparatus charge explicitly.

## Physical initialization of the complete ground instrument

The source energy preparation in the owning proof starts with `I320/320`.
This source marginal has charged sectors. It is physically admissible with
an explicitly retained dual-representation apparatus, rather than being an
uncompensated source-only Gauss state. Use a second 320-dimensional register
with currents `-Q_source^T`, the conjugate source parity, and the entangled
state `Omega_D=sum_i |i>source |i>dual / sqrt(320)`. Its marginal is exactly
`I320/320`. In vectorized coefficient form the twelve total Gauss actions
are `Q_source I-I Q_source=0`; the code checks each actual source current.
Its total parity is even because each source basis state is paired with
its equal-parity dual state. The full purification dimension is 102400.

A declared gauge- and parity-preserving preparation exists: the empty source
and empty dual state and `Omega_D` are vectors in the same even total Gauss
kernel. Extend their normalized-vector map to a unitary within that kernel
and use identity outside it. This is a finite apparatus specification,
without an efficiency claim or native idle vacuum-selection claim.

Each source energy projector commutes with its gauge currents and parity.
The energy PVM tensored with the dual identity therefore preserves the total
constraints and retains all fourteen outcomes. Their probabilities are
exactly the original ranks divided by 320. The ground probability remains
`2/320=1/160`, and its conditional source marginal is the same rank-two
ground projector divided by two. Because that projector is annihilated by
every source current, both the source and the dual register are separately
gauge singlets in this branch. They are both odd in parity; their combined
parity stays even. The dual quantum state is retained, including on all
failure branches. It remains a spectator during the charged readout below.
Tracing that spectator for the computation of local source outcomes is
valid; it is not an assertion that the preparation discarded its record.

## Source and boundary registers

Retain the six lepton CAR modes, the five-state quantum scalar, one fermionic
reservoir r of charge plus one, and a compensator b of charge minus one.
The compensator is a static apparatus/boundary mode, not a field identified
with a measured particle. The active dimension is
`2^8*5=1280`, or 409600 with the retained dual preparation register.
Its two retained initial settings are

`(n_r,n_b)=(1,1)` or `(0,0)`, each with recorded probability one half.

The source reference is the normalized rank-two gauge-singlet ground of
the same 320-state lepton/scalar Hamiltonian, with empty quark spectators.
Both settings have total hypercharge zero. Every other local gauge generator
annihilates the source reference, and both apparatus modes are weak/color
singlets. The combined initial density is therefore supported in the
complete source-plus-apparatus Gauss kernel. The source reference is not
replaced by a different charge ensemble during the measurement.

## A full, charge-conserving fermionic exchange

Order the source modes first, then r and b. The unitary is the actual
fermionic permutation exchanging source mode zero with r and fixing b.
Every one of its 256 fermion basis columns, including doubly occupied and
compensator-occupied columns, is retained with the permutation sign. The
code checks all eight annihilator intertwiners. The identity is tensored
on the entire scalar memory.

For the hypercharge generator the conserved operator is

`G_total=G_source+n_r-n_b`.

For every other current it is the source current alone. The source mode
being exchanged is an e^c component, so its charge equals that of r and it
is a weak/color singlet. Actual full-Fock matrix commutators verify that the
permutation commutes with all twelve total currents. The scalar generators
commute because the operation acts as identity on their register. The
compensator number also commutes with the permutation in every sector.
Thus the combined evolution stays in its initial Gauss kernel.

The source Hamiltonian commutes with its twelve currents. Its complete
energy projectors therefore commute with the total currents when extended
by the identity on r and b. The reservoir and compensator occupation PVMs
also commute. Consequently each complete retained measurement effect is
compatible with the total Gauss constraint. The pointer can copy these
classical outcomes without copying an arbitrary quantum field state.

## Exact correspondence to the complete energy records

The code evolves every nonzero density-matrix entry of the declared source
reference and both apparatus preparations through this actual permutation.
It then reads the reservoir and compensator and evaluates the same fourteen
source energy effects. There are 56 records: two initial settings, two final
reservoir values and fourteen energy labels. The compensator outcome equals
the retained setting with probability one; its value is included as boundary
metadata rather than silently omitted. Impossible compensator changes have
zero effects on the admitted initialized apparatus domain. The full unitary
and all apparatus PVMs exist on the entire 1280-state space.

All 56 source energy probabilities agree exactly with the source instrument
in the owning proof, including zeros and failed addition/removal. These 56
records are conditional on successful ground preparation. The fourteen
initial energy-preparation outcomes form a separate complete stage. Retaining
its thirteen non-ground outcomes as failure records and multiplying each
subsequent readout weight by the original success probability `1/160` gives
one complete 69-record preparation-and-readout distribution. The code
exports those unconditional probabilities and checks that their sum is one.
The fermionic
permutation can induce source-parity signs in an ordinary tensor ordering;
the executed density equality uses the specified odd-parity ground. No
unqualified equality of cross-parity source CP maps is claimed. Classical
mixtures of other parity sectors require their own stated convention.

This supplies the charge-compensation part of that energy readout. The exchange
is a driven apparatus operation; its work supply is separate. It does not
turn the onsite source energy into a relativistic rest mass, select a
Higgs background, or establish an interacting continuum.
