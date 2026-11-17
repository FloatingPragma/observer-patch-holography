# The common matter, field and response construction

## What is being built

The official executable realization is [CCF-LT-2](SPEC.md). Its creator
specifies a finite oriented carrier, the full native matter algebra, one
Hamiltonian, initial states, idle references, and complete observer
instruments before execution. The simulated inhabitants are not required
to design or install their own laws. The scientific question is which
compatible records this construction produces and how they compare with
physical observations.

An instance consists of

```text
(oriented complex, native CAR algebra, rotor/reference carriers,
 source coefficients, local energy allocation, auxiliary configurations,
 initial joint state, complete instruments, observation horizon).
```

Each vertex retains 96 CAR modes. The source is the spatially replicated
114-factor matter snapshot plus signed charged/neutral Dirac transport,
complete leptonic and quark-current contacts, quark-family splitting,
electric squares and magnetic face terms. All normal-order constants,
neutral mixing, complementary occupations and calibrated coefficients
remain part of the law. The exact source and generator are in
[ccf_lt2](ccf_lt2/README.md).

## The joint law

Write the matter/field Hamiltonian as `H_S=sum_v h_v`. The exact energy
allocation shares words over their support, electric energy over endpoints,
magnetic energy over face vertices, and pure scalars uniformly. The law is

```text
g_v = local response energy + half of each incident gradient energy
k_v = g_v + N_v (h_v-e_v I) + e_v I
H_total = sum_v [A_v + M_v (k_v-F I) + F I].
```

All source and auxiliary terms are installed before free evolution
`exp(-i tau H_total)`. Every vertex has an independent three-state response
register and 16-state capacity/memory register. Both local lapses are strictly
positive. The capacity lapse acts on the complete vertex energy, including
its response-energy shard. The outer reference `F` is one common number.
Uniform source shifts with matched inner and outer reference shifts change
only a common phase. Nonuniform shifts cannot generally be absorbed into
that common outer reference. Disconnected components factorize exactly.

For `G=BE-n-r`, the exact native charge is `n=3(T3+Y)`. Dressed matter
words, electric squares, closed face loops and both auxiliary layers
preserve every Gauss operator. Hard-wall shifts retain their actual ordered
matrix action, including the wall. Charge-flow identities never assume a
finite raising matrix has an inverse.

## Preparation, repair and readout

Preparation is a declared operation on the complete source/reference state.
[Charge-preserving Gauss preparation](ccf_lt2/FIELD_PREPARATION.md) acts on
flux consistency and retains its
syndrome history, capacity and work. It is separate from free evolution.
No repair, reset or conditional successful-arrival replacement is inserted
in that interval. A Gauss-compatible initial state remains compatible by
the same Hamiltonian.

A complete instrument specifies every Kraus operator, including failed,
no-arrival, complement and capacity outcomes. Retain reference and memory
registers through the stated protocol. Public copies are copies of orthogonal
record labels. They do not clone arbitrary private states. Premeasurement
operator diagnostics and detector outcomes have distinct labels. Terminal
meter work is the change in mean energy induced by the nonselective
instrument; total energy is conserved during the free interval.

## Completion target and acceptance criteria

The finite acceptance calculation must regenerate the complete native
polynomial, verify its exact source recovery and local Gauss identities,
check local energy allocation, and apply the complete local operator to full native
basis columns. The finite matrix diagnostics must retain complete records,
measurement work and the declared negative controls. The
[dependency map](particle_joint_reconstruction/DEPENDENCIES.md) binds these
checks and the associated proofs.

A small matrix is an execution of the full native law only when an exact
reducing-subspace proof connects it to that law. Otherwise it is a separately
declared finite diagnostic. Full native column action evaluates every term
on the specified occupation and flux labels; it is not arbitrary many-body
state propagation. The interacting quark and local-energy terms do not
permit an unproved colourful-vacuum or neutral-spectator reduction.

The physical target is compatible complete observation laws with common
units, states, clocks, geometry and calibrated uncertainties. It includes
interacting spectra, outgoing decay channels, electromagnetic and matter
propagation on one metric, and the actual prerequisites of the retained
Einstein theorem. Algebraic compatibility and a finite response establish
only their stated results. Observed masses, couplings, physical time and
continuum identification require source-specific evidence.

## Reusable mathematical foundations

Finite Hermitian dynamics, positive state valuations, complete instruments,
CAR identities and conditional tensor reconstruction retain their stated
mathematical validity. Their hypotheses must be instantiated with this
source and its actual observations. Source-specific spectra, ground states,
geometry, clock frequencies and refinement estimates are recalculated or
connected by explicit maps. An independent component's successful replay
is not itself such a map.

[Premise bindings](PREMISES.md), [physical scope](PHYSICAL_AUDIT.md) and the
[source contract](SOURCE_SELECTION.md) give the corresponding domains.
Architectural choice is legitimate; uniqueness is a separate mathematical
claim, and empirical calibration remains visible as input ancestry.
