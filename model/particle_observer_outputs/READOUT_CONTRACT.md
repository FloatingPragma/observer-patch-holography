# Complete observer records from one finite interacting source

A specified finite Hamiltonian, preparation and apparatus determine a complete
observer record law. The calculation does not require a theorem that this
architecture is unique among possible computers. It does require that each
reported probability use the same Hamiltonian, state, clock and instrument.
This contract connects the actual finite Higgs sector in
`higgs_records.py` to the existing complete-instrument compiler and propagates
its errors through an entire recorded experiment.

There are two admitted implementation routes. The first declares the finite
unitary preparation and pointer coupling as primitive source services. Their
complete matrices define the experiment. The second encodes those matrices
into the coherent spin-control interface and consumes the existing compiler
theorem. Both routes compute observer outputs. Neither route identifies a
finite detector label with a laboratory particle without a field and
experimental comparison map.

## One Hamiltonian and a complete spectral apparatus

Use the owner's scalar tags `0,1,2`, canonical fermion masks and conserved
sector `N_F-2N_H=0`. With three family coefficients its dimension is 307.
The fixture specifies `Y=diag(1,2,3)` and common pair cost `d=1`; these integer
inputs contain no particle reference values. The actual source interaction
gives `g²=4 Tr(Y†Y)=56` and gap `Delta=sqrt(d²+4g²)=15`.
Its complete spectral decomposition is

```math
H=-7\Pi_-+8\Pi_++D,\qquad
\operatorname{rank}\Pi_-=\operatorname{rank}\Pi_+=1,
\quad\operatorname{rank}D=305,
\quad\Pi_-+\Pi_++D=I_{307}.
```

`energy_menu` checks the projections, their orthogonality, every equation
`H Pi=E Pi`, and the full reconstruction of H. The dark block is one
305-dimensional quantum outcome. It is not replaced by a finer basis
measurement that would erase coherences within that eigenspace. The
two-outcome scalar-number detector is independently
`{P_empty=|Omega><Omega|, P_occupied=I-P_empty}` on exactly this space.
These projectors need not commute with H, so sequential number records
include their real measurement backaction. The complete sequential laws are
computed in `higgs_records.py`.

For any complete projective menu `{Pi_r}`, choose a pointer with M orthogonal
labels, where M is a power of two at least as large as the menu. The pointer
starts at zero and its shift is `S|j>=|j+1 mod M>`. Define

```math
U_{\rm pointer}=\sum_r\Pi_r\otimes S^r.
```

Orthogonality gives `U_pointer† U_pointer=I`, including every input pointer
column. Its blank-pointer action is the isometry
`V psi=sum_r Pi_r psi tensor |r>`. Unused pointer labels have zero ideal
maps and remain in the alphabet. Committing the pointer gives the complete
instrument

```math
\mathcal I(\rho)=\bigoplus_r\Pi_r\rho\Pi_r,
\qquad p(r)=\operatorname{Tr}(\rho\Pi_r).
```

Every identity persists after tensoring an arbitrary untouched reference.
The code checks the full operator identities, rather than only testing a
ground-state distribution. It also tests a degenerate-block coherence and
an entangled reference control. A ground input has spectral probabilities
`(1,0,0)`; the zero-weight upper and dark outcomes remain in its experiment.
Other preparations may occupy all three outcomes.

The returned label names a spectral value in the chosen Hamiltonian unit.
It does not name an electron, a relativistic mass, or a value in GeV.
With the same source ground as the energy reference, the three reported
values are respectively `0,15,8`; the unshifted values `-7,8,1` use the
Hamiltonian's declared empty-state reference. This additive shift changes
no probabilities or sequential gap.
For a finite unresolved or numerically approximated spectrum, the existing
[finite-resolution spectral readout](../spectral_readout/README.md) supplies
the complete positive Fourier-clock POVM and explicit alias/resolution
bounds. Sharp spectral projectors with computable entries are available in
the exact three-block source here; no general equality oracle for arbitrary
computable eigenvalues is assumed.

## Preparation and one operational clock

The owner constructs a gauge-invariant phase-shifted vertex preparation
unitary from the same coupling tensor. On the vacuum/bright sector it sends
the empty configuration to the exact lower-energy vector and fixes every
dark vector. Its inverse and its complete matrix are retained. The empty
configuration is an admitted initial preparation. This is a coherent
preparation protocol, with number terms paused or refocused as specified by
the protocol, rather than an assertion of autonomous relaxation.

Choose one macro tick `delta_t=pi/45` in the fixture's Hamiltonian unit and
one fixed finite word approximating `U=exp(-i delta_t H)`. Repeating that
same word n times represents time `n delta_t`. Preparation and detector
services occur between these ticks; their resource counts are separate.
The declared schedule pauses or refocuses autonomous H during those
services. A design that lets H act during them must include that evolution
in its apparatus matrix and error budget.

This supplies a common operational program clock for the one-tick and
two-tick sequential experiments. The lower-level pulse word and its bound
are counted explicitly. The construction does not equate compilation cost,
repair-attempt count, modular parameter or laboratory seconds. In
particular, compiling a different arbitrary-time unitary at each requested
time and ignoring the different overhead would not establish this fixed
macro-clock contract.

## Reuse of the finite instrument compiler

The constructive theorem in
[spin_operational_position, section 3](../spin_operational_position/PROOF.md#3-every-computable-finite-instrument-has-an-explicit-source-word)
covers every specified complete computable finite instrument. It
uses admitted coherent addressed load/export, blank pointers, refined
integrated-port controls and classical pointer commitment. It retains fine
Kraus/environment labels, unused pointer labels and source export failures.
Its approximation bound is uniform over all inputs and spectator references.
`readout_contract.py` imports its existing `generic_instrument_plan`;
it introduces no competing compiler or synthesis bound.

An explicit isometry J sends the ordered 307-element source basis to the
first 307 computational states of nine logical bits. The other 205 states
are retained. Every measurement projection is extended as `J Pi_r J†`,
with an additional `I-JJ†` outcome labeled outside the field encoding.
Every source unitary is extended by identity on that complement. Thus the
detector admits all 512 input columns, with no code postselection. The
energy menu has four nonzero projectors after padding; the number menu has
three, plus one unused zero pointer label.

This encoding is a chosen computational realization of the declared finite
field and its apparatus. It is not an unproved inference that bare
occupation repair implements arbitrary coherent operations. For the wider
original scalar/CAR carrier of dimension `3*2^18=786432`, the same observable
can be extended with its `786125`-dimensional outside-conserved-sector
projection. That outcome is retained without assigning it an energy
calculated from the 307-state restriction. The entire ambient spectrum is
not computed here. A binary implementation of that wider domain must use
its larger register and compiler resource count.

For q logical spins and m fine Kraus branches, the inherited theorem uses
`D=2^(q+ceil(log2 m))` and at most `20D²+4D+10q` primitive pulses. With the
declared refined controls, it constructs a complete instrument within any
requested positive half-diamond error. In the 307-state encoding, the
four-output energy detector has q=9, m=4 and a pulse upper bound of
`83,894,362`. One extended unitary tick, or the preparation unitary, has
q=9, m=1 and upper bound `5,245,018`. These are finite existence/resource
bounds, not a claim that those long tapes were executed or that the
construction is efficient.

An apparatus chosen as a source primitive instead has its declared exact
matrix and exact complete record law. Proving that this primitive follows
from a narrower instruction set is required only when that narrower-source
claim is made. Configuration is permitted; using measured particle masses
to set an input does not make those same masses independent predictions.

## Complete histories and their error bound

Let the prepared input have complete half trace-distance error e0. Let the
preparation, one fixed macro tick and each detector have complete
half-diamond errors `ep`, `es`, `ed` from their declared ideal operations.
Retain a common first-failure alphabet and make termination an explicit
absorbing record; ideal failure outputs have zero probability. For n ticks
and m detector uses, hybrid replacement and CPTP contractivity give

```math
\tfrac12\|\widetilde\rho_{\rm complete}-\rho_{\rm complete}\|_1
\le\min(1,e_0+e_p+n e_s+m e_d).
```

The terminal object includes every quantum output, classical record and
untouched reference. The same bound therefore controls total variation of
the full history and the absolute probability error of any event in that
history. The argument applies to common record-adaptive continuations;
each conditional implementation must satisfy its uniform bound. It uses
no independence between successive quantum states and divides by no
success probability. Comparing different successful subensembles would
require another contract.

At `e0=0` and `ep=es=ed=10^-6`, two ticks and two number readouts have
complete event error at most `5*10^-6`. This error belongs to a declared
compiler precision, not to a sample-size calculation. Finite-shot uncertainty
must be added from the chosen sampling model. The gap inversion in the
owning record calculation uses unconditional events from the one-tick and
two-tick experiments and keeps its nonzero-denominator and unaliased-phase
conditions explicit.

## Born weights, state updates and consensus

The finite theorem
[`finite_busch_gleason`](../../Lean/EventAlgebra/FiniteBuschGleason.lean#L1133)
represents normalized nonnegative valuations additive on coexistent complex
matrix effects by a unique density matrix. The
[operational bridge](../../Lean/EventAlgebra/OperationalBornBridge.lean#L214)
adds equal-effect probability agreement and common disjoint refinements.
The chosen finite quantum-source semantics satisfies these requirements,
and the complete instrument calculates its concrete record probabilities.
This is a same-source use of the finite quantum formalism, not a derivation
of arbitrary complex operations from the narrow classical occupation
repair alone.

The pointer coupling also fixes the state update. Born weights alone do
not select it. The
[certain-state theorem](../../Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean)
characterizes the Lüders map under its stated Kraus/effect and invariance
premises, and the
[finite calibration protocol](../finite_instrument_calibration/PROOF.md)
provides a separate experimental identification and reference-stable error
bound if the apparatus is only known through trusted calibration probes.
No extra calibration assumption is needed for an apparatus whose exact
complete coupling is part of the chosen source specification.

Consensus copies and compares orthogonal detector records. It preserves the
calculated outcome weights; it does not select outcomes by agreement or
manufacture a particle identification. The finite interacting record law,
its chosen clock and a controlled apparatus are available on this domain.
A spatial relativistic limit and matching to laboratory species are distinct
claims requiring their own observable maps.

## Reproduction

```sh
python3 -m pytest -q model/particle_observer_outputs/test_readout_contract.py
```

The executable certificate binds the same Y, d, scalar tags, Hamiltonian and
projections as the owning Higgs records. It reports complete projector ranks,
padding and pointer menus, inherited compiler bounds and the multi-time
error allowance. It exports no physical particle or mass label.
