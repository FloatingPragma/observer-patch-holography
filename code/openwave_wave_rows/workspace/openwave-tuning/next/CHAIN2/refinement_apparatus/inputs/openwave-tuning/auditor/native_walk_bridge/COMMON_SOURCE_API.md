# Native local twelve-event coherent walk

The proposed common source identifier is
`native-local12-current-v0-proposed`. Its development member is exactly
`model.native_spatial_refinement.source.graph(2)`: eight cubic cells,
ninety-six ports and 288 native seams. The C10 source and the global-event
walk are comparison constructions. No source equivalence to them is
inferred from matching labels.

## State and interfaces

The prepared data array has shape `(128,16)`. Position address
`(cell<<4)|port` has `cell=4*x+2*y+z` and `port<12`; local ports 12 through
15 are source padding. Event address `2*direction+sigma` has directions
zero through five and a move/hold bit sigma. Events 12 through 15 are
event padding. All addresses are retained.

`common_walk.py` exposes the shared interfaces:

| Function | Contract |
|---|---|
| `source_descriptor()` | Complete geometry, encoding, scheduler and clock declaration |
| `native_transfer()` | The 96 by 96 local-attempt population matrix |
| `transfer_numerators()` | Its integer numerator matrix, denominator twelve |
| `initial_cell_state()` | Source-word preparation of the twelve-port cell and uniform event register, after separately retained blank resets |
| `apply_walk(state)` | One source-compiled `W=S R`, on every padded source/event amplitude |
| `complete_position_records(state)` | All joint source/event probabilities and separate source/event padding weights, without normalization |
| `raw_record(position,event,law)` | Original native seam, coin and current for the terminal inverse-recoding intervention |
| `compile_walk.pulse_word()` | Explicit ordered M6 pulse records for one W, with scalar phases |
| `compile_walk.preparation_word()` | Full event, local-port and law-bit preparation word |

Other workers attach measurements to this W. A retained current register
requires its declared joint preparation and a proved lift: a definite
which-path register can change reduced coherent propagation. An unproved
trace over that register is not an identity source map.

## Local direction and source event

The thirty internal icosahedral seams split into the five perfect matchings
listed by `source_descriptor()`. The inherited external seams form a sixth
matching. Each port belongs to exactly one edge of each colour, so reversal
preserves the direction label. The controlled involution S holds sigma-zero
columns and exchanges the endpoints of the indicated matching on sigma-one
columns. It fixes source and event padding. Relabeling the six neighbour
directions by these colours changes no walk or position observable.

For a uniform twelve-event isometry A,

\[
A^*SA=T_{\rm local}=I-\mathcal L/12.
\]

Every port has self numerator six and one numerator for each neighbour.
This matrix is reconstructed from the executed reversible circuit and
compared with the actual native endpoint columns. Its self coefficient
retains the native coin hold.

The pulse program reuses the full native eight-bit C/E event and its
reversible record recoder. Endpoint extraction copies only orthogonal
computational information into three work bits. A source-controlled coin
flip sets the original native placement coin, the actual native event runs,
and the recoder exposes sigma. The geometry swap updates the source address;
endpoint work is then uncomputed. All 4,096 source/event/law columns,
including both paddings, execute this circuit and return clean work and
scratch. Its action extends to arbitrary untouched quantum references.

The inverse record intervention is emitted too. For a valid terminal
source position y and sigma, the preceding source position is
`x=matching(y)` when sigma is one and x=y otherwise. If `(a,b)` is the
canonical oriented seam, the original native coin and current are

\[
c=\sigma\mathbin{\mathrm{xor}}1_{x=a},\qquad
j=1_{x=a}-1_{y=a}.
\]

The inverse recoder restores the actual raw C/E labels. All 2,304 terminal
source/event/law combinations replay the literal native endpoint rule.
This is a terminal measurement convention. Coherently reused event
amplitudes do not constitute fresh uniform coins at each later W step.

## Reflection, spectrum and padding

The uniform event preparation factors as a three-component state on two
bits and two Hadamards, since twelve equals three times four. Every factor
is emitted through the existing M6 scalar/Pauli compiler. The reflection
word is the inverse preparation, followed by
`diag(1,−1,...,−1)`, followed by preparation. Its scalar phase is retained:

\[
R=2|u_{12}\rangle\langle u_{12}|-I_{16},\qquad W=SR.
\]

All sixteen reflection columns are source-executed and compared in full
matrix phase. A global minus sign would shift controlled-W quasienergies
by pi and fails this comparison. Event padding has eigenvalue minus one.
Source padding remains a separately flagged sector on which W acts as R.

For `psi=Aa+SAb`, direct multiplication gives

\[
(a,b)\mapsto(-b,a+2T_{\rm local}b).
\]

Thus a Laplacian eigenvalue lambda gives unit-circle phases with
`cos(omega)=1−lambda/12`. The graph is connected, regular and contains
triangles, so its transfer spectrum lies in `(0,1]`, with one constant
eigenvector at one. All ninety-six modes are checked against the same
compiled W. On the 1,152-dimensional valid source/event space the full
multiplicities are 194 at plus one, 768 at minus one and 190 dispersive
modes. Event padding adds 384 minus-one modes. Source padding adds 32
plus-one and 480 minus-one modes. None is removed by a fitted dispersion
or species assignment.

## Scheduler and three clocks

The local scheduler is a declared amendment. With M=288 global seams,

\[
T_{\rm global}=(1-6/M)I+(6/M)T_{\rm local}.
\]

A global attempt contacts the one-load source with probability 1/48.
The complete classical thinning recipe retains every nonincident seam,
coin, law and zero current. The number of such letters before contact has
probability `(47/48)^k/48`. A finite wait cap retains its timeout; it never
silently replaces the wait by one local tick. Corresponding Poisson
population generators agree when the local rate is the global rate divided
by 48. This is a diagonal classical-process identity. It supplies no
equivalence of raw quantum channels or of the two coherent walks.

The native attempted-event clock, ordered M6 pulse/service count and
coherent W-step clock are distinct. One W contains 432,425 source pulses
and 864,850 proper-code load/export service events. These are logical
operation counts, not physical duration. The source supplies no measured
light speed, Lorentz limit or massive scalar equation through this count.

## Complete apparatus and accuracy

The coherent program uses twenty-three work bits: eleven data bits, three
endpoint bits, three raw-record bits, one law bit and five scratch bits.
The retained native full-sector reset supplies each blank. Native pointer
readout has its own auxiliary reset before each measured work bit. Any
reset timeout retains its flag, final work and fine native tape. The pulse
compiler retains the complete M6 export channel, including every invalid
processor outcome. Source padding, event padding, dirty work and processor
export failure are different flags.

The receipt executes full M6 primitive matrices and the complete four-bit
reflection word. The reversible twenty-three-bit circuit executes through
its exact Boolean factorization, with complete clean source/event/law
columns and complete local dirty-scratch truth tables. It does not claim a
dense floating multiplication of all 432,425 pulses on the entire register.
The ideal circuit identities follow from the stated primitive and proper-code
services. Floating residuals diagnose their implementation; they do not
certify accumulated physical pulse error. A physical error bound requires
certified primitive errors composed over the emitted word.

`COMMON_WALK_PROGRAM.json` contains the routine calls and source pulses.
`COMMON_WALK_RECEIPT.json` binds the program, source ancestry, complete
reset/readout recipes, transfer matrix and spectral checks. Coherent
memory, code load/export, routing and phase-preserving control are explicit
architecture choices. No official OpenWave status or common physical
completion is claimed.
