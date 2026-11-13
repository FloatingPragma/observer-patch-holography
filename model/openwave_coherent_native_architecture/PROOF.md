# A complete coherent native event

The architecture uses the native C/E endpoint permutations, a fair coherent
law bit and a fair coherent placement bit. Coherent storage and proper-code
load/export are declared services. The construction does not require these
choices to be unique.

## Complete event and recorded marginal

There are three source-address bits, three raw-record bits, one placement
bit and one law bit. On all 256 addresses the event applies the exact C
repair swaps from the source compiler; E additionally copies the output
parity into the high record bit. Each factor is a permutation. Source
addresses 0 through 5 encode all local occupations of total at most two;
addresses 6 and 7 are retained padding. This is the local factor of the
existing nine-bit source encoding, with spectator bits untouched.

Prepare the law and placement bits in their equal superpositions and the
record in zero. The resulting columns have four amplitudes 1/2 with disjoint
joint support. Thus V*V=I, including source padding and arbitrary references.
Measuring the original record, placement and law bits gives exactly
K_r=sqrt(w_r)P_r. Direct independent floor/ceiling enumeration checks all
36 local physical matrix units. All 16 binary matrix units are compared to
all 25 slots of the native process. The zero hold slot and every zero raw
slot are retained. The full encoded register has 32 labels; its additional
C high-bit padding labels have probability zero on physical inputs.
Kraus completeness and tensoring identity establish the arbitrary-reference
identity and local no-signalling. The finite Kraus and reference lemmas in
`Lean/Dynamics/ChoiCPTP.lean` and
`Lean/ObserverPatchHolography/Locality/NoSignalling.lean` apply directly.

## Reversible record re-encoding

For binary endpoints set sigma=coin xor(left-only). Native repair then equals
a swap of the endpoints when sigma=1 and identity when sigma=0. The joint
output, original coin and raw record recover the input. A reversible
permutation therefore sends the sixteen binary event columns to the source
after that controlled swap, the same law bit, the sigma bit and blank raw
record. The full 256-address extension is explicitly specified and executed.
Its inverse recovers the original readable record. No classical committed
record is erased: this operation requires uncommitted coherent registers.

For a uniform seam register the re-encoded event is the controlled-swap
isometry S A, with a factorized equal law bit. Original native records can
be recovered and measured by the inverse recoding. Tracing the re-encoded
coin instead is a different quantum channel: the source matrix unit
|10><01| maps to zero under unread C/E, and to
(|10><01|+|01><10|)/2 under the traced swap coin. Both channels have the same
one-load population transition. Repeated coherent reflection is a declared
continuation; it is not equated to repeated committed repair.

## Source compilation and errors

The event and recoder are compiled into Gray-path adjacent transpositions.
An adjacent transposition is a multiply controlled bit flip, implemented
by reversible Toffoli ladders with five scratch bits. Exact integer
execution checks every one of the 8192 clean and dirty scratch inputs for
each service. All 256 clean inputs return clean scratch with the specified
permutation. Dirty inputs are retained by the complete circuit.

The Toffoli, flip and Hadamard routines execute the original M6 source
exponentials on all their input columns. Their checked Frobenius errors
bound operator errors. If G routine calls each have error at most epsilon,
then the full product error is at most s/(1-s), s=G epsilon<1. This follows
from product expansion and includes arbitrary reference systems. The
reported epsilon=10^-10 is a checked floating guard, not an outward
physical service certificate. Each pulse requires a proper-code load and
export. The exact block-preserving generators make their complement
probabilities zero on valid inputs; the complete export instruments retain
all invalid-processor failures. Address padding and processor failure are
different outcomes. No polar correction or normalization is applied.

The finite circuit provides an implemented coherent architecture under the
stated proper-code services. Initialization, storage, reset energy and
laboratory error bounds are not inferred from numerical gates. The native
clock, particle interpretation and continuum field transfer require their
own observation and unit maps. No OpenWave physical criterion is awarded.
