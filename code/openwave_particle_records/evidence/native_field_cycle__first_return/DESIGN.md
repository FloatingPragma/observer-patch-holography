# Coherent first-return field cycles

Candidate CAR-first-return-field-clock-v1 extends the sealed native_field_cycle
source on exactly the same486 CAR/Gauss columns, q=0,1, L=1, all memory and
electric addresses. It replaces only the loop operation in the complete
common clock by the first-return map of that actual native word. The invalid
source flag is identity and returns0. No empirical value enters.

The theorem sought is a complete same-source operation F that preserves
every final tagged position and spin, changes only its memory/electric
records, and has a common coherent execution horizon including all boundary
cases. Native U itself is the original three-seam permutation. For every
valid state s define r+(s) as the first t>=1 with tag(U^t s)=tag(s), and
r-(s) using U*. F(s)=U^r+(s)s. A finite permutation returns to its initial
state within486 steps, so both lengths exist with1<=r<=486. The backward
first return from F(s) equals r+(s). Signed CAR phases are retained.

## Complete reversible implementation

A nine-bit length register has512 values. C+ XORs r+(s) into that register,
and C- XORs r-(s). Each oracle has the following declared reversible circuit:
copy the computational source address into a scratch address by modular
addition; execute486 controlled actual U or U* words on scratch; record all
486 tag-equality bits; XOR the first-hit priority encoder into the output;
reverse all486 words and erase their equality bits; undo the initial modular
copy. All scratch returns to zero. The copy entangles basis labels and is
not an arbitrary quantum-state clone. Comparators, modular addition and
priority encoding are declared classical reversible services.

V executes U exactly c times using511 fixed time slots, applying U in slot t
iff t<=c. Its action on every counter value is defined. The full operation
C- V C+ maps input (s,0) to (F(s),0), since r-(F(s))=r+(s). On arbitrary
nonzero counter inputs the same three reversible gates define the full
complement. Enumerate all487*512 columns and their inverse, including invalid.
The full scratch circuit is a composition of permutations on all ancillary
inputs; its specified preparation is scratch0. The audit retains all hit
bits and intermediate scratch source addresses before their uncomputation.

Every amplitude uses2455 native loop-word slots:972 for each return oracle
and511 for V, hence7365 seam slots, plus fixed reversible arithmetic. Gates
that are not selected are identity. This is a declared macroclock, not a
claim that all route lengths have equal native-hop count. No orthogonal
return-time history is left in the coherent output. Keeping that history
defines the explicit dephasing control.

## Records and comparisons

Execute every actual native loop and inverse column, compute their finite
first-return paths, and retain every path length, intermediate address,
realized current, electric increment and sign. F must commute with the tag
and every charge-density pointer; every field increment must have zero
divergence. Its inverse must equal the backward first-return map.

The same positive lazy amplitude is (6I+F+F*)/8. Replace the loop dilation
in the sealed common clock by [[0,F*],[F,0]], with unchanged eight-address
signal/reflection/current service. Regenerate every complete reference law
and all Gauss0/uniform-memory spin-resolved stationary states. There is no
nonreturn term at this macrolevel, so exact stationarity reads
<selected charge current>=<first-return field electric change>. Underlying
charged loop microcurrents remain visible and equal minus that change.
Calling this macrofield operation an autonomous physical magnetic field
requires its spatial measurement/clock attachment; that is not assumed.

The coherence experiment prepares the uniform Gauss0 source basis, with
every memory value. A complete pointer for the ideal F output and its
orthogonal complement detects whether return length leaked. If length r is
retained, the overlap probability is sum_r p_r^2; full uncomputation gives1.
Both source-pointer outcomes and all512 counter values are retained. A
two-length plant must show the predicted1/2 dephased overlap. This control
tests exactly the instrument difference relevant to the first-return claim.

The physical row5 clauses are unchanged: stationary closed-cycle component
balance is assessed; no far-zone dipole, orbital g or Larmor result is
inferred from a triangle or a Pauli identity. The exact finite theorem stack
and original documentary snapshots are inherited with their hashes. The
source identity is proved through actual native columns and the reversible
controller, not through coincident candidate names.
