# An executed weak detector on the encoded family chart

Platform row 11 remains untested. This construction implements the finite,
unclocked family detector on a new explicitly encoded apparatus. It leaves
the coherent source services and their load/export semantics declared. It
does not identify a charged lepton or a physical rest-energy threshold.

## 1. The chart and the actual native family tick

Each original seam object has the source index

    i=4 load+2 record+orientation,       i=0,...,11.

Use its literal four-bit encoding. The two-object isometry is

    E|i_C,i_E> = |i_C+16 i_E>.

It selects 144 addresses of the existing 512-address nine-bit work register.
All other 368 work addresses are padding; choose address 511 as the one
ideal destructive output. The original H66 wrapper, codebook and validator
are unchanged. This is a separately specified source chart and export.
The family uses two additional logical bits, with codes 0,1,2 carrying its
original cyclic register and code 3 retained as an extra output. Every object
filter acts as the identity on the family and any external reference.

The single-object source operator is B=V Q R2. Its Q is the actual coherent
`Observer.repair`, compared at all 360 admitted oriented local columns.
Let the four object bits be orientation o, record r, low-load l and high-load z.
For z=0, Q is the reflection about the normalized vector

    u_o = X_l^o (|00>-|01>-|10>-|11>)/2

on l,r. For z=1, Q is the identity. A source circuit prepares u_o with
Hadamards, Z phases, one controlled Z and the controlled flip o→l. Conjugating
the phase on l=r=z=0 by this circuit gives Q. Explicit source words then
exchange the orientation and apply V. A padding-only inverse of V R2 makes
the load-code-3 extension exactly identity. The initial implementation omitted
this last padding correction. Its failed full-16-address comparisons are
preserved; all original 12-state comparisons already agreed.

The controlled family shift is the exact permutation

    |i_C,f> -> |i_C, f+step(i_C) mod3>,

for i_C<12 and f<3, with the other addresses fixed. The existing Gray-path
compiler emits its full six-bit permutation using three clean scratch bits.
Every one of the 512 dirty-scratch addresses has a specified bijective action.
The certificate checks that full truth table and executes the actual native
flip/Toffoli blocks on all 64 clean input columns, including their phases and
dirty-scratch outputs. A clean-scratch export is part of this new apparatus;
any nonzero dirty output stops the run with a retained flag.

The inverse shift is independently compiled from the inverse permutation.
The inverse object word consists of actual reversed negative source requests.
No returned floating array adjoint is substituted for either inverse service.
Together these words agree with the original 432-state unclocked U0 on all
family columns in the stated chart. They also specify the padding extension.

This uses the source's coherent proper-code services. Classical native-byte
equality alone does not select those coherent controls, and this packet does
not claim to derive them from the bare repair interpreter.

## 2. Why the family ray is now admitted

The original analyzer Phi is the normalized sum of its three orthogonal flags.
Its actual source coefficients are positive. Eight nonzero squared amplitudes
are 1/24 and eight are 1/12. The original product preparation p likewise has
positive coefficients and rational squared amplitudes. Both statements are
read from `finite_family_clock.primitives`, and normalization is checked
exactly. No magnitude or sign is chosen from a measured mass.

The restricted H66 interface rejects their dimension and the relative sqrt(2)
amplitude. The underlying source `_state_word` has a different input contract:
all 512 exact nonnegative squared amplitudes summing to one. It therefore
admits E Phi, E p, and their mixtures with one padding ray. Its native word
implements every conditional binary split, including zero-weight prefixes.
No diagonal preparation phase word is required because every actual phase
is zero. Omitting that identically zero word is an explicit circuit choice;
the retained amplitude-tree schedule still contains all its zero requests.

The ideal split angles are approximated on the existing grid Q=2^41, with
error at most 1/Q per split. At one tree level the conditional rotations act
on orthogonal prefixes, so its operator error is bounded by their maximum,
not their number. Nine levels give full unitary target error at most 9/Q.
Every inverse-Walsh equality is checked exactly by the owning compiler.
This is a finite approximation, not an assertion of exact transcendental
angles or arbitrary precision at fixed hardware.

The nine-bit tree has 12,291 native pulses. A complete diagonal phase word has
11,270. Each rank-one phase below therefore has

    2*12291+11270 = 35852

pulses. Every word is below the unchanged 65,536-pulse admission ceiling.
The source evaluator's generic six-qubit resource stop is also unchanged;
the direct source-word array executor is used and named explicitly.

## 3. The actual weak Phi instrument

Write P=|Phi><Phi|, s=3/5, epsilon=1-s=2/5 and theta=arccos(s).
Let W prepare w=(E Phi+|511>)/sqrt(2) from literal work address zero.
Execute its actual inverse word, the source phase 2 theta at address zero,
and W. The ideal full work operator is

    U_w = I+(exp(2i theta)-1)|w><w|.

Coherent chart export gives the successful object update

    E* U_w E = I+(s exp(i theta)-1)P.

The ideal padding row is

    <511|U_w E = (exp(2i theta)-1)<Phi|/2.

Its effect is (1-s²)P. Every padding row belongs to the export's one observed
flag; the implementation does not introduce a measurement that resolves all
padding addresses. Those rows are retained as a Kraus representation of that
flag and preserve the untouched family/reference state.

The second actual rank-one phase word uses E Phi without padding and phase
-theta. Its success operator is I+(exp(-i theta)-1)P. Multiplication gives

    L = I+(s-1)P = I-epsilon P.

The correcting word has zero padding target, with its actual numerical padding
and every processor first-failure output retained. The successful update is
derived from these words. It is not inferred from an effect and replaced by
a chosen Lueders operation. The ideal complete effect identity is

    (1-s²)P + L*L = I.

The new certificate executes those words on the actual family state and on
fixed full-chart/reference probes. A uniform source-word error bound supplies
the all-input statement, including entangled references; the finite probe
menu is not substituted for that operator argument.

## 4. The unclocked family return is an actual branch

Prepare the original p using its source amplitude-tree word and chart export.
Keep preparation-padding and processor-failure outcomes. On success, execute
U0, the two-word Phi filter, and the actual U0 inverse. The final object
detector is a destructive test of p using the same rank-one phase construction
at s=0. Its padding row is -<p| and its residual is I-|p><p|.
The padding flag returns the family while destroying the object ancilla.
Its complete source residual is the nonreturn outcome and stops this trial.

Let M0=<Phi|U0|p> and A0=M0* M0. Independent exact contraction of the original
432-state tick gives

    M0 = I/(3sqrt(3)) + (R+R*)/(3sqrt(6)),      R³=I.

The actual returned family Kraus map is

    K_return = -<p| U0* L U0 |p> = -(I-epsilon A0).

The minus sign comes from the actual destructive padding row and is retained
in the source comparison. It cancels in the returned branch effect. The
complete ideal coarse effects are

    first Phi click:       (1-s²)A0 = (16/25)A0,
    ancillary nonreturn:   epsilon²(A0-A0²) = (4/25)(A0-A0²),
    destructive return:    (I-epsilon A0)².

Their sum is I. The actual finite source data retain preparation padding,
forward and inverse scratch failures, correction padding, all 368 object
padding rows, every family-code-3 output and the final nonreturn residual.
No successful path is normalized. Tensoring these actual branch maps with
the identity proves the arbitrary-reference interpretation.

The first-click and nonreturn state updates are their emitted maps. Their
effects alone are not used to invent terminal family states. The saved arrays
retain their full source/family/reference amplitudes for every family input
column, so their original effects can be reproduced independently.

## 5. Controls on the same family observable

Omitting the compensating phase preserves the first click effect but changes
the surviving coefficient to c=-16/25+12i/25. The returned family map becomes
-(I+c A0). The executor detects its nonzero imaginary part and its difference
from the compensated return. This control tests the actual state update,
rather than merely testing the same first-click effect twice.

The second fixed control replaces Phi by its normalized hold flag and uses
the same source compiler and complete sequence. Its first-click family effect
differs from (16/25)A0. This is an apparatus control with no mass comparison.
The actual unclocked native family compression is the positive planted signal.
Normalization, padding counts and projector identities are consistency checks;
they discharge no physical row clause by themselves.

## 6. Full numerical and record error accounting

Every emitted primitive uses its actual full six-level source array. The
interval checker compares every distinct kind/angle request with its ideal
full source matrix, retaining negative and zero requests. The bound
g=2*10^-11 includes the stored primitive envelope and the fixed pi-payload
allowance of the owning source. Logical placement and tensoring a reference
preserve operator norm. All proper-code first-failure outputs have their
original prefix maps and are exposed by a lazy query.

For a block of n pulses and independent target-angle error t/Q, put

    x=n g,
    e=x/(1-x)+t/Q,
    f=n g²/(1-x)²,
    delta=(e(2+e)+f)/2,
    gain=1/(1-x)²+f.

Then e bounds the stored good-code product, f bounds all first processor
failures, and delta bounds half diamond distance of the complete exported
instrument from its ideal target. The chart and clean-scratch exports are
declared trace-preserving instruments, so they do not increase the bound.
For the object preparation use t=9. For the Phi weak word use t=20, and for
its correction t=19. The destructive preparation detector uses t=20, a
conservative bound covering the separately approximated pi phase. The small
fixed object/permutation circuits carry their fixed pi error inside g and an
additional conservative grid allowance stated in the receipt.

Extend every block as the identity on earlier terminal flags. Sequential
replacement by ideal instruments gives

    total_halfdiamond <= sum_j delta_j product_(k<j) gain_k.

The receipt recomputes this rational bound after every circuit edit. It is
distinct from the 3e-8 numerical array-comparison gate. Floating propagation
is a diagnostic of the stored products; the uniform bound applies to their
exact mathematical products. No hardware fidelity follows from either one.

The full source run exceeds the small-script target and its runtime is retained.
A separate bounded reader verifies the recorded source data and reproduces
the family maps; it does not claim to recompile the entire source experiment.

## 7. What this implements, and what remains open

The new explicit chart and emitted native words implement an actual weak Phi
filter and destructive p-return for the unclocked finite family register.
This removes the earlier H66 coordinate-interface barrier for this particular
encoded apparatus under its declared coherent source services.

The same local Phi detector tensor identity can act with an arbitrary retained
count archive as reference. A return test of a full count preparation requires
its separate encoding, word and export. It is not supplied by this nine-bit
object detector. No N>0 count apparatus is claimed executed here, nor is a
growing-register compatibility or resource theorem silently inherited.

At fixed angle grid and primitive error, a limit epsilon→0 is not certified by
repeating these finite error bounds. A uniform refinement of pulse precision,
source capacity and detector resources would be a further task. The supplied
family architecture, phase preparation, observer calibration and coherent
source contracts retain their empirical ancestry. A common physical clock,
electric charge and three lepton rest minima are still unresolved.
