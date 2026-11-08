# Complete native weak-measurement menus

The carrier is the weighted binary two-load source of dimension 66. Its
embedding C uses the existing nine-bit source codes in a 512-dimensional
work register. All coherent requests are the X, Y, Z, I, H and J pulses in
`observer/source_scalar/compiler.py`. Proper-code load, coherent gluing,
processor blanks, complete export and Kraus trace semantics are declared
source services. The repair/hold interpreter does not itself issue these
coherent words. No gate, register or host limit is enlarged here.

`native_effect_synthesis` constructs individual rank-one effects. The
construction below supplies an instrument with several jointly accessible,
possibly noncommuting effects and a retained quantum residual. It uses the
same coherent source export. It does not insert the destructive native tree
detector between stages.

## One weak filter on the source

Let v be a nonzero Gaussian-rational 66-vector, u=v/||v|| and P=|u><u|.
Let s be rational with 0<=s<=1. Every supplied numerator and denominator has
at most 64 bits. Choose one of the 446 existing padding addresses p and set

\[
 w=(Cu+|p\rangle)/\sqrt2,\qquad
 \theta=\arccos s\in[0,\pi/2],\qquad \alpha=e^{2i\theta}.
                                                               \tag{1}
\]

The amplitude-tree and diagonal-phase construction in
`native_effect_synthesis/synthesis.py` supplies a unitary W with W|0>=w.
Its squared amplitudes are rational even when ||v|| is irrational. The
existing complete diagonal phase word supplies a phase only at literal work
address 0. Executing W inverse, that phase, and W gives

\[
 U=W\bigl[I+(\alpha-1)|0\rangle\langle0|\bigr]W^*
   =I+(\alpha-1)|w\rangle\langle w|.                            \tag{2}
\]

The source-success part of the complete address export is

\[
 K=C^*UC=I+\frac{\alpha-1}{2}P
        =I+(s e^{i\theta}-1)P.                                 \tag{3}
\]

The only nonzero ideal padding row is
<p|UC=(alpha-1)<u|/2. Its effect is (1-s^2)P. All 446 padding rows belong to
the export's one observed abort flag; address-resolved padding measurements
are not introduced. The abort retains the untouched reference and stops the
experiment. Its source state is not passed to a later stage.

Conditional on the unnormalized source-success branch, execute a second
rank-one phase word with source vector Cu, no padding amplitude, and phase
-theta. Its ideal source operator is

\[
 V=I+(e^{-i\theta}-1)P,\qquad
 L=VK=I+(s-1)P.                                                 \tag{4}
\]

The second word's ideal padding rows vanish. Their actual numerical rows
and all processor failures remain in the complete output. Equations (3) and
(4) prove the successful state update; it is not inferred from the outcome
effect. In the ideal instrument, the abort effect and surviving effect obey

\[
 (1-s^2)P+L^*L=I.                                               \tag{5}
\]

At s=0 the construction is a destructive test of P with complementary
quantum residual. At s=1 it is the identity target with a zero abort effect.
Both endpoints retain the entire emitted word and all zero failure labels.

## A complete coexistent menu

For 1<=m<=8 choose the stated inputs (v_j,s_j). Put A_0=I and
A_j=L_j A_{j-1}. Stop at the first abort and otherwise proceed to the next
stage. The ideal menu has effects

\[
 E_j=A_{j-1}^*(1-s_j^2)P_j A_{j-1},\qquad
 E_{m+1}=A_m^*A_m.                                             \tag{6}
\]

The last outcome retains the quantum source operator A_m rho A_m^*. An
abort retains the reference operator obtained by contracting the source
row <u_j|A_{j-1}, multiplied by 1-s_j^2. These formulas extend by linearity
to arbitrary operators and by tensoring with the identity to any untouched
reference. They include every earlier measurement's disturbance.

Each E_j is positive. Equation (5) gives
E_j=A_{j-1}^*A_{j-1}-A_j^*A_j, hence

\[
 \sum_{j=1}^{m+1}E_j=I.                                       \tag{7}
\]

All entries of the effects, successful contractions and residual in (6) are
Gaussian rational. `ideal_menu` evaluates them exactly, without square roots
or floating matrix diagonalization. Any two subsets of this one outcome
set define coexistent effects. Their joint four-event refinement is obtained
by the intersection, the two differences and the complement of the subsets.
For disjoint subsets, the sum of their effects is bounded by I. These are
coarse-grainings of one executed instrument, rather than separate tests
whose probabilities are added without a common experiment.

An explicit noncommuting example takes v_1=e_0, s_1=3/5,
v_2=e_0+e_65 and s_2=4/5. On the indicated two-dimensional corner,

\[
 E_1=\begin{pmatrix}16/25&0\\0&0\end{pmatrix},\qquad
 E_2=\begin{pmatrix}81/1250&27/250\\27/250&9/50\end{pmatrix}.
                                                               \tag{8}
\]

The (0,65) entry of [E_1,E_2] is 216/3125, which is nonzero. Nonetheless
E_1+E_2<=I by (7). Ignoring the first success contraction would instead use
(9/25)P_2 for E_2 and would produce the wrong complete measurement.
The residual effect includes the identity on the other 64 source directions.

## Native words and execution domains

The source amplitude tree uses 12291 pulses and its complete diagonal phase
uses 11270. A rank-one phase (2), including preparation and its actual
reversed negative word, therefore uses

\[
 2(12291+11270)+11270=58392                                   \tag{9}
\]

pulses. Each of the two words in a stage is individually accepted by
`native_program.interpreter.source_word`: nine bits, 58392 pulses, unchanged
65536-pulse admission ceiling, exact binary64 angle payloads and the existing
description-size limit. At eight stages there are 16 separately admitted
words and at most 934272 pulses along the surviving path. This is a bounded
composition of existing services, not an enlarged single admitted word.
The stage controller advances only on the coherent source-success branch;
every first failure terminates it. The stage and word labels are retained.

Generic `execute_word` returns its unchanged `RESOURCE_STOP` before a nine-bit
word: that evaluator accepts six bits and 256 pulses. The direct source
compiler's nine-bit tape interface supplies the coherent products. The
complete instrument is given lazily by the actual prefix products, processor
first-failure queries and source-address export. No expanded complete-history
run is claimed. Host exhaustion outside these evaluated diagnostic arrays is
not silently converted into success. The native record-tree readout is not
used in this experiment, so no tree host bound is extended.

At each pulse, the full six-level source array is used. A processor-complement
row resets only the addressed proper code, retains every other work bit and
reference, records its pulse and complement, and stops the entire menu.
After a word with no processor failure, coherent source export keeps C^*
of its work state on success. Padding export traces the padding work address
while retaining its reference, with all rows summed into one flag. The
success source is loaded through C before the next word. There is no
source-chart projection inside any coherent word.

## Finite-precision complete-channel bound

Let Q=2^41. The directed 120-bit integer interval construction in
`native_effect_synthesis/intervals.py` emits theta with error at most 1/Q
from the exact rational split (s^2,1-s^2). The vector preparation has full
unitary error at most 11/Q. Consequently the first rank-one phase has target
error at most 24/Q and the correcting phase at most 23/Q: two preparation
errors and respectively 2/Q or 1/Q for the central phase. Exact inverse
Walsh checks retain every zero coefficient. All phase angles have magnitude
at most pi, so the existing binary64 representability proof applies.

The ideal target and pulse-resource construction hold for every input in the
stated rational domain. A numerical certificate is conditional on
`certify_word` passing its exact stored-array checks for every emitted
primitive, including negative and zero requests. A failed primitive rejects
the certificate. There is no assumed all-input SciPy rounding theorem.
The fixed basis/parity angles carry the separate pi-payload allowance from
the same owner. Thus g=2*10^-11 bounds each stored primitive against its
required ideal gate. This concerns exact products of the stored arrays;
floating array propagation is a diagnostic, and hardware error is outside
the certificate.

For one word of N=58392 pulses and a=24 or 23, define

\[
 x=Ng,\quad e_a=\frac{x}{1-x}+\frac aQ,\quad
 f=\frac{Ng^2}{(1-x)^2},\quad
 \delta_a=\frac{e_a(2+e_a)+f}{2},\quad
 G=\frac1{(1-x)^2}+f.                                         \tag{10}
\]

Telescoping its stored good-code products gives operator error e_a. Summing
all first processor failures gives mass at most f. The complete address
export is an exact completely positive trace-preserving map under its
declared semantics. Therefore delta_a bounds half the diamond norm of the
whole word instrument error, including padding and processor failures, and
G bounds the actual complete map's gain. Ideal processor-failure blocks are
zero; they occupy the same complete alphabet in the comparison.

Extend each word map as the identity on all earlier terminal output flags.
Its gain is then at most G>=1. If word j has error delta_j and gain G_j,
successive replacement by ideal channels yields

\[
 \delta_{\rm menu}\le
 \sum_{j=1}^{2m}\delta_j\prod_{k<j}G_k,
 \qquad G_{\rm menu}\le\prod_{j=1}^{2m}G_j.                    \tag{11}
\]

This follows from the triangle inequality, complete bounded gain of each
actual prefix, and contraction by each ideal suffix. It applies to arbitrary
entangled inputs, with no independent-source preparation assumption. For
m<=8 the exact rational upper bound is below 1/50000. Rounded maps need not
be exactly trace preserving. No postselection or survivor normalization is
used to obtain (11). A score bounded in absolute value by one changes by at
most twice this half-diamond error. Coarse-graining the observed flags cannot
increase it.

## Operational theorem boundary

This constructs a bounded family of complete positive-operator-valued
measurements and explicit instruments on the same source. It supplies common
noncommuting coexistent effects beyond a diagonal context. It does not
compile every supplied positive effect or every coexistent pair. The fixed
64-bit input and 41-bit angle domains do not give arbitrary precision.
`Lean/EventAlgebra/OperationalBornBridge.lean` has universal operational
premises; this finite family does not discharge those quantifiers or derive
Kraus trace semantics from the classical repair interpreter.

The carrier, work register, coherent services and preparation input are fixed.
No growing-capacity quantum embedding or compatible measurement refinement
is supplied. Any use on a growing source must produce those maps and preserve
the complete records, states, clocks and error bounds. No physical Hilbert
carrier, relativistic field, physical energy or particle identification is
inferred from this finite menu.
