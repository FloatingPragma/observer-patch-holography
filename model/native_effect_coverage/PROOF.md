# Native measurement and preparation coverage

The carrier is the original binary total-two sector
\(\mathcal H=\operatorname{span}\{|n\rangle:n\in\{0,1\}^{12},\sum n_i=2\}\),
of dimension66. The source weights are exactly
\(\theta=(25,25,8,\ldots,8)\). Every native step below is a fixed-seam E
command of `native_program/interpreter.py`, with its complete25-slot raw
instrument, one-slot FIFO and `export_oldest` rule. Coherent commands use
the declared nine-bit chart and proper-code services of
`observer/source_scalar/compiler.py`. Reading these specifications is
allowed. None of the constructions selects a field action or physical clock.

## A complete detector from native records

Use the eleven existing seams

```
(0,1) (0,5) (0,7) (0,10) (0,11) (1,8)
(1,9) (4,5) (6,7) (2,10) (3,8).
```

They form a spanning tree. An E record contains its coin, incoming contrast
\(d=n_a-n_b\) and current \(j\). Its actual update is
\(n_a\mapsto n_a-j\), \(n_b\mapsto n_b+j\). Start a displacement vector
\(s=0\). At each recorded step calculate

\[
 d^{\rm initial}_{ab}=d-s_a+s_b,
 \qquad s_a\leftarrow s_a-j,\quad s_b\leftarrow s_b+j. \tag{1}
\]

These are differences of the original source configuration. Set an offset
at vertex0 to zero and propagate the eleven differences along the tree.
If the resulting offsets are \(o_i\), total occupation fixes the common
constant:

\[
 n_i^{\rm initial}=\frac{2-\sum_k o_k}{12}+o_i,
 \qquad n_i^{\rm terminal}=n_i^{\rm initial}+s_i. \tag{2}
\]

The decoder reads only records. It rejects nonbinary or inconsistent
words. `run()` uses the initial state solely to verify the independent
record decoder against literal native execution.

Each primitive raw Kraus map is a positive weighted partial injection.
Composition preserves this property. Equations(1)–(2) show that a complete
word supports at most one input column. Consequently every nonzero fine map
has the exact form

\[
 K_h=\sqrt{w_h}\,|y_h\rangle\langle x_h|,
 \qquad w_h=\prod_{r=1}^{11}p_{e_r}(c_r). \tag{3}
\]

Every source column has one outcome for each coin at every step and the
two coin probabilities sum to one. Thus
\(\sum_{h:x_h=x}w_h=1\), and coarse-graining by decoded initial label
gives exactly the66 effects \(|x\rangle\langle x|\). It does not give
the Lüders update: the terminal source is \(y_h\). All zero native addresses
are retained. There are \(25^{11}\) formal raw words and
\(66\,2^{11}=135168\) nonzero rank-one words. The actual source kernel
check compares all1452 endpoint/coin columns with `observer.native.native`.

For any source/reference operator, the fine output is

\[
 (K_h\otimes I)\rho(K_h^*\otimes I)
 =w_h|y_h\rangle\langle y_h|\otimes\rho^R_{x_hx_h}. \tag{4}
\]

This matrix-unit identity establishes arbitrary-reference coverage. The
66-by66 rational product of native population kernels computes the exact
initial/terminal coarse channel; it does not discard the fine-word API.

The admitted program declares eleven attempts, ten exported records and
one FIFO record, with480 record-storage bits including terminal status.
The decoder can use64 signed16-bit working cells in addition to the
immutable source chart and retained archive. The generic expanded
`interpreter.execute` has a65536-branch host ceiling: it returns
`RESOURCE_STOP` after three steps with15625 blocks. That returned stop is
preserved and tested. `continue_raw` is a separate exact lazy query of one
complete unnormalized word, calling the existing continuation kernel at
each step with the actual FIFO/archive. The existing
`native_causal_dataflow.replay_trace` independently authenticates concrete
eleven-step runs. This factorization represents the complete mathematical
instrument without pretending that the expanded evaluator completed it.
It is not a normalized outcome sampler.

## Record-controlled restoration and preparation

Let \(c(x)\) be the existing nine-bit code of the source configuration.
After the native detector, both \(x_h\) and \(y_h\) are public record
functions. Apply the existing X word on the set bits of
\(c(y_h)\mathbin{\mathrm{xor}}c(z)\). Each X is the actual two-pulse
word \(I(-\pi/2),X(\pi/2)\). At most nine X blocks, or eighteen source
pulses, are needed. The identity branch invokes no empty service word.
All4356 source-address pairs are checked by the exact permutation code.
Intermediate work addresses may lie outside the source chart; the full512
work states are retained until the terminal export.

For the exact ideal service target, selecting \(z=x_h\) turns(3) into
\(\sqrt{w_h}|x_h\rangle\langle x_h|\). This supplies the full basis
Lüders instrument with every coin word retained. Selecting a fixed
\(z=a\) instead gives \(\sqrt{w_h}|a\rangle\langle x_h|\), hence

\[
 \rho_{SR}\longmapsto |a\rangle\langle a|\otimes\operatorname{Tr}_S\rho_{SR}
 \tag{5}
\]

after tracing the archive. With the archive retained, the prepared source
factors from archive and reference, while archive and reference can remain
correlated. No pure source input or postselection is required. The coherent
code, code gluing, work initialization and six-level processor blanks are
declared services. This is a bounded observer composition of native
commands and `source_word` services; the bare repair/hold grammar does not
issue coherent words. The generic `execute_word` also has a six-qubit host
ceiling. The existing source compiler, used here on the admitted nine-qubit
word, executes that word directly.

## Full finite measurement and preparation span

`native_interface/tomography.py` supplies4291 settings and the exact
4356-dimensional Hermitian effect design. For each pair \(a<b\), its
reversible X/CNOT route takes the two source codes to local addresses2,1;
the existing mixer or imaginary analyzer follows, then inverse routing.
The route is checked on all512 addresses for all2145 source pairs. The
maximum analyzer uses31 certified blocks and113 source pulses. This proof
consumes that coherent design, replacing its fair current-only detector
and2520-attempt timeout contract with(1)–(4) on the weighted E source.

If \(U\) denotes the full512 ideal analyzer, the successful fine map is
\(K_h C^*UC\), where \(C\) embeds the66 source codes, accompanied by
the complete terminal padding export. Some other source inputs can enter
padding even for an ideal analyzer. Its structural abort probability is
part of the ideal experiment, not numerical error. Only the two selected
output rows are used for the pair coordinate; no claim that other source
states are fixed is needed. These rows give the exact quadrature
observables

\[
 X_{ab}=|a\rangle\langle b|+|b\rangle\langle a|,
 \qquad Y_{ab}=-i|a\rangle\langle b|+i|b\rangle\langle a|. \tag{6}
\]

Together with the66 basis effects, they span every Hermitian matrix and
generate \(M_{66}(\mathbb C)\). The native record instrument remains
destructive unless the record-controlled restoration is requested.

After preparing \(|a\rangle\), execute the reverse-sign inverse
analyzer word. Its target pair is reducing, so the resulting ideal pure
states are \(|a\rangle\), \((|a\rangle+|b\rangle)/\sqrt2\) and
\((|a\rangle+i|b\rangle)/\sqrt2\). The inverse word is actually
executed and certified; an adjoint of a returned floating array is not
substituted for execution. The128 local projector entries are checked
exactly in \(\mathbb Q(\sqrt2)+i\mathbb Q(\sqrt2)\). These4356
preparations separate all Hermitian effects: their expectations determine
the diagonal entries, then the real and imaginary entries using
\((E_{aa}+E_{bb})/2+\Re E_{ab}\) and
\((E_{aa}+E_{bb})/2-\Im E_{ab}\). Equality on this family is therefore
exact effect equality within the supplied matrix semantics. If each such
expectation differs by at most \(\eta\), each off-diagonal real or
imaginary error is at most \(2\eta\), so

\[
 \|E-F\|_F\le\sqrt{34386}\,\eta<186\eta. \tag{7}
\]

## Complete service errors and stable actual span

Every original returned small-block array is checked on its entire2,4 or8
dimensional carrier against its exact target, with fixed operator gate
\(10^{-11}\). Both inverse analyzer words are executed independently.
The serialized gate is fixed; floating residual digits are not receipt
data. Every primitive's full six-level complement and code gain are
checked by `modular_clock._processor_bounds` against \(g=2\,10^{-11}\).

For local binary64 multiplication let \(u=2^{-53}\) and
\(\gamma_{32}=32u/(1-32u)\). Dot lengths are at most four. A local
gate's entrywise absolute norm is bounded by its Frobenius norm
\(2(1+g)\). The full-eight-address running block has Frobenius norm
at most4: its initial norm is \(\sqrt8\), and the checked envelope
\(35[g+8\gamma_{32}(1+g)]+2^{-900}<1/4\) closes this bound.
Each multiplication therefore contributes at most
\(8\gamma_{32}(1+g)\) in Frobenius norm; the subsequent propagation
gain is below2. For at most35 pulses the exact product of stored primitive
blocks differs from its returned array by at most
\(16\cdot35\gamma_{32}(1+g)+2^{-900}<10^{-11}\).
Thus every exact stored-pulse block has the uniform target error
\(b=2\,10^{-11}\). This covers ordinary IEEE multiplication with the
stated tiny-underflow allowance; analytic matrix-exponential and hardware
accuracy are outside the certificate.

For a word containing \(B\) certified blocks and \(N\) pulses, set
\(e=(1+b)^B-1\). The owning processor theorem bounds the sum of all
first-failure output masses by

\[
 f_N=\frac{N g^2(1+\texttt{NATIVE\_GATE})^6}{(1-Ng)^2}.
\]

Retaining every first failure and padding export gives the complete
half-diamond bound \(\delta=[e(2+e)+f_N]/2\) and map gain
\((1+e)^2+f_N\), with arbitrary untouched references. Conditional reset
uses the maximum bound over branches. Analyzer, reset and full coherent
preparation have respectively \((B,N)=(31,113),(9,18),(40,131)\), giving
\(\delta<6.21\,10^{-10},1.81\,10^{-10},8.01\,10^{-10}\).
These are complete-output comparisons, not a normalization of numerical
success blocks. Literal array propagation is a diagnostic with separately
stated rounding slack. It retains structural ideal padding abort.

The actual stored-primitive menu has full rank as well. In orthonormal
real Hermitian coordinates, the ideal analysis map given by the diagonal
and(6) has smallest singular value at least1. Each of its4290 compiled
quadrature rows changes by at most \(2\sqrt{66}\delta\) in Euclidean
norm. The total perturbation is below \(1066\delta<1\); hence the
actual complete-record linear functionals also span4356 dimensions. The
rounded maps have the stated gain and need not be exactly trace preserving;
their numerical positive functionals are not asserted to be normalized
positive-operator-valued measurements. Likewise, for
one fixed identical source input, each actual preparation's successful
source block is within complete half trace distance \(\delta_p\) of
its ideal projector. Failures remain recorded and are scored zero. If a
Hermitian \(D\) is annihilated by all actual preparation tests, all ideal
means are bounded by \(2\delta_p\|D\|_F\). Equation(7) gives
\(\|D\|_F\le372\delta_p\|D\|_F\), forcing \(D=0\).
These successful blocks have bounded normalization error; no exact pure replacement
is asserted for rounded processor services.

## Operational scope and the finite Born premises

The finite ideal menu consisting of one pair analyzer and basis readout,
with arbitrary classical coarse-grainings and mixtures, does not realize
every positive effect. Its rank-one atoms have either one coordinate or
two equal-modulus coordinates with relative phase \(\pm1,\pm i\).
Padding abort is diagonal. A positive sum equal to a rank-one operator
can use only summands with that same one-dimensional range. Therefore the
projector onto \((3|a\rangle+4|b\rangle)/5\), with pair matrix
\(\frac1{25}\begin{pmatrix}9&12\\12&16\end{pmatrix}\), is outside
this positive menu. This statement is about the specified ideal menu;
it does not rule out other sequences of admitted coherent services.
Linear span and algebra generation do not provide universal positive
effect realization or common refinements of all coexistent effects.

There is an exact positive refinement on the diagonal algebra. After the
tree detector, append at most53 E repairs on seam(0,1). Its two source
weights are25, so every placement bit is an independent fair bit even
conditional on the preceding source state and records. With \(m\) bits,
the integer \(U\) is uniform on \(0,\ldots,2^m-1\). For diagonal
dyadic effects \(f_x,g_x\ge0\), \(f_x+g_x\le1\), accept the disjoint
intervals of lengths \(2^mf_x\) and \(2^mg_x\) after decoding \(x\).
This produces their joint refinement by actual record counting. Every
coin, current and rejected interval is retained. The suffix changes the
source; no nondemolition suffix claim is made. Its scope is the diagonal
context with at most64 native attempts, not arbitrary real effects or a
coherent abort-reset context.

The bindings to `OperationalBornBridge.lean`, `FiniteBuschGleason.lean`
and `OperationalAdditivityBoundary.lean` distinguish complete record
weights from universal effect realizability and coexistence. This package
supplies a complete finite instrument menu, operationally separating
preparations and diagonal refinements. It consumes the runtime's Kraus
trace probability semantics. It does not independently derive that
semantics from counting on all complex effects or discharge the universal
Born premises. The nine-bit services have no proved growing-capacity
embedding here, and these results identify no physical position-momentum
commutator, relativistic field or measured particle.
