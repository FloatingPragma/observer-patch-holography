# Native contact transport has an ordinary spatial action

## The microscopic lift fixes the contact phase

On blank memory a native seam has generator

\[
h=(I-S)/2+\kappa D^2,
\]

where S exchanges its ternary loads and D is their signed difference.
The two terms commute. Since D² has integer eigenvalues 0, 1 and 4,
at κ=8

\[
e^{-i\pi h}=e^{-i8\pi D^2}e^{-i\pi(I-S)/2}=S.
\]

The blank code is reducing in the full 99-state interaction. The identity
therefore concerns the actual native contact with its collar, rather than
a chosen phase convention for a reduced swap. All nine load columns are
independently checked against the full 99-state exponential. Isolating this
contact parks every other interaction and reference term during the pulse.

A product of these contacts is exactly the ordinary occupation-space
permutation U_p associated with the product p of endpoint transpositions.
Thus U_p U_q=U_pq and U_identity=I. Any two contact words implementing the
same full port map have the same complete action. This holds on arbitrary
ternary load superpositions with blank collars, including entanglement with
an untouched reference. It assumes no intermediate measurement, memory
write, extra phase pulse or additional Hamiltonian during transport.

The two-cell source uses independently prepared identical vectors ψ in
properly matched local frames. The desired exchange maps each of its
twelve support addresses onto the corresponding address of the other
cell. Its action on the physical occupation wave is

\[
U_{\rm exchange}(\psi_A\otimes\psi_B)
 =\psi_A\otimes\psi_B.
\]

No ordering sign is introduced when the indistinguishable loads are sorted
into a physical word. The microscopic gate identity supplies the lift.
This excludes a negative exchange sign for this source in this isolated
π-contact class, independently of the compiler's choice of paths.

## The selected negative character is spatial parity

The native guarded trimer has a simple energy line at 100.8894585203. Its
six arm-permuting graph symmetries commute with the guarded Hamiltonian.
The signed average K of these permutations is an exact projector. The
source's certified spectral gap and K overlap put the exact line in its
image. The finite timing-filter error is inherited from the owning
[source certificate](../../../../candidate_laws/interval_echo/openwave/row_13_baryon/arm_symmetry/README.md).

For these six physical cell symmetries the sign of the permutation of
the three arms equals the determinant of the spatial icosahedral map.
Three are proper rotations and have character +1. Three are improper
reflections and have character −1. Independent coordinate fits and the
complete proper-rotation census verify this identification.

Consequently, exchanging the two identically framed cores is positive,
while additionally reflecting one returned core has phase −1 on the exact
energy line. Finite preparation gives a small bright-port allowance. That negative
character is selected by native energy; it is not assigned to an internal
label. Its spatial operation is improper, so it supplies no sign for a
proper 2π rotation. A calibrated 72-degree proper template permutation
repeated five times has identity port map and exact positive lift by the
contact theorem above. Calibration against an observer-reconstructed
continuous spacetime rotation is a separate requirement.

## An address read distinguishes every proper template half-turn

The half-turn census uses the same four-filter source and the same blank
π contacts with the other interactions parked. Exactly 15 proper
icosahedral permutations have order two. A coordinate fit with determinant
one and trace −1 verifies all 15. For each turn and the idle setting, one
apparatus reads the three occupied addresses. An independent bank compares
the turn with idle on the two arms of the phase reader. For the native
contact the negative-port probability is
$\|(I-U)\psi\|^2/4$ on an admitted source. Each involution has
$U^2=I$, so repetition implements a full template turn with positive phase.
At changed transport contrast 7.5 the address operator is the same.
The additional phases on unequal contacts enter the phase-port banks.

A whole-state return, up to either phase, implies equality of every
address frequency with the idle frequency. Each terminal address bin
has a binomial marginal conditional on its bank's source admissions.
Apply Hoeffding's two-sided radius $\sqrt{\log(2M/0.01)/(2n)}$ to each
frequency with $n$ admissions and $M=2\times16\times220$ bins, including
the sharp polar source. Enlarge each radius by the recorded allowance
for discretizing the sampling probabilities. For a given axis, subtract
both radii from the largest absolute difference of corresponding idle
and turned frequencies. The result lower-bounds total variation with
simultaneous 99% coverage, also when the distinguishing bin or axis is
chosen from these records. Every axis has a positive lower bound on both
sources. The comparison tests full address identification on this finite
template. A coarser reader may erase the residual directional information.
[Detector and statistics](half_turns.py) · [All address and phase counts](half_turns.json.gz)

## The public component audit covers the continuous native pulses

Let a terminal diagonal projector C ask whether the six occupied ports
form exactly two connected components of size three. It uses only the
observed ports and carrier adjacency, without material-origin labels.
There are 160 six-vertex simple graphs satisfying this predicate. The
independent test enumerates all 2^15 such graphs and checks the decoder.

For any real incoming load wave ψ and an isolated native contact, the
probability of any diagonal load predicate during that pulse satisfies

\[
p_C(t)=\cos^2(t/2)p_C(0)+\sin^2(t/2)p_C(\pi),
\qquad 0\le t\le\pi.
\]

Equal-load configurations are fixed. Each unequal-load pair evolves as
e^{-it(κ(a−b)²+1/2)}[cos(t/2)I+i sin(t/2)S] for the two endpoint loads a,b. The two incoming coefficients are real,
so their probability cross term vanishes. Summing the selected diagonal
entries gives the formula. The contrast phase cancels inside each pair.
After a complete native κ8 pulse the state is a real permutation of the
incoming real state, so the argument applies successively to the whole
word. The minimum of each ideal probability occurs at a pulse endpoint.

If each finite source differs from its real ideal line, up to a common
phase, by vector norm at most δ, the two-source vector differs by at most
2δ. Unitary transport preserves that bound. The difference of any event
probability is at most the pure-state trace distance, itself at most 2δ.
Here the inherited bound gives 2δ≤2.17824×10⁻⁸. The program uses the actual
finite source for each sampled endpoint probability, and the real ideal
line for this continuous-pulse audit.

Both exchanged-source routes have an ideal minimum component probability
0.07635598 at stage 36, against 0.99867133 at the source. The finite-source
allowance is negligible compared with this loss. The pulse compilation
therefore fails to keep two spatially separate connected cores along the
route. It cannot be promoted to a free or intact-core exchange merely
because its endpoint port map swaps the source supports.

## Changed contact law

At κ4 the same π contact is the same positive swap. At κ7.5 it exchanges
an unequal binary pair with phase +i, while equal pairs have phase one.
For a fixed sharp input its total phase is i raised to the number of
unequal-load contacts along that word. All positions are exchanged as
before, but different paths can give different phase interference. The
full ternary and collar calibration retains these phases. The κ7.5
comparison keeps the native κ8-prepared source and all routing settings.
It demonstrates a changed contact phase; it supplies no fermionic lift.

These statements do not cover adiabatic band transport, nonblank memory
dressing, intervening projections or an independently reconstructed
spacetime rotation. Such constructions require their own microscopic
transport and phase reference.
