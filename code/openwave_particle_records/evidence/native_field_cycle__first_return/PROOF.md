# A complete coherent native first-return operation

Let U be the signed permutation obtained by the three actual native seams
on the full486-state domain. The invalid flag is fixed. Its unsigned
permutation is a bijection, and every orbit has length at most486. A tag
fiber is the set of source addresses with one fixed tagged position.

For valid s, define r+(s) as the least positive t such that U^t s has the
same tag as s, using the unsigned source address for this test. Define r-
with U*. Both lengths are at most486. Define F(s)=U^r+(s)s, including the
product of actual native CAR signs along the word. On each orbit of U,
list the addresses in a fixed tag fiber in cyclic order. F takes each
listed address to the next one. Thus F is a signed permutation preserving
every tag fiber, and backward first return is its exact inverse. In
particular

    r-(F(s))=r+(s).

This equality is a source-permutation fact, not a variable-time instrument
claim. The coherent implementation below supplies the missing instrument.

## Fixed-horizon return-length oracles

The computational address has487 values, including the invalid flag, and
the length register has512. Initialize a scratch source address to0 and
486 scratch hit bits to0. Modular addition copies the orthogonal source
address into scratch. On a superposition this creates an entangled pair
of address registers; it does not copy an arbitrary unknown quantum state.

Execute486 U words on scratch, recording at each step whether its tag
equals the unchanged input's tag. Invalid input records false throughout.
A reversible priority encoder XORs the index of the first true bit into
the length target, with index0 for invalid. Then reverse the486-word
computation, erasing each equality bit before its inverse word, and undo
the modular addition. The scratch source, every hit bit and every phase
are restored. Native phases from the work-source scan cancel because the
reverse word is the exact adjoint. The target length is the only output.

Every constituent operation is a permutation or signed permutation on
its full register domain. The encoder is of the form target XOR f(bits),
so it also acts reversibly on dirty target registers. On prepared scratch0
the resulting C+ acts as |s,c> -> |s,c XOR r+(s)>. Using U* in the scan
gives C-. The code stores all fixed-horizon traces and checks the complete
forward/reverse scratch action on every source address. It separately
checks the direct first-return path against the oracle result. Arithmetic
and comparators are declared coherent services; no full pulse compiler is
claimed for those gates.

## The complete field operation

V acts on every source and counter address by |s,c> -> |U^c s,c>, with
the signed source action. It has a fixed511-slot program: in slot t apply
U if c>=t, and identity otherwise. It is unitary for every c, including
the nonphysical dirty-counter complement. The composition

    C- V C+

maps prepared counter0 as

    |s,0> -> |s,r+(s)> -> |F(s),r+(s)> -> |F(s),0>.

The last equality uses r-(F(s))=r+(s). Every ancilla is clean at the end,
so coherences between different return lengths and arbitrary untouched
references are preserved. No return-time label is discarded or measured.
The inverse is C+ V* C-, and every one of487*512 source/counter columns
is checked with its sign. Bijectivity on the complete domain establishes
the operator identity on every matrix unit. The original source domain
and the full dirty-counter complement have explicit actions.

Each oracle uses972 native loop-word slots. The middle uses511. The common
macroinstruction therefore has2455 native loop-word slots, or7365 native
seam slots, plus fixed reversible arithmetic. A nonselected gate is the
identity, not a shorter physical time for that amplitude. The clock is
declared in these macroinstructions; conversion to laboratory duration or
another reconstructed source clock is not supplied.

## Boundary completion and charge continuity

Every component of each word uses the existing whole-event boundary hold.
The first-return construction has no assumed electric translation at the
boundary. Since F preserves the final tag and compensator, it preserves
the charge density. Gauss preservation on every microstep gives

    divergence(E_after-E_before)=0

on every source address. The integrated native charged microcurrent along
the complete return path equals minus this electric increment. That path
may have a different number of actual moving hops on different inputs;
the common coherent schedule does not erase its recorded ancestry.

The first-return source operation therefore changes field/memory records
while commuting with every final charge-position pointer. This is an exact
macrolevel field-cycle operation on the same source. Its implementation
contains charged intermediate native routes. A physical field interpretation
must state which intermediate interventions are accessible at the claimed
clock scale; this theorem does not remove those charged routes from the
full source specification.

For q=0, F fixes every electric link and final tag. Native position updates
and boundary decisions are independent of initial memory, so the whole
first-return route and its length are the same on all six memory inputs
with fixed tag and flux. The active-memory update is a fixed translation
modulo6 on that route. The CAR crossing signs telescope through the
signed tagged creation map; return to the same tagged pair gives product+1.
Consequently F acts as identity on the uniform-memory vector of each of
the81 tag/flux fibers. This is an exact81-dimensional subspace statement,
tensored with both spins, rather than a zero-charge expectation statement.
The current queries also preserve this subspace because their memory
actions are translations. Hence the complete q=0 reference law is
identical with the loop omitted on this specified preparation domain.
The original frozen requirement of a nonzero difference therefore fails;
the explicit amended control checks this null equality and retains that
failed original gate. The q=1 positive discrimination gates are unchanged.

## Coherence control

For a normalized source preparation psi, let p_r be its weight on inputs
of return length r. Without the final C-, the joint output is
sum_s psi_s |F(s),r+(s)>. The pointer projecting the source onto ideal
F psi, with all counter outcomes retained, has success probability
sum_r p_r². With the full uncomputation, it has probability1 in counter0.
The source-pointer complement is also retained. The control therefore
detects the exact loss that would occur if orthogonal history records
were left in the environment. A two-length equal-weight plant gives1/2.

This calculation checks coherence of the implemented first-return source
map. The equality of basis maps alone would not establish it for an
instrument with an unobserved retained history. The code evaluates both
complete instruments on the fixed uniform Gauss-zero preparation.

## Common stationary clock

Insert F into the established loop dilation [[0,F*],[F,0]] while keeping
the same current queries, reflection and signal preparation. SELECT and
W=SELECT R remain complete unitaries. The signal compression is
(4I+sum_e sigma.n_e J_e+(F+F*)/2)/8. Its positivity follows from the same
operator-norm bound. The lazy loop amplitude (6I+F+F*)/8 lies in[1/2,1].
Neither amplitude kernel is silently substituted for W as the clock.

For each electric link the exact stationary identity is

    <selected realized charge current>
      = <first-return electric increment after reflection>.

The nonreturn term vanishes on every column, rather than merely in a
selected expectation. Both sides have zero divergence on stationary
states. The observer can measure them on fresh identically prepared
states with the complete query reference kept until the terminal read.
All Gauss0/uniform-memory signal bands and both walk phases are retained,
including spin-blind and omitted-loop controls.

This supplies a finite coherent, boundary-complete field-cycle source and
its same-clock current balance. It does not supply a spatial magnetic
dipole, an orbital comparison fixing g, a magnetic Larmor intervention,
or an electron interpretation of the neutral CAR pair. Those observables
are not replaced by a commutator definition of B or by assigning D²/(2m)
as an energy.
