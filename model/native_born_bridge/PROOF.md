# Exact native record measure and its classical effect closure

## Source and probability inputs

Let X be the 78 nonnegative integer configurations of total load two on the
twelve source ports. The source function `model.observer.native.native`
repairs a selected ordered seam (a,b) using a placement bit c. If s=x_a+x_b,
the new endpoints are (floor(s/2),ceil(s/2)) for c=0 and the reverse for c=1.
The committed record includes c and j=x_a-x'_a. Our record also retains the
selected seam, an action already known to the public controller. Only the
thirty committed carrier edges are admitted. Other loads
are unchanged. Conservation makes X invariant. `endpoint` executes the source
and independently checks this floor/ceiling description.

Supply a probability law mu on X and independent fair placement bits. Neither
mu nor fair independence is derived here. Each specified input and bit tape
has a deterministic execution. Probability describes the supplied ensemble
over those tapes and inputs; it is not produced by declaring deterministic
execution random. The code accepts any exact rational mu, including arbitrary
zero entries. The proof works for every real probability law on finite X.

A finite policy chooses the next seam or an explicit halt from the committed
record prefix alone, before the next bit. It does not read hidden loads or
future bits. The executable snapshots the policy as a finite prefix table,
enumerates it on all 78 configurations independently of mu, rejects missing
rules and retains every leaf. The primitive has no abort branch; explicit
halts, all bit outcomes and zero-probability records remain present. The
analytic result covers arbitrary finite horizons; the executable has at most
ten attempts and 100,000 enumerated leaves. These are resource limits, not
physical assumptions or an infinite-history claim.

## Complete execution, coarse graining and continuation

For a leaf h started at x, let k(h) be its number of actual attempts, r(h) its
record and y(h) its terminal configuration. Its conditional weight is exactly
2^(-k(h)). At every internal node the two child weights sum to its weight;
an early halt retains the whole node weight. Induction on the remaining
horizon gives sum_{h from x} 2^(-k(h))=1. Thus

    P_mu(r,y) = sum_{x,h: r(h)=r,y(h)=y} mu(x) 2^(-k(h))

is nonnegative and normalized. Its public marginal p_mu(r)=sum_y P_mu(r,y)
comes from finite classical counting. No trace rule is used. Reversing the
record restores x by adding j to the first endpoint and subtracting it from
the second, in reverse order. Therefore terminal (y,r) identifies its initial
configuration and bit history; the executable independently checks this
custody property. The record marginal can still forget that initial state.

For a public event A, define

    e_A(x) = sum_{h from x, r(h) in A} 2^(-k(h)).

Then 0<=e_A<=1 and P_mu(A)=sum_x mu(x)e_A(x). Disjoint union gives
e_(A union B)=e_A+e_B and P_mu(A union B)=P_mu(A)+P_mu(B).
The sure event pulls back to the constant function one. This is an operational
coarse-graining result from complete record sets, rather than matrix-effect
additivity inferred from a displayed probability formula.

If an extended policy preserves every earlier decision and every explicit
halt, summing all future descendants of an earlier leaf returns that leaf's
weight. Hence pushforward along prefix restriction gives exactly the earlier
record law. `refinement_certificate` checks this equality separately at all
78 initial configurations, so it holds for every initial law. Replacing mu(x)
by any joint classical law mu(x,z) proves the same statement with an arbitrary
untouched classical reference z and its conditional distributions. This is
not a quantum-reference or physical no-signalling theorem.

## Unique diagonal density, without choosing coherences

Let R be the complete finite set of public outcomes. The native weights give
the functional L(f)=sum_r p_mu(r) f(r) on the commutative algebra C^R. It is
complex linear, positive on nonnegative real functions, and normalized.
Conversely, any such functional is fixed by its values on the point indicators:
f=sum_r f(r)1_r. Thus the observed atomic weights determine L uniquely.
Under the canonical diagonal embedding of C^R in matrices this functional is

    L(f) = Tr(diag(p_mu) diag(f)).

This trace expression is a consequence of the finite sum, not its probability
postulate. It supplies the unique DIAGONAL density. It does not select a unique
state of the full complex matrix algebra. For example, the two positive
rank-one matrices [[1/2,1/2],[1/2,1/2]] and
[[1/2,-1/2],[-1/2,1/2]] agree on every diagonal observable. The same argument
gives the diagonal representation of the joint (record,terminal configuration)
law and of the supplied input law; it does not make hidden configurations
publicly accessible or infer mu from an insufficient record experiment.

## Native implementation of dyadic effects and real extension

Fix a nonnegative integer m and integers 0<=k_r<=2^m for every original
record. Append m actual native attempts on a fixed seam. The new recorded
placement bits form the uniform integer U in {0,...,2^m-1}, independently of
the original record. Their native currents and terminal loads may change;
this does not alter the retained earlier record or the bit law. Accept iff
U<k_r. All reject outcomes remain. The native acceptance probability is
exactly sum_r p_mu(r) k_r/2^m. Thus every dyadic classical effect has an actual
finite native record-event implementation, without a separately supplied
biased coin, amplitude or detector.

For any real f:R->[0,1], the dyadic lower/upper approximants

    f_m^- = floor(2^m f)/2^m,   f_m^+ = ceil(2^m f)/2^m

are native-implementable and satisfy f_m^-<=f<=f_m^+ and
||f_m^+-f_m^-||_infinity<=2^-m. Their finite-sum weights have the same bound.
They therefore determine exactly one order-preserving extension to all real
classical effects, namely L(f). No additional continuity axiom is necessary.
This is a limit characterization; no finite program exactly implements an
arbitrary irrational effect, and the bounded executable does not execute
arbitrarily large m. It executes dyadic cases and certifies exact rational
instances of the order enclosure.

If f+g<=1, L(f+g)=L(f)+L(g). For dyadic f and g with a common denominator,
one suffix integer realizes them as disjoint intervals [0,k_f) and
[k_f,k_f+k_g); their union implements f+g. Equal classical effects have equal
weights in every initial law, independently of the chosen acceptance
implementation. These are actual classical instances of effect equivalence
and coexistent-effect realization. They do not assert those properties for
noncommuting matrix effects.

## Relation to the existing quantum and Born results

For either existing quantum repair lift C or E, a fixed recorded branch maps
each allowed basis vector to its repaired basis vector with amplitude
1/sqrt(2). The record makes the map a partial injection. On a diagonal input,
K diag(mu) K* therefore moves each allowed diagonal entry with weight 1/2;
there are no interference terms. Summing the complete branches reproduces the
classical endpoint law. Composing along a public adaptive policy proves the
same equality for every complete finite record history. This proves equality
of the diagonal restrictions; it does not select one coherent lift. The
executable independently compares all 78 input columns and 4,680 source
seam/coin branches with the existing exact C/E diagonal blocks, after deriving
the classical weights.

Full Born representation is proved by
`EventAlgebra.finite_busch_gleason` and `isEffectValuation_iff_born` in
`Lean/EventAlgebra/FiniteBuschGleason.lean`: normalized nonnegative valuations
additive on every coexistent complex matrix-effect sum have a unique density
representation. Real affinity follows in that theorem; continuity is not an
extra premise. The native argument here instantiates the classical diagonal
part of that operational architecture, rather than deriving the whole complex
effect algebra from scalar repair.

The existing `model/native_interface/tomography.py` is real executable positive
content: compiled real and imaginary analyzers give all 4,356 Hermitian
coordinates of the 66-state binary sector, with export abort and timeout
outcomes retained. This packet rechecks its exact design and all 156 source
pins. It does not rerun the entire compiler/quantum execution. Coherent
proper-code services, quantum state/channel semantics and repeated preparations
remain premises of that tomography result.

Tomographic completeness identifies a state within a represented matrix model;
it does not prove that arbitrary operational assignments on all effects are
matrix states. We reuse the kernel-checked countermodel in
`Lean/EventAlgebra/OperationalAdditivityBoundary.lean`:
`producedCubicValuation_not_additive` and `producedCubicValuation_not_born`.
On its committed tomographically complete qubit fixture, a normalized positive
complement-respecting valuation matches every displayed value yet fails a
coexistent sum, with 35/64 versus 143/256, a defect of -3/256. This package
rechecks that rational arithmetic and binds the existing proof source; it
does not present the countermodel as a native law or recompile Lean.

The remaining complex bridge is operational: supply a source implementation
family covering the required complex effects, prove that equal represented
effects have equal record probabilities across implementations, and realize
coexistent sums by compatible disjoint event unions. The complete fair-coin
tree proves those clauses for its classical record effects. Neither the
classical theorem nor the conditional complex tomography proves them for all
complex effects, selects the quantum representation, or establishes physical
laboratory identification.
