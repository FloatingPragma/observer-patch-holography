# Native local observers for the interacting spin limit

The protected single-occupancy spin algebra has a direct native gate and
measurement implementation. This package executes finite preparations,
Heisenberg exchange, Pauli probes and the complete recorded C/E field repair.
It also supplies an explicit finite-resource approximation to local observer
histories in the infinite, fixed-spacing spin dynamics of
[FA90](../protected_spin_limit/README.md). The construction uses supplied
local input marginals, or a specifically prepared Néel product. It never
queries or prepares an unspecified infinite ground state.

There are two distinct precision statements. The finite examples execute
the existing binary64-controlled six-level source and certify its actual
literal coefficients. Arbitrary accuracy of the resource construction uses
a new, explicit **dyadic integrated-port command algorithm** for that same
mathematical source law. It requires the declared controller to supply those
coefficients at the requested accuracy. A fixed binary64 interface is not
asserted to implement an arbitrarily accurate sequence.

## The exact protected encoding and exchange gate

There are two occupation bits per site, with up first. Encode a logical bit
as `0 -> |down>` and `1 -> |up>`. On the full four-dimensional onsite space put

```
Xbar=Xup Xdown,       Ybar=Yup Xdown,       Zbar=Zup.
```

These obey the Pauli algebra and each commutes with the single-occupancy
projector. In this basis `(Xbar,Ybar,Zbar)=(2Sx,-2Sy,-2Sz)` on the protected
space: the same fixed spin-frame rotation is made at every site. Heisenberg
exchange is unchanged by this rotation. Starting with a logical bit on the
up register and a blank down register, `Xdown` followed by `CNOT(up,down)`
is the encoding isometry. The phase-exact X gate uses two native pulses;
CNOT uses the existing three-pulse word. Encoding costs five pulses per site.
A specified Néel product costs at most seven per site, from blank registers.

For two sites the extended bond is

```
hxy=(J/4)(Xbar_x Xbar_y+Ybar_x Ybar_y+Zbar_x Zbar_y-I).
```

The three nonconstant Pauli products commute, also on the full occupation
carrier. Their restrictions give precisely `J(Sx.Sy-1/4)` of FA90. Thus
`exp(-it hxy)` is exactly three Pauli rotations and its scalar phase. The
source compiler uses 27, 27, 7 and 1 pulses respectively, for **62 pulses**.
No large-U virtual process is used to implement this gate. An arbitrary
single-spin X/Y rotation uses eleven pulses; Z uses one. A local Pauli probe
rotates the chosen encoded axis, copies the up bit to a blank pointer, and
undoes the basis change: three pulses for Z and 25 for X or Y. Every probe
has both outcomes and the reserved six-level pointer slots.

The ideal encoded gates preserve single occupancy. Intermediate pulses in
their decomposition need not do so. Complete field repairs are scheduled
**between completed gate words**, rather than assumed harmless during each
individual pulse. Controlled execution and that schedule are apparatus.

## Complete field repair, including its off-code action

We use the full-domain two-flavor C/E map in `defect_field_bridge/full_domain.py`.
For each fair coin and every local occupation mask, the exact rule supplies
the output mask, raw record and CAR sign. Collisions are held with raw 8;
they are not removed. The 32 input/coin assignments on adjacent two-site
fields inject into the nine-bit field/coin/four-record register. Complete
them deterministically to a permutation of all 512 addresses.

Each permutation transposition is decomposed along a Gray path. A
Gray-neighbor swap is a nine-bit controlled flip. Its source word is the
complete 512-term commuting parity expansion of a rank-one reflection,
including the scalar term, the final Gray-cycle CNOTs and the two target
basis rotations. There are 2044 pulses before fixed-pattern conjugations.
All padding addresses are part of the permutation. With the fair-coin pulse,
the adjacent C and E dilations cost **69,635** and **94,207** pulses. The code
streams their actual six-level pulse matrices; it does not insert an ideal
permutation in place of the native word.

For separated endpoints in a CAR ordering, ordinary wire swaps would lose
the intervening fermion sign. Move the two endpoint modes past the intervening
sites using adjacent **fermionic SWAP = SWAP times CZ**, perform the adjacent
repair, and reverse the routing. Each fermionic SWAP has eleven native
pulses, including the scalar phase because `J(-pi)=-CZ`. A two-mode block
passes another two-mode block in four swaps. Forward and reverse routing
therefore cost 88 pulses per intervening site. Conjugating the canonical
CAR modes proves this is the original signed, full-domain repair for every
input and reference, including number-changing superpositions. The crossing
execution explicitly includes both even and odd intervening occupation.

On single occupancy, the exact result is
`M_(coin,raw) P = delta_(raw,0) P/sqrt(2)`. Thus ideal repairs contribute
independent fair coin records and do not disturb the encoded spin. Their
complete off-code action is still required: native gate errors are not
postselected back into P. All ideal/native comparisons are comparisons of
complete dilations, so the usual telescoping bound covers these errors and
their later repair outcomes.

## Source-read precision without a hidden binary64 limit

The retained source has twelve anti-Hermitian generators `D_p` with exact
entries in `Q(sqrt(5))`. Its independently checked identities express each
X/Y/Z/I/J target generator with coefficients in that field; the H target
uses the checked `(X+Z)/sqrt(2)` combination. The exact absolute-row norm
majorants in `precision.py` sum to

```
B=sum_p ||D_p||_upper=29.
```

Anti-Hermiticity makes the absolute-entry matrix symmetric, so its maximum
row sum also bounds the operator norm; no dimension factor is omitted.

For a rational requested angle, or rational multiple of pi, the algorithm
encloses the target integrated coefficients using integer square roots and
`pi=16 atan(1/5)-4 atan(1/239)`. The arctangent series has explicit signed
alternating remainder. Radicals have dyadic enclosures obtained by integer
square root. Interval multiplication and division by the positive sqrt(2)
enclosure preserve containment. Increase the internal precision until each
nearest p-bit dyadic is proved within `2^-p` of its target coefficient.
This terminates because all twelve target coefficients are computable real
numbers and interval widths tend to zero. No float decides termination.

If q is the emitted dyadic vector, anti-Hermitian Duhamel gives

```
||exp(sum q_p D_p)-exp(theta T)|| <= 29*2^-p.
```

Complete coherent load/export, including complement records, is an
isometry and preserves this bound. For N pulses it is at most `29 N 2^-p`
in complete-dilation operator norm, and **58 N 2^-p** in unhalved trace norm
on arbitrary entangled inputs. Choose the smallest p satisfying the desired
budget. The program emits the twelve actual rational commands and their
interval errors at any finite requested precision; the source law is unchanged.

The finite executions below instead use the existing binary64 pulse API.
Their coefficient errors are calculated from the actual emitted binary64
rationals, the exact source matrices, and rational enclosures for the desired
angles. These fixed errors are not sent to zero. Implementing the scalable
family requires the finer integrated-port interface, coherent memory and
load/export services at their specified precision. No physical noise floor,
hardware latency or fault-tolerance theorem is inferred from the ideal source
specification.

## A complete local-history cutoff, uniform in references

Use FA90's infinite cubic spin dynamics, with fixed lattice spacing and
bond norm J. Let X contain s sites in the radius-r cube. A finite protocol
has m interventions on X plus finite observer memory, at fixed times between
0 and T. Measurements, classically adaptive operation choices and their
records are represented by local Stinespring unitaries of norm one with fresh
ancillas. No external ancilla couples to the field away from these local
interventions. Final observables of norm at most one can include X and all
observer records and an arbitrary spectator reference.

The constructive calculator takes J, T and the intervention times as exact
rationals. More general computable inputs need their own certified input
intervals. Ideal field repairs act as identity on spin and only append
independent fair-coin records, so these particular interventions have no
additional spin support; their complete native off-code errors are charged
to the global finite-circuit budget below.

FA90's commutator-chain proof applies unchanged after tensoring spectator
systems: the first chain edge has at most 6s choices; no Hamiltonian bond
touches an ancilla. Therefore its cutoff estimate is completely uniform in
the reference dimension. For `R>r`, write

```
B_R(T)=min(2,6 J s (2R+1)^2 T 2^(-(R-r)) (exp(48JT)-1)).
```

In the interaction picture the protocol is the product of m evolved local
unitaries. Replacing each `alpha_(t_j)(V_j)` by its radius-R counterpart
costs at most B_R(T). A final evolved local observable costs one more such
term. Telescoping the product and its adjoint gives

```
||rho_observer,infinite-rho_observer,R||_1 <= (2m+1) B_R(T).       (1)
```

This is unhalved trace norm, clipped at two. The initial states are the same
on the compared field/reference algebra; the finite experiment only needs
their marginal on the R cube plus its reference. An approximate supplied
marginal adds its trace-norm error. Taking the supremum over final
reference-valued observables proves (1); no false Markov assumption for the
reduced field evolution is used. The finite input may be entangled with the
outside. A specified infinite Néel product has an explicit native finite
marginal, while an unspecified ground marginal has not been supplied by this
argument. Variable stopping times or an unbounded number of interventions
need a separate truncation; (1) concerns the stated finite protocol.

## Finite-box splitting and a constructive resource algorithm

The radius-R cube has `V=(2R+1)^3` sites and `E=3(2R)(2R+1)^2` bonds.
Disjoint bonds commute. The number A of unordered overlapping bond pairs is
at most `15V`. The first-order product formula for one interval of length d
has operator error at most `A J^2 d^2`: the general one-step bound is
`(d^2/2) sum_(e<f)||[h_e,h_f]||`, and each commutator is at most `2J^2`.
Splitting every free interval into steps of size at most delta and summing
their durations gives complete-channel trace error

```
2 A J^2 T delta.                                             (2)
```

Local interventions and complete repairs are contractions in trace norm,
so this estimate telescopes through all records. Each step uses E exact
62-pulse exchange words. Choose R making (1) at most epsilon/3, a step count
making (2) at most epsilon/3, and p making `58 N 2^-p<=epsilon/3`. Add the
native counts for the specified preparations, repairs, Pauli probes and
rotations. Their sum N is finite. `bounds.resource_plan` carries out those
integer/rational choices, using `exp(x)<=2^ceil(2x)` for an explicit cutoff
without a fixed-time limitation.

The returned pulse count is for the specified fixed schedule of those
words. The general history estimate also permits adaptive interventions,
but their actual controlled-choice or controller implementation costs
must be included in N; the fixed-schedule example does not budget arbitrary
adaptive Stinespring synthesis.

For example, J=1, T=1/32, one repair and one Pauli probe on a one-site support,
with requested trace error 0.1, gives R=19, V=59,319, at most
560,522,294,853 pulses and p=50. This conservative finite resource upper
bound is not a proposed practical experiment. Smaller exact finite examples
are executed separately. Taking epsilon to zero supplies a finite native
source program for each local-history accuracy under the refined control
interface. It does not prepare the infinite ground, identify relativistic
particles, or turn the fixed-spacing Heisenberg dynamics into a continuum QFT.

## Execution, complete outcomes and custody

The six executions are the two-site encoding Choi state, the exchange Choi
state, complete C and E repair Choi states on all sixteen field inputs, a
crossing-sign/reference test, and a Néel preparation/exchange/repair/Z-probe
history. Matrix columns are a retained quantum reference, not a classical
mixture substituted for it. All field and record rows survive execution.
Every native pulse is evaluated on the full six-level source, with all
proper-code complement outputs. First-failure weights carry their
unnormalized prefixes; the literal word and replay algorithm determine every
failure vector. There is no success renormalization or occupation projection.

The repair records retain all 32 coin/raw combinations. The physical pointer
also retains all reserved six-level outcomes; those outside the binary
encoding have zero ideal/final-code probability, while intermediate export
failures have their own absorbing labels. The observer-history experiment
retains all 64 coin/raw/probe combinations. The lossless archive and literal
program digests authenticate the resulting states and failure tables.
Host-arithmetic deviations and tiny cancellation diagnostics are separate
from the rational source-error certificates.

Run `python -W error -m pytest -q model/native_spin_observer`, then
`python -W error -m model.native_spin_observer --write` and `--check`.
Cold replay reconstructs every executed pulse. Current source bindings and
inherited receipt custody are reported separately; old receipts and their
historical status are not rewritten.
