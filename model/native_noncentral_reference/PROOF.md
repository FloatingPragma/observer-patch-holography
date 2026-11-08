# Exact faithful local references from retained repair records

## Native preparation with a known branch reference

Use the fair twelve-port cell, its thirty original internal seams and the
rooted eleven-edge E tree from
[the complete local preparation](../native_fair_local_preparation/PROOF.md).
Run that tree, then two ordered thirty-seam E sweeps. This is one admitted
word with 71 attempts and 71 literal nodes. Each attempt retains all 25
original raw addresses. One-slot export-oldest memory leaves 70 archived
records and one final record.

For every nonconstant tree word `h`, the actual branch operator is
`K_h=2^(-11/2)|y_h><x_h|`. The records determine both configurations and
`q=|x_h|=|y_h|` without an extra occupation query. For arbitrary input on the
cell and any exterior/reference system, its output is

```
2^-11 |y_h><y_h| tensor rho_R,(x_h,x_h).                         (1)
```

The source state in (1) is pure conditional on the actual tree flag; its
probability is not replaced by one. Zero words are retained. The exceptional
constant words act on
`span{|empty>,|full>}` without resolving that coherent two-dimensional code.
They are explicit inconclusive outputs and are not assigned a local total.

Write `A_q` for the owning exact integer population matrix of the following
sixty E attempts. Its entries count fair coin words, its row and column sums
are `2^60`, and its conditional source reference is

```
sigma_(q,y) = diag(A_q[:,y]/2^60).                              (2)
```

The certificate checks every entry and every column of all eleven
nonconstant sectors. Every entry is strictly positive, every column is
nonuniform, and every eigenvalue is at least `1/8192`. Thus (2) is an exact
faithful noncentral reference on its `binom(12,q)`-dimensional sector.
It is produced after a finite word, without an equilibrium approximation or
an imported Hamiltonian. This observation retains `h` and marginalizes the
sixty mixing-record values. The actual source archives those values; keeping
them reveals the final basis configuration instead of (2). No arbitrary
archive/source factorization is asserted.

A full nonconstant tree branch factors the local source from the untouched
reference as in (1). Its normalized reference density can be arbitrary,
including an entangled exterior or an earlier archive. This is conditional
on `h` and preserves its actual weight. It makes no claim about factorizing
all tree flags or the exceptional constant code.

## A compatible regional modular reference

The references (2) generally have exterior-context-dependent local ratios.
A second exact construction supplies a common regional modular action.
Its port weights are a declared observer protocol computed from the same
native rooted tree:

```
a_i = 1+depth_tree(i) = (1,2,3,4,3,2,3,2,3,3,2,2)_i,
W(x) = product_(i occupied in x) a_i,
Z_q = sum_(|x|=q) W(x),
L_q = min_(x,y: |x|=|y|=q) floor(A_q[x,y]/W(x)).                (3)
```

These integer choices introduce no fitted parameter or auxiliary random
coin. The observer's choice of this tree-weight protocol is explicit; the
source does not force it to be a unique vacuum or physical thermal state.
The finite source matrices give `L_q>=1` and

```
alpha_q = L_q Z_q/2^60 >= 1/12,       1<=q<=11.                 (4)
```

For fixed tree endpoint `y`, each sixty-coin word has probability `2^-60`
and a known endpoint `x`. Order the coin words lexicographically within that
endpoint. Accept the first `L_q W(x)` words and reject the rest. This is a
classical function of the actual retained tree and mixing records. The
records determine the coins and endpoint; the rule does not inspect a
hidden state or change the native transition law. Formula (3) guarantees
that there are enough real words for every endpoint. It follows exactly
that

```
Pr(accept and x | h) = L_q W(x)/2^60,
Pr(accept | h) = alpha_q,
sigma_q^accepted = diag(W(x)/Z_q).                             (5)
```

Acceptance does not depend on `y` within the same local sector. Equation
(5) follows from counting equally weighted original paths, not from a
normalized fit or a postulated rejection coin. The complete output retains
accept, reject, constant and zero outcomes with their actual probabilities.
For example `alpha_2=7771/65536`; rejection has probability `57765/65536`.
The unconditional acceptance probability is
`sum_(q=1)^11 alpha_q Pr_input(N_cell=q)`, bounded below by the actual
nonconstant-sector weight divided by twelve. It can be zero on a constant
input. There is no unconditional preparation claim for the accepted density.

On the comparison observation retaining `h` and the acceptance flag, each
accepted source/reference block is exactly

```
2^-11 alpha_q sigma_q^accepted tensor rho_R,(x_h,x_h).          (6)
```

The rejected source block is also retained, with coefficients
`(A_q[x,y]-L_q W(x))/2^60`; these are nonnegative and sum to `1-alpha_q`.
They can be combined with any subsequent complete instrument. Equation
(6) applies only while the sixty marginalized record values are not fed
back into that instrument. Conditioning the full archive instead would
prepare known basis states and gives a different reference.

The accepted density is faithful on each nonconstant sector, with

```
lambda_min(sigma_q^accepted) = min_x W(x)/Z_q >= 1/8192.        (7)
```

All eleven inequalities in (4) and (7), including every endpoint column in
(3), are exact integer/rational checks. They do not require increasing the
number of attempts with exterior capacity. A conditional flag can have zero
probability for a particular initial input, so no faithful joint density
on all flags and arbitrary reference states is claimed.

## Bounded record processing

The acceptance predicate does not enumerate `2^60` coin words. For a known
terminal `x`, compute backward suffix counts

```
C_60(z)=1_(z=x),
C_t(z)=C_(t+1)(f_(t,0)(z))+C_(t+1)(f_(t,1)(z)),               (8)
```

using the original native endpoint maps. Starting at `y`, each recorded
coin one contributes `C_(t+1)(f_(t,0)(z_t))` to the lexicographic rank.
Update `z_t` by the actual recorded coin and continue. This partitions all
smaller words by their first differing bit, proving that the computed rank
ranges from zero to `A_q[x,y]-1` exactly once. The inverse construction
selects zero or one by the same suffix counts. Boundary tests check the
last accepted and first rejected word against the unchanged native kernel.

At most `61*924=56364` integers suffice for a suffix table. Every count is
at most `2^60` and fits in 61 bits; no native attempts or random bits are
added. The implementation is a bounded classical observer/analyst readout
of the existing tape. It is not a native repair opcode or a compiled online
feedback instruction. Programming an apparatus to condition later quantum
controls on this flag has that additional implementation obligation.
The original source retains its FIFO, archive and expanded-evaluator stops.

## One modular action on local overlaps

Let `N_i` be occupation on port `i`. On the accepted global local-number
sector,

```
K_q=-log sigma_q^accepted = (log Z_q)I - sum_i(log a_i) N_i.    (9)
```

For a subset `B` of the twelve ports, use its natural number-preserving
algebra
`M_B=direct_sum_k B(H_(B,k)) tensor I_(Bc,q-k)`.
The generator's exterior terms commute with this algebra, so global modular
flow restricts to the onsite generator `-sum_(i in B)(log a_i)N_i`.
The actual regional partial trace has weights

```
R_B(x)=W_B(x) e_(q-|x|)(a_Bc)/Z_q,                            (10)
```

where `e_j` is the elementary symmetric sum over exterior port weights.
Its extra term in `-log R_B` is scalar in each local number sector.
Consequently its intrinsic modular group agrees exactly with the restricted
ambient group. Nested restrictions compose, overlap actions agree on every
shared observable, and disjoint number-preserving algebras commute by the
occupation tensor embedding. No extra boundary-context enlargement is
needed for (9). The source conditional reference (6) is common to these
regions; separate independently prepared regional states are not substituted.

For completeness, within sector `k` of `B`, the conditional exterior
reference has weights `W_Bc(z)/e_(q-k)(a_Bc)`. Taking its expectation gives
an explicit completely positive unital, bimodular, state-preserving
retraction onto `M_B`. Dephasing different local-number sectors precedes the
expectation. The accepted state is diagonal, so this dephasing preserves it.
These are algebraic regional maps; an implementation of each map by a
bounded observer is a separate operation claim.

A coherence that moves one occupation from native port zero to port one
has a reference eigenvalue ratio two, independent of every exterior
configuration. Thus its modular frequency has magnitude `log 2`. This
phase can be probed using the separately declared proper-code occupation
coupling and coherent services, with their own implementation bounds. No
compiled full modular-flow experiment is asserted here. Every occupation
observable commutes with (9): this flow gives local charge phases, not
spatial transport. Noncentrality alone supplies neither a Lorentz boost
nor geometric modular stress.

## Growing attachment and gravitational boundary

Only internal seams are used, so the same 71-command word attaches to any
admitted federation cell with the original record-width and storage
charges. The exact source-word controls are at evaluated source levels;
integer address admission at larger levels is distinct from expanding the
global carrier or its complete histories. The reference is indexed by the
actual nonconstant tree sector, and its acceptance rate remains explicit.
Neither the fair global density nor the distribution of local sectors is
selected by this protocol.

Equations (5)--(10) provide a source-prepared faithful noncentral local
reference with a consistent modular action on its natural regional net.
They bypass the centrality of a sector-uniform reference on that declared
observation. Native repairs do not preserve (9) as autonomous stationary
quantum time. The acceptance observation selects a prepared ensemble,
while a universal physical vacuum and its continuing preparation require
a different argument. No map from modular parameter to native attempt
count, proper time or energy units is supplied. The twelve-port construction
provides no physical small-ball tower, conserved Lorentz stress, responsive
area, fixed-proper-volume generalized stationarity or positive Newton
normalization. These geometric requirements are distinct from the finite
faithfulness, preparation and overlap results proved here.
