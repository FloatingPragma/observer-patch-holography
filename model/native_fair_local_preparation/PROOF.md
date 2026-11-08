# Local preparation by fair native repair

This construction uses the existing fair binary federation, a specified
root cell, the original fixed-seam E commands and the export-oldest memory
policy. It uses no coherent gate, reset, source encoding or hidden occupation
read. The input can be any state of that cell and an arbitrary untouched
reference, including the rest of a fixed-global-number source. Local number
is conserved during the controlled preparation word. Its distribution is an
input property, not a selected density.

## Complete word on every local sector

The cell has the twelve native ports and thirty internal seams, in their
literal source order. Run the eleven-edge spanning E tree from
[the all-number detector](../native_sector_readout/PROOF.md), followed by
b blocks. Each block runs two ordered sweeps through all thirty seams.
Every command has its original 25 raw addresses, including hold and zero
C addresses. The program has exactly 11+60b attempts and literal nodes.
There is no local-uniform selector or nested-loop primitive in this word.

The tree theorem gives, for every nonzero record h outside the constant
corner, K_h=2^(-11/2)|y_h><x_h|. Both configurations are decoded from the
records. Without a supplied local total, the only exceptional words have
K_h=2^(-11/2)P_c, where
P_c=|empty><empty|+|full><full|. The empty/full coherence is preserved.
These statements concern all 4096 local configurations and all references,
not a conditional total-two cut. The fine source/reference output is
K_h rho K_h* with the tree record retained.

Every later complete E map has a diagonal record effect and is a weighted
partial injection on the binary basis. On diagonal inputs its unread action
for seam e is exactly Q_e=(I+Swap_e)/2. This follows by evaluating both
actual native placement coins on all four endpoint inputs. These are
population matrices; they are not the complete coherent E instrument.

## Exact two-sweep contraction

For each q in {0,...,12}, let d_q=binomial(12,q) and let P_q be the product
of the sixty population maps, in chronological source order. The owning
integer computation starts from I, and for each seam replaces A by
A+Swap_e A. The result A_q satisfies P_q=A_q/2^60. Every entry is a
nonnegative integer, and both its row and column sums are 2^60.
At intermediate step j, each entry and each sum is at most 2^j. Thus signed
64-bit addition and summation are exact throughout; the final exponent is
60, below 63. Sums across the selected row minima are Python integers.
The certificate enumerates every q and every source configuration in that
sector. It checks

```
mu_q(y) = min_x A_q(y,x)/2^60,
beta_q = sum_y mu_q(y) >= 1/2.                                      (1)
```

All thirteen exact numerator matrices are bound by canonical little-endian
integer hashes. The minimum occurs at q=6 and is greater than 0.55;
no rounded value is used in the acceptance test. The complementary sectors
are checked directly rather than omitted. This is a finite exhaustive
source calculation. It asserts no unenumerated-graph minorization.

Subtract the common column mu_q from P_q. The remainder is nonnegative
and every column sums to 1-beta_q. Hence, for any two probability columns
u,v, their half-l1 distance contracts by at most 1-beta_q. The uniform
column tau_q is stationary because P_q is doubly stochastic. It follows
for every starting configuration y that

```
D(P_q^b delta_y, tau_q) <= (1-beta_q)^b <= 2^(-b).                   (2)
```

This is the usual finite minorization argument, with the exact native word
supplying its premise. It does not infer quantum convergence from population
irreducibility alone: the preceding complete E tree supplies the required
diagonalization, including its exceptional two-dimensional corner.

## Arbitrary references and retained records

Define the channel R on the cell and reference by

```
R(rho_AR) = P_c rho_AR P_c
          + sum_(q=1)^11 tau_q tensor sum_(x:|x|=q) rho_R,xx.       (3)
```

Here tau_q is the normalized diagonal identity in sector q. The first term
includes both constant populations and their off-diagonal reference blocks.
The remaining cross-configuration blocks are removed. Equation (3) is a
completely positive trace-preserving channel and an idempotent conditional
reference map. It preserves every local-number probability and is not a
replacement by a single state.

A stronger comparison retains every tree record. Let T be the complete
tree instrument and J=(R tensor id_archive) T. Let J_b be the full native
word with only the subsequent sixty-b mixing-record values traced. Then
for every joint input, even when R contains an entangled exterior,

```
D(J_b(rho_AR), J(rho_AR)) <= 2^(-b),
(1/2)||J_b-J||_diamond <= 2^(-b).                                  (4)
```

Proof: on each nonconstant tree branch the positive reference block is
rho_R,xx times its nonnegative branch weight. Apply (2), multiply by the
trace of that block, and sum. The weights sum to at most one. Every constant
branch is unchanged by all later E commands, including its off-diagonal
reference entries. Its error is exactly zero. This argument works for any
reference dimension and proves the complete channel bound. Tracing the
tree record too yields (3), because local number is conserved and the tree
kills precisely the stated cross-configuration terms.

The tensor entry identity `liftLeft_conj_apply` in
[NoSignalling.lean](../../Lean/ObserverPatchHolography/Locality/NoSignalling.lean)
supplies the local-operator/reference algebra used above. Its theorem
`ptraceFst_local_kraus` supplies the complementary untouched-reference
marginal after an unread complete Kraus sum. The native tree attachment,
integer minorization and archive comparison are proved here; they are not
asserted to be formalized by those two existing lemmas.

The physical output of the declared instrument retains **all** records;
(4) specifies a comparison observation, not erasure by the source. Every
later common complete instrument contracts this state error if it consumes
the comparison state and does not consult the marginalized mixing records.
A bounded-gain numerical map multiplies the input error by its certified
gain and adds its own implementation error. Full archive-driven feedback
cannot use (4) without a separate argument.

The distinction is quantitative. Starting from a known nonconstant basis
configuration, the tree and every retained current determine the final
configuration y. A comparison which replaces that source by tau_q while
keeping the same full archive has half-trace distance exactly
1-1/d_q. For q=6 it is 923/924, regardless of the number of mixing blocks.
Thus full archive/source factorization is false even when (4) is small.
The retained tree archive can remain correlated with the untouched reference.

## Growing-source attachment and resources

A declared cell c is addressed by ports 12c,...,12c+11. The original
federation grammar admits precisely these unchanged internal seams at every
level where the source's integer resource checks pass. At evaluated levels
0,1,2 the owning controls query the actual native kernel, retaining each
raw slot and the original FIFO/archive. At higher declared levels the
original host resource stop precedes graph expansion; algebraic finite-word
attachment does not assert a completed host run or an infinite routing law.
Intercell routing is unnecessary for this single-cell preparation.

The unchanged node limit gives 0<=b<=68. The choice b=32 uses 1931 attempts,
1931 nodes and error at most 2^-32. The largest b=68 uses 4091 attempts and
leaves five literal nodes for any appended protocol. Each combined
preparation/readout/body word must pass the same total-resource checks.
The preparation is a fixed-seam reservation of the selected source cell;
it is not an interval of the autonomous globally uniform scheduler.
Its 11+60b attempts must be charged to the common counted clock.

One FIFO slot and export-oldest retain L-1 exported records and one final
record for L=11+60b. Record width depends on the global address level, so
only the local command count and error are capacity-independent. All
program-description, storage and attempt checks remain active. The generic
expanded evaluator reaches its unchanged branch ceiling and returns
RESOURCE_STOP. Lazy individual raw-word queries provide exact finite
instrument controls, not a normalized sampling algorithm or an assertion
that the exponential archive was enumerated.

## Compatible references and the operational boundary

A canonical global restriction is diagonal and constant within each local
number sector. A Bernoulli product with any declared p is also a mixture
of the same tau_q. Equation (3) fixes these local states exactly. Conversely,
(3) cannot select their sector weights from an arbitrary input, nor replace
exterior correlations by an independent state. In particular an independently
fixed two-particle cell is not the full local restriction of a growing
canonical source. Multiple prepared cells preserve their actual joint local
number distribution. This theorem does not prepare a globally canonical
state or a pure transfer ground.

The complete native branch operators are weighted partial injections, so
K_h* K_h is diagonal for every fixed or record-adaptive repair/hold history.
Two input states with the same occupation-basis diagonal therefore have
identical complete classical record laws. Feeding only those records to
an independent quantum service cannot distinguish their phases: its output
is the same mixture over the same classical labels. The executable control
checks all 1500 raw seam/law addresses for two orthogonal relative-phase
states. The induction over histories supplies the all-word statement.
This does not deny coherent source output; it identifies the additional
source-to-service coupling needed to read it.

The apparatus of [the finite phase experiment](../native_phase_experiment/README.md)
prepares a declared twelve-bit register from blank inputs using coherent
services and measures phase after an explicit occupation-register embedding.
Its fair state P/Tr(P) is not a sector-uniform state or the transfer ground.
That experiment supplies no automatic load from an arbitrary growing
federation into its apparatus. The separate nine-bit source-word grammar
also supplies no such map. This package derives the conditional mixed
preparation (4) from native repair alone. A growing coherent preparation,
selected density, physical vacuum, relativistic action and particle masses
are outside its conclusion.
