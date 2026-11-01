# The complete binary repair algebra and its protection boundary

The existing [all-occupation instrument](../exchange_sector/MANY_BODY.md)
has a complete source-only conservation catalogue. For either exchange lift,
the native C/E branches generate

```text
A = (direct_sum N=1,...,11 M_{binomial(12,N)}) direct_sum C I_E,
E = span{|empty>, |full>}.

A' = (direct_sum N=1,...,11 C I_N) direct_sum M_2(E).
dim_C A = binomial(24,12)-1 = 2,704,155;   dim_C A' = 15.
```

The exceptional two-dimensional multiplicity is essential: empty and fully
occupied configurations give identical complete raw records and both remain
unchanged. Their coherent span is an exact noiseless qubit for the declared
native repair. It has distance one against local errors and requires
coherence between global occupation numbers zero and twelve. It is not a
localized topological defect.

This is the full 4096-state component domain already used by
`many_body.py`. The frozen [class V2 admission](../completion/program_class_v2.json)
instead fixes global total number; this result does not change that contract.
Within any single number sector the fixed observable algebra is scalar.
If all numbers are allowed only as a classical superselection sum, the fixed
observable algebra is `span{P_0,...,P_12}`, of dimension thirteen. The
off-diagonal extremal memory requires a separately admitted coherent apparatus.

No proper-port-local observable is exactly fixed under the full admitted
repair. This excludes a nontrivial **fixed local source charge** in this
finite class. It does not exclude particles that move while carrying a
global conserved charge: total occupation is exactly such a global charge.
It also does not exclude controlled lifetimes, active correction, charges
in a joint source/controller/history algebra, a different carrier or a
refining limit.

## Actual branches and a finite constructive proof

[`occupation.py`](occupation.py) uses the actual signed partial injections
returned by `many_body.branches`. For each law and character it retains every
coin, raw, zero and reserved abort slot: 240 slots for C or 480 for E, with
245,760 input columns in either law. All twelve ports and all thirty native
seams are retained. The fair Kraus factor is `1/sqrt(2)` before seam
selection. Positive selector factors do not alter the generated algebra.

1. **Effects identify every nonconstant occupation.** On seam `(a,b)`, the
   moving raw effects are `n_a(1-n_b)` and `(1-n_a)n_b`; their difference is
   `n_a-n_b`. If two bit strings have the same edge differences on a
   connected seam graph, their difference is constant at every port. Binary
   values leave only identical strings or the pair empty/full. Twenty-two
   actual moving effects on an eleven-edge native spanning tree therefore
   give 4095 distinct signatures. Products of effects and their complements
   construct every nonextremal rank-one diagonal projector, and the remaining
   projector `P_E` onto the two extremes. The verifier checks every signature
   directly against the actual branch support.

2. **Transport generates each full fixed-number matrix algebra.** For a
   nonzero branch column, multiplying on both sides by the constructed
   endpoint projectors gives `P_y M_h P_x = s_h(x)|y><x|`, with
   `s_h(x)=+1` or `-1`. Scalar multiplication removes this nonzero sign inside
   the algebra. A checked forest supplies 4083 actual directed transport
   witnesses, rooted once in each of the eleven intermediate number sectors.
   Every forest edge is independently checked against its recorded native
   branch and sign. Multiplication and adjoints along the tree produce all
   matrix units in that sector. This is an algebra proof; it is not a gauge
   equivalence between positive and negative complete instruments.

3. **There are no additional algebra blocks.** Every native branch preserves
   number. On `E` its restriction is exactly the identity for the appropriate
   raw-zero slot or exactly zero otherwise, including both coins. The full
   inventory checks this equality, rather than checking only branch norms.
   Thus the generated algebra is contained in the displayed direct sum.
   The first two steps generate that entire direct sum, proving equality.
   Its commutant is consequently the displayed eleven scalar blocks and
   the full extremal `M_2` multiplicity.

`algebra_certificate` produces the witness; `verify_algebra` replays supplied
effect labels, all signatures and the signed forest without rerunning the
producer's search. It rejects incomplete inventories, false dimensions,
changed signatures and invalid or incomplete transport witnesses. No dense
4096-square matrix, eigenvalue threshold or floating rank is used.

The retained `certificate()` is compact: it executes and verifies all four
full witnesses, then records their SHA-256 digests, dimensions, complete raw
inventory counts and transport census. It similarly checks and summarizes
the faithful-reference and peripheral witnesses. It does not retain tens
of thousands of regenerable signature, forest and diagonal entries in a
second public JSON file. A full raw witness can be inspected directly:

```sh
python -c "import json; from model.native_fixed_algebra import occupation as o; print(json.dumps(o.algebra_certificate('C', -1), indent=2))"
```

`peripheral_certificate`, `faithful_certificate`, `locality_certificate`
and `extremal_certificate` expose the corresponding full witnesses. The
digests use sorted-key JSON with separators `(',', ':')`, UTF-8 encoding,
and no nonfinite values. Changing a raw witness changes its retained digest.

Both C and E have the same moving effects and transition support. E also
refines the C raw partition. The proof therefore holds with any per-seam
C/E mixture, including pure C or pure E on different seams, provided every
native seam and both endpoint coins have positive support. An identity hold
and a positive repair probability do not change the algebra. The argument
also tolerates either supported sign choice: scalar phases cannot create
an internal multiplicity once the effects separate a sector's vertices.

The underlying algebra argument generalizes to any finite connected seam
graph with this binary endpoint law. Edge differences have the same
separation property and the exclusion-token graph is connected in each
fixed-number sector. One elementary connectivity proof restricts to a
spanning tree and fills or empties a leaf to match a target configuration:
move the nearest needed token or vacancy along its path, fix that leaf and
continue on the smaller tree. This is a conditional family theorem; no
particular OPH refinement or limiting physical locality is selected by it.
The executed certificate here is on the actual twelve-port native graph.

## Faithful references, conservation and stationary states

For fixed strictly positive endpoint weights `theta_p`, the endpoint coin
puts a single seam load at `a` with probability
`theta_a/(theta_a+theta_b)`. Define the faithful diagonal proof reference

```text
pi(S) = product_{p in S} theta_p / product_p(1+theta_p).
```

Conditioned on the occupations outside a seam, the two single-load weights
are in ratio `theta_a:theta_b`. The endpoint heat bath preserves that ratio;
equal endpoints hold. Signed columns do not alter diagonal probabilities.
Every seam's complete C and E maps therefore preserve the same `pi`.
`faithful_certificate` checks both completeness and this identity exactly
on every column, independently for all thirty seams. The retained fair and
alternating `(1,2)` controls use rational arithmetic. General positive
weights follow from the displayed two-endpoint identity. This is a
mathematical weighted extension; a biased full-carrier source service has
not been executed by this module.

This proof reference is not a dynamically selected preparation. Multiplying
all `theta_p` by a common factor does not change the endpoint coins, while
it does change the full-reference number populations. More generally any
positive choice of number-sector weights supplies a faithful stationary
reference when combined with the corresponding conditional product laws.

For the source Heisenberg channel `Phi(X)=sum_h K_h* X K_h`, completeness and
faithful invariance imply

```text
Re Tr[pi X*(X-Phi(X))] = (1/2) sum_h ||X K_h-K_h X||_pi^2.
```

If `Phi(X)=X`, every commutator vanishes. Applying the same identity to
`X*` gives commutation with the adjoints, so the fixed space is precisely
`A'`. Conversely every member of `A'` is fixed by completeness. This extends
the branchwise classification to **all** source Heisenberg conserved
observables, not merely conserved expectation in one chosen state.

Positive source-independent selector weights, a positive repair probability
and an arbitrary C/E mixture on each seam preserve this conclusion. For a
classical controller state `h`, universal conservation of a source-only
observable requires it to be fixed by the conditional source map `Phi_h`
on every admitted source/controller input. The same classification applies
at each such state. This assertion does not supply a stationary density or
a mixing theorem for a growing, correlated record-feedback process.

For a **fixed** selector and endpoint weights, the full stationary density
catalogue is

```text
rho_stationary = B_E + sum_{N=1}^{11} p_N pi_N,
B_E >= 0 on E, p_N >= 0, Tr(B_E)+sum_N p_N = 1,
pi_N(S) = product_{p in S} theta_p / sum_{|T|=N} product_{p in T} theta_p.
```

The sector populations are initial data. Within any intermediate sector,
the sole normalized stationary state is `pi_N`, which is faithful on that
sector and is mixed when the sector has dimension greater than one. The
only nontrivial passive quantum information multiplicity is the extremal
two-dimensional block. Scalar multiplicity in the other sectors means no
protected quantum information; it does not assert the existence of a
stationary pure one-dimensional code there.

## A convergence statement without an inferred clock or rate

There is also an exact obstruction to nontrivial peripheral eigenvalues.
Sum the unweighted holding Kraus maps across all seams and coins, including
both E holding raw labels when present. The result is a strictly positive
diagonal operator `D`. The actual certificate checks its 4096 diagonal
entries; they lie between 40 and 60. For any supported per-seam C/E choice,
the same `D` lies in the linear span of its Kraus maps by rescaling their
positive coefficients.

If `Phi(X)=lambda X` with `|lambda|=1`, faithful invariance gives

```text
sum_h ||X K_h-lambda K_h X||_pi^2 = 0.
```

Thus `X D=lambda D X`. Any nonzero matrix entry satisfies
`D_jj=lambda D_ii`. Both diagonal entries are strictly positive, so
`lambda=1`. A unital completely positive map is power bounded and hence
has no nontrivial Jordan block at eigenvalue one. Finite dimensionality
therefore gives convergence to the fixed space, with all other eigenvalues
strictly inside the unit circle. No quantitative uniform gap or physical
time has been inferred.

The dual Schrödinger limit is consequently

```text
lim_t Phi_*^t(rho) = P_E rho P_E
                    + sum_{N=1}^{11} Tr(P_N rho) pi_N.
```

To identify this limit, note that the extremal matrix units and every
`P_N` are conserved observables. The displayed stationary family has
dimension fifteen, equal to the Heisenberg fixed-space dimension, and
contains all states with those conserved data. This proves both completeness
of the stationary catalogue and the asserted limit. It is not a claim of
adaptive-history convergence, an energy gap or a particle lifetime.

## Locality and the extremal memory control

Let `R` be a proper subset of ports and suppose a fixed observable has the
form `A_R tensor I_outside` in the occupation tensor product. The general
fixed observable consists of scalar number blocks and an arbitrary
extremal two-by-two block. Its extremal off-diagonals change every outside
bit, so they cannot occur in a proper-local operator. The remaining
diagonal depends only on total number, say `f(k+l)`, where `k` is the local
count and `l` the outside count. Locality requires this value to be
independent of `l`. The intervals `[k,k+|outside|]` overlap for successive
`k`, forcing every value of `f` to coincide. Thus the proper-local fixed
algebra is scalar.

The same conclusion holds for the proper CAR-local observable algebra.
Every invariant in the complete catalogue is parity even. The only
off-diagonals are proportional to the product of all twelve creation
operators or its adjoint; these are not supported on a proper set. Its
remaining diagonal part reduces to the preceding count argument. Fermionic
Jordan-Wigner strings are not silently identified with tensor-local odd
observables.

On the extremal code every complete native branch is a scalar multiple of
the identity or zero. Its quantum information, including an untouched
reference, therefore survives every finite classical record-adaptive native
schedule. However, for **every** individual port,

```text
P_E n_p P_E = |full><full|,
P_E (I-2n_p) P_E = logical Z.
```

A one-port phase flip swaps the cat states
`(|empty> +/- |full>)/sqrt(2)`. Independent phase flips of probability `p`
per added service reduce logical-X expectation by `(1-2p)^t`. Computational
measurement of one port removes their coherence. A one-port bit flip sends
the two codewords completely into the singleton and one-hole sectors.
These are explicitly added local-error controls, not additional native C/E
branches. They show code distance one and the failure of the local
Knill-Laflamme condition: compressed local occupation is not scalar.

The two cat states alone have identical reduced densities on every proper
subset. That fact does not imply local indistinguishability of their whole
two-dimensional code: the empty/full codewords are distinguished by one
port. The native protection comes from the specified endpoint law's
identical action on two global constant configurations. It supplies no
reference-like exterior for a localized charged core and no topological
robustness against admitted arbitrary local perturbations.

Finally, the exchange homology catalogue concerns holonomy of paths under
declared face and square deformations. Those two-cell relations do not
change the effect-separation or transport proof above. A choice of scalar
exchange character supplies phases on transitions; it does not by itself
produce a protected local state algebra. The explicit orientation-memory
construction remains a different enlarged source/observer system. No
species, Lorentz spin, physical rest mass, calibrated energy or completed
CP-03/CP-05 gate follows from this finite classification.

## Complete source histories on the twelve occupation bits

[`occupation_service.py`](occupation_service.py) compiles the existing
six-level proper-code pulses into the complete native endpoint operation on
seam `(0,5)`. It uses all twelve occupation bits, three C record bits or
four E record bits, and both exchange characters. This is an encoded finite
component, not an execution of the V2 federation interpreter. Its records
include the coin, signed transfer flags and E parity flag; every terminal
address and every first failed proper-code export is retained.

The general-input guard checks the ideal compiled word against the native
branch rule on every one of the 4096 input occupations. Only six occupation
bits participate in this seam's rule, including the four sign spectators;
an independent 64-column matrix comparison and an explicit check of every
full-carrier column establish the tensor extension. This prevents an
identity replacement from passing just because the extremal code is fixed.
Each actual primitive's full success and complementary output columns are
then compared with its independent ideal columns. A fresh numerical
operator-norm admission of `epsilon = 10^-11` per pulse gives, for a word
of length `L`, the complete-output bound `2 L epsilon/(1-L epsilon)`.
The general-input repair certificate also includes a separate `10^-11`
ideal-word allowance. These are fixed numerical admission ceilings, not
outward-rounded interval proofs or measured hardware accuracies.

Sixteen preparation/repair/readout histories cover four cases for each C/E
law and exchange character:

| Preparation and intervention | Complete ideal readout |
| --- | --- |
| Extremal cat, native repair | Inverse preparation recovers the logical input and its untouched reference |
| Extremal cat, one-port phase flip | Logical phase is reversed; the inverse-preparation output changes |
| Extremal cat, retained one-port pointer | The pointer distinguishes empty/full and removes their coherence when ignored |
| Complementary six-occupation code, native repair | Complete native records distinguish the components; no passive qubit survives |

The six-occupation code uses
`S={0,1,2,3,4,6}` and its complement. A two-mode beam operation and five
controlled swaps prepare it from two declared six-occupation inputs. All
four beam columns and all eight controlled-swap columns independently
verify number conservation of those ideal macros. Intermediate proper-code
pulses need not conserve the encoded occupation number. The extremal probe
instead uses supplied zero/one-occupation inputs and a number-phase
preparation; its encoded codewords have occupations zero and twelve. It
does not amend V2's prohibition on cross-global-number coherence.

Each history executes two logical columns, so its complete amplitude bound
also controls arbitrary superpositions and an untouched reference. The
largest register has seventeen work bits and 131072 terminal addresses.
Preparation, native repair, intervention and inverse preparation all count
toward the pulse, load/export and failure-slot resources in the
[receipt](occupation_results.json). Blank apparatus bits and supplied initial
occupations are explicit. No successful-output conditioning, restart or
renormalization is used. Literal immutable primitive and input archives
reconstruct every failed prefix, including complex off-diagonal reference
information; aggregate failure weight is not substituted for those outputs.
The archive keeps immutable payloads and shape/type descriptions; public
arrays are fresh views, so changing their metadata cannot change later
reconstruction. Input and target columns must convert to complex128 without
value loss, and reference Gram matrices must remain finite. Unsupported
precision or overflow is refused before it can become a successful admission.

The owned replay starts a fresh isolated Python process with an empty bytecode
cache, binds its complete source and data closure before and after execution,
and admits every complex output before producing compact exact ideal laws
and fixed budgets. Unlisted *ideal* addresses are zero; actual small leakage
is still retained and checked. The retained JSON contains no platform-specific
residual floors. Full numerical outputs remain available through
`experiment(..., return_run=True)`.

```sh
python -m model.native_fixed_algebra.occupation_packet --check
python -m pytest -q -W error model/native_fixed_algebra
```
