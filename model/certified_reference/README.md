# Exact error of the literal plaquette source controls

The complete plaquette preparation in `plaquette_reference` no longer needs
an assumed `2e-11` mathematical complete-output tolerance. This package
proves an error below `1.700e-12` for its **existing literal source controls**,
including angle rounding, twelve-port coefficient rounding, all 256
classical records, and every first export failure. No source pulse is changed.

The mathematical source is the retained exact six-level response with
entries in `Q(sqrt(5))`. Each actual `integrated_controls` binary64 output
is interpreted as its exact rational value. This is neither an assertion
that the floating `expm` matrices are correctly rounded nor a theorem that
physical hardware realizes the specified source. Initialization, coherent
load/export, complete readout, conditional scheduling and record memory
remain explicit simulator services.

## Independent finite certificates

An independent integer Horner recurrence checks all 65,536 entries of
`Q8=sum_{k=0}^8[-(H-I)/64]^k/k!` against the original full-Fock action.
It also checks all exact column norms and their probabilities. The producer
uses summed powers; this checker does not repeat that computation.

The checker parses every native word against an independently written
Pauli/CNOT/frame template. It reads the **emitted** binary64 Walsh angles,
forms each inverse Walsh sum in rational arithmetic, and encloses its
sine and cosine. Each resulting two-vector is compared with the exact
integer subtree norms, including every signed final leaf. It does not
trust the producer's `atan2`, floating norm or angle-rounding accuracy.
The 65,535 conditional nodes cover the probability word and all 256
conditional column words. Maximum prefix errors add across eight tree
levels. The proof concerns the specified initialized input; it does not
claim arbitrary-input accuracy for a state-preparation word.

The outward dyadic kernel has 160 bits. Square roots use integer square
roots; sine/cosine use 64 Taylor terms with an explicit geometric absolute
tail for `|angle|<=4`; pi uses Machin's identity and 80 alternating inverse
tangent terms. No floating transcendental evaluation or decision enters
these bounds. All fixed basis, CNOT, blanking-phase and inverse-frame
angles are separately compared with their exact multiples of pi.

The twelve original source generators and all ten control identities are
checked exactly. For each of 7,612 distinct literal pulse-control vectors,
the checker encloses the full `6x6` difference
`G_actual - angle*T_kind`. Anti-Hermiticity makes its absolute row bound
an operator-norm bound. Duhamel then gives
`||exp(G_actual)-exp(angle*T_kind)|| <= ||G_actual-angle*T_kind||`.
This includes rounding-induced coupling outside the intended proper code.
Both coherent and dephasing export Kraus families are checked exactly on
all six inputs; their common Stinespring isometry preserves the error.

Chronological error telescoping includes every abort output. Orthogonal
classical records give the complete preparation bound

```
epsilon_prep <= epsilon_probability + max_j epsilon_conditional_j
             < 1.700e-12.
```

The maximum is over all 256 branches; it is not an omitted-outcome average
or successful-sample normalization. The full exact source instrument is
trace preserving, and its abort probability is at most
`epsilon_prep^2 < 3e-24`. Actual source controls may have nonzero abort
probability even though their ideal exact-generator counterparts do not.

## Regional modular consequence and continued action

Combine the proved preparation error with the preceding rational polynomial
Gibbs design error `delta`. Put `d_exec=(2+epsilon_prep)epsilon_prep`.
If `m` is the preceding rational regional eigenvalue floor, then

```
lambda_min(rho_A,success) >= m-d_exec > 1/20,
||log rho_A,success - log rho_A,Gibbs||
  <= (delta+d_exec)/(m-d_exec).
```

The source success field remains **subnormalized**, with its abort
complement retained separately. Its logarithm is compared directly, with
no success conditioning. Adding the exact regional eighth-order remainder
gives a total bound below `2e-9`, now proved for the literal mathematical
source rather than conditional on an unproved execution tolerance. The
normalized Gibbs scalar `log Tr_A(E exp(-beta(H-I))) I` is retained in the
absolute comparison and drops only for trace-zero entropy tangents.

An independent route also checks the original 4,320-pulse field continuation
at modular parameter `s=1/2`, including every FSWAP, phase and contact word.
Its literal-source error is below `5.380e-13`, plus the existing exact
Strang bound `1463/805306368`. This continuation has an all-input operator
bound, including arbitrary untouched references and all export aborts.
The full-reference logarithm bound connects it to the actual success
reference modular action, up to the stated Gibbs normalization scalar.

This closes one finite precision premise. It does not remove the declared
source/service specification, prepare a scalable thermodynamic reference,
derive a physical clock or establish any geometric Einstein premise.
The negative controls reject corrupted tree phases, missing routing or
blanking, wrong polynomial coefficients, changed port controls and omitted
export outcomes, including zero-weight complements.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/certified_reference
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.certified_reference --check model/certified_reference/results.json
```
