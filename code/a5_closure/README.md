# Icosahedral closure certificate (C19)

Exact representation-theoretic certificate for the icosahedral closure theorem:
the twelve-port screen module, the compact-Lie trichotomy, the rank-five
noncentral corollary, and the oriented face-phase multiplicities.

Referenced as **[C19]** in `extra/compact_proof_of_oph.tex` and as the
closure certificate in `paper/observers_are_all_you_need.tex` and
`paper/recovering_relativity_..._compact.tex`.

## Run

```bash
python3 a5_screen_sm_closure.py          # module identity, trichotomy, face phases, kinetic relation
python3 exterior_sm_completion.py        # exterior matter, anomalies, weak load, deck control
python3 a5_selection_certificate.py      # sharpness data for Cohn-Kumar universal optimality
python3 a5_compact_lie_classifier.py     # compact-Lie enumeration
python3 a5_harmonic_decomposition.py     # angular multiplet sequence
python3 bh_log_correction.py             # conditional horizon log-coefficient decision tree
python3 independent_trichotomy_check.py  # independent re-derivation, trusts nothing above
python3 test_audit.py                    # regression suite
```

Requires Python 3.11+ and SymPy. The suite exits 0.

## What is certified (exact; no floating-point fit, no measured number)

| Object | Statement |
|---|---|
| Vertex module | `chi_P12 = (12,0,0,2,2)`, so `P12 = 1 + 3 + 3' + 5`, multiplicity-free |
| Adjacency spectrum | `det(xI-A) = (x-5)(x+1)^5(x^2-5)^3` → canonical ranks `1,3,3,5` |
| SM adjoint restriction | `ad su(3) = End_0(3') = 3' + 5`; with `su(2) → 3`, `u(1) → 1`, the total is `1 + 3 + 3' + 5` |
| Icosahedral selection | 3 distinct inner products `{-1, ±1/sqrt5}`, spherical 5-design ⇒ **sharp** (`m=3`, strength `2m-1=5`). By Cohn–Kumar (JAMS 20, 2007) it uniquely minimizes every strictly completely monotonic pair cost of squared distance, up to `O(3)`. |
| D-optimal selector | Maximizing `det(F1) det(F2)` at fixed vector/quadrupole trace gives `F1=2 I3`, `F2=(4/5) I5`; the six centered projectors form a regular simplex, hence the unique real ETF(3,6), the icosahedral axes. The objective remains a source-side branch premise. |
| Compact-Lie trichotomy | Exactly three algebras: `u(1)^12`, `su(2)^2 + u(1)^6`, `su(3) + su(2) + u(1)` |
| **Inner-action closure** | If the `A5` action is **inner**: `dim Z(g) <= 1` (inner autos fix the center pointwise); `Z(g)=0` forces `su(2)^4` whose fixed-space dimension is a multiple of 3, contradicting `dim g^{A5}=1`; hence `dim Z(g)=1`, semisimple dim 11, and `11=8+3` uniquely ⇒ `su(3)+su(2)+u(1)`. **Needs no `W5-NONCENTRAL` receipt.** |
| Angular multiplets | `l=2: 5` (irreducible); `l=3: 3'+4`; `l=4: 4+5`; `l=5: 3+3'+5`; `l=6: 1+3+4+5`. First nonconstant invariant at `l=6`. |
| Horizon log coefficient | **Conditional decision tree, not a universal law.** Unconstrained `q`-state cells give `S=K log q` exactly, i.e. `c=0`. Exact 12-port balance: `c=11/2`. 24-slot balance: `c=23/2`. Nonabelian `SU(3)xSU(2)` singlet saddle: `c=11/2`. Full 12-dim gauge singlet incl. `U(1)`: `c=6`. The horizon measure must be derived before `c` is frozen. |
| Rank-five corollary | `W_5` noncentral ⇒ `su(3) + su(2) + u(1)` uniquely |
| Face-phase multiplicities | `m_w = (dim V - chi_V(3A))/3 = (0,1,1,1,2)` on `(1,3,3',4,5)`; minimal nontrivial extension has dimension 3 |
| Missing-four/Higgs no-go | **Rejected.** `End_{A5}(4) = R`, so no `A5`-invariant complex structure exists for a commuting hypercharge `U(1)`. The identity `1+3+3+4+5=16` is dimensional only. |
| Exterior matter witness | On a declared trace-balanced carrier `V=C+W`, the selected non-vacuum even package `Lambda^2 V + Lambda^4 V` is one 15-state chiral Standard Model generation. Its hypercharges, three one-Higgs invariant lines, perturbative anomalies, and weak-doublet parity check exactly. |
| Weak multiplicity | The exterior package contains three colored weak-doublet copies and one lepton doublet, so the per-generation multiplicity is exactly `4`. Its identification with a physical port load requires `PORT-LOAD-TRACE`. |
| Free-deck control | Any declared free `Z6` action on 24 slots has four orbits and a four-dimensional invariant-function space. A physical result requires the deck action and `PORT-WEAK-INTERTWINER`. |

## Why the trichotomy holds (proof sketch)

A compact Lie algebra is reductive, `g = z + [g,g]`, and both summands are
characteristic, hence `A5`-stable. The identity component of the center is a
torus whose exponential kernel is an integral lattice preserved by every group
automorphism, so the `A5`-action on `z` is defined over `Q`. Because `3` and `3'`
are Galois conjugate over `Q(sqrt5)`, a rational submodule contains them with
equal multiplicity. Enumerating rational centers against compact semisimple
dimensions leaves exactly three cases; the other five die because no compact
semisimple algebra has dimension 1, 5, or 7, or because an `A1^4` / `A1^2`
adjoint cannot supply `5` (any homomorphism `A5 → S_4` is trivial, `A5` being
simple of order 60).

## Claim boundary

The exact chain is a conditional finite `A5` module, a compact coefficient
bracket, trace-balanced block integration under coefficient `1/2`, and a
six-axis lattice quotient `Z6`. Central record projectors commute; physical
gauge currents and global form require the full-rank compact skew-adjoint
receipt **PORT-CURRENT-INNER**,
**PORT-SPIN-LIFT**, block/axis descent, and the MAR matter/tensor-kernel premises.
The coefficient construction supplies no observational discriminator.

The exterior calculation supplies a second conditional matter route. It does
not turn coefficient directions into physical currents or select the global
form. It is also not the full even Clifford module: that module contains
`Lambda^0 V` as well as `Lambda^2 V + Lambda^4 V`. Removing the singlet and
excluding other anomaly-free light sectors requires the explicit MGFC/MAR
premise or a source-derived observer-visible discriminator. Selecting `H=W`
as the physical scalar is a separate gate. If the face-phase `A5` action remains
exact on the family multiplicity space, it restricts Yukawa tensors to
`A5`-invariant pairings; general family matrices require a source-derived
breaking, hiding, or forgetting mechanism.

### Open receipts

- **UD12**: emit a feasible integer charge domain and additive strict
  unit-splitting cost.
- **RP-A5**: emit either the D-optimal tomography objective or a strictly
  completely monotonic pair kernel without an icosahedral template.
- **PORT-CURRENT-INNER**: construct a full-rank `A5`-equivariant map onto a
  twelve-dimensional compact commutator-closed current algebra, prove that the
  induced action lies in `Int(g)`, and intertwine refinement.
- *(weaker group-level branch only)* **A5-COMMON-ACTION** + **W5-NONCENTRAL**:
  one group-level action shared by ports and gauge reconstruction, plus one
  source-derived repair composition with nonzero projected `W_5` commutator.
- **BLOCK-DETERMINANT-BALANCE**, **PORT-SPIN-LIFT**, and
  **AXIS-CENTER-DESCENT**: select total trace zero, land the mathematical
  Clifford lift physically, and identify port-loop winding with the deck class.
- **A5-FAMILY-ATTACHMENT** (family corollary only): prove the chiral family
  fiber's local face-corner phase is the restriction of a global `3` or `3'`
  action, independent of port labels, worker IDs, chart choices, and refinement
  presentation.
- **A5-FAMILY-DESCENT**: derive how the selection symmetry is broken, hidden,
  or forgotten before general family Yukawa matrices are admitted.
- **EXTERIOR-PACKAGE-SELECTION**: select the non-vacuum even exterior package
  and `H=W`, then exclude the vacuum singlet and other anomaly-free light
  generators.
- **PORT-WEAK-INTERTWINER** and **PORT-LOAD-TRACE**: identify a physical
  four-dimensional screen invariant with the four weak-doublet copies and
  prove that its normalized additive load is the physical `4P` readout.

## Conditional coupling relation

If the physical quadratic gauge kinetic operator on the port module is an
adjacency polynomial `K = f(A)`, color coherence (the color triplet and quintet
sharing one kinetic coefficient) forces, at degree two,

```text
k1 = 3*k2 - 2*k3
```

independent of which inequivalent triplet carries color (both assignments
verified). At degree one it collapses to bare coupling unification. A physical
discriminator requires a source-derived polynomial degree, kinetic
normalization, carrier scale, and complete threshold/RG map. The relation is
not a forward test.

## Formal status

`Lean/ObserverPatchHolography/Screen/Compact12.lean` formalizes the abstract
`u(3) + so(3)` matrix commutator, Lie laws, dimension, and a matrix
noncentrality witness. `Z6Exact.lean` formalizes the six-axis lattice quotient,
`S2DesignSignature.lean` formalizes the `11/25` arithmetic, and
`UnitSplit12.lean` proves only that twelve positive integer weights summing to
twelve are all one. The slot count, Euler/cost law, source selector,
icosahedral frame, `Phi`, face orientation, `Theta`, `A5` equivariance,
compact-Lie trichotomy, trace-balanced group, and physical descent are paper
proofs or open receipts.

## Novelty boundary

Finite-group flavor models use triplet representations, and Standard-Model
global quotients and line sectors are established subjects. This certificate
makes no priority claim. Its scoped result is the full gauge-adjoint module on
the icosahedral vertices plus the conditional compact-Lie classification.
