# OPH Selection Ledger

OPH builds public reality from finite observers that read patches of a shared
screen, compare overlaps, and repair disagreement until every observer agrees.
On the declared carrier lineage that consensus requirement forces the twelve-port
screen geometry outright, with machine-checked proofs, and every remaining
structural choice is finite: a declared premise with a stated menu, or an open
item with a named owner. There is no tunable layer anywhere in the list. This
ledger records every discrete structural selection in the framework exactly
once, classified as FORCED (derived from the axioms and source objects), EXPOSED
PREMISE (chosen and declared in theorem antecedents), or OPEN (neither, with the
owning GitHub issue). It shows exactly what is proved and exactly what is left:
pick an open row and work on it. It is the canonical free-versus-forced surface
for
[issue #554](https://github.com/FloatingPragma/observer-patch-holography/issues/554).

## Forced by consensus

| # | Selection | Where it is forced |
| --- | --- | --- |
| 1 | Twelve-unit port splitting (all twelve port weights equal one) | Issue #565 closure packet (`code/a5_closure`): unique all-one total-12 allocation with exact floor 12, next floor 14, gap 2. Lean: `Lean/ObserverPatchHolography/Screen/UnitSplit12.lean`, `unit_split_of_positive_sum`. Conditional on row 4. |
| 2 | Inverse port pairing and six axes | Issue #565 closure packet: unique distance-3 antipodal pairing and six axes on the port incidence, with negative controls. Conditional on row 4. |
| 3 | Icosahedral screen selector and port frame ($A_5$ action, rank-3 Gram frame) | Issue #565 closure packet: $\mathrm{Aut}^+\cong A_5$, full $\mathrm{Aut}\cong A_5\times C_2$, Gram identity $G^2=4G$ with trace 12. Lean: `PortFrameGram.lean`, `A5PortAction.lean` (kernel `decide`). Conditional on row 4. |

## Exposed premises

| # | Selection | Where it is declared | Menu |
| --- | --- | --- | --- |
| 4 | Echosahedral carrier lineage (icosahedral geometry for the screen carrier) | Declared quotient-visible carrier lineage of the #565 selector; the packet does not derive that every OPH carrier is echosahedral. Lean `PhysicalA5ForcingNoGo.lean` (`noSourceOnlyChargeReconstruction`) proves the bare total-charge reduct cannot force the choice, so a declared premise is the exact boundary. | 3 (tetrahedral, octahedral, icosahedral) |
| 5 | Compact-gauge refinement receipt | Declared in the D7 antecedent of the compact paper and the synthesis paper: cofinal tail with coherent surjective charge pullbacks, finite block embeddings, tensor realizations, $3{+}1$D symmetry, and compatible forgetful fibers. | declared receipt |
| 6 | MAR matter class with the one-generation/one-Higgs witness ($N_g=3$ as least economy value) | Minimal admissible realization is axiom five; the low-energy packet with the one-generation/one-Higgs chiral witness is declared in the D8/D9 antecedents. The packet admits $3\le N_g\le 5$ and selects $N_g=3$ as the least value of the declared economy class. | 3 ($N_g\in\{3,4,5\}$) |
| 7 | Z6 reserve pricing input $\ell_{\text{shared}}=P/4$ | Declared branch input to the reserve-trace theorem (`thm:z6-reserve-trace`, screen microphysics paper); every downstream reserve number inherits it. | declared input |
| 8 | BR-0: cell product reading of the screen count | `code/capacity_readback/F_CONSTRUCTION_2026-07-14.md`, Step 2: reading (i) is the only executable reading and is carried with its warning as an inherited hypothesis; it reproduces the declared screen dimension exactly. | 2 |
| 9 | Publicness policy | Frozen publicness policy declared in the finite closure construction of the synthesis paper; deriving it from source objects belongs to the issue #505 program. | frozen policy |
| 10 | Scheduler selection (fair-block contraction branch) | Declared in the consensus-paper fair-block antecedent, where the implementation supplies $(\lambda,\varepsilon,A,\beta,L)$ for the contraction statement. | declared branch |
| 11 | Quantitative closure map declaration (pixel closure) | Each candidate closure map is declared before evaluation; the interval certificate proves exactly one fixed point per declared map on the physical domain (compact proof, certificate C1). The physical Thomson identification is a separate receipt. | per declared map |

## Open: the work-wanted list

| # | Selection | Owner and state |
| --- | --- | --- |
| 12 | Physical port-current realization | Lean `PhysicalA5ForcingNoGo.lean` (`noSourceOnlyCurrentReconstruction`) proves port count and total charge cannot select it; the full-rank commutator-closed current lift from source response fields is the outstanding construction. Issue [#567](https://github.com/FloatingPragma/observer-patch-holography/issues/567), with the closed #566 packet as declared input. |
| 13 | Determinant, spin, loop, axis, and refinement descent to the physical $\mathbb Z_6$ global form | Abstract six-axis quotient arithmetic is proved in Lean `Z6Exact.lean`; the physical cocharacter-lattice descent and exact line spectrum are the target of issue [#567](https://github.com/FloatingPragma/observer-patch-holography/issues/567). |
| 14 | Physical matter spin lift | The cover-level matter tensors and spin lift are the target of issue [#314](https://github.com/FloatingPragma/observer-patch-holography/issues/314) (Super-Tannakian matter lift), feeding row 13. |
| 15 | Screen-to-family attachment (rank-45 map to three rank-15 families) | The canonical lowest rank-three screen band exists; its source-derived complex-linear attachment to the physical matter-pole residue space is the target of issue [#569](https://github.com/FloatingPragma/observer-patch-holography/issues/569). |
| 16 | Capacity-closure selector (exact finite-size slack law with one physical zero) | Source-only correctable public-record closure and the unique slack zero are the targets of issues [#505](https://github.com/FloatingPragma/observer-patch-holography/issues/505) and [#551](https://github.com/FloatingPragma/observer-patch-holography/issues/551). |
| 17 | Capacity-readback branch axes BR-1 through BR-6 (reserve semantics, reserve attachment, readback count effect, subfederation marking, symmetry quotient, readback family; menus 3, 2, 4, 3, 3, 4) | Declared menus in `F_CONSTRUCTION_2026-07-14.md`, Steps 3-7, with the fifth bridge reading and the CAP-B family barred before evaluation; the branch axes are selected by the slack law of row 16. Owner: issue [#551](https://github.com/FloatingPragma/observer-patch-holography/issues/551). |

Totals: 17 selections. 3 FORCED, 8 EXPOSED PREMISE, 6 OPEN. Every open row
belongs to the particle-realization and capacity programs; the recovery of
relativity and quantum structure from observer consensus consumes only the
forced rows and declared premises above.

## How to read a row

A FORCED row is a theorem: given the axioms and the cited premises, the
selection has no alternative, and the citation names the machine-checked or
certified proof. An EXPOSED PREMISE row is a declared choice: the theorem
downstream of it states the choice in its antecedent, and the menu column
records how many alternatives the declaration passed over. An OPEN row is a
selection that is neither derived nor stated as a clean antecedent; the owning
issue tracks the work that moves it into one of the other two classes. A row
can move up (OPEN to EXPOSED PREMISE to FORCED) and never silently disappears:
the consensus argument is exactly as strong as the FORCED column plus the
declared premises it names, and every one of those premises is a finite menu.
Every open row is a well-posed problem with an issue link; solving one moves
it up the ledger.
