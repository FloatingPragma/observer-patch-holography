import Tower.ConsensusTower
import Tower.PublicWorldQuotient
import Tower.FixedPointEndpoint
import Tower.EventGeometryReadout
import Tower.FullPremiseInhabitant
import Tower.OperationalObserver

/-!
# OPH construction tower umbrella

`Tower.ConsensusTower` defines the timeless finite refinement root in the
`OPH.Tower` namespace. `Tower.PublicWorldQuotient` and
`Tower.FixedPointEndpoint` add the bounded A4 kernel quotient, the descended
OPH repair endpoint, and the conditional schedule-independent endpoint
theorem. The full endpoint theorem retains explicit termination, confluence,
semantic-completeness, and quotient-congruence premises.
`Tower.OperationalObserver` packages clauses 1 through 6 and the own-observer
half of clause 7 over an arbitrary tower, with a nontrivial finite witness and
four negative controls. The separate fixed-regulator cross-observer half is
implemented over the QFT access cut in `QFT.OperationalOverlapEvidence`; it is
kept out of this lower-layer umbrella to avoid reversing the import boundary.
Neither receipt selects a unique observer or attaches an instrument,
consciousness claim, or source realization. This
umbrella supplies no source realization, physical world, clock, causal net,
geometry, or continuum limit. `Tower.EventGeometryReadout` adds the selected
settled branch, the all-quotient event and geometry readout fragment with
its fine-to-coarse naturality squares, and the boundary-fibre theorem from
an independently constructed complete event signature, with the
boundary-as-signature circularity made definitional.
`Tower.FullPremiseInhabitant` exhibits one explicit branch
satisfying every premise of the conditional endpoint theorem
simultaneously, with genuine choice points and Newman-lemma confluence,
so schedule-independent representative-independent unique public
endpoints hold on it unconditionally.
-/
