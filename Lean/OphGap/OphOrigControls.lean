import OphGap.OphOrigWitness

/-!
# Ending B measured, and negative controls for the transport on the REAL data

`transportCheck` is a conjunction; each control makes one conjunct false by kernel evaluation
and `transportCheck_false_of_*` turns that into refusal of the whole certificate. Evaluating
only the failing conjunct keeps each control to a single bounded kernel run (see the report's
BUILD section for why the whole certificate cannot be re-run in one process).
-/

namespace GapScout.Oph

/-! ## ENDING B, measured: the existing checker cannot take the original ids at all.
`check` needs vertices `0..n-1` and tree edges first; on OPH's raw data it is `false`
(kernel-evaluated), which is why a transport theorem is needed rather than a re-run. -/

set_option maxRecDepth 100000 in
theorem endingB_direct_check_false : check origNodes.length origEdges cycle = false := by
  decide +kernel

/-- Collide two carriers: send 14860 (the witness's second node) to the label of 5923. -/
def collidedNewT : NTree := newT.set 64 14860 (newT.get 5923)

set_option maxRecDepth 100000 in
theorem collision_nodesOk_false : nodesOk n collidedNewT oldT origNodes = false := by decide +kernel

/-- A non-injective vertex "relabelling" is REFUSED (no left inverse at 14860). -/
theorem collision_refused :
    transportCheck n collidedNewT oldT origNodes origEdges edges permT invPermT edgesT = false :=
  transportCheck_false_of_nodes collision_nodesOk_false

/-- Edge reordering that is not a bijection: send original edge 1 where edge 0 goes. -/
def collidedPermT : NTree := permT.set 64 1 (permT.get 0)

set_option maxRecDepth 100000 in
theorem edge_collision_permOk_false :
    permOk (newF newT) collidedPermT invPermT edgesT edges.length 0 origEdges = false := by
  decide +kernel

/-- A non-bijective edge reordering is REFUSED (`invPermT` is no longer a left inverse). -/
theorem edge_collision_refused :
    transportCheck n newT oldT origNodes origEdges edges collidedPermT invPermT edgesT = false :=
  transportCheck_false_of_perm edge_collision_permOk_false

set_option maxRecDepth 100000 in
theorem sign_mismatch_permOk_false :
    permOk (newF newT) permT invPermT edgesT edges.length 0 (flipAt 0 origEdges) = false := by
  decide +kernel

/-- A sign that does not match (OPH's edge 0 flipped to +1) is REFUSED. -/
theorem sign_mismatch_refused :
    transportCheck n newT oldT origNodes (flipAt 0 origEdges) edges permT invPermT edgesT = false :=
  transportCheck_false_of_perm sign_mismatch_permOk_false

theorem wrong_vertex_set_len_false : decide (origNodes.tail.length = n) = false := by decide +kernel

/-- Drop one carrier from the node list (a bijection of the WRONG vertex set): REFUSED. -/
theorem wrong_vertex_set_refused :
    transportCheck n newT oldT origNodes.tail origEdges edges permT invPermT edgesT = false :=
  transportCheck_false_of_len wrong_vertex_set_len_false

#print axioms endingB_direct_check_false
#print axioms collision_refused
#print axioms edge_collision_refused
#print axioms sign_mismatch_refused
#print axioms wrong_vertex_set_refused

end GapScout.Oph
