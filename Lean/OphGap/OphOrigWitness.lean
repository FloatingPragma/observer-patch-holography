import OphGap.Slices.P4

/-!
# The gap theorem for OPH's graph WITH ITS ORIGINAL CARRIER IDS

`OphWitness.lean` certifies the BFS-relabelled graph (`check_true`). The kernel checks
`transportCheck` on OPH's own node list and edge list (`OphOrigData.lean`, freeze sha256
`a0be6fc6…`) in thirteen slices (`OphGap/Slices/*`, one kernel evaluation per module), reassembled
here by the proved `_split` lemmas; `origGram_posDef_of_transportCheck` then carries positive
definiteness back to the matrix whose columns are OPH's carrier ids and whose rows are OPH's
edges in OPH's order.

Everything is `decide +kernel`; no shortcuts or unproved placeholders.
-/

open Matrix

namespace GapScout.Oph

theorem len_nodes_take : (origNodes.take 3000).length = 3000 := by
  rw [List.length_take, origNodes_length]; decide
theorem len_nodes_drop_take : ((origNodes.drop 3000).take 3000).length = 3000 := by
  rw [List.length_take, List.length_drop, origNodes_length]; decide
theorem len_orig_take : (origEdges.take 4000).length = 4000 := by
  rw [List.length_take, origEdges_length]; decide
theorem len_orig_drop_take : ((origEdges.drop 4000).take 4000).length = 4000 := by
  rw [List.length_take, List.length_drop, origEdges_length]; decide
theorem len_orig_take3 : (origEdges.take 3000).length = 3000 := by
  rw [List.length_take, origEdges_length]; decide
theorem len_orig_drop_take3 : ((origEdges.drop 3000).take 3000).length = 3000 := by
  rw [List.length_take, List.length_drop, origEdges_length]; decide
theorem len_orig_drop2_take3 : (((origEdges.drop 3000).drop 3000).take 3000).length = 3000 := by
  rw [List.length_take, List.length_drop, List.length_drop, origEdges_length]; decide
theorem len_edges_take : (edges.take 4000).length = 4000 := by
  rw [List.length_take, edges_length]; decide
theorem len_edges_drop_take : ((edges.drop 4000).take 4000).length = 4000 := by
  rw [List.length_take, List.length_drop, edges_length]; decide

theorem part_nodes : nodesOk n newT oldT origNodes = true :=
  nodesOk_split n newT oldT 3000 origNodes nodes_1
    (nodesOk_split n newT oldT 3000 (origNodes.drop 3000) nodes_2 nodes_3)

theorem part_endpoints : endpointsOk n newT oldT origEdges = true :=
  endpointsOk_split n newT oldT 4000 origEdges endpoints_1
    (endpointsOk_split n newT oldT 4000 (origEdges.drop 4000) endpoints_2 endpoints_3)

theorem part_edgesT : edgesTOk edgesT 0 edges = true :=
  edgesTOk_split edgesT 0 4000 4000 edges len_edges_take rfl edgesT_1
    (edgesTOk_split edgesT 4000 4000 8000 (edges.drop 4000) len_edges_drop_take rfl edgesT_2 edgesT_3)

theorem part_perm : permOk (newF newT) permT invPermT edgesT edges.length 0 origEdges = true :=
  permOk_split _ _ _ _ _ 0 3000 3000 origEdges len_orig_take3 rfl perm_1
    (permOk_split _ _ _ _ _ 3000 3000 6000 (origEdges.drop 3000) len_orig_drop_take3 rfl perm_2
      (permOk_split _ _ _ _ _ 6000 3000 9000 ((origEdges.drop 3000).drop 3000) len_orig_drop2_take3 rfl
        perm_3 perm_4))

/-- THE TRANSPORT CERTIFICATE PASSES, kernel-checked (in slices): the relabelling is injective
on the 8,662 carriers, every edge endpoint is a carrier, and the edge permutation carries OPH's
edge list onto the certified list up to orientation. -/
theorem transport_true :
    transportCheck n newT oldT origNodes origEdges edges permT invPermT edgesT = true := by
  unfold transportCheck
  rw [part_len, part_sorted, part_nodes, part_endpoints, part_lenEq, part_edgesT, part_perm]
  rfl

/-- Every OPH edge joins two distinct carriers of OPH's node list. -/
theorem oph_endpoints : EndpointsIn origNodes origEdges := endpointsIn_of_check transport_true

/-- MAIN RESULT, about the ORIGINAL graph: `Dᵀ D` built from OPH's edge list on OPH's
carrier ids is positive definite (`unfold OrigGramPosDef origD` for the literal content). -/
theorem oph_original_gram_posDef : OrigGramPosDef origNodes origEdges oph_endpoints :=
  origGram_posDef_of_transportCheck check_true transport_true

/-- Every eigenvalue of that matrix is strictly positive: the gap is strictly positive. -/
theorem oph_original_eigenvalues_pos :
    OrigEigenvaluesPos origNodes origEdges oph_endpoints oph_original_gram_posDef :=
  origEigenvaluesPos_of_transportCheck check_true transport_true

/-- The original signed derivative has trivial kernel. -/
theorem oph_original_kernel_trivial : OrigKernelTrivial origNodes origEdges oph_endpoints :=
  origKernelTrivial_of_transportCheck check_true transport_true

#print axioms part_nodes
#print axioms part_endpoints
#print axioms part_edgesT
#print axioms part_perm
#print axioms transport_true
#print axioms oph_endpoints
#print axioms oph_original_gram_posDef
#print axioms oph_original_eigenvalues_pos
#print axioms oph_original_kernel_trivial
#print axioms origD_injective
#print axioms origGram_posDef_of_transport
#print axioms origGram_posDef_of_transportCheck
#print axioms mem_of_labelled
#print axioms index_equiv_of_check
#print axioms triBadPerm_refused
#print axioms sqRelabBij_refused
#print axioms sqRelabCollapse_refused
#print axioms sqToTri_refused
#print axioms triPerm_accepted
#print axioms square_not_posDef
#print axioms check_true
#print axioms tampered_refused

end GapScout.Oph
