# Row 6: theorems behind spin ℏ/2 (field-carried J on the electron state)

Producer TB, 2026-11-14. Read-only collection under `openwave/docs/THEOREM_BRIEF.md`. RER = `reverse-engineering-reality/`.
Criterion: "field-carried J = ℏ/2 on the electron state".

## Same-carrier spin, charge and field transport

[The CAR/Gauss join](../shared/spin_gauge_join/PROOF.md) now puts local field spin, signed charge, integer links and coherent loop observations on one carrier. All 332,640 coherent and 665,280 coin columns pass; full overflow and admission complements remain. The independent two-path reader realizes the finite Clifford-square precursor and retains its nonzero zero-charge geometry/memory term. [Root replay](../shared/spin_gauge_reader_amendment/REVIEW.md) passes through a separately registered receipt-schema correction; the original failed comparison is preserved. The neutral pair carries two spin-half factors, so its total spin is not identified with one electron. Physical spatial-rotation/electron attachment and magnetic dipole, g and Larmor remain open.

## 1. Chain

1. The sixty committed twelve-port rotations form a group isomorphic to PSL(2,5); composing with the central
   quotient gives a surjection SL(2,5) → PortGroup with kernel exactly {±I}.
   `OPH.SL2F5PortCover.pslToPort_injective`, `slToPort_surjective`, `slToPort_ker`, `slToPort_neg_one`,
   `cover_square`, RER `Lean/Screen/SL2F5PortCover.lean`. Kind: `lean`.
2. A5 has two distinct involutions, so it embeds in no group with a unique involution; in particular A5 ⊄ SU(2),
   and a two-dimensional spin action needs the double cover. `OPHInner.no_embedding_of_unique_involution`,
   `OPHInner.two_involutions_of_embedding`, RER `Lean/Screen/A5OPH.lean` 844-861.
   Kind: `lean`.
3. The face-field representations factor through the quotient: −1 ∈ SL(2,5) acts trivially on both rank-three
   golden pieces. `OPH.GoldenSectorPSL2F5Representations.goldenSL_center_trivial`, `goldenSL_neg_one_trivial`,
   `goldenPlusSL_neg_one_trivial`, `goldenMinusSL_neg_one_trivial`, RER
   `Lean/Screen/GoldenSectorPSL2F5Representations.lean`. Kind: `lean` (face records carry integer J only).
4. The sixty proper implementers lift to SU(2) with lift group of order 120 and a unique involution; the lift is
   non-split; an explicit isomorphism with SL(2,5) is checked on all 120² products.
   RER `code/a5_closure/super_tannakian_matter_lift_certificate.md` § 2.2 "PORT-SPIN-LIFT",
   `code/a5_closure/sl2f5_port_spin_bridge_certificate.py`. Kind: `exec` (exact over Q(√5); target σ·v supplied).
5. On the five-level code V = C³ ⊕ C², the central involution z of the rotation cover acts as
   diag(1,1,1,−1,−1) = e^{iπ·6Y}; on the Gauss sector the 2π rotation equals (−1)^N; exactly one grading
   agrees with it. RER `flagship/from_observer_consensus_to_standard_physics.tex` 2340-2346 and
   `thm:defect-classes` (2307-2318). Kind: `proof` (paper), on the declared gauge carrier (PR-59).
6. Read through the area counter modulo 20 (decoder A+k), the slowest non-stationary mode of the
   (port, counter) record is a doublet, m = +0.500003 ± 0.000011 and −0.499994 ± 0.000009, rate 0.987554
   (next 0.976967). For 300 of 300 arbitrary laws the label is half-odd; the doublet leads in 119 of 300.
   `openwave-push/mf/S_rows_05_06/ROW_06/STATUS.md`. Kind: `exec` (class N+; the half-odd value is an identity
   of the decoder; which sector is slowest depends on the law).
7. None of eleven tried record quantities of value 1/2 is conserved by the native law. Same packet. Kind: `exec`.
8. Declared source: the lowest band is an isolated spin-half doublet with Casimir 3/4 (min-max, Schur).
   RER `model/openwave_infinite_photon/PROOF.md` § "Exact isolated spin-half band". Kind: `proof` (installed source).
9. In the charge −1 sector the lowest band is exactly two spin-half doublets of the light Dirac block; on the
   prepared line native Q = −1 and the full dressed matter-plus-Maxwell Casimir is 3/4.
   RER `model/openwave_native_quantum_charge/PROOF.md` § "Complete charged band and observations". Kind: `proof`
   (installed source; Clifford matrices are architectural inputs).
10. Paper statement: `paper/tex_fragments/OBSERVED_STRUCTURE_CONNECTIONS.tex` 61-70 ("physical species
   identification is an additional comparison"). Kind: `proof` (paper).

Links 1-4 are finite group facts of the carrier. Link 6 is on the native law read through a supplied counter.
Links8–9 and the earlier reversible integer-field candidate (−0.5374 ℏ at the event offset, `openwave-push/rowfull/ROW_06`, classF) are installed-source results. The later CAR-spin-current construction below supplies conserved and locally transported spin on an actual declared CAR extension of the native source.

## 2. Premises

- unregistered: area counter modulo 20 and the lifted rotations (A+k). Realized only on an installed source (reader).
- unregistered: ℏ per full turn (unit of the counter). Open; no record fixes it.
- unregistered: the reading "electron" of the doublet. Open (PR-15 for the energy anchor).
- PR-20 uniform seam schedule. Constructed on declared inputs.
- PR-47 spacetime Spin and fermion-chirality attachment. Open (declared input): link 4 uses σ·v as the target.
- PR-59 declared five-mode exterior grammar (link 5). Open (declared input).
- PR-54 source gauge field (field-carried J needs a field record). Open (obstructed).

## 3. Historical gaps of the earlier counter candidates

1. **Law-selected sector, not reader-selected.** Theorem to prove (finite, exact): on the 240-state chain
   (12 ports × Z20 counter) of the committed fair-coin law, the decay factors γ_k of the counter characters
   exp(2πi k a/20) satisfy γ_1 = γ_19 > γ_k for every k ∉ {0, 1, 19}, i.e. the character-1 doublet is the unique
   slowest non-stationary mode, with exact value 0.987554234. Discharge: **theorem** (exact algebraic
   eigenvalue computation over the 240-state kernel). It turns the reader's choice of k = 1 into a property of
   the law; it does not make the counter a record of the law.
2. **The counter as a register of the law.** On candidate A the counter is not written by any event.
   Discharge: **owner decision** (pin a law whose face events write the counter, candidate J) followed by
   **execution** of the same doublet measurement on that candidate.
3. **A conserved J of value 1/2.** Link 7 finds none. Theorem to prove: on the pinned candidate, the
   counter-character-1 projection commutes with the event kernel on the doublet sector (an exact invariant
   subspace with J_z = ±1/2). Discharge: **theorem** if the candidate law commutes; else **execution** shows absence.
4. **Field-carried.** Read angular momentum carried by the field of the same electron object. A remote face register was an earlier candidate-specific proposal, not a universal requirement; a local CAR field observable can qualify if its physical rotation and electron identification are established.

## 4. Shared premises

A+k decoder: rows 5, 6, 7. PR-47: rows 5, 6, 7. PR-20: rows 5, 6, 7, 21. PR-59: rows 6, 8, 9. PR-54: rows 5, 6, 8, 9.


## Follow-up correction and exact certificate (Codex CRT, 2026-11-14)

See [the shared proof and certificate](../shared/rows05_07_record_spin/PROOF.md).
The full 240-state A+k chain has leading nonstationary multiplicity four, not
two, at 0.987554234029336661. This is now certified by its exact characteristic
polynomial and an independent trace/transitivity calculation. The previous
sector-specific doublet description must not be read as a full-chain claim.

A+k is irreducible, so its pathwise conserved state functions are constant.
Commutation of a mode projector with K is therefore insufficient to establish
a conserved physical J. The history quotient theorem includes adaptive
controllers only when all memory, interventions and references factor through
the port quotient. It does not exclude a different, path-sensitive operation.

The literal row-7 static criterion does not require a relative central phase or
exchange statistics; those are stronger questions, not additional acceptance gates.


### Complex-doublet clarification

The [separate exact addendum](../shared/rows05_07_complex_structure/PROOF.md)
shows that the full leading real four-dimensional space carries the invariant
complex structure `J=(T-T^-1)/(2 sin(pi/10))`. Its complex dimension is two, and
`T^10=-I` there. Exact first/conjugate-character polynomials each have the top
root with multiplicity two. Thus real dimension four is compatible with one
complex doublet and does not refute a spinor interpretation. The supplied
counter connection, physical-rotation interpretation and field-carried physical
angular momentum remain open.

## Native CAR field-spin/current join

The [new exact source proof](../shared/magnetic_source_join/PROOF.md) uses the
existing full-domain two-flavour CAR repair, including all collision flags
and fermionic crossing signs. Every branch commutes with total field spin;
the same recorded hop transports local spin through an exact continuity
equation. The one-particle sector has Casimir3/4 and projection±1/2. A spatially
uniform polarized marginal is stationary for each fixed-seam channel.

The exhaustive audit passed153,216 local continuity checks. Full native source
dilations and8192-address pointer laws read the spin. A source-compiled flavour
swap preserves the charge-record effects but reverses spin readings, providing
a same-object changed-law control. Exact integer identities and binary64 source
diagnostics are separately stated. The supplied Spin2 representation has the
existing finite port-spin lift; selection of that architecture is explicit.

This is a field-carried, locally transported half-unit on the declared source.
The remaining physical join is the electron object and its spatial-rotation
interpretation, including units. Full electromagnetism, Poynting angular momentum
or a compulsory remote-face detector are not extra clauses of row6. Row5 still
needs its own dipole, g and Larmor construction. Root's [same-group review](../shared/source_joins_root_review/REVIEW.md)
is not G6 or a platform award.

## Same-clock current and loop continuation

The [current-query clock](../shared/native_current_clock/PROOF.md) has an exact complete stationary electric-continuity identity, with coherent current quadratures distinguished from realized transport. Its separately certified geometric error replaces an uncovered conditional allowance without erasing the original receipt. The [full native loop extension](../shared/native_field_cycle/PROOF.md) gives nonzero selected stationary circulation with complete return/nonreturn balance. All loop microcurrents still contribute to total transport; spin-blind circulation prevents a magnetic-moment inference. No supplied Pauli square is identified with a measured field.

## Complete first-return field-cycle instrument

The [first-return proof](../shared/native_field_cycle/first_return/PROOF.md) uses fixed-length reversible scans and full counter uncomputation to preserve references across distinct native return lengths. Every final charge position returns and field increments have zero divergence, with all intermediate charged routes retained. The same-clock stationary component balance is exact; the spin-blind spectrum and larger circulation prevent a magnetic g inference. The [precise Maxwell interaction interface](../shared/native_field_cycle/first_return/NEXT_INTERFACE.md) distinguishes native occupation-field positives from an unproved product-source interaction.
