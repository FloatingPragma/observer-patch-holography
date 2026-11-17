# Row 5: theorems behind the magnetic moment (g ≈ 2, magnetostatics)

Producer TB, 2026-11-14. Read-only collection under `openwave/docs/THEOREM_BRIEF.md`. RER = `reverse-engineering-reality/`.
Criterion: "a stationary state sources a static dipole B from internal circulation, g ≈ 2; dipole Larmor response".

## Same-carrier spin, charge and field transport

[The CAR/Gauss join](../shared/spin_gauge_join/PROOF.md) now puts local field spin, signed charge, integer links and coherent loop observations on one carrier. All 332,640 coherent and 665,280 coin columns pass; full overflow and admission complements remain. The independent two-path reader realizes the finite Clifford-square precursor and retains its nonzero zero-charge geometry/memory term. [Root replay](../shared/spin_gauge_reader_amendment/REVIEW.md) passes through a separately registered receipt-schema correction; the original failed comparison is preserved. The neutral pair carries two spin-half factors, so its total spin is not identified with one electron. Physical spatial-rotation/electron attachment and magnetic dipole, g and Larmor remain open.

## 1. Chain

Native process (one unit on the twelve-port carrier, uniform seam schedule, fair coin):

1. Every seam event obeys the discrete continuity law x' − x = −b_e j; conditional on the past and the selected
   seam, E[j] = d/2 and Var(j) = (|d| mod 2)/4. RER `paper/executable_observer_construction.tex`,
   `prop:fields` (i), (ii), lines 4096-4101. Kind: `proof` + `exec` (all 78 states).
2. Mean load relaxes as μ_{k+1} − μ_k = −Lμ_k/60 by KX = XT on the 78-state kernel; slow multiplier
   q = (55+√5)/60, no oscillating mode. Same proposition (iii), (iv). Kind: `proof` + `exec`.
3. On the committed law the mean circulation record and the moment along an axis are exactly zero.
   `openwave-push/mf/S_rows_05_06/ROW_05/STATUS.md`, gate B2. Kind: `exec` (class N, exact 0).
4. Any record that is a function of port labels carries an SL(2,5) action through the port quotient, and the
   central −I acts as the identity rotation: `OPH.SL2F5PortCover.slToPort_ker` (kernel = center {±I}),
   `slToPort_neg_one`, `slToPort_surjective`. RER `Lean/Screen/SL2F5PortCover.lean`. Kind: `lean`.
   On the face-field golden pieces −1 acts trivially: `OPH.GoldenSectorPSL2F5Representations.goldenSL_neg_one_trivial`,
   RER `Lean/Screen/GoldenSectorPSL2F5Representations.lean`. Kind: `lean`.
5. Through the area counter modulo 20 (decoder A+k) and a swept-area weight exp(iBa) (decoder A+k+a), the
   record has a dipole-form circulation, moment 0.39260 ± 0.0012 (exact 0.391857), and
   g = (2/3)/cos(β/2) = 0.783714; response 4e-15 when seams enclose no face (orbital); under flux the record
   aligns and does not precess. Same packet, gates B3-B7. Kind: `exec` (class N+, decoder supplied).
6. Declared source (complete matter shells n ≤ 2, mass 16, field coefficient 1/8, trap, full photon Fock
   space): the two lowest states form an isolated spin-half band, Casimir 3/4, exactly degenerate (Schur).
   RER `model/openwave_infinite_photon/PROOF.md` § "Exact isolated spin-half band". Kind: `proof` (installed source).
7. On that source PXP = μσ_z on the ground band (rotation covariance, π rotation), g_static = 4E0|μ|, and the
   enclosed static coefficient is g ∈ [1.99827530667, 1.99827540121]. RER `model/openwave_infinite_current/PROOF.md`
   §§ "Coordinate-matrix derivative bound" (‖K'‖ ≤ √5, ‖K''‖ ≤ 141/400), "Exact line selection and enclosed current".
   Kind: `proof` + `exec` (interval certificate), installed source.
8. Finite-time precession: ‖U_b(t)P − exp[−it(E0+bμσ_z)]P‖ ≤ 2√2 v/(γ−2v) + |t|D with v = √5|b|.
   RER `model/openwave_infinite_current/PROOF.md` § "Full-reference finite-time precession";
   `model/openwave_quantum_larmor/PROOF.md` § "Uniform finite-time comparison". Kind: `proof`, installed source.
9. Paper statement of the result: `paper/tex_fragments/OBSERVED_STRUCTURE_CONNECTIONS.tex` 50-59 and
   `NATIVE_COMPLETIONS_AND_COMPOSITE_RESPONSE.tex` 41-53 ("the deviation from two does not reproduce the
   observed electron anomaly; field normalization, trap and regulator are declared inputs"). Kind: `proof` (paper).

Links 1-3 are on the native process. Links 6-8 are on an installed Dirac/Clifford source (class F). No theorem
or execution connects link 5 (g = 0.78 on records) to link 7 (g ≈ 2 on the installed source).

## 2. Premises

- unregistered: area counter modulo 20 with lifted rotations (decoder A+k). Realized only on an installed source
  (a reader's register, outside the pinned law).
- unregistered: swept-area weight exp(iBa) as the applied flux (decoder A+k+a). Realized only on an installed reader.
- PR-20 equal source-counting measure on the sixty directed seams (uniform schedule, weight 1/60). Constructed on
  declared inputs (register status); the 1/60 of link 2 depends on it.
- PR-47 spacetime Spin and fermion-chirality attachment. Open (declared input): the Clifford matrices and the
  spin-half band of link 6 are installed.
- PR-54 source gauge-field, current and action attachment. Open (obstructed): the photon field coupled through
  K(b,x) = sqrt(16² + (B(x)+bE_z)²) is declared.
- unregistered: mass 16, field coefficient 1/8, trap, radial form factor, matter-shell cap n ≤ 2 of
  `openwave_quantum_rotational`. Realized only on an installed source.
- PR-15 clock and energy calibration anchors. Open: no Tesla conversion, no anomaly comparison in units.
- unregistered: rule that ties a frame precession rate ω_s to an orbital rate ω_c of the same object
  (g = 2ω_s/ω_c). Not found in RER; open.

## 3. Historical gaps of the earlier counter candidates

1. **A spinor-valued record on the native law.** Every port-label record carries SL(2,5) with −I trivial
   (link 4). Theorem to prove (finite, a no-go that fixes the requirement): for every finite record alphabet R
   and every record map r: (port histories of the native 78-state kernel) → R that is equivariant under the
   sixty port rotations, the pulled-back SL(2,5) action on functions of R has −I acting as +1. A spin-half
   (−I = −1) response therefore needs a register outside the port histories. Discharge: **theorem** (provable
   from `slToPort_neg_one`), then **owner decision**: which register the pinned law writes (the area counter
   written by face events on candidate J, or the unread choice register of retention schedule D′).
2. **An intrinsic response.** The old A+k circulation of link5 vanishes when its tested seams enclose no face.
   This is a boundary of that packet. A nonzero moment on a spanning-tree schedule is not a universal extra
   clause of OpenWave's test. The new candidate must exhibit stationary internal circulation and its sourced
   dipole on the same object, with a changed-law control and a functioning planted dipole detector.
3. **g as a ratio of two record rates.** Required: precession rate of the frame record and orbital rate of the
   circulation record in one flux sector on one object, g = 2ω_s/ω_c, gate [1.9, 2.1]. Discharge: **execution**;
   the rule linking the two rates is an **owner decision** (not found in RER).
4. **Sourced static field.** A record held away from the unit (face register of J, far zone). **Execution**.

## 4. Shared premises

A+k decoder: rows 5, 6, 7 (also 26, 27 by the R5 review). PR-47: rows 5, 6, 7. PR-54: rows 5, 8, 9, 20, 29.
PR-15: rows 5, 10, 11. PR-20: rows 5, 6, 7, 21 (uniform schedule).


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

### Exact finite Pauli precursor

The actual `Lean/EventAlgebra/PauliMinimalCoupling.lean` supplies a missing
positive link in this register. `sigmaDot_mul_self` proves the Clifford square
for any triple in a complex algebra; `pauliSquare` proves
`(sigma dot X)^2=(X dot X)I-q(sigma dot B)` under the explicit source
commutator hypothesis. The coefficient is exactly q. With the separately
identified orbital coefficient q/(2m) and spin sigma/2, this is the algebraic
factor two. The commuting control has no spin term, and the nonzero Pauli
inhabitant discriminates its coefficient.

This theorem does not install an electron, a physical magnetic field, or a
native kinetic generator. The constructive use is to derive the kinetic and
curvature operators from the same observer process and verify the theorem's
commutator hypothesis on them, then measure the stationary dipole and Larmor
response. Merely defining B by a commutator proves an algebraic identity and
does not discharge those observable identifications. The module was read;
no new Lean build or physical moment measurement is claimed here.

### Actual circulation and coherent loop records

`model/native_circulation/README.md` constructs a nonzero fluctuating record
circulation on the continuing binary source. For a triangular cycle its exact
stationary variance is `N/66`, even though its mean is zero. In the recovered
three-dimensional frame, the axial weights obey `W B^T=0`, and the accumulated
current has covariance `(4N/99)I`. Its native branch implementation and complete
finite counter laws are explicit. A finite mean of zero therefore does not
mean that every circulation record vanishes. This positive construction alone
does not source a stationary polarized magnetic dipole.

The zero-loop result in `model/observer/native_loop_readout.py` concerns the
specific additive paired-coin contrast `j0+j1=n_a-n_b`, which is a gradient.
The same library contains coherent native-history loop comparisons and
`holonomy_reconstruction.py`:40 oriented comparisons retain320 route/coin
histories, and exact integral identities recover19 gauge coordinates from19
face phases under the stated phase/readout contract. Coherent tape access,
relative phases and dynamic-phase calibration remain explicit. These actual
finite results supply a usable curvature instrument; they do not claim that
the ordinary current records select a unique connection.

The existing `model/native_observer_program/README.md` also joins a declared
integer/Gauss extension to its actual source/current maps. Before proposing a
new such extension, consume its full implementation and distinguish its memory
and field content from the CAR spin source. Declared extensions can be tested;
joining different sources by their names is insufficient. This is a source
inspection, not a new replay or magnetic-moment measurement.
