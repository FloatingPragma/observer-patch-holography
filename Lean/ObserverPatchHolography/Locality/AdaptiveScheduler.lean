import ObserverPatchHolography.Locality.DependencyCone

namespace OPH.Locality

open OPH

/-!
# Adaptive-scheduler locality: dependency cones with a consultation region

`DependencyCone` bounds the influence cone of every fixed exogenous repair
word: after a word of `n` moves, the readout on a region `S` is a function
of the initial records on `ball S n`.  The B4 boundary accepts that result
as the fixed-word upper bound only; the adaptive form belongs to E2
obligation 5 (issue `#693`).  This module supplies that adaptive form.

An adaptive scheduler chooses the next firing site from the step index and
the current records rather than from a fixed word.  When the scheduler
consults only a declared region `R` (`ConsultsOnly`), the dependency cone
of a probe region `S` after `n` steps is bounded by
`ball (S ∪ R) n = ball S n ∪ ball R n`: the consultation region enters the
cone alongside the accumulated move radii (`adaptiveRun_agree_on`,
`adaptive_no_influence`).  The single-site move has declared locality
radius one (it reads one closed edge neighborhood and writes one site), and
`ball` is `n`-fold graph growth, so `ball S n ∪ ball R n` is exactly "at
most `n` radii from the probe region or from the consultation region".
Sites outside both accumulations are untouched by every run of the
scheduler.

The enlargement over the fixed-word cone is genuine.  A two-cell
countermodel exhibits an adaptive scheduler that reads a remote cell and
thereby decides whether the probe cell is repaired: the probe readouts of
two runs differ although the two initial records agree on `ball S n` for
every `n`, an influence that `applyWord_agree_on` forbids for every fixed
word of every length (`adaptivity_enlarges_cone`).  Hence the consultation
region cannot be dropped from the adaptive bound
(`consultation_region_not_droppable`).

Refinement naturality is stated as the law of a declared structure with an
identity witness.  A `ConeRefinement` declares a patch map compatible with
closed neighborhoods and a records map intertwining one adaptive step;
those declared laws force cones to map into cones
(`ConeRefinement.ball_image`), runs to commute with the readback
(`ConeRefinement.run_natural`), and the coarse cone bound to govern the
readback of the fine run (`ConeRefinement.readback_cone_bound`).  The
identity refinement inhabits the structure (`idRefinement`).

## Claim boundary

Algebraic locality of the declared system only.  No physical causality,
speed limit, clock, probability, or continuum limit enters; distance is
graph growth of the declared patch graph, and the bipartition into probe
and consultation regions is a named input.  No theorem here produces the
scheduler or the refinement maps from the OPH source: source production of
the scheduler stays open.
-/

variable {C : OPHCarrier}

/-! ## The adaptive scheduler and its runs -/

/-- An adaptive scheduler: the next firing site is a function of the step
index and the current records, rather than a letter of a fixed word. -/
def AdaptiveScheduler (C : OPHCarrier) : Type :=
  ℕ → Records C → Site C

/-- The scheduler consults only the region `R`: its choice at every step is
a function of the records restricted to `R`. -/
def ConsultsOnly (σ : AdaptiveScheduler C) (R : Set C.Patch) : Prop :=
  ∀ (n : ℕ) (x y : Records C), (∀ p ∈ R, x p = y p) → σ n x = σ n y

/-- `n` steps of scheduler-driven repair: at each step the scheduler reads
the current records and the accumulated step index, and the chosen site
fires one `localRepair` move. -/
noncomputable def adaptiveRun : ℕ → AdaptiveScheduler C → Records C → Records C
  | 0, _, x => x
  | n + 1, σ, x => adaptiveRun n (fun k => σ (k + 1)) (localRepair C (σ 0 x) x)

@[simp] theorem adaptiveRun_zero (σ : AdaptiveScheduler C) (x : Records C) :
    adaptiveRun 0 σ x = x := rfl

theorem adaptiveRun_succ (n : ℕ) (σ : AdaptiveScheduler C) (x : Records C) :
    adaptiveRun (n + 1) σ x
      = adaptiveRun n (fun k => σ (k + 1)) (localRepair C (σ 0 x) x) := rfl

/-- A fixed word packaged as a state-blind scheduler (with a default site
past the end of the word). -/
def wordScheduler (w : List (Site C)) (d : Site C) : AdaptiveScheduler C :=
  fun n _ => w.getD n d

/-- A fixed word consults nothing: the empty region suffices. -/
theorem wordScheduler_consultsOnly_empty (w : List (Site C)) (d : Site C) :
    ConsultsOnly (wordScheduler w d) (∅ : Set C.Patch) :=
  fun _ _ _ _ => rfl

/-- On its own length, the word scheduler replays the fixed word: the
fixed-word runs of `DependencyCone` are the state-blind instances of
adaptive runs. -/
theorem adaptiveRun_wordScheduler (w : List (Site C)) (d : Site C)
    (x : Records C) :
    adaptiveRun w.length (wordScheduler w d) x = applyWord w x := by
  induction w generalizing x with
  | nil => rfl
  | cons k w' ih =>
      rw [List.length_cons, adaptiveRun_succ, applyWord_cons]
      exact ih (localRepair C k x)

/-! ## Cone algebra: regions grow independently -/

/-- A region is contained in each of its grown cones. -/
theorem subset_ball (S : Set C.Patch) (n : ℕ) : S ⊆ ball (C := C) S n := by
  induction n with
  | zero => exact fun _p hp => hp
  | succ n ih => exact fun p hp => subset_grow _ (ih hp)

/-- Growth is monotone in the region. -/
theorem grow_mono {A B : Set C.Patch} (hAB : A ⊆ B) :
    grow (C := C) A ⊆ grow (C := C) B := by
  rintro p (hp | hp)
  · exact subset_grow B (hAB hp)
  · obtain ⟨i, hi, hpi⟩ := Set.mem_iUnion₂.mp hp
    exact closedNbhd_subset_grow (hAB hi) hpi

/-- Growth distributes over unions of regions. -/
theorem grow_union (A B : Set C.Patch) :
    grow (C := C) (A ∪ B) = grow A ∪ grow B := by
  unfold grow
  rw [Set.biUnion_union, Set.union_union_union_comm]

/-- The cone of a union is the union of the cones: the adaptive bound
below is literally "`n` radii from the probe region or `n` radii from the
consultation region". -/
theorem ball_union (A B : Set C.Patch) (n : ℕ) :
    ball (C := C) (A ∪ B) n = ball A n ∪ ball B n := by
  induction n with
  | zero => rfl
  | succ n ih =>
      show grow (ball (C := C) (A ∪ B) n) = _
      rw [ih, grow_union]
      rfl

/-! ## The adaptive dependency-cone bound -/

/-- **Adaptive dependency-cone bound.** For a scheduler consulting only
`R`, agreement of two initial records on `ball (S ∪ R) n` forces the two
`n`-step adaptive runs to agree on `S ∪ R`.  The consultation region is the
content over the fixed-word case: initial agreement on the grown cone of
`R` keeps the two runs choosing the same site at every step, and initial
agreement on the grown cone of `S` then transports through those common
moves exactly as in the fixed-word bound. -/
theorem adaptiveRun_agree_on (σ : AdaptiveScheduler C) {R : Set C.Patch}
    (hσ : ConsultsOnly σ R) (S : Set C.Patch) {n : ℕ} {x y : Records C}
    (hagr : ∀ p ∈ ball (C := C) (S ∪ R) n, x p = y p) :
    ∀ j ∈ S ∪ R, adaptiveRun n σ x j = adaptiveRun n σ y j := by
  induction n generalizing σ x y with
  | zero => exact fun j hj => hagr j hj
  | succ n ih =>
      have hsched : σ 0 x = σ 0 y :=
        hσ 0 x y fun p hp =>
          hagr p (subset_ball (S ∪ R) (n + 1) (Set.mem_union_right S hp))
      intro j hj
      rw [adaptiveRun_succ, adaptiveRun_succ, ← hsched]
      exact ih (fun k => σ (k + 1)) (fun m a b h => hσ (m + 1) a b h)
        (localRepair_agree_on_grow (σ 0 x) (ball (S ∪ R) n)
          (by simpa only [ball] using hagr))
        j hj

/-- **No adaptive influence outside both accumulations.** Changing the
initial record at a site outside the `n`-fold grown cone of the probe
region and outside the `n`-fold grown cone of the consultation region
cannot change the probe readout of any `n`-step run of a scheduler
consulting only `R`. -/
theorem adaptive_no_influence (σ : AdaptiveScheduler C) {R : Set C.Patch}
    (hσ : ConsultsOnly σ R) (S : Set C.Patch) (n : ℕ) (x : Records C)
    {p : C.Patch} (hpS : p ∉ ball (C := C) S n)
    (hpR : p ∉ ball (C := C) R n) (s : C.State p) :
    ∀ j ∈ S, adaptiveRun n σ (Function.update x p s) j
      = adaptiveRun n σ x j := by
  have hp : p ∉ ball (C := C) (S ∪ R) n := by
    rw [ball_union]
    rintro (h | h)
    · exact hpS h
    · exact hpR h
  intro j hj
  exact adaptiveRun_agree_on σ hσ S
    (fun q hq =>
      Function.update_of_ne (fun h => hp (by rw [← h]; exact hq)) s x)
    j (Set.mem_union_left R hj)

/-! ## Countermodel: adaptivity genuinely enlarges the cone

A two-cell system.  The probe cell (`false`) carries one self-loop overlap
that is satisfied exactly when the probe records `true`; the remote cell
(`true`) carries no edge at all, so the probe's cone `ball {false} n` is
`{false}` at every depth and the remote cell lies outside it.  The
scheduler reads the remote cell and fires at the probe exactly when the
remote flag is down.  Two initial records agreeing on the probe cone and
differing only at the remote cell then produce different probe readouts
after one adaptive step, an influence that the fixed-word bound
`applyWord_agree_on` forbids for every word of every length. -/

namespace TwoCell

/-- Two cells, one self-loop overlap at the probe cell `false` demanding
the probe record `true`; the remote cell `true` is edge-free. -/
def twoCellCarrier : OPHCarrier where
  Patch := Bool
  State := fun _ => Bool
  Edge := Unit
  src := fun _ => false
  tgt := fun _ => false
  Iface := fun _ => Bool
  projSrc := fun _ s => s
  projTgt := fun _ _ => true
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases h : a = b
    · rw [if_pos h]
      exact ⟨fun _ => h, fun _ => rfl⟩
    · rw [if_neg h]
      exact ⟨fun h1 => absurd h1 one_ne_zero, fun h2 => absurd h2 h⟩

/-- The probe's closed neighborhood is the probe alone: the self-loop has
both endpoints there. -/
theorem closedNbhd_probe :
    closedNbhd (C := twoCellCarrier) false = {false} := by
  ext j
  simp only [closedNbhd, Set.mem_setOf_eq]
  constructor
  · rintro (rfl | ⟨e, _, h | h⟩)
    · rfl
    · exact h.symm
    · exact h.symm
  · rintro rfl
    exact Or.inl rfl

/-- The probe's cone is the probe alone at every depth: the remote cell is
outside every fixed-word dependency cone of the probe. -/
theorem ball_probe (n : ℕ) :
    ball (C := twoCellCarrier) {false} n = {false} := by
  induction n with
  | zero => rfl
  | succ n ih =>
      show grow (ball (C := twoCellCarrier) {false} n) = {false}
      rw [ih]
      ext p
      constructor
      · rintro (hp | hp)
        · exact hp
        · obtain ⟨i, hi, hpi⟩ := Set.mem_iUnion₂.mp hp
          have hi' : i = false := hi
          subst hi'
          rw [closedNbhd_probe] at hpi
          exact hpi
      · exact fun hp => Or.inl hp

theorem remote_outside_probe_cone (n : ℕ) :
    (true : Bool) ∉ ball (C := twoCellCarrier) {false} n := by
  rw [ball_probe]
  intro h
  exact Bool.noConfusion (Set.mem_singleton_iff.mp h)

/-- A single-site replacement repairs the probe's overlaps exactly when it
is the demanded flag `true`. -/
theorem solvesAt_probe_iff (x : Records twoCellCarrier) (s : Bool) :
    SolvesAt (C := twoCellCarrier) false x s ↔ s = true := by
  constructor
  · intro h
    have hcons := h () (Or.inl rfl)
    have hval : Function.update x false s false = true := hcons
    rwa [Function.update_self] at hval
  · rintro rfl
    intro e _
    show Function.update x false true false = true
    rw [Function.update_self]

/-- A broken probe triggers its own repair move. -/
theorem trigger_probe_of_broken (x : Records twoCellCarrier)
    (hx : x false = false) :
    LocalTrigger (C := twoCellCarrier) false x := by
  refine ⟨(), Or.inl rfl, ?_⟩
  show ¬ (x false = true)
  rw [hx]
  exact Bool.false_ne_true

/-- Firing at a broken probe installs the demanded flag. -/
theorem localRepair_probe_broken (x : Records twoCellCarrier)
    (hx : x false = false) :
    localRepair twoCellCarrier false x false = true := by
  have h : LocalTrigger (C := twoCellCarrier) false x
      ∧ LocallySolvable (C := twoCellCarrier) false x :=
    ⟨trigger_probe_of_broken x hx, ⟨true, (solvesAt_probe_iff x true).mpr rfl⟩⟩
  rw [localRepair_of_fire twoCellCarrier false x h, Function.update_self]
  exact (solvesAt_probe_iff x _).mp (Classical.choose_spec h.2)

/-- The remote cell has no incident edge, so its move is the identity. -/
theorem localRepair_remote_quiescent (x : Records twoCellCarrier) :
    localRepair twoCellCarrier true x = x := by
  apply localRepair_of_quiescent
  rintro ⟨⟨e, hinc, _⟩, _⟩
  rcases hinc with h | h
  · exact Bool.noConfusion h
  · exact Bool.noConfusion h

/-- The adaptive scheduler of the countermodel: it reads the remote flag
and fires at the remote cell when the flag is up, at the probe when the
flag is down. -/
def remoteReader : AdaptiveScheduler twoCellCarrier :=
  fun _ z => z true

/-- The scheduler consults only the remote cell. -/
theorem remoteReader_consultsOnly :
    ConsultsOnly remoteReader ({true} : Set Bool) :=
  fun _ _ _ h => h true rfl

/-- Both start records leave the probe broken; they differ only in the
remote flag. -/
def flagDown : Records twoCellCarrier := fun _ => false

def flagUp : Records twoCellCarrier := fun i => i

theorem flagDown_ne_flagUp : flagDown ≠ flagUp := by
  intro h
  exact Bool.false_ne_true (congrFun h true)

/-- The two start records agree on the probe's cone at every depth. -/
theorem agree_on_probe_cone (n : ℕ) :
    ∀ p ∈ ball (C := twoCellCarrier) {false} n, flagDown p = flagUp p := by
  intro p hp
  rw [ball_probe] at hp
  have hp' : p = false := hp
  subst hp'
  rfl

/-- One adaptive step: with the remote flag down the scheduler fires the
probe and repairs it; with the remote flag up it fires the quiescent
remote cell and the probe stays broken.  The remote flag decides the probe
readout. -/
theorem remote_flag_moves_probe :
    adaptiveRun 1 remoteReader flagDown false = true
      ∧ adaptiveRun 1 remoteReader flagUp false = false := by
  constructor
  · rw [adaptiveRun_succ, adaptiveRun_zero]
    exact localRepair_probe_broken flagDown rfl
  · rw [adaptiveRun_succ, adaptiveRun_zero]
    have hq : localRepair twoCellCarrier (remoteReader 0 flagUp) flagUp
        = flagUp := localRepair_remote_quiescent flagUp
    rw [hq]
    rfl

/-- **Adaptivity genuinely enlarges the cone.** Every fixed word of every
length leaves the two probe readouts equal, because the start records agree
on the probe's fixed-word cone; the one-step adaptive runs of the
remote-reading scheduler make them differ.  The influence travels through
the scheduler's consultation of a cell outside every `ball {false} n`. -/
theorem adaptivity_enlarges_cone :
    (∀ w : List (Site twoCellCarrier),
        applyWord w flagDown false = applyWord w flagUp false)
      ∧ adaptiveRun 1 remoteReader flagDown false
          ≠ adaptiveRun 1 remoteReader flagUp false := by
  constructor
  · intro w
    exact applyWord_agree_on w {false} (agree_on_probe_cone w.length) false rfl
  · rw [remote_flag_moves_probe.1, remote_flag_moves_probe.2]
    exact fun h => Bool.noConfusion h

end TwoCell

/-- **The consultation region cannot be dropped.** There is a carrier, a
scheduler consulting a declared region, a probe region, a depth, and two
initial records agreeing on the probe's own `n`-fold cone whose adaptive
runs differ on the probe.  The bound of `adaptiveRun_agree_on` is therefore
tight in its shape: agreement on `ball S n` alone, the fixed-word cone,
does not bound adaptive runs, and the `ball R n` term carries content. -/
theorem consultation_region_not_droppable :
    ∃ (Cw : OPHCarrier) (σ : AdaptiveScheduler Cw) (R S : Set Cw.Patch)
      (n : ℕ) (x y : Records Cw),
      ConsultsOnly σ R
        ∧ (∀ p ∈ ball (C := Cw) S n, x p = y p)
        ∧ ∃ j ∈ S, adaptiveRun n σ x j ≠ adaptiveRun n σ y j := by
  refine ⟨TwoCell.twoCellCarrier, TwoCell.remoteReader, {true}, {false}, 1,
    TwoCell.flagDown, TwoCell.flagUp, ?_, ?_, ?_⟩
  · exact TwoCell.remoteReader_consultsOnly
  · exact TwoCell.agree_on_probe_cone 1
  · exact ⟨false, rfl, TwoCell.adaptivity_enlarges_cone.2⟩

/-! ## Refinement naturality as a declared structure -/

/-- A declared refinement between two adaptive repair systems: a patch map
and a records map from the fine carrier to the coarse carrier, one
scheduler on each side, and the two laws that make the cone bound commute.
`nbhd_law` declares that the patch map respects the moves' locality radii
(closed neighborhoods map into closed neighborhoods); `step_law` declares
that the records map intertwines one scheduled repair step at every step
index.  Both are declared data of the refinement, matching the claim
boundary: naturality is a law of the declared system, and no theorem
produces these maps from the OPH source. -/
structure ConeRefinement (Cf Cc : OPHCarrier) where
  /-- Fine patches to coarse patches. -/
  patchMap : Cf.Patch → Cc.Patch
  /-- Fine records to coarse records (the declared readback). -/
  recordsMap : Records Cf → Records Cc
  /-- The fine-side scheduler. -/
  fineSched : AdaptiveScheduler Cf
  /-- The coarse-side scheduler. -/
  coarseSched : AdaptiveScheduler Cc
  /-- Declared radius compatibility: closed neighborhoods map into closed
  neighborhoods. -/
  nbhd_law : ∀ i j : Cf.Patch, j ∈ closedNbhd (C := Cf) i →
    patchMap j ∈ closedNbhd (C := Cc) (patchMap i)
  /-- Declared step naturality: the readback of one fine scheduled step is
  one coarse scheduled step of the readback, at every step index. -/
  step_law : ∀ (n : ℕ) (x : Records Cf),
    recordsMap (localRepair Cf (fineSched n x) x)
      = localRepair Cc (coarseSched n (recordsMap x)) (recordsMap x)

namespace ConeRefinement

variable {Cf Cc : OPHCarrier}

/-- One growth step commutes with the patch map: the image of a grown
region lies in the grown image. -/
theorem grow_image (F : ConeRefinement Cf Cc) (S : Set Cf.Patch) :
    F.patchMap '' grow (C := Cf) S ⊆ grow (C := Cc) (F.patchMap '' S) := by
  rintro _q ⟨j, hj, rfl⟩
  rcases hj with hj | hj
  · exact subset_grow _ ⟨j, hj, rfl⟩
  · obtain ⟨i, hi, hji⟩ := Set.mem_iUnion₂.mp hj
    have hk : F.patchMap i ∈ F.patchMap '' S := ⟨i, hi, rfl⟩
    exact closedNbhd_subset_grow hk (F.nbhd_law i j hji)

/-- **The cone commutes with the refinement map.** The image of the fine
`n`-fold cone lies in the coarse `n`-fold cone of the image region, for
every depth. -/
theorem ball_image (F : ConeRefinement Cf Cc) (S : Set Cf.Patch) (n : ℕ) :
    F.patchMap '' ball (C := Cf) S n
      ⊆ ball (C := Cc) (F.patchMap '' S) n := by
  induction n with
  | zero => exact fun q hq => hq
  | succ n ih =>
      exact fun q hq =>
        grow_mono ih (F.grow_image (ball (C := Cf) S n) hq)

/-- The refinement with both schedulers advanced by one step; carrier of
the induction behind `run_natural`. -/
def shift (F : ConeRefinement Cf Cc) : ConeRefinement Cf Cc where
  patchMap := F.patchMap
  recordsMap := F.recordsMap
  fineSched := fun k => F.fineSched (k + 1)
  coarseSched := fun k => F.coarseSched (k + 1)
  nbhd_law := F.nbhd_law
  step_law := fun n x => F.step_law (n + 1) x

/-- **Run naturality.** The declared step law propagates to every depth:
the readback of the `n`-step fine adaptive run is the `n`-step coarse
adaptive run of the readback. -/
theorem run_natural (F : ConeRefinement Cf Cc) (n : ℕ) :
    ∀ x : Records Cf,
      F.recordsMap (adaptiveRun n F.fineSched x)
        = adaptiveRun n F.coarseSched (F.recordsMap x) := by
  induction n generalizing F with
  | zero => exact fun _x => rfl
  | succ n ih =>
      intro x
      rw [adaptiveRun_succ, adaptiveRun_succ, ← F.step_law 0 x]
      exact ih F.shift (localRepair Cf (F.fineSched 0 x) x)

/-- **The cone bound commutes with the refinement.** When the coarse
scheduler consults only `R` and two fine records have readbacks agreeing on
the coarse cone `ball (S ∪ R) n`, the readbacks of their `n`-step fine
adaptive runs agree on `S`.  The fine dynamics is bounded through the
declared refinement by the coarse adaptive cone. -/
theorem readback_cone_bound (F : ConeRefinement Cf Cc)
    {R : Set Cc.Patch} (hσ : ConsultsOnly F.coarseSched R)
    (S : Set Cc.Patch) (n : ℕ) {x y : Records Cf}
    (hagr : ∀ p ∈ ball (C := Cc) (S ∪ R) n,
      F.recordsMap x p = F.recordsMap y p) :
    ∀ j ∈ S, F.recordsMap (adaptiveRun n F.fineSched x) j
      = F.recordsMap (adaptiveRun n F.fineSched y) j := by
  intro j hj
  rw [F.run_natural n x, F.run_natural n y]
  exact adaptiveRun_agree_on F.coarseSched hσ S hagr j
    (Set.mem_union_left R hj)

end ConeRefinement

/-- The identity witness: every carrier with any scheduler refines itself
by identity maps, and both structure laws hold definitionally. -/
def idRefinement (C : OPHCarrier) (σ : AdaptiveScheduler C) :
    ConeRefinement C C where
  patchMap := id
  recordsMap := id
  fineSched := σ
  coarseSched := σ
  nbhd_law := fun _i _j h => h
  step_law := fun _n _x => rfl

/-- On the identity witness, run naturality is the literal identity. -/
theorem idRefinement_run_natural (C : OPHCarrier)
    (σ : AdaptiveScheduler C) (n : ℕ) (x : Records C) :
    (idRefinement C σ).recordsMap (adaptiveRun n (idRefinement C σ).fineSched x)
      = adaptiveRun n σ x :=
  (idRefinement C σ).run_natural n x

end OPH.Locality

#print axioms OPH.Locality.adaptiveRun_agree_on
#print axioms OPH.Locality.adaptive_no_influence
#print axioms OPH.Locality.adaptiveRun_wordScheduler
#print axioms OPH.Locality.TwoCell.adaptivity_enlarges_cone
#print axioms OPH.Locality.consultation_region_not_droppable
#print axioms OPH.Locality.ConeRefinement.ball_image
#print axioms OPH.Locality.ConeRefinement.run_natural
#print axioms OPH.Locality.ConeRefinement.readback_cone_bound
