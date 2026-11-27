import ObserverPatchHolography.EinsteinBranch.Consensus

/-!
# Append-only cohorts and the independent clock-pair outcomes

At a fixed blank preparation, appointment and sensor calibration, each of
two independent runs of one fully enrolled pair emits one of four sign
words with positive probability. Each added cohort keeps the actual old
signs, timestamps, masks and weighing tickets. The public quotient below
describes completed records on this subexperiment. A second presentation
tag permits two serializations of exactly the same public record.

Prefix restriction preserves every old row readout. The stronger
Einstein-branch arrows would require the entropy and geometry readings to
have identical fibres. This fails for the raw sampled entropy and the raw
independent geometry on this blank. A transcript accompanying both
readings would be a different output type. The module makes no claim
about deterministic ensemble entropies or continuum physics.
-/

namespace OPH.EinsteinBranch.ClockPairRecords

abbrev Word := Bool × Bool
abbrev Bank := Word × Word
abbrev Prefix (A : ℕ → Type) (r : ℕ) := (i : Fin r) → A i.val
abbrev Banks := Prefix (fun _ => Bank)

/-- Precision and sensor values belong to each retained appointment. -/
structure Sensors (bits : ℕ) where
  appointment : Fin (2^bits)
  localClock : Fin 2 → Fin (2^bits)
  beaconStamps : Fin 4 → Fin (2^bits)
  mask : Fin 2 → Bool
  rodStop : Fin (2^bits)
  boundaryArea : Fin (2^bits)
  shellMass : Fin 2 → Fin (2^bits)
  shellRadius : Fin 2 → Fin (2^bits)
  ordinaryRelease : Fin (2^bits)

structure Program where
  bits : ℕ → ℕ
  atSetting : (i : ℕ) → Sensors (bits i)

def truncate {A : ℕ → Type} {r s : ℕ} (h : r ≤ s) (x : Prefix A s) : Prefix A r :=
  fun i => x ⟨i.val, Nat.lt_of_lt_of_le i.isLt h⟩

def readRows {A B : ℕ → Type} (f : (i : ℕ) → A i → B i) {r : ℕ}
    (x : Prefix A r) : Prefix B r := fun i => f i.val (x i)

theorem truncate_readRows {A B : ℕ → Type} (f : (i : ℕ) → A i → B i)
    {r s : ℕ} (h : r ≤ s) (x : Prefix A s) :
    truncate h (readRows f x) = readRows f (truncate h x) := rfl

/-- The old bit width and its sensor settings survive appending finer cohorts. -/
def joined (P : Program) {r : ℕ} (x : Banks r) :
    Prefix (fun i => Sensors (P.bits i) × Bank) r :=
  readRows (fun i word => (P.atSetting i, word)) x

theorem truncate_joined (P : Program) {r s : ℕ} (h : r ≤ s) (x : Banks s) :
    truncate h (joined P x) = joined P (truncate h x) := rfl

/-- Clock/area, work, modular, entropy and scale decoders on any one retained
row obey this identity even when their output types differ by setting. -/
theorem truncate_instrument (P : Program) {B : ℕ → Type}
    (read : (i : ℕ) → Sensors (P.bits i) × Bank → B i)
    {r s : ℕ} (h : r ≤ s) (x : Banks s) :
    truncate h (readRows read (joined P x)) =
      readRows read (joined P (truncate h x)) := rfl

/-- Conditional append law. Normalized fresh finite records preserve each
old record's probability, including when the new protocol depends on it. -/
theorem append_marginal {New : Type} [Fintype New]
    (oldProb : ℚ) (conditional : New → ℚ) (normalized : ∑ y, conditional y = 1) :
    (∑ y, oldProb * conditional y) = oldProb := by
  rw [← Finset.mul_sum, normalized, mul_one]

/-- The protected sensors contain the engineering settings and public material
and geometry tickets. Completing an allowed independent two-shot bank
creates no inconsistency on this fixed source. -/
def recordTower (P : Program) : BareConsensusTower ℕ where
  State := fun r => Banks r × Bool
  stateFinite := fun _ => inferInstance
  Quot := Banks
  quotFinite := fun _ => inferInstance
  Boundary := Prefix (fun i => Sensors (P.bits i))
  quotient := fun _ x => x.1
  quotient_surjective := fun _ q => ⟨(q, false), rfl⟩
  mismatch := fun _ _ => 0
  step := fun _ _ _ => False
  normalForm := fun _ x => x
  consistent := fun _ _ => True
  boundary := fun _ _ i => P.atSetting i.val
  consistent_iff := by intros; simp
  step_descends := by intros; contradiction
  normalForm_consistent := by intros; trivial
  normalForm_idempotent := by intros; rfl
  boundary_normalForm := by intros; rfl
  coarse := fun h x => truncate h x
  coarse_refl := by intros; rfl
  coarse_trans := by intros; rfl
  normalForm_natural := by intros; rfl

theorem gauge_iff_public (P : Program) {r : ℕ} (x y : Banks r × Bool) :
    GaugeEquivalent (recordTower P) r x y ↔ joined P x.1 = joined P y.1 := by
  change x.1 = y.1 ↔ _
  constructor
  · exact congrArg (joined P)
  · intro h
    funext i
    exact congrArg Prod.snd (congrFun h i)

/-- The blank correlation is 1/2 and each endpoint is fair. -/
def blankProb (w : Word) : ℚ := if w.1 = w.2 then 3/8 else 1/8
def bankProb (b : Bank) : ℚ := blankProb b.1 * blankProb b.2
def repeats : Bank := ((false, false), (false, false))
def splits : Bank := ((false, false), (false, true))

theorem blankWord_normalized : ∑ w : Word, blankProb w = 1 := by
  norm_num [Fintype.sum_prod_type, Fintype.sum_bool, blankProb]
theorem actualHistories : bankProb repeats = 9/64 ∧ bankProb splits = 3/64 := by
  norm_num [bankProb, blankProb, repeats, splits]

/-- Empirical Shannon entropy of the four-word histogram from two launches:
one populated bin or two equally populated bins. -/
noncomputable def sampleEntropy (b : Bank) : ℝ :=
  if b.1 = b.2 then 0 else Real.log 2

theorem entropy_separates : sampleEntropy repeats ≠ sampleEntropy splits := by
  have pos : 0 < Real.log 2 := Real.log_pos (by norm_num)
  simpa [sampleEntropy, repeats, splits] using ne_of_lt pos

def one (b : Bank) : Banks 1 := fun _ => b

theorem sensor_boundary_fails (P : Program) : ¬ BoundaryFiber (recordTower P) := by
  apply not_boundaryFiber_of_witness (recordTower P) 1 (one repeats) (one splits)
      True.intro True.intro rfl
  intro h
  have same : repeats = splits := congrFun h 0
  exact entropy_separates (congrArg sampleEntropy same)

/-- The published typed arrows contain a genuine cycle. Equality of the
independent geometry and entropy readings must induce the same partition. -/
theorem geometry_eq_iff_entropy_eq {ι : Type} [Preorder ι] {T : BareConsensusTower ι}
    {R : CommonReadoutTower T} (arrows : TypedArrowCommutation R)
    (r : ι) (x y : T.Quot r) :
    R.geometryRead r x = R.geometryRead r y ↔ R.entropyRead r x = R.entropyRead r y := by
  obtain ⟨gm, _, _, me, eg, _⟩ := arrows
  constructor
  · intro h
    rw [← me r x, ← me r y, ← gm r x, ← gm r y, h]
  · intro h
    rw [← eg r x, ← eg r y, h]

/-- Any readout whose geometry comes from the fixed public sensors and whose
entropy contains the independently sampled histogram entropy fails EB2
on two positive-probability histories of this same blank. -/
theorem no_raw_typed_arrows (P : Program) (R : CommonReadoutTower (recordTower P))
    (fromSensors : (recordTower P).Boundary 1 → R.Geometry 1)
    (geometryMeasured : ∀ q, R.geometryRead 1 q =
        fromSensors ((recordTower P).boundary 1 q))
    (entropyValue : R.Entropy 1 → ℝ)
    (containsSample : ∀ q, entropyValue (R.entropyRead 1 q) = sampleEntropy (q 0)) :
    ¬ TypedArrowCommutation R := by
  intro arrows
  have geom : R.geometryRead 1 (one repeats) = R.geometryRead 1 (one splits) := by
    rw [geometryMeasured, geometryMeasured]
    rfl
  have ent := (geometry_eq_iff_entropy_eq arrows 1 (one repeats) (one splits)).mp geom
  have h := congrArg entropyValue ent
  rw [containsSample, containsSample] at h
  exact entropy_separates h

#print axioms truncate_instrument
#print axioms append_marginal
#print axioms gauge_iff_public
#print axioms actualHistories
#print axioms sensor_boundary_fails
#print axioms no_raw_typed_arrows

/-! ## Deterministic, completed population censuses

The protected boundary is one binned exterior pulse measurement. The independent
calorimeter and four-word census distinguish the same three prepared stocks.
An exhaustive finite census is a different instrument from an individual outcome
or a finite random sample. Every view keeps its own raw anchor and only its own
public data from each additional cohort.
-/

inductive Kind | geometry | modular | event | stress | entropy | scale

structure View (S : Type) where
  Anchor : Type
  Row : ℕ → Type
  anchor : S → Anchor
  emit : (i : ℕ) → S → Row i

def View.Value {S} (v : View S) (r : ℕ) := v.Anchor × Prefix v.Row r

def View.read {S} (v : View S) (r : ℕ) (x : S) : v.Value r :=
  (v.anchor x, fun i => v.emit i.val x)

def View.coarse {S} (v : View S) {r s : ℕ} (h : r ≤ s) (x : v.Value s) : v.Value r :=
  (x.1, truncate h x.2)

noncomputable def View.translate {S} [Nonempty S] (v w : View S)
    (r : ℕ) (x : v.Value r) : w.Value r :=
  w.read r (Function.invFun v.anchor x.1)

theorem View.translate_read {S} [Nonempty S] (v w : View S)
    (h : Function.Injective v.anchor) (r : ℕ) (s : S) :
    v.translate w r (v.read r s) = w.read r s :=
  congrArg (w.read r) (Function.leftInverse_invFun h s)

structure Catalogue (S : Type) where
  view : Kind → View S
  geometryKey : Function.Injective (view .geometry).anchor
  modularKey : Function.Injective (view .modular).anchor
  eventKey : Function.Injective (view .event).anchor
  entropyKey : Function.Injective (view .entropy).anchor

/-- A stock is a canonical coordinate for its reachable complete boundary
census. The independently emitted geometry anchor is injective on this menu. -/
def censusTower {S : Type} [Fintype S] (c : Catalogue S) : BareConsensusTower ℕ where
  State := fun _ => S × Bool
  stateFinite := fun _ => inferInstance
  Quot := fun _ => S
  quotFinite := fun _ => inferInstance
  Boundary := fun _ => (c.view .geometry).Anchor
  quotient := fun _ p => p.1
  quotient_surjective := fun _ s => ⟨(s, false), rfl⟩
  mismatch := fun _ _ => 0
  step := fun _ _ _ => False
  normalForm := fun _ s => s
  consistent := fun _ _ => True
  boundary := fun _ => (c.view .geometry).anchor
  consistent_iff := by intros; simp
  step_descends := by intros; contradiction
  normalForm_consistent := by intros; trivial
  normalForm_idempotent := by intros; rfl
  boundary_normalForm := by intros; rfl
  coarse := fun _ => id
  coarse_refl := by intros; rfl
  coarse_trans := by intros; rfl
  normalForm_natural := by intros; rfl

theorem census_gauge_iff {S : Type} [Fintype S] (c : Catalogue S)
    (r : ℕ) (p q : S × Bool) :
    GaugeEquivalent (censusTower c) r p q ↔
      ∀ k, (c.view k).read r p.1 = (c.view k).read r q.1 := by
  constructor
  · intro h k
    exact congrArg ((c.view k).read r) h
  · intro h
    exact c.geometryKey (congrArg Prod.fst (h .geometry))

noncomputable def censusReadouts {S : Type} [Fintype S] [Nonempty S]
    (c : Catalogue S) : CommonReadoutTower (censusTower c) where
  Geometry := (c.view .geometry).Value
  Modular := (c.view .modular).Value
  Event := (c.view .event).Value
  Stress := (c.view .stress).Value
  Entropy := (c.view .entropy).Value
  Scale := (c.view .scale).Value
  geometryRead := (c.view .geometry).read
  modularRead := (c.view .modular).read
  eventRead := (c.view .event).read
  stressRead := (c.view .stress).read
  entropyRead := (c.view .entropy).read
  scaleRead := (c.view .scale).read
  geometryToModular := (c.view .geometry).translate (c.view .modular)
  modularToEvent := (c.view .modular).translate (c.view .event)
  eventToStress := (c.view .event).translate (c.view .stress)
  modularToEntropy := (c.view .modular).translate (c.view .entropy)
  entropyToGeometry := (c.view .entropy).translate (c.view .geometry)
  geometryToScale := (c.view .geometry).translate (c.view .scale)
  geometryCoarse := (c.view .geometry).coarse
  modularCoarse := (c.view .modular).coarse
  eventCoarse := (c.view .event).coarse
  stressCoarse := (c.view .stress).coarse
  entropyCoarse := (c.view .entropy).coarse
  scaleCoarse := (c.view .scale).coarse

/-- Each arrow decodes the independent stock key and emits the named meter's
own retained prefix. On unenrolled records the total extension is arbitrary.
None of the continuum, scale-limit or balance hypotheses is proved here. -/
noncomputable def admissibleCensus {S : Type} [Fintype S] [Nonempty S]
    (c : Catalogue S) : EinsteinAdmissibleTower ℕ where
  consensus := censusTower c
  readouts := censusReadouts c
  boundaryFiber := fun _ _ _ _ _ h => c.geometryKey h
  typedArrows := ⟨View.translate_read _ _ c.geometryKey, View.translate_read _ _ c.modularKey,
    View.translate_read _ _ c.eventKey, View.translate_read _ _ c.modularKey,
    View.translate_read _ _ c.entropyKey, View.translate_read _ _ c.geometryKey⟩
  refinementNaturality := by
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;> intros <;> rfl

inductive Stock | blank | ordinary | held deriving DecidableEq, Fintype, Inhabited

/-- Exterior receiver-spacing divided into ten-billion-tick bins. -/
def clockBin : Stock → ℕ
  | .blank => 1759 | .ordinary => 1798 | .held => 1831

/-- Counts in the 00 cell from the simultaneous, always-enrolled calibration
pair in clock_pairs/census.py. The 11 cell has the same count. -/
def count00 : Stock → ℕ
  | .blank => 8979466059761067529443439612208842872766049544044544
  | .ordinary => 8979466059084433769322628515964704813945367275665563
  | .held => 8979466058592336482626388080561005611989660595829918

def censusPeriod : ℕ := 23945242826029513411849172299223580994042798784118784

def mass : Stock → ℚ
  | .blank => 0 | .ordinary => 1/3 | .held => 19/33

def jointWords (s : Stock) : (ℕ × ℕ) × (ℕ × ℕ) :=
  let a := count00 s
  let b := censusPeriod/2 - a
  ((a,b),(b,a))

theorem complete_word_counts (s : Stock) :
    (jointWords s).1.1 + (jointWords s).1.2 +
      (jointWords s).2.1 + (jointWords s).2.2 = censusPeriod := by
  cases s <;> norm_num [jointWords, count00, censusPeriod]

noncomputable def agreementRate (s : Stock) : ℝ := 2*(count00 s : ℝ)/(censusPeriod : ℝ)

theorem agreementRate_range (s : Stock) :
    agreementRate s ∈ Set.Icc (2 : ℝ)⁻¹ 1 := by
  cases s <;> norm_num [Set.mem_Icc, agreementRate, count00, censusPeriod]

theorem agreementRate_injective : Function.Injective agreementRate := by
  intro x y h
  cases x <;> cases y <;> first | rfl | norm_num [agreementRate, count00, censusPeriod] at h

/-- H of a symmetric four-word distribution. The equality flag has entropy
binEntropy p; conditioning on that flag leaves one independent fair sign. -/
noncomputable def jointEntropy (s : Stock) : ℝ := Real.log 2 + Real.binEntropy (agreementRate s)

theorem fourWord_entropy (p : ℝ) (h₀ : 0 < p) (h₁ : p < 1) :
    -2*(p/2)*Real.log (p/2) - 2*((1-p)/2)*Real.log ((1-p)/2) =
      Real.log 2 + Real.binEntropy p := by
  rw [Real.log_div h₀.ne' (by norm_num), Real.log_div (by linarith : 1-p ≠ 0) (by norm_num)]
  simp only [Real.binEntropy, Real.log_inv]
  ring

theorem jointEntropy_injective : Function.Injective jointEntropy := by
  intro x y h
  apply agreementRate_injective
  apply Real.binEntropy_strictAntiOn.injOn (agreementRate_range x) (agreementRate_range y)
  exact add_left_cancel (by simpa only [jointEntropy] using h)

theorem clockBin_injective : Function.Injective clockBin := by
  intro x y h
  cases x <;> cases y <;> first | rfl | norm_num [clockBin] at h

theorem jointWords_injective : Function.Injective jointWords := by
  intro x y h
  have hh : count00 x = count00 y := congrArg (fun p : (ℕ × ℕ) × (ℕ × ℕ) => p.1.1) h
  cases x <;> cases y <;> first | rfl | norm_num [count00] at hh

/-- Finite quantized sensor and census registers from each appended cohort.
Only completed deterministic population acquisitions belong in this table. -/
structure CensusRow where
  sensor : Fin 9 → ℕ
  mask : Fin 2 → Bool
  modularCounts : Fin 4 → ℕ
  eventStamps : ℕ × ℕ
  weighedShells : ℚ × ℚ
  entropyCounts : Fin 4 → ℕ
  releaseScale : ℕ × ℕ

noncomputable def stockCatalogue (rows : ℕ → Stock → CensusRow) : Catalogue Stock where
  view
    | .geometry => ⟨ℕ, fun _ => (Fin 9 → ℕ) × (Fin 2 → Bool), clockBin,
      fun i s => ((rows i s).sensor, (rows i s).mask)⟩
    | .modular => ⟨(ℕ × ℕ) × (ℕ × ℕ), fun _ => Fin 4 → ℕ,
      jointWords, fun i s => (rows i s).modularCounts⟩
    | .event => ⟨ℕ × ℕ, fun _ => ℕ × ℕ, fun s => (1,clockBin s), fun i s => (rows i s).eventStamps⟩
    | .stress => ⟨ℚ, fun _ => ℚ × ℚ, mass, fun i s => (rows i s).weighedShells⟩
    | .entropy => ⟨ℝ, fun _ => Fin 4 → ℕ, jointEntropy, fun i s => (rows i s).entropyCounts⟩
    | .scale => ⟨ℕ × ℕ, fun _ => ℕ × ℕ,
      fun _ => (2475880078570760549798248448,2475880029379442694124883014),
      fun i s => (rows i s).releaseScale⟩
  geometryKey := clockBin_injective
  modularKey := jointWords_injective
  eventKey := fun _ _ h => clockBin_injective (congrArg Prod.snd h)
  entropyKey := jointEntropy_injective

/-- The finite EB structure on the three exhaustive stock censuses. All
additional deterministic cohorts may have different clocks, masks, units
and word counts; those rows are retained componentwise. -/
noncomputable def clockPairCensusTower (rows : ℕ → Stock → CensusRow) : EinsteinAdmissibleTower ℕ :=
  admissibleCensus (stockCatalogue rows)

#print axioms admissibleCensus
#print axioms fourWord_entropy
#print axioms jointEntropy_injective
#print axioms clockPairCensusTower

end OPH.EinsteinBranch.ClockPairRecords
