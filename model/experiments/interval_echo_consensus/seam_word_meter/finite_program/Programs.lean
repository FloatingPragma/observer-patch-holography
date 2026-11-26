import Mathlib.Tactic

namespace OPH.MaterialMeanPrograms

/-- Low place first, with native digits zero, one or two on a write. -/
def places : List ℕ → ℕ
  | [] => 0
  | a :: rest => a + 2 * places rest

def carry : ℕ → List ℕ → ℕ
  | c, [] => c
  | c, a :: rest => (c + a) % 2 + 2 * carry ((c + a) / 2) rest

theorem carry_exact (c : ℕ) (ts : List ℕ) : carry c ts = c + places ts := by
  induction ts generalizing c with
  | nil => simp [carry, places]
  | cons a xs ih => rw [carry, ih, places]; omega

theorem observed_sum (a b : List ℕ) (h : a.length = b.length) :
    places (a.zipWith (· + ·) b) = places a + places b := by
  induction a generalizing b with
  | nil => cases b <;> simp_all [places]
  | cons x as ih =>
    cases b with
    | nil => simp at h
    | cons y bs =>
      simp only [List.length_cons, Nat.succ.injEq] at h
      simp only [List.zipWith_cons_cons, places, ih _ h]
      omega

theorem exact_even_midpoint (a b : List ℕ) (h : a.length = b.length)
    (even : (places a + places b) % 2 = 0) :
    2 * (carry 0 (a.zipWith (· + ·) b) / 2) = places a + places b := by
  rw [carry_exact, Nat.zero_add, observed_sum a b h]
  omega

/-- Any supplied finite real calibration decodes componentwise means linearly. -/
def decode {n : ℕ} (standard weight : Fin n → ℝ) : ℝ := ∑ i, weight i * standard i

theorem decode_midpoint {n : ℕ} (s x y : Fin n → ℝ) :
    decode s (fun i => (x i + y i) / 2) = (decode s x + decode s y) / 2 := by
  unfold decode
  rw [← Finset.sum_add_distrib, Finset.sum_div]
  exact Finset.sum_congr rfl (fun i _ => by ring)

/-- A finite decision tree for finite instruction and observation alphabets. -/
inductive Program (Move Obs Result : Type) where
  | stop : Result → Program Move Obs Result
  | act : Move → Program Move Obs Result → Program Move Obs Result
  | inspect : (Obs → Program Move Obs Result) → Program Move Obs Result

def run {M O R S : Type} (step : M → S → S) (read : S → O) :
    Program M O R → S → R
  | .stop result, _ => result
  | .act m next, s => run step read next (step m s)
  | .inspect next, s => run step read (next (read s)) s

/-- Any supported interpretation with commuting primitive actions and material
readouts has the same complete result for every finite adaptive program. The
result can itself include the entire public transcript. -/
theorem adaptive_execution {M O R S T : Type} (realize : S → T)
    (step : M → S → S) (machine : M → T → T) (read : S → O) (station : T → O)
    (move_ok : ∀ m s, machine m (realize s) = realize (step m s))
    (read_ok : ∀ s, station (realize s) = read s)
    (prog : Program M O R) (s : S) :
    run machine station prog (realize s) = run step read prog s := by
  induction prog generalizing s with
  | stop => rfl
  | act m next ih => simp only [run, move_ok, ih]
  | inspect next ih => simp only [run, read_ok, ih]

def mass (p : ℝ) : List Bool → ℝ
  | [] => 1
  | b :: rest => (if b then p else 1-p) * mass p rest

theorem interior_full_support (p : ℝ) (h : 0 < p ∧ p < 1) (word : List Bool) :
    0 < mass p word := by
  induction word with
  | nil => norm_num [mass]
  | cons b rest ih =>
    cases b
    · exact mul_pos (sub_pos.mpr h.2) ih
    · exact mul_pos h.1 ih

/-- An apparatus may condition its response on the entire finite launch tape.
Every attainable terminal observation is shared by both interior sources. -/
theorem same_pooled_observations {T : Type} (protocol : List Bool → T → Prop)
    (p q : ℝ) (hp : 0 < p ∧ p < 1) (hq : 0 < q ∧ q < 1) (t : T) :
    (∃ w, 0 < mass p w ∧ protocol w t) ↔ (∃ w, 0 < mass q w ∧ protocol w t) := by
  constructor
  · rintro ⟨w, _, ht⟩; exact ⟨w, interior_full_support q hq w, ht⟩
  · rintro ⟨w, _, ht⟩; exact ⟨w, interior_full_support p hp w, ht⟩

end OPH.MaterialMeanPrograms

namespace OPH.MaterialMeanPrograms

def bits : ℕ → ℕ → List ℕ
  | 0, _ => []
  | w+1, n => n % 2 :: bits w (n / 2)

theorem bits_length (w n : ℕ) : (bits w n).length = w := by
  induction w generalizing n <;> simp [bits, *]

theorem places_bits (w n : ℕ) (hn : n < 2^w) : places (bits w n) = n := by
  induction w generalizing n with
  | zero => simp at hn; simp [bits, places, hn]
  | succ w ih =>
    rw [pow_succ] at hn
    simp only [bits, places]
    rw [ih _ (by omega)]
    omega

abbrev Code (w : ℕ) := Fin (2^w)

/-- A physical digit word is read at each of the two endpoints. Its carries
are normalized locally, exactly as in the finite-program machine. -/
def physicalMean {w : ℕ} (a b : Code w) : Code w :=
  ⟨carry 0 ((bits w a).zipWith (· + ·) (bits w b)) / 2, by
    have eqLen : (bits w a).length = (bits w b).length := by simp [bits_length]
    rw [carry_exact, Nat.zero_add, observed_sum _ _ eqLen, places_bits _ _ a.isLt,
      places_bits _ _ b.isLt]
    have ha := a.isLt
    have hb := b.isLt
    omega⟩

theorem physicalMean_val {w : ℕ} (a b : Code w) :
    (physicalMean a b).val = (a.val+b.val)/2 := by
  change carry 0 _ / 2 = _
  rw [carry_exact, Nat.zero_add, observed_sum _ _ (by simp [bits_length]),
    places_bits _ _ a.isLt, places_bits _ _ b.isLt]

/-- The two real readings agree exactly when a spare low place is reserved. -/
theorem two_mul_physicalMean {w : ℕ} (a b : Code w) (ha : 2 ∣ a.val) (hb : 2 ∣ b.val) :
    2 * (physicalMean a b).val = a.val + b.val := by
  rw [physicalMean_val]
  obtain ⟨u, hu⟩ := ha
  obtain ⟨v, hv⟩ := hb
  omega

abbrev Shelf (Site : Type) (n w : ℕ) := Site → Fin n → Code w

def Ready {Site : Type} {n w : ℕ} (d : ℕ) (s : Shelf Site n w) :=
  ∀ i j, 2^d ∣ (s i j).val

def act {Site : Type} [DecidableEq Site] {n w : ℕ} (a b : Site) (s : Shelf Site n w) :
    Shelf Site n w := fun i k => if i = a ∨ i = b then physicalMean (s a k) (s b k) else s i k

theorem reserve_weaken {Site : Type} {n w d : ℕ} (s : Shelf Site n w) (h : Ready (d+1) s) :
    Ready d s := by
  intro i j
  exact dvd_trans (by rw [pow_succ]; exact dvd_mul_right _ _) (h i j)

theorem reserve_after_mean {Site : Type} [DecidableEq Site] {n w d : ℕ}
    (s : Shelf Site n w) (a b : Site) (h : Ready (d+1) s) : Ready d (act a b s) := by
  intro i j
  unfold act
  split
  · obtain ⟨x, hx⟩ := h a j
    obtain ⟨y, hy⟩ := h b j
    refine ⟨x+y, ?_⟩
    rw [physicalMean_val, hx, hy, pow_succ]
    have : 2^d*2*x+2^d*2*y = (2^d*(x+y))*2 := by ring
    rw [this, Nat.mul_div_cancel _ (by decide : 0 < 2)]
  · exact reserve_weaken s h i j

noncomputable def amount {n w : ℕ} (standard : Fin n → ℝ) (den : ℝ) (key : Fin n → Code w) : ℝ :=
  decode standard fun j => (key j).val / den

noncomputable def realAct {Site : Type} [DecidableEq Site] (a b : Site) (s : Site → ℝ) : Site → ℝ :=
  fun i => if i = a ∨ i = b then (s a + s b)/2 else s i

theorem material_step {Site : Type} [DecidableEq Site] {n w d : ℕ}
    (standard : Fin n → ℝ) (den : ℝ) (s : Shelf Site n w) (a b : Site) (h : Ready (d+1) s) :
    (fun i => amount standard den (act a b s i)) =
      realAct a b (fun i => amount standard den (s i)) := by
  have ev (i : Site) (j : Fin n) : 2 ∣ (s i j).val := by
    exact dvd_trans (dvd_pow_self 2 (Nat.succ_ne_zero d)) (h i j)
  funext i
  by_cases hi : i = a ∨ i = b
  · simp only [act, realAct, if_pos hi, amount]
    rw [← decode_midpoint]
    apply congrArg (decode standard)
    funext j
    have eqNat := two_mul_physicalMean (s a j) (s b j) (ev a j) (ev b j)
    have eqReal : 2 * ((physicalMean (s a j) (s b j)).val : ℝ) =
        ((s a j).val : ℝ) + ((s b j).val : ℝ) := by exact_mod_cast eqNat
    rw [show ((physicalMean (s a j) (s b j)).val : ℝ) =
        (((s a j).val : ℝ) + ((s b j).val : ℝ))/2 by linarith]
    ring
  · have hsame : act a b s i = s i := by funext k; simp [act, hi]
    rw [hsame]; simp [realAct, hi]

/-- Compile the supplied finite read menu in advance. At run time the station
indexes an ordinary finite table with its measured coefficient word. -/
noncomputable def materialTable {K O : Type} [Fintype K] (f : K → O) :
    Vector O (Fintype.card K) := Vector.ofFn fun i => f ((Fintype.equivFin K).symm i)

theorem table_at_key {K O : Type} [Fintype K] (f : K → O) (k : K) :
    (materialTable f).get (Fintype.equivFin K k) = f k := by
  change (Vector.ofFn fun i => f ((Fintype.equivFin K).symm i))[(Fintype.equivFin K k).val] = f k
  simp

inductive Tree (Site Obs Result : Type) where
  | stop : Result → Tree Site Obs Result
  | mean : Site → Site → Tree Site Obs Result → Tree Site Obs Result
  | read : Site → (Obs → Tree Site Obs Result) → Tree Site Obs Result

def execute {S I O R : Type} (step : I → I → S → S) (inspect : I → S → O) :
    ℕ → Tree I O R → S → Option R
  | _, .stop r, _ => some r
  | 0, _, _ => none
  | d+1, .mean a b next, s => execute step inspect d next (step a b s)
  | d+1, .read i next, s => execute step inspect d (next (inspect i s)) s

/-- There are no local correctness hypotheses in this instance. The same
width bound, calibration and finite read table apply to the entire program.
Fuel also bounds the required spare low bits. -/
theorem supported_program {I O R : Type} [DecidableEq I] {n w : ℕ}
    (standard : Fin n → ℝ) (den : ℝ) (policy : I → ℝ → O) (d : ℕ)
    (tree : Tree I O R) (s : Shelf I n w) (ready : Ready d s) :
    execute (act (n:=n) (w:=w)) (fun i t =>
      (materialTable fun key : Fin n → Code w => policy i (amount standard den key)).get
        (Fintype.equivFin _ (t i))) d tree s =
    execute realAct (fun i t => policy i (t i)) d tree (fun i => amount standard den (s i)) := by
  induction d generalizing tree s with
  | zero => cases tree <;> rfl
  | succ d ih =>
    cases tree with
    | stop => rfl
    | mean a b next =>
      simp only [execute]
      rw [ih next (act a b s) (reserve_after_mean s a b ready), material_step standard den s a b ready]
    | read i next =>
      simp only [execute, table_at_key]
      simpa only [table_at_key] using ih (next (policy i (amount standard den (s i)))) s (reserve_weaken s ready)

end OPH.MaterialMeanPrograms

namespace OPH.MaterialMeanPrograms

/-- An unknown member of any finite real preparation menu is represented by
its material axis marker. The calibration and width are independent of the
particular choice of marker. -/
def enrolled {I : Type} {n : ℕ} (tag : I → Fin n) (d : ℕ) : Shelf I n (d+1) :=
  fun i j => ⟨if j = tag i then 2^d else 0, by
    have hp : 0 < 2^d := by positivity
    rw [pow_succ]
    split <;> omega⟩

theorem enrollment_ready {I : Type} {n : ℕ} (tag : I → Fin n) (d : ℕ) :
    Ready d (enrolled tag d) := by
  intro i j
  change 2^d ∣ (if j = tag i then 2^d else 0)
  split <;> simp

theorem enrollment_value {I : Type} {n : ℕ} (s : Fin n → ℝ) (tag : I → Fin n) (d : ℕ) (i : I) :
    amount s (2^d) (enrolled tag d i) = s (tag i) := by
  classical
  unfold amount decode
  rw [Finset.sum_eq_single (tag i)]
  · simp [enrolled]
  · intro j _ h; simp [enrolled, h]
  · simp

theorem complete_enrolled_run {I O R : Type} [DecidableEq I] {n : ℕ}
    (standard : Fin n → ℝ) (tag : I → Fin n) (policy : I → ℝ → O)
    (d : ℕ) (tree : Tree I O R) :
    execute (act (n:=n) (w:=d+1)) (fun i t =>
      (materialTable fun key : Fin n → Code (d+1) => policy i (amount standard (2^d) key)).get
        (Fintype.equivFin _ (t i))) d tree (enrolled tag d) =
    execute realAct (fun i t => policy i (t i)) d tree (fun i => standard (tag i)) := by
  simpa only [enrollment_value] using supported_program standard (2^d) policy d tree
    (enrolled tag d) (enrollment_ready tag d)

end OPH.MaterialMeanPrograms
