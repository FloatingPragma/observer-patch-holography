import Std

namespace PatchClosure

variable {A : Type} [DecidableEq A]

/-- Addresses that must be present for every address of a patch to run locally. -/
inductive Needed (seeds : List A) (parents : A → List A) : A → Prop where
  | seed (a) : a ∈ seeds → Needed seeds parents a
  | parent {a b} : Needed seeds parents a → b ∈ parents a → Needed seeds parents b

/-- Scan an execution order backwards, retaining a seed or an actual
    prerequisite of a retained operation. The result is in execution order. -/
def generate (seeds : List A) (parents : A → List A) : List A → List A
  | [] => []
  | a :: tail =>
    let chosen := generate seeds parents tail
    if seeds.contains a || chosen.any (fun b => (parents b).contains a)
    then a :: chosen else chosen

/-- Every predecessor belongs to the execution menu and precedes its dependent. -/
structure Installation (parents : A → List A) (order : List A) : Prop where
  unique : order.Nodup
  noSelf : ∀ a ∈ order, a ∉ parents a
  present : ∀ a ∈ order, ∀ b ∈ parents a, b ∈ order
  earlier : order.Pairwise (fun a b => b ∉ parents a)

private theorem tail_sub (seeds : List A) (parents : A → List A) (a : A) (s : List A) :
    generate seeds parents s ⊆ generate seeds parents (a :: s) := by
  intro b hb
  simp only [generate]
  split
  · exact List.mem_cons_of_mem _ hb
  · exact hb

theorem generate_subset (seeds : List A) (parents : A → List A) (s : List A) :
    generate seeds parents s ⊆ s := by
  induction s with
  | nil => exact fun _ h => h
  | cons a s ih =>
    intro b hb
    simp only [generate] at hb
    split at hb
    · rcases List.mem_cons.mp hb with e | h
      · exact List.mem_cons.mpr (.inl e)
      · exact List.mem_cons_of_mem _ (ih h)
    · exact List.mem_cons_of_mem _ (ih hb)

theorem generate_only_needed (seeds : List A) (parents : A → List A) (s : List A) :
    ∀ a ∈ generate seeds parents s, Needed seeds parents a := by
  induction s with
  | nil => simp [generate]
  | cons b s ih =>
    intro a ha
    simp only [generate] at ha
    split at ha
    next h =>
      rcases List.mem_cons.mp ha with hEq | hRest
      · subst a
        rcases Bool.or_eq_true_iff.mp h with hs | hp
        · exact .seed b (by simpa using hs)
        · obtain ⟨j,hj,hp⟩ := List.any_eq_true.mp hp
          exact .parent (ih j hj) (by simpa using hp)
      · exact ih a hRest
    next => exact ih a ha

private theorem generated_seed (seeds : List A) (parents : A → List A) {a : A}
    (seed : a ∈ seeds) : ∀ s : List A, a ∈ s → a ∈ generate seeds parents s := by
  intro s
  induction s with
  | nil => simp
  | cons b s ih =>
    intro ha
    rcases List.mem_cons.mp ha with e | h
    · subst b
      simp [generate,seed]
    · exact tail_sub seeds parents b s (ih h)

private theorem generated_parent (seeds : List A) (parents : A → List A)
    {a b : A} (hp : b ∈ parents a) (s : List A)
    (hs : s.Pairwise (fun x y => y ∉ parents x))
    (hself : ∀ x ∈ s, x ∉ parents x)
    (ha : a ∈ generate seeds parents s) (hb : b ∈ s) :
    b ∈ generate seeds parents s := by
  induction s with
  | nil => simp at hb
  | cons x xs ih =>
    obtain ⟨hfront,htail⟩ := List.pairwise_cons.mp hs
    have ha' : a ∈ generate seeds parents xs := by
      by_cases e : a = x
      · subst x
        rcases List.mem_cons.mp hb with eb | hbt
        · subst b
          exact False.elim ((hself a (by simp)) hp)
        · exact False.elim (hfront b hbt hp)
      · simp only [generate] at ha
        split at ha
        · exact (List.mem_cons.mp ha).resolve_left e
        · exact ha
    rcases List.mem_cons.mp hb with e | hbt
    · subst x
      have h : (generate seeds parents xs).any (fun j => (parents j).contains b) = true :=
        List.any_eq_true.mpr ⟨a,ha', by simpa using hp⟩
      simp only [generate,h,Bool.or_true,if_true,List.mem_cons_self]
    · exact tail_sub seeds parents x xs
        (ih htail (fun j hj => hself j (List.mem_cons_of_mem _ hj)) ha' hbt)

theorem generate_exact {seeds order : List A} {parents : A → List A}
    (install : Installation parents order) (seedIn : seeds ⊆ order) (a : A) :
    a ∈ generate seeds parents order ↔ Needed seeds parents a := by
  refine ⟨generate_only_needed _ _ _ _, ?_⟩
  intro h
  induction h with
  | seed a ha => exact generated_seed _ _ ha _ (seedIn ha)
  | parent prev hp ih =>
    exact generated_parent _ _ hp order install.earlier install.noSelf ih
      (install.present _ (generate_subset _ _ _ ih) _ hp)

/-- Independent certificate predicate. Every listed address has a justification,
    read from the physical menu, in the suffix closer to a requested address. -/
def justified (seeds : List A) (parents : A → List A) : List A → Bool
  | [] => true
  | a :: tail =>
      (seeds.contains a || tail.any (fun b => (parents b).contains a)) &&
      justified seeds parents tail

def check (seeds : List A) (parents : A → List A) (cert : List A) : Bool :=
  seeds.all cert.contains &&
  cert.all (fun a => (parents a).all cert.contains) &&
  justified seeds parents cert

private theorem justified_sound (seeds : List A) (parents : A → List A) (s : List A)
    (hj : justified seeds parents s = true) : ∀ x ∈ s, Needed seeds parents x := by
  induction s with
  | nil => simp
  | cons a s ih =>
    obtain ⟨ha,ht⟩ := Bool.and_eq_true_iff.mp hj
    intro x hx
    rcases List.mem_cons.mp hx with rfl | h
    · rcases Bool.or_eq_true_iff.mp ha with r | p
      · exact .seed _ (by simpa using r)
      · obtain ⟨b,hb,hp⟩ := List.any_eq_true.mp p
        exact .parent (ih ht _ hb) (by simpa using hp)
    · exact ih ht _ h

theorem check_sound {seeds cert : List A} {parents : A → List A}
    (ok : check seeds parents cert = true) (x : A) :
    x ∈ cert ↔ Needed seeds parents x := by
  obtain ⟨h,hj⟩ := Bool.and_eq_true_iff.mp ok
  obtain ⟨hs,hp⟩ := Bool.and_eq_true_iff.mp h
  constructor
  · exact justified_sound seeds parents cert hj x
  · intro hx
    induction hx with
    | seed a ha => simpa using List.all_eq_true.mp hs a ha
    | parent h hpar ih =>
      exact (by simpa using List.all_eq_true.mp (List.all_eq_true.mp hp _ ih) _ hpar)

private theorem generated_justified (seeds : List A) (parents : A → List A) (s : List A) :
    justified seeds parents (generate seeds parents s) = true := by
  induction s with
  | nil => rfl
  | cons a xs ih =>
    simp only [generate]
    split
    next h => exact Bool.and_eq_true_iff.mpr ⟨h,ih⟩
    next => exact ih

theorem generated_certificate_passes {seeds order : List A} {parents : A → List A}
    (install : Installation parents order) (seedIn : seeds ⊆ order) :
    check seeds parents (generate seeds parents order) = true := by
  apply Bool.and_eq_true_iff.mpr
  refine ⟨?_,generated_justified _ _ _⟩
  apply Bool.and_eq_true_iff.mpr
  constructor
  · exact List.all_eq_true.mpr (fun a ha => by
      simpa using (generate_exact install seedIn a).mpr (.seed _ ha))
  · exact List.all_eq_true.mpr (fun a ha => List.all_eq_true.mpr (fun b hb => by
      simpa using (generate_exact install seedIn b).mpr
        (.parent ((generate_exact install seedIn a).mp ha) hb)))

theorem least_extension {seeds order : List A} {parents : A → List A}
    (install : Installation parents order) (seedIn : seeds ⊆ order)
    (patch : A → Prop) (contains : ∀ a ∈ seeds, patch a)
    (closed : ∀ a, patch a → ∀ b ∈ parents a, patch b) :
    ∀ a ∈ generate seeds parents order, patch a := by
  intro a ha
  have h := (generate_exact install seedIn a).mp ha
  clear ha
  induction h with
  | seed a h => exact contains a h
  | parent h hp ih => exact closed _ ih _ hp

/-- Reject a permutation that names the right vertices in an illegal order. -/
def orderedCheck (parents : A → List A) : List A → Bool
  | [] => true
  | a :: tail => !tail.contains a && !(parents a).contains a &&
      tail.all (fun later => !(parents a).contains later) && orderedCheck parents tail

theorem orderedCheck_spec (parents : A → List A) (xs : List A) :
    orderedCheck parents xs = true ↔ xs.Nodup ∧ (∀ a ∈ xs, a ∉ parents a) ∧
      xs.Pairwise (fun a b => b ∉ parents a) := by
  induction xs with
  | nil => simp [orderedCheck]
  | cons a xs ih =>
    simp [orderedCheck,ih,List.nodup_cons,List.pairwise_cons]
    grind

def acceptCertificate (seeds : List A) (parents : A → List A) (cert : List A) : Bool :=
  check seeds parents cert && orderedCheck parents cert

theorem acceptCertificate_sound {seeds cert : List A} {parents : A → List A}
    (h : acceptCertificate seeds parents cert = true) :
    Installation parents cert ∧ (∀ a, a ∈ cert ↔ Needed seeds parents a) := by
  obtain ⟨hc,ho⟩ := Bool.and_eq_true_iff.mp h
  obtain ⟨hn,hs,he⟩ := (orderedCheck_spec _ _).mp ho
  have equiv := check_sound hc
  exact ⟨⟨hn,hs,fun a ha b hb => (equiv b).mpr (.parent ((equiv a).mp ha) hb),he⟩,equiv⟩

private theorem generate_sublist (seeds : List A) (parents : A → List A) (s : List A) :
    (generate seeds parents s).Sublist s := by
  induction s with
  | nil => exact .slnil
  | cons a s ih =>
    simp only [generate]
    split
    · exact .cons₂ a ih
    · exact .cons a ih

theorem producer_accepted {seeds order : List A} {parents : A → List A}
    (install : Installation parents order) (seedIn : seeds ⊆ order) :
    acceptCertificate seeds parents (generate seeds parents order) = true := by
  have sub := generate_sublist seeds parents order
  exact Bool.and_eq_true_iff.mpr ⟨generated_certificate_passes install seedIn,
    (orderedCheck_spec _ _).mpr ⟨install.unique.sublist sub,
      fun a ha => install.noSelf a (generate_subset _ _ _ ha), install.earlier.sublist sub⟩⟩

/-- Records on admissible patches restrict contravariantly. -/
def restrict {P Q : A → Prop} (i : ∀ a, P a → Q a) {R : Type} (r : Subtype Q → R) :
    Subtype P → R := fun a => r ⟨a.val,i a.val a.property⟩

omit [DecidableEq A] in
theorem restrict_id {P : A → Prop} {R : Type} (r : Subtype P → R) :
    restrict (fun _ h => h) r = r := rfl

omit [DecidableEq A] in
theorem restrict_comp {P Q U : A → Prop} {R : Type} (i : ∀ a, P a → Q a)
    (j : ∀ a, Q a → U a) (r : Subtype U → R) :
    restrict i (restrict j r) = restrict (fun a h => j a (i a h)) r := rfl

#print axioms generate_exact
#print axioms least_extension
#print axioms acceptCertificate_sound
#print axioms producer_accepted
#print axioms restrict_comp
end PatchClosure
