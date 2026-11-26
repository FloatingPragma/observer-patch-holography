import Mathlib

namespace OPH.ScreenCheckpoint

abbrev Word := Fin 128
abbrev Label := Fin 8

def syndrome (x : Word) : BitVec 3 :=
  (List.finRange 7).foldl (fun s i =>
    if (BitVec.ofNat 7 x.val).getLsbD i.val then s ^^^ BitVec.ofNat 3 (i.val+1) else s) 0

def flip (x : Word) (i : Label) : Word :=
  if i.val = 0 then x else ((BitVec.ofNat 7 x.val) ^^^ ((1 : BitVec 7) <<< (i.val-1))).toFin

def correct (x : Word) : Word := flip x (syndrome x).toFin

set_option maxRecDepth 8192 in
set_option maxHeartbeats 4000000 in
theorem correct_read : ∀ (x : Word) (i : Label),
    correct (flip (correct x) i) = correct x := by decide +kernel

set_option maxRecDepth 8192 in
theorem code_size : (Finset.univ.filter (fun x : Word => correct x = x)).card = 16 := by
  decide +kernel

theorem packing_upper {M : Type*} [Fintype M] (publish : M → Word)
    (decode : Word → M) (right : ∀ m e, decode (flip (publish m) e) = m) :
    Fintype.card M ≤ 16 := by
  have errors : ∀ (x : Word) (i j : Label), flip x i = flip x j → i = j := by
    decide +kernel
  have inj : Function.Injective (fun p : M × Label => flip (publish p.1) p.2) := by
    rintro ⟨m,i⟩ ⟨n,j⟩ eq
    have mn : m = n := (right m i).symm.trans ((congrArg decode eq).trans (right n j))
    subst n
    exact Prod.ext rfl (errors _ _ _ eq)
  have bound := Fintype.card_le_of_injective _ inj
  rw [Fintype.card_prod] at bound
  simp only [Label, Word, Fintype.card_fin] at bound
  omega

end OPH.ScreenCheckpoint
