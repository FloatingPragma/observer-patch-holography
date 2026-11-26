import Mathlib.Data.List.Basic

namespace OPH.SeamFlags

/-- The first load is supplied separately. Each subsequent load gives one seam read. -/
def flags (first : Bool) : List Bool → List Bool
  | [] => []
  | b :: rest => (first.xor b) :: flags b rest

/-- Recover the other loads along an observed path once a root load is specified. -/
def recover (first : Bool) : List Bool → List Bool
  | [] => []
  | b :: rest => (first.xor b) :: recover (first.xor b) rest

def arrival (first : Bool) : List Bool → Bool
  | [] => first
  | b :: rest => arrival (first.xor b) rest

theorem flags_recover (first : Bool) (bits : List Bool) :
    flags first (recover first bits) = bits := by
  induction bits generalizing first with
  | nil => rfl
  | cons b xs ih => cases first <;> cases b <;> simp [flags, recover, ih]

theorem recover_flags (first : Bool) (tail : List Bool) :
    recover first (flags first tail) = tail := by
  induction tail generalizing first with
  | nil => rfl
  | cons b xs ih => cases first <;> cases b <;> simp [flags, recover, ih]

/-- Replaying the complete measured path returns to its first load on a closed loop. -/
theorem arrival_close (first last : Bool) (tail : List Bool) :
    arrival first (flags first (tail ++ [last])) = last := by
  induction tail generalizing first with
  | nil => cases first <;> cases last <;> rfl
  | cons b xs ih => cases first <;> cases b <;> simp [flags, arrival, ih]

end OPH.SeamFlags
