import Mathlib

/-!
# Physical A5-to-SM forcing: finite non-identifiability core

These theorems formalize the exact logical boundary used by the paper.  The
completion constructors are abstract finite tags, not physical OPH models.
They prove that a forgetful source interface with two distinct completions has
no source-only function that reconstructs every completion.  A richer
observer-like operational packet may still select a completion.
-/

namespace OPHPhysicalA5NoGo

inductive ChargeProfile where
  | tetrahedral
  | octahedral
  | icosahedral
  deriving DecidableEq, Repr

def charges : ChargeProfile → List Nat
  | .tetrahedral => [3, 3, 3, 3]
  | .octahedral => [2, 2, 2, 2, 2, 2]
  | .icosahedral => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

theorem everyProfileHasTotalTwelve (profile : ChargeProfile) :
    (charges profile).sum = 12 := by
  cases profile <;> decide

structure EulerSource where
  totalCharge : Nat
  deriving DecidableEq

def eulerTwelve : EulerSource := ⟨12⟩

def forgetCharge (_profile : ChargeProfile) : EulerSource := eulerTwelve

/-- Euler total charge alone cannot reconstruct every admissible profile. -/
theorem noSourceOnlyChargeReconstruction :
    ¬ ∃ reconstruct : EulerSource → ChargeProfile,
        ∀ profile, reconstruct (forgetCharge profile) = profile := by
  rintro ⟨reconstruct, h⟩
  have ht := h ChargeProfile.tetrahedral
  have hi := h ChargeProfile.icosahedral
  change reconstruct eulerTwelve = ChargeProfile.tetrahedral at ht
  change reconstruct eulerTwelve = ChargeProfile.icosahedral at hi
  rw [ht] at hi
  cases hi

inductive CurrentCompletion where
  | abelianTwelve
  | compactStandardModelType
  deriving DecidableEq, Repr

structure PortSource where
  portCount : Nat
  totalCharge : Nat
  deriving DecidableEq

def twelvePortSource : PortSource := ⟨12, 12⟩

def forgetCurrent (_completion : CurrentCompletion) : PortSource :=
  twelvePortSource

/-- Port count and total charge do not select a physical current bracket. -/
theorem noSourceOnlyCurrentReconstruction :
    ¬ ∃ reconstruct : PortSource → CurrentCompletion,
        ∀ completion, reconstruct (forgetCurrent completion) = completion := by
  rintro ⟨reconstruct, h⟩
  have ha := h CurrentCompletion.abelianTwelve
  have hs := h CurrentCompletion.compactStandardModelType
  change reconstruct twelvePortSource = CurrentCompletion.abelianTwelve at ha
  change reconstruct twelvePortSource =
    CurrentCompletion.compactStandardModelType at hs
  rw [ha] at hs
  cases hs

inductive MatterCompletion where
  | exteriorPacket
  | exteriorPacketPlusSterile
  deriving DecidableEq, Repr

def forgetMatter (_completion : MatterCompletion) : CurrentCompletion :=
  .compactStandardModelType

/-- A bosonic current reduct does not exclude an invisible sterile extension. -/
theorem noCurrentOnlyMatterReconstruction :
    ¬ ∃ reconstruct : CurrentCompletion → MatterCompletion,
        ∀ completion, reconstruct (forgetMatter completion) = completion := by
  rintro ⟨reconstruct, h⟩
  have he := h MatterCompletion.exteriorPacket
  have hs := h MatterCompletion.exteriorPacketPlusSterile
  change reconstruct CurrentCompletion.compactStandardModelType =
    MatterCompletion.exteriorPacket at he
  change reconstruct CurrentCompletion.compactStandardModelType =
    MatterCompletion.exteriorPacketPlusSterile at hs
  rw [he] at hs
  cases hs

def settle4 : Fin 4 → Fin 4 := ![0, 0, 1, 2]

theorem settle4NotInjective : ¬ Function.Injective settle4 := by
  intro h
  have h01 : settle4 0 = settle4 1 := by decide
  have : (0 : Fin 4) = 1 := h h01
  omega

theorem fixedPointFreePairingCountArithmetic :
    Nat.factorial 12 / (2 ^ 6 * Nat.factorial 6) = 10395 := by
  norm_num [Nat.factorial]

#print axioms everyProfileHasTotalTwelve
#print axioms noSourceOnlyChargeReconstruction
#print axioms noSourceOnlyCurrentReconstruction
#print axioms noCurrentOnlyMatterReconstruction
#print axioms settle4NotInjective
#print axioms fixedPointFreePairingCountArithmetic

end OPHPhysicalA5NoGo
