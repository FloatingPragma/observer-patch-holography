import ObserverPatchHolography.Locality.NoSignalling

open scoped BigOperators
open Matrix OPH.Locality

namespace OPH.MaterialSeparator

/-- Tracing out the material separator and the source side of each
conditional product flight preserves the remote record. The complete source
operation can vary by sector, including on correlated initial states. -/
theorem pin_sectors_preserve_remote
    {B L R I : Type*} [Fintype B] [Fintype L] [Fintype R]
    [DecidableEq L] [DecidableEq R]
    (s : Finset I) (K : B → I → Matrix L L ℂ)
    (hc : ∀ b, ∑ i ∈ s, (K b i)ᴴ * K b i = 1)
    (M : B → Matrix (L × R) (L × R) ℂ) :
    (∑ b, ptraceFst (∑ i ∈ s, liftLeft (β := R) (K b i) * M b
               * (liftLeft (β := R) (K b i))ᴴ)) =
    ∑ b, ptraceFst (M b) := by
  classical
  apply Finset.sum_congr rfl
  intro b _
  exact ptraceFst_local_kraus s (K b) (hc b) (M b)

#print axioms pin_sectors_preserve_remote

end OPH.MaterialSeparator
