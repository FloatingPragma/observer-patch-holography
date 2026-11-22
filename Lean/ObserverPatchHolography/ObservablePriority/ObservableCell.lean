import ObserverPatchHolography.ObservablePriority.ObservablePriority
import ObserverPatchHolography.ObservablePriority.CellFlagReadout

namespace ObservableCell
open ObserverPatchSharp ObservablePriority CellBits

/-- A checked cell decoder consumes a closure selected from observed winning
    addresses. Rank only enters the correctness theorem, not recovery. -/
def FitsMenu (c : CellChart) (cert : List Edge) : Prop :=
  cert.Perm (menu c)

private theorem recovered_closed {world cert seeds order : List Edge} {law : Candidate world}
    {o : HitOracle} (real : Realizes law o)
    (built : checkedClosure world o seeds order = some cert) :
    PredecessorClosed law.rank law.schedule cert.contains := by
  have checked := (recovered_closure_exact real built).1
  have installed := (PatchClosure.acceptCertificate_sound checked).1
  intro e he hk f hf touch higher
  have rootMem := law.complete.mem_iff.mp he
  have selected : e ∈ cert := by simpa using hk
  have hp : f ∈ recoverParents world o e :=
    (recoverParents_member real rootMem).mpr ⟨law.complete.mem_iff.mp hf,touch,higher⟩
  simpa using installed.present e selected f hp

private theorem take_cert {world cert seeds order : List Edge} {law : Candidate world}
    {o : HitOracle} (real : Realizes law o)
    (built : checkedClosure world o seeds order = some cert) :
    (law.schedule.filter cert.contains).Perm cert := by
  obtain ⟨check,exactness⟩ := recovered_closure_exact real built
  have structureOK := (PatchClosure.acceptCertificate_sound check).1
  have subset : cert ⊆ world := by
    intro e he
    dsimp only [checkedClosure] at built
    split at built
    next ok => cases built; exact ok.1.mem_iff.mp (PatchClosure.generate_subset _ _ _ he)
    next => contradiction
  apply (List.perm_ext_iff_of_nodup ((law.legal.imp (fun h => h.1)).filter _) structureOK.unique).mpr
  intro e
  simp only [List.mem_filter,List.contains_iff_mem]
  exact ⟨And.right,fun he => ⟨law.complete.mem_iff.mpr (subset he),he⟩⟩

theorem recovered_cell_cost {world cert order : List Edge} (c : CellChart)
    {law : Candidate world} {o : HitOracle} (real : Realizes law o)
    (built : checkedClosure world o c.inner order = some cert)
    (fit : FitsMenu c cert) (hm : Ordered law.rank (menu c)) (nonemptyPatterns : cellCheck c = true) :
    Communicates c (fun L => localWord c (sharpRun law.rank law.schedule L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun law.rank law.schedule L)) b → 3 ≤ b :=
  exactly_three_full c law.rank law.schedule cert.contains law.legal hm (recovered_closed real built)
    ((take_cert real built).trans fit) nonemptyPatterns

/-- A physical menu that is closed and required by every closed installation
    containing the seeds is precisely the one produced from the hit oracle. -/
theorem fits_of_minimum {world cert order : List Edge} (c : CellChart)
    {law : Candidate world} {o : HitOracle} (real : Realizes law o)
    (built : checkedClosure world o c.inner order = some cert)
    (onBoard : menu c ⊆ world) (noDup : (menu c).Nodup)
    (closed : PredecessorClosed law.rank world (menu c).contains)
    (necessary : ∀ other : Edge → Bool, PredecessorClosed law.rank law.schedule other →
      (∀ a ∈ c.inner, other a = true) → ∀ e ∈ menu c, other e = true) :
    FitsMenu c cert := by
  obtain ⟨ok,exactness⟩ := recovered_closure_exact real built
  have certNoDup := (PatchClosure.acceptCertificate_sound ok).1.unique
  apply (List.perm_ext_iff_of_nodup certNoDup noDup).mpr
  intro e
  rw [exactness e]
  constructor
  · intro need
    induction need with
    | seed a ha => exact List.mem_append_left _ ha
    | parent prev par ih =>
      obtain ⟨hb,meet,higher⟩ := PatchComposition.physical_member.mp par
      simpa using closed _ (onBoard ih) (by simpa using ih) _ hb meet higher
  · intro he
    have compulsory := necessary cert.contains (recovered_closed real built) (fun a ha => by
      simpa using (exactness a).mpr (.seed a ha)) e he
    exact (exactness e).mp (by simpa using compulsory)

#print axioms fits_of_minimum
#print axioms recovered_cell_cost
end ObservableCell
