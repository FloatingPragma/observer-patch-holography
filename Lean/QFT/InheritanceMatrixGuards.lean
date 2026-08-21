import Mathlib
import QFT.InheritanceMatrix

/-!
# Stale-absence guards for the E4 status matrix

Every typed non-evaluable exit of `QFT/InheritanceMatrix.lean` cites an
absence: a declaration or structure that exists neither in the tree nor
in pinned Mathlib.  Such a citation holds against the pin and goes false
silently when the pin is bumped to a Mathlib that has grown the missing
structure, or when a module lands in the tree that supplies it.  Commit 4012dea5 records a hand correction of a stale row-4 citation
on its merge day.  This module makes the staleness loud: it re-checks the cited
absences on every build and **fails to elaborate** the moment one of
the named probes matches, so a stale exit turns the build red instead of
lying quietly in the matrix.

Two check layers:

* the probe tables and their row links are ordinary definitions with
  `rfl` receipts, kernel-checked against `e4StatusMatrix`;
* the scans themselves run at elaboration time (`run_cmd`) over the
  loaded environment (the full pinned Mathlib plus the import closure
  of the matrix) and `throwError` on any hit, which fails the build.

Fail-closed clauses: the scan asserts it actually saw pinned Mathlib
(module- and constant-count floors), so an import restructure that
emptied the scan would itself be a red, never a silent green.  A cited
module or declaration that stops resolving is likewise a red: presence
citations are asserted by name, and `import Mathlib` fails outright if
the library root disappears.

## What each guard covers (and its companion)

The tree-file-level companion `code/audit/test_e4_absence_guards.py`
covers claims this module structurally cannot see: files outside the
import closure of `QFT.InheritanceMatrix` (e.g. a future
`QFT/TimeSlice.lean` carrying a stronger Cauchy construction, or the row-4
artifact paths in `code/`), and comment-level text such as the cited
Tomita TODO in `Mathlib/Analysis/InnerProductSpace/StandardSubspace.lean`.

| Row | Guarded absence | Layer |
|-----|-----------------|-------|
| 1 | no Tomita/modular-conjugation/KMS declaration in pinned Mathlib; no Tomita declaration in the closure | env scan |
| 1 | `Mathlib.Analysis.InnerProductSpace.StandardSubspace` exists (cited TODO home) | presence |
| 2 | `Mathlib.LinearAlgebra.CliffordAlgebra.SpinGroup` and `spinGroup` exist (cited as present-but-unconnected) | presence |
| 2 | no spin-structure declaration in the closure | env scan |
| 3 | no DHR quasi-local net, selected vacuum representation, or Doplicher-Roberts declaration in the closure | env scan |
| 6 | no Lorentz/globally-hyperbolic/Cauchy-surface/Cauchy-embedding declaration or module in pinned Mathlib | env scan |
| 6 | no Cauchy-embedding declaration in the closure | env scan |
| 6 | the finite `TimeIndexedNetEvolution` interface and its row anchor remain present | tree companion |
| 7 | no Haag/Ruelle/scattering declaration or module in pinned Mathlib | env scan |
| 7 | no energy-momentum, spectrum-condition, or spectral-adapter declaration in the closure | env scan |

## What these guards do NOT cover

Stated here so a green build is never quoted as more than it is:

* A guard checks that **named probes are absent from the pinned
  library and the loaded closure**.  It does not prove the row's
  mathematical claim, that the target is unstatable by any route,
  and it cannot: absence of a name is not absence of a concept.
* A structure landing under a name none of the probes match passes
  the guard.  The probes are substrings chosen against the citation
  vocabulary of the matrix, not a semantic net.
* The closure scans see only modules imported (transitively) by this
  file.  Tree files outside that closure are the Python companion's
  job, and row 4's gate, "the full field/action reconstruction packet is
  not supplied", names no single declaration, so no name probe can watch
  it; only its concrete precursor citations are checked in the companion.
* Row 5's horizon exit cites no external structure (it is the typed
  exit of `QFT/StructuralInheritance.lean` itself, pinned by
  the row-5 anchors), so there is nothing external for a guard to
  watch: it cannot go stale through Mathlib growth.
* A guard firing is a *finding*, not automatically a discharge: a hit
  means the cited absence needs re-review, and the correct response
  may be to reword the row, not to claim the target statable.
-/

namespace OPH.QFT.E4AbsenceGuards

open Lean

/-- One stale-absence probe set for a matrix row: lowercase substrings
matched against the lowercased full names of declarations and modules.
A match means the cited absence has (probably) ended, and the build
fails so a human re-reviews the row. -/
structure AbsenceProbe where
  row : Nat
  claim : String
  probes : List String
  deriving Repr

/-- Probes over pinned-Mathlib-origin names (module root `Mathlib`).
These watch the "pinned Mathlib lacks X" halves of rows 1, 6, 7. -/
def mathlibAbsenceProbes : List AbsenceProbe :=
  [{ row := 1
     claim := "Tomita conjugation / modular theory / KMS is an explicit \
       TODO in pinned Mathlib (StandardSubspace.lean)"
     probes := ["tomita", "modularconjugation", "modularoperator", "kms"] },
   { row := 6
     claim := "no Lorentzian or globally-hyperbolic structure of any \
       kind in pinned Mathlib"
     probes := ["lorentz", "globallyhyperbolic", "globalhyperbolic",
       "cauchysurface", "cauchyembedding", "pseudoriemannian",
       "semiriemannian"] },
   { row := 7
     claim := "zero Haag/scattering material in pinned Mathlib"
     probes := ["haag", "ruelle", "scattering", "waveoperator"] }]

/-- Probes over project-origin names (module root not a dependency
root) in the import closure of `QFT.InheritanceMatrix`.  These watch
the "no X exists in the tree" halves of rows 1, 2, 3, 6, 7 as far as
this closure reaches; the whole-tree file scan is the Python
companion's job. -/
def projectAbsenceProbes : List AbsenceProbe :=
  [{ row := 1
     claim := "no Tomita-type conjugation on the record structures in tree"
     probes := ["tomita"] },
   { row := 2
     claim := "no spin structure exists in the tree"
     probes := ["spinstructure"] },
   { row := 3
     claim := "no DHR quasi-local net, selected vacuum representation, or \
       Doplicher-Roberts reconstruction in tree (full-DHR exit)"
     probes := ["quasilocal", "vacuumrepresentation", "doplicherroberts"] },
   { row := 6
     claim := "no relative-Cauchy evolution, Lorentzian metric \
       perturbation, or stress-response derivative in tree; the E3 Cauchy \
       embedding class and the typed time-indexed net-evolution interface \
       are cited as present separately"
     probes := ["relativecauchy", "metricperturbation", "stressresponse"] },
   { row := 7
     claim := "no energy-momentum or spectrum-condition object and no \
       B9 spectral adapter in tree"
     probes := ["energymomentum", "spectrumcondition", "spectraladapter"] }]

/-- A citation of something the matrix asserts is PRESENT (row 1's TODO
home; row 2's present-but-unconnected `spinGroup`).  If it stops
resolving the citation is stale in the other direction, and that is
equally a red. -/
structure PresenceCitation where
  row : Nat
  claim : String
  module? : Option Name := none
  decl? : Option Name := none

/-- The presence citations of the matrix that live in pinned Mathlib. -/
def presenceCitations : List PresenceCitation :=
  [{ row := 1
     claim := "the cited Tomita-TODO module exists"
     module? := some `Mathlib.Analysis.InnerProductSpace.StandardSubspace },
   { row := 2
     claim := "SpinGroup.lean exists (cited as unconnected)"
     module? := some `Mathlib.LinearAlgebra.CliffordAlgebra.SpinGroup
     decl? := some `spinGroup },
   { row := 6
     claim := "the order-theoretic E3 Cauchy embedding class exists \
       (cited with no Lorentzian reading attached)"
     module? := some `QFT.LocallyCovariantLimit
     decl? := some `OPH.QFT.IsCauchyEmbedding }]

/-- Module roots that belong to dependencies, not to this project.
Anything else in the loaded environment is treated as project origin. -/
def dependencyRoots : List Name :=
  [`Mathlib, `Std, `Lean, `Init, `Lake, `Batteries, `Aesop, `Qq,
   `ProofWidgets, `Plausible, `Cli, `ImportGraph, `LeanSearchClient]

/-- Fail-closed floor: the scan must have seen at least this many
pinned-Mathlib modules, or it did not scan the library it claims to
have scanned.  Pinned Mathlib has 7872. -/
def mathlibModuleFloor : Nat := 7000

/-- Fail-closed floor: the scan must have seen at least this many
Mathlib-origin constants.  Pinned Mathlib contributes ~530000. -/
def mathlibConstFloor : Nat := 400000

/-! ## Kernel-checked links between the guard tables and the matrix -/

/-- Every guarded row number is a row of `e4StatusMatrix`. -/
theorem guarded_rows_exist_in_matrix :
    ((mathlibAbsenceProbes ++ projectAbsenceProbes).map (·.row) ++
      presenceCitations.map (·.row)).all
        (fun i => (e4StatusMatrix.map (·.index)).contains i) = true := rfl

/-- Every absence-probed row carries the typed non-evaluable exit
verdict in the live matrix, computed from `e4StatusMatrix` itself, so a
shape change on a probed row breaks this receipt.  Row 5 also carries
that exit and is deliberately unprobed: its exit is internal, with
nothing external to watch. -/
theorem absence_probe_rows :
    ((mathlibAbsenceProbes ++ projectAbsenceProbes).map (·.row)).all
      (fun i =>
        ((e4StatusMatrix.filter (·.index == i)).map (·.shape)).all
          (fun s => s.any (fun v => v matches .typedNonEvaluableExit)))
      = true := rfl

/-- The exit rows whose citations are prose-only (1, 6, 7; a declared
scope list of this guard, since prose-only-ness is a property of the
citation text, not of the matrix data) each carry at least one
pinned-Mathlib absence probe. -/
theorem prose_only_rows_probed :
    ([1, 6, 7].map fun i =>
      (mathlibAbsenceProbes.filter (·.row == i)).length != 0) =
      [true, true, true] := rfl

/-- No probe list is empty: a guard with nothing to check would be a
silent green, which is the defect this module exists to kill. -/
theorem no_empty_probe_list :
    ((mathlibAbsenceProbes ++ projectAbsenceProbes).map
      (·.probes.isEmpty)).all (· == false) = true := rfl

/-- Every probe is nonempty lowercase alphabetic.  The scans lowercase
their haystacks, so an uppercase probe would be a dead probe, and the
Python sync regex recognizes exactly this character class. -/
theorem probes_lowercase_alphabetic :
    ((mathlibAbsenceProbes ++ projectAbsenceProbes).map (·.probes)).all
      (fun ps => ps.all (fun p =>
        p.length != 0 && p.toList.all (fun c => c.isLower)))
      = true := rfl

/-! ## The elaboration-time scan

Everything below runs on every build of this module.  A hit, a missing
presence citation, or a scan that saw too little of Mathlib throws,
and the build fails with the full list of findings. -/

open Elab Command in
run_cmd do
  let env ← getEnv
  let mods := env.header.moduleNames
  let contains (hay probe : String) : Bool := (hay.splitOn probe).length > 1
  let mut mathlibModules := 0
  let mut mathlibConsts := 0
  let mut findings : Array String := #[]
  for (md, i) in env.header.moduleData.zipIdx do
    let mname := mods[i]!
    let isMathlib := mname.getRoot == `Mathlib
    let isDependency := dependencyRoots.contains mname.getRoot
    if isMathlib then
      mathlibModules := mathlibModules + 1
      let mlow := mname.toString.toLower
      for p in mathlibAbsenceProbes do
        for probe in p.probes do
          if contains mlow probe then
            findings := findings.push
              s!"row {p.row}: Mathlib MODULE name '{mname}' matches absence probe '{probe}': cited absence has ended ({p.claim})"
    for c in md.constNames do
      let clow := c.toString.toLower
      if isMathlib then
        mathlibConsts := mathlibConsts + 1
        for p in mathlibAbsenceProbes do
          for probe in p.probes do
            if contains clow probe then
              findings := findings.push
                s!"row {p.row}: Mathlib declaration '{c}' (module {mname}) matches absence probe '{probe}': cited absence has ended ({p.claim})"
      else if !isDependency then
        for p in projectAbsenceProbes do
          for probe in p.probes do
            if contains clow probe then
              findings := findings.push
                s!"row {p.row}: project declaration '{c}' (module {mname}) matches absence probe '{probe}': cited tree absence has ended ({p.claim})"
  -- Presence citations: a stale citation path is a finding, never a pass.
  for pc in presenceCitations do
    if let some m := pc.module? then
      unless mods.contains m do
        findings := findings.push
          s!"row {pc.row}: cited module '{m}' is GONE from the loaded environment: citation path is stale ({pc.claim})"
    if let some d := pc.decl? then
      unless env.contains d do
        findings := findings.push
          s!"row {pc.row}: cited declaration '{d}' is GONE: citation is stale ({pc.claim})"
  -- Fail closed: a scan that did not actually see pinned Mathlib is a
  -- finding, not a green.
  if mathlibModules < mathlibModuleFloor then
    findings := findings.push
      s!"scan integrity: only {mathlibModules} Mathlib modules seen (floor {mathlibModuleFloor}): the guard did not scan the library it claims to guard"
  if mathlibConsts < mathlibConstFloor then
    findings := findings.push
      s!"scan integrity: only {mathlibConsts} Mathlib-origin constants seen (floor {mathlibConstFloor}): the guard did not scan the library it claims to guard"
  if findings.isEmpty then
    logInfo m!"E4 stale-absence guards GREEN: {mathlibModules} Mathlib modules, {mathlibConsts} Mathlib constants scanned; {mathlibAbsenceProbes.length} Mathlib probe sets, {projectAbsenceProbes.length} project probe sets, {presenceCitations.length} presence citations checked."
  else
    throwError ("E4 stale-absence guard FAILURES "
      ++ s!"({findings.size}):\n"
      ++ String.intercalate "\n" findings.toList
      ++ "\nA hit means a cited absence in QFT/InheritanceMatrix.lean has (probably) stopped being true. Re-review the row; do not delete the probe to make the build green.")

end OPH.QFT.E4AbsenceGuards

-- Axiom audit: the table receipts are ordinary kernel-checked rfl proofs.
#print axioms OPH.QFT.E4AbsenceGuards.guarded_rows_exist_in_matrix
#print axioms OPH.QFT.E4AbsenceGuards.absence_probe_rows
#print axioms OPH.QFT.E4AbsenceGuards.prose_only_rows_probed
#print axioms OPH.QFT.E4AbsenceGuards.no_empty_probe_list
#print axioms OPH.QFT.E4AbsenceGuards.probes_lowercase_alphabetic
