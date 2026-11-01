# Five completion issues: source results and remaining joins

This pass advances **#1138, #1139, #1142, #1143 and #1145**, connecting
operational clocks, particle protection and spectra, interacting fields,
and the modular premises of gravity. The [completion program](COMPLETION_PROGRAM.md)
records the revised priorities and source-bound evidence.

The readable simulator specification, observer and scheduler implementation,
repair law and memory law are admitted model inputs. Their consequences
must be proved; implementation uniqueness is not a prerequisite. Agreement
with nature and compatibility of the particle, field and gravity
constructions still require their own evidence.

| Issue | Completed result in this pass | Remaining physical join |
| --- | --- | --- |
| [#1138: references and clocks](https://github.com/FloatingPragma/observer-patch-holography/issues/1138) | [Regional clock comparisons](../model/modular_clock_comparison/README.md) are source-operated with all inputs and failures. Lean proves cocycle and frequency-intertwiner implications; exact source arithmetic excludes nonstationary uncorrected common-parameter intertwiners on the default dimer. | One continuing field reference, scalable comparisons and common physical energy/time normalization; active correction disturbance stays counted. |
| [#1139: protected defects](https://github.com/FloatingPragma/observer-patch-holography/issues/1139) | [Charged doublet protection](../model/charged_defect/README.md) gives a maximal internal qutrit under complete block repairs, with a full-Fock outside outcome and reference-compatible exterior. | The uniform interacting field breaks this local protection. A localized invariant under common propagation/repair and its physical boundary charge remain missing. |
| [#1142: interacting fields](https://github.com/FloatingPragma/observer-patch-holography/issues/1142) | [On-shell two-body lattice scattering](../model/contact_scattering/ONSHELL.md) has complete ac wave operators, a unitary scattering operator and an exact eight-node contact spectrum reaching zero energy. | Common physical vacuum and asymptotic detectors, chiral/taste control, dynamical gauge fields and physical interacting refinement. |
| [#1143: particle spectrum](https://github.com/FloatingPragma/observer-patch-holography/issues/1143) | The same protected carrier has an exact charge/gap census. Its triplet states have identical gap and charge, excluding three distinct mass generations under any common scale or chemical potential. | Nondegenerate masses and spin/current/kinetic matching on the same protected excitations, followed by independent physical comparisons. |
| [#1145: Einstein premises](https://github.com/FloatingPragma/observer-patch-holography/issues/1145) | [Regional modular expansion](../model/regional_modular_expansion/README.md) constructs an exact eighth-order operator with error below 1.917e-9 and a uniform entropy-tangent bound, including the noncommuting boundary correction. | Plaquette reference preparation, geometric modular weights, Lorentz stress, entropy-area stationarity, positive Newton normalization and controlled coarse-graining. |

## Exact results and domains

The charged carrier retains 264 internal/spatial states. Its complete branch
commutant is `M1 + M3`, proved through all 720 raw repair slots, effect
separation, a connected movement tree and an odd loop. The protected
triplet is internal spin one. Its supplied decoupled-dimer gap is
`2*(sqrt(17)/8-1/16)` for every internal state. A full-Fock completion retains
the outside-sector outcome. For the **four local charged-pair columns on
two adjacent dimers**, uniform four-site field leakage has squared operator
`(3 P_singlet + P_triplet)/512`, with squared Frobenius norm `3/256`.
Every nonzero state in that local carrier leaks at first order. This is
not an exclusion of arbitrary spatial superpositions or other geometries.

The scattering theorem applies the trace-class theorem to the actual
antisymmetric fixed-total-momentum fiber. Completeness concerns absolutely
continuous subspaces; dark point states and the bound pair remain separate.
At zero total momentum the contact density is
`E^2 A((E/2c)^2)/(4*pi^2*c^3)`, with `0<A_n<=1` and a rational geometric
tail. The interacting density has the same quadratic threshold. This is
an exact infinite fixed-spacing two-body result, not a finite-volume
extrapolation. Zero model energy is not a derived physical vacuum energy,
and no photon or physical massless particle is identified.

Regional clock comparisons use twenty pulses. Global-to-regional correction
at modular parameter one half uses 2,826 pulses with ideal operator error
at most `11/32768`, retaining all sixteen input columns and all failures.
Exact source characteristic-polynomial and regional moment checks
authenticate the spectral formulas. All 243 frequency comparisons have
separation greater than `1/10000`; an uncorrected intertwiner reaches only
the six-dimensional stationary regional algebra. Active cocycle correction
supplies a different valid operation and retains its state disturbance.

The regional modular construction retains all 256 plaquette Fock states and
all sixteen regional states. At `beta=1/32`, exact coefficients through
degree eight approximate `-log X_A`, the logarithm of the unnormalized
reduced Gibbs weight, with operator error at most
`9582815241821/5000000000000000000000`. For the normalized generator `K_A`,
the approximation includes the exact scalar `log(Tr_A X_A) I`.
That scalar cancels in trace-zero entropy variations. The same bound multiplies the trace
norm of every trace-zero entropy variation. The cubic boundary term does
not commute with bare local energy. A rational witness proves distance
greater than `2.7e-8` from that bare generator plus any scalar. The finite
temperature/region are specified; the norm bound may grow with volume.

## Claim disposition

The [premise map](../model/PREMISES.md) records **38 realized scoped
attachments among 43 reviewed domains**, covering all 84 canonical rows.
The four additions are FA-40 through FA-43. This is a map of explicit
proofs and domains, not a completeness score. The generic clock implications
have new Lean proofs; the concrete operator and rational source attachments
have separate exact and numerical checks.

No canonical physical premise, full completion gate or measured particle
mass is newly established. No eligible empirical comparison is added.
The accepted CP-00/01/13 foundations and frozen executable class remain
unchanged. Block protection cannot be combined with uniform field dynamics,
and dimer preparation cannot be relabeled as plaquette preparation.

The flagship and seven companion paper sources share the stronger results
and retain the failed identifications. Historical receipts and frozen
predictions keep their original scope. This is local work; GitHub issue
closure, publishing and CI repair are deferred.

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error \
  model/modular_clock_comparison model/charged_defect \
  model/contact_scattering model/regional_modular_expansion
python -m model.completion.finite_attachment_map --check
python tools/build_claim_views.py --check
python tools/build_completion_program.py --check
python tools/build_model_theorem_map.py --check
```

Each package documents its retained receipt replay. The new Lean module is
`Lean/QFT/ModularClockComparison.lean`. Targeted validation does not claim
that every repository or CI suite passed.

## Validation completed

- The four scientific packages pass **145 tests** with warnings treated as
  errors. Their retained receipts replay successfully.
- The broader selected regression batch covered **587 cases and 22 subtests**.
  Its only initial failure was the stale frozen-review binding; the audited
  repair passed the affected test and nineteen finite-map tests. A separate
  forty-two-test class/map review also passed without weakening a validator.
- `lake build QFT.ModularClockComparison QFT` completed successfully
  (**8,443 jobs**). The added module contains no `sorry` or new axioms.
- All **1,414 model artifact bindings**, claim and premise views, theorem-map
  inventory, public claim surfaces, and completion-program checks pass.
- The complete postdiction-ledger builder reran successfully against the
  final reviewed sources, including its existing numerical and Lean checks.
- The flagship and seven active companion PDFs rebuild successfully. All
  eight pass the existing warning gate and manifest checks; their changed
  pages were rendered and visually inspected. The release manifest is a
  local preview, not a published release.

The frozen class digest and observation exposure remain unchanged. The
dated independent review refresh binds the four new finite attachments
without accepting a new physical gate.
