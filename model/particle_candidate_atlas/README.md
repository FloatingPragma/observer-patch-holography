# Particle candidates and their source connections

The [large SVG atlas](../../docs/research/OPH_PARTICLE_CANDIDATE_ATLAS.svg)
shows every elementary Standard Model target, its antiparticle convention,
the exact internal-current match, the charged action calculations and their
comparison with experimental references. The supplied Standard Model benchmark
and measured GeV inputs are separate
from the calculated source actions. Physical source masses and defect-to-particle
identifications are unproved.

The graphic is generated from the authenticated
[charged-field receipt](../charged_lepton_response/results.json) and its
[candidate map](../charged_lepton_response/CANDIDATE_MAP.md). Experimental
reference editions and the discrepancy audit are separate inputs. The SVG
contains clickable proof and measurement references and a machine-readable
input-hash manifest. Measured inputs are confined to labeled calibrations and
the conditional benchmark; no primitive source dynamics are inferred from them.

The [reference table](reference_values.json) preserves CODATA2022 and
PDG2026 editions, quark mass schemes and confidence levels. The generated
[comparison audit](audit.json) retains exact rational source intervals and
residuals. It distinguishes numerical enclosure width from unspecified
physical coordinate and observable-matching errors. The atlas authenticates
the child receipt and its complete input pins; it does not repeat the child
scientific experiments. It executes bounded interval and closure calculations
and labels the larger operator results as retained theorem evidence. All owning
Python and proof files and their direct source inputs are hashed. The
[compatibility analysis](COMPATIBILITY.md)
states the complete-controller, full-CAR, finite gauge/Higgs and compact-defect
joins, their exact exclusions and their physical assumptions.

```sh
python -m model.particle_candidate_atlas --write
python -m model.particle_candidate_atlas --check
```
