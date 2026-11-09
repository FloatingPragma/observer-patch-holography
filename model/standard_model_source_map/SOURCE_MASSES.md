# Source action mass coefficients

These are the static mean-feedback coefficients with canonical kinetic form and unit neutral Higgs norm.
They are not measured pole masses. No GeV conversion or physical family assignment is supplied.
Every complete duty policy is retained; no policy is selected by comparison with experiment.

| Duty | Internal channel | Ascending action coefficients | Ratios to smallest | Koide quotient |
| --- | --- | --- | --- | --- |
| always | up | 0.0001576421455, 0.04279127122, 0.9570510866 | 1 : 271.4456282 : 6071.035659 | 0.6971072253 |
| always | down | 0.0001401415555, 0.04369077654, 0.9561690819 | 1 : 311.7617497 : 6822.880469 | 0.6959509268 |
| always | charged_lepton | 3.17943645e-05, 0.02148604675, 0.9784821589 | 1 : 675.7816076 : 30775.33312 | 0.7675778354 |
| C | up | 0.0001530816675, 0.009877422387, 0.05687365707 | 1 : 64.52387506 : 371.5249383 | 0.5454087188 |
| C | down | 0.0001361768715, 0.009038222896, 0.05772976136 | 1 : 66.37120384 : 423.9322047 | 0.5556112847 |
| C | charged_lepton | 3.13009615e-05, 0.003753003469, 0.02806182667 | 1 : 119.900581 : 896.5164431 | 0.5797503556 |
| E | up | 0, 0.0004332398375, 0.1503181764 | undefined (exact zero) | 0.9032911354 |
| E | down | 0, 0.0003842220765, 0.1303982284 | undefined (exact zero) | 0.9023276182 |
| E | charged_lepton | 0, 8.97218395e-05, 0.1306927287 | undefined (exact zero) | 0.9502391008 |
| CE | up | 0, 0.0004284587455, 0.009657438293 | undefined (exact zero) | 0.7125717065 |
| CE | down | 0, 0.0003798778175, 0.008370012324 | undefined (exact zero) | 0.7104394085 |
| CE | charged_lepton | 0, 8.92321815e-05, 0.004075682882 | undefined (exact zero) | 0.7754397616 |

Exact interval endpoints and authenticated source provenance are in [particle_map.json](particle_map.json).
All policies have three zero coefficients in the minimal neutrino action and zero charged-current CP invariant.
These do not identify physical neutrino masses or reproduce observed quark CP violation.
The [feedback theorem](FEEDBACK_LIMIT.md) controls the declared weak-pulse limit without changing these coefficients.
The [Koide test](KOIDE_ATTACHMENT.md) excludes balanced positive-root Koide for all four charged-lepton rows.
A [different coherent program](KOIDE_SOURCE.md) produces a balanced family generator with an unselected phase.
See the [mass interpretation](MASS_READOUT.md) for the physical attachment and clock requirements.
