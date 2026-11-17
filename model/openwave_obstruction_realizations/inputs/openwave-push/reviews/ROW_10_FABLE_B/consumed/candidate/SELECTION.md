# Does the theory repository select the neutral basis on the carrier? Search record

Checkout `reverse-engineering-reality` at fd461f1d (working tree with uncommitted changes in other files; every
file used is pinned by its own SHA256 in `SOURCE_CONTRACT.json`). Rule: `PREREGISTRATION.md` section 6.

## Searches (case-insensitive, text files `*.md`, `*.tex`; `*.yaml` in `claims/`)

Pattern: `neutral register|klein|majorana|tri-?bimaximal|golden[- ]ratio|pmns|residual (symmetr|subgroup)|half-turn|
fifth-turn|third-turn|neutrino|icosahedral flavou?r|flavou?r symmetr|mixing angle` with `grep -rIihoE`.

| Term | `model/` | `paper/`, `flagship/`, `claims/` |
| --- | --- | --- |
| neutrino | 337 | 155 |
| Majorana | 226 | 30 |
| Klein | 29 | 48 |
| half-turn | 23 | 7 |
| PMNS | 17 | 5 |
| mixing angle | 14 | 4 |
| third-turn | 4 | 7 |
| neutral register | 4 | 4 |
| flavour symmetry | 2 | 2 |
| golden ratio | 0 | 23 |
| tri-bimaximal, residual symmetry or subgroup, fifth-turn, icosahedral flavour | 0 | 0 |

A first search that included `*.json` and `*.py` did not finish under the machine's load and was stopped at 16:35
(one term completed: "neutral register", 10 files in `model/` with those types). The bodies of `*.py` and `*.json`
were read only for the packages named below.

## What the hits are

| Hit | File | Reading |
| --- | --- | --- |
| "Klein" | Klein–Gordon in several spellings, Kaluza–Klein, a bibliography key, a German title word; "Klein-four section obstruction" (spin statistics artifact, `claims/assumption_dictionary.md:272` and its historical copies) and "conjugation on the five Klein-four subgroups identifies the latter faithfully" (`paper/screen_microphysics_and_observer_synchronization.tex:308`). Lines that hold "Klein" together with a neutrino, flavour, mixing or Majorana term: 0 | group facts; no neutral basis |
| "neutral register" | `model/lepton_register_round/README.md:3`: "Two seam objects on one face of the carrier and a neutral register carrying the face third-turn form one 432-state tick"; `model/finite_family_clock/PROOF.md:16` | the charged-lepton generation register is called neutral; it is the charged side |
| "golden ratio" | one line mentions mixing: the bibliography entry Everett and Stuart (2009), "Icosahedral (A5) family symmetry and the golden ratio prediction for solar neutrino mixing" | credited prior work of the particle paper: "A5 flavor potentials and residual-stabilizer selection", "Residual-symmetry flavor mixing and icosahedral axis geometry" (`paper/tex_fragments/PRIOR_WORK_DERIVING_THE_PARTICLE_ZOO_FROM_OBSERVER_CONSISTENCY.tex:8, 107`) |
| claim OPH-A5-REAL-AXIS-CABIBBO-NOGO | `claims/claim_registry.yaml:3762` | the 31 residual axes of the same carrier as a menu for the Cabibbo angle: "the narrowly defined hypothesis that the Cabibbo angle is directly an acute Euclidean angle between two axes in this menu is excluded". Quark side. Its numbers (6 + 10 + 15, sixteen rows, 20.905157°, 0.356822) are recomputed in `s5_supplement.py` |
| claim OPH-W5-STABILISER-POTENTIAL-BOUNDARY | `claims/claim_registry.yaml:3800` | "the C3- and C5-fixed loci are one-dimensional and every nonzero fixed point has a doubly degenerate spectrum … Thus stabilizer symmetry alone does not select two charged-family mass ratios". Charged families; the same fact limits a cyclic neutral group (PROOF.md, Statement 3) |
| PMNS | `model/openwave_neutrino_oscillation` (uncalibrated), `model/neutrino_mass_observables`, `model/majorana_defect_bridge`, `model/conditional_particle_spectrum`, `model/particle_candidate_atlas` (calibrated or conditional); `model/common_source/chiral_candidate/ACTION.md:169`: "their degenerate sector supplies no physical PMNS prediction" | see below |
| claims ledger, OPH-OBSERVER-MATTER-INTERACTION-WITNESS (`claims/claim_registry.yaml:18852`; packages per this side's sheet: `neutrino_mass_observables`, `neutrino_mass_response`, `native_flavor_orbit`, `native_flavour_transport`, `charged_lepton_response`) | `claims/claim_registry.yaml` | "Physical source identity, defect assignment, mass and neutrino mass-eigenstate assignment are null; the physical particle map is empty." |
| particle paper | `paper/deriving_the_particle_zoo_from_observer_consistency.tex:1405` | the packages "do not supply measured PMNS parameters, a physical neutrino mass or a muon lifetime" |
| flagship | `flagship/from_observer_consensus_to_standard_physics.tex:2216` | "A separate 288-mode Dirac source calibrates two neutral mass splittings, three mixing angles and one phase to the specified NuFIT 6.0 data edition." |

## The packages the task names

| Package | What it fixes for the neutral states | Enumerated subgroup? |
| --- | --- | --- |
| `model/openwave_neutrino_oscillation` | neutral coefficient Mν = V* Y V†: Y of `native_repair_mass` (a five-state memory on five port lines of one fifth-turn orbit, weights 12⁴ : 12³ : 12² : 12 : 1), V a quarter-turn pulse on the line of port 2; "The pulse, unit Higgs component h=(0,1), and equality of charged/neutral singular values are architecture choices made without a neutrino mass, angle, or mass-splitting input. This source is distinct from the analyzer-effect charged-lepton mass architecture" | no: Y commutes with none of the 59 non-trivial rotations (smallest commutator 0.0044); Mν is invariant under none (0.0035); its frame is 0.25 to 0.28 from the nearest subgroup frame |
| `model/native_flavour_transport` | the pulses on ports 2, 3, 10 and the 40 words; "Specifying its ports, duration and order is a declared flavour transport choice. Gauge covariance does not force it."; "It does not show that the unchanged repair memory selects this transport" | no; quark side; "the three minimal neutral directions remain massless" |
| `model/native_flavor_orbit` | "different charged spectra and nontrivial real quark mixing"; "All four source duty policies have exactly zero charged-current CP invariant" | no neutral basis |
| `model/neutrino_mass_response` | a propagating neutral composite and a sterile extension; "Lepton number forbids a Majorana term in that original law"; "These finite gaps have no proved relativistic rest-mass attachment" | no family frame |
| `model/lepton_register_round` | the charged frame: lines of the face third-turn (the product of the two seam half-turns of seams (2, 3) and (2, 10); faces (4, 5, 11) and (6, 7, 8)) | the charged side of every pair of menu H |
| `model/openwave_neutrino_calibration`, `model/openwave_neutrino_matter`, `model/neutrino_mass_observables` | the six measured parameters as inputs | calibration |

## Classification

S3: no statement selects a subgroup or a set of lines of the carrier for the neutral states. The one neutral frame on
file (first row of the table) is carried to the end in `s3_source_join.py`: on its own charged frame the angles are
(1.89 × 10⁻⁸, 1.764 × 10⁻⁵, 2.138 × 10⁻³); on the register's rows (0.0338, 0.4802, 0.9304).
