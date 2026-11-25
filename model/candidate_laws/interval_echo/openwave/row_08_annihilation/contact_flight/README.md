# What leaves the ring when its two walls disappear?

The source and both guides read their port loads once, at the end of free Interval Echo evolution. The observers compare the ten dimer records from the same launch. Four contacts report whether the source ring recovered its order. Two report whether energy is parked at the attachments. Four distinguish the near and far contacts inside the two guides. The fixed dimer guards make each mismatched contact an energy reading of 32.5 exchange units at native contrast 8. The initially reversed pair of dimers adds two such contacts, costing 65 units. The ten dimer orientations can flip freely. Their exchange energy is a separate observable. [Preparation and readout](PROTOCOL.md)

At the published exposure of 32 exchange units, 2,065,048 of 67,108,864 wall launches gave an ordered ring. Among these, **2,061,173 (99.81%) showed exactly two mismatches elsewhere on the fixture**. In 963,259 launches the attachments were clear and each wing held exactly one mismatch. Both mismatches landed at the far contacts in 646,526 launches. These simultaneous material records account for a pair of 32.5-unit contact energies at the remote contacts, with the source free throughout flight. They do not track trajectories of individual walls or measure an outgoing radiation mode. [Complete local word counts](results.json)

| Same exposure and initial guide preparation | Ordered ring | One mismatch in each wing, both attachments clear | Both mismatches at far contacts | Far pair among ordered rings |
| --- | ---: | ---: | ---: | ---: |
| Registered, two initial walls | 2,065,048 | 963,259 | 646,526 | 31.31% |
| Registered, ordered sham | 67,074,812 | 159 | 72 | 0.000107% |
| Disconnected guides, two walls | 3,431 | 0 | 0 | 0% |
| Zero contrast, two walls | 22,089,847 | 3,452,219 | 316,129 | 1.43% |
| Zero contrast, ordered sham | 33,270,740 | 3,588,444 | 170,204 | 0.512% |

All denominators are 67,108,864 original launches per setting. The simultaneous 99% interval for the far pair among ordered wall launches is 31.18%–31.44% on the registered law, 1.42%–1.45% at zero contrast and 0%–0.294% with both attachments cut. At zero contrast each mismatched contact costs 0.5 unit. A far pair itself occurs under both laws. The contrast changes the energy, its concentration conditional on source recovery, and the background. At exposure 8 the registered wall preparation gave 51,224 far pairs among 114,346 ordered rings. These finite guides can return energy to the source.

## An all-exposure limit for complete disappearance of the material walls

With both wings connected, zero mismatches on the entire fixture fixes exactly two opposite orientation words. Let P project onto them and E₀ denote their diagonal energy, 325 in the simulator's offset convention. The guarded Hamiltonian is a contact-count diagonal plus ten single-dimer flips of amplitude −1/2. Each of the two words has ten distinct one-flip neighbors, and the two neighbor lists are disjoint. Thus P H P = E₀ P and the norm of P H (1−P) is √10/2.

For the wall launch let E be its mean energy. The four ring flips have mutually orthogonal images of amplitude 1/2. The guides begin in their boundary-conditioned eigenstates, so the initial energy standard deviation is exactly 1 and is conserved in free flight. The two initial ring walls raise the diagonal by 65 relative to E₀. The six wing flips together can lower the guide energy by at most 3, since each has norm 1/2 and each contact cost is nonnegative. Hence E−E₀ ≥ 62. At any exposure, projecting (H−E)ψ onto the two ordered words and using the triangle inequality bounds their total probability by

$$
\|P\psi\|^2 \le \left(\frac{1+\sqrt{10}/2}{62}\right)^2 < 0.001734.
$$

This bound covers all original wall launches at every readout exposure. It concerns complete ordering of this guarded ring and both three-dimer guides in a joint load reading. On the registered run, 2,738 of 67,108,864 wall launches gave that record at exposure 32. The measurements support disappearance of the source walls and displacement of their contact energy into neutral material. The full OpenWave criterion requires conversion to radiation. Two- and three-photon positronium channels are follow-ups in the platform's wording. [OpenWave criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L110)

## Replay

From the RER root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_08_annihilation.contact_flight.measure
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_08_annihilation/contact_flight
```

The replay checks all sampled words and the matrix identities used in the bound. The test decodes the recorded words independently against the physical contact addresses. The observer receives the port orientations at readout. Exposure labels are settings of the apparatus in simulator units.
