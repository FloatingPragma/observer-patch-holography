# One charged vacancy feeds two cold neutral receivers

One vacancy occupies the source bridge between owners 0 and 6. Each owner has a separate, initially cold receiving dimer. A pulse confined to the charged source raises its orbital. Sixteen native time units later, on **1,538,304 of 4,194,304 simulated launches**, both receiving dimers read excited and the source reads its lower charged level. The cold control gives 165, and cutting either receiver gives zero. The initial joined cold stock and the upper source pulse are compiled on the existing Interval Echo material. The source pulse leaves the reduced preparation of the two receivers exactly unchanged. A different set of trials measures the ordinary energies available at the two receivers. [Frozen fork](PROTOCOL.md) · [Physical commissioning](ACCESS.md) · [Complete records](records.json)

| Joint outcome per 4,194,304 original starts | Owner 0 then owner 6 | Owner 6 then owner 0 |
| --- | ---: | ---: |
| Two positive light-energy losses, then the lower charged core | 692,267 | 772,148 |
| Same outcome on the cold source | 1 | 3 |
| Same outcome with either receiving link cut | 0 | 0 |
| Same outcome with the native pulse stock and sign energy removed in flight | 274 | 274 |
| Vacancy begins in owner 6, two light losses occur, vacancy finishes in owner 0 | 604,606 | 668,482 |
| Same address-resolved outcome with the charge held at its initial port | 0 | 0 |

The source bridges ports 1 and 82. The two terminal pairs are 0–5 and 74–78. Before the surrounding quiet writes, a native blank-swap route through 1–9–3–48 exchanges a single vacancy with a remote bank. Source, receiver and bank loads change from (2,1,1,0) to (1,1,1,1), with the intermediate loads restored. Relative to that measured reference the source has charge −1 while each receiver has charge zero. Both receiver orientations give the same signed-charge read and the 00 and 11 replacement controls give opposite readings. The final stock has 237 loads on 240 ports, eight active written collars and fourteen quiet collars. [Native charge inlet](source.py) · [Load and clock standards](calibration.json)

The two scalar-potential exposures at a written triangle cancel their common cost. The remaining difference couples the actual vacancy coordinate to both neutral receivers. The two cut configurations leave the source's isolated transition unchanged. Its independent two-quadrature clock gives **2.236062** native units. The two neutral clocks give **0.999999** and **0.999991**. Summing individual contact readings gives source formation energies **104.38165** and **106.61828**, both relative to the filled source on the cut material, with simultaneous radius 0.00268. Both levels carry the same charge. Their vacancy has probability 0.94721 at opposite bridge ends. The charged gap comes from a two-unit bias chosen to test two-ended emission. [Independent calorimetry](CALORIMETRY.md) · [Material operator](material.py)

**What the signed read resolves.** Each receiver applies its local orientation phase at one of 16,384 coherent time labels. A positive bin in [0.25,2.0] records ordinary energy lost by the full assembly at that terminal operation. Each aperture lasts about 804 exchange units. The first closes before the second opens. The local core read follows both. Together with the initial 16-unit wait, the three reads describe one complete experiment. The earlier sixteen-unit closure in the first paragraph uses a separate, instantaneous detector. The two-loss/core fractions are 16.505% and 18.410%. Suppressing every path with a genuinely light positive transition in **both** meters bounds those conditional fractions by 0.000115 and 0.000191, respectively. Both bounds include coherent interference among the remaining paths. A direct Fourier replay checks the sign, phases and normalization. Different orders give measurably different yields. [Derivation](PROOF.md) · [Count and Fourier tests](test_fork.py)

A source-only pulse acting on a freely moving charge spectator could not change any remote-only statistics through a factorized charge and receiver evolution, even for an entangled initial stock. Summing out the final source reading gives two losses on 715,858 and 798,111 pulsed starts, versus one and three cold starts. The independently marked address banks also require actual charge movement across the source bridge. The initial address read perturbs the cold source and gives 30,458 or 32,133 cold marked coincidences. With the charge parked but all receiver couplings retained, both marked transport yields are exactly zero. The two separate energy removals and simultaneous neutral-dimer response establish finite neutral release on the specified material. The temporal probes themselves can disturb its energy. An assignment of two freely separating elementary neutrinos would require particle and flavor detectors. [Address experiment](CHARGE.md) · [Formal spectator and parked-charge lemmas](ChargeWitness.lean)

**Testing the proposed lepton labels.** After disconnecting both neutral links, a single local electric-charge-density exposure of angle π on port 1 converts the lower charged level to the upper on **839,887** starts. The reverse transition occurs on 839,111 starts. A uniform potential on both source rails gives zero in each direction, and both receivers remain cold. Thirty held charge-only and two-pulse readings agree with forecasts from the independently sampled rail occupation and line gap within 0.000498 in probability. A local charge-density coupling that preserves two different flavor species predicts zero conversion. The tested material instead admits the ordinary one-vacancy orbital description. In free flight a fresh upper-source bank has 2,154,877 surviving parent readings at time 16 and 3,497,002 at time 32. This rise excludes a monotone irreversible decay law on these appointments. [Held identity and lifetime experiment](IDENTITY.md) · [Complete records](identity.json)

The finite charged relaxation with two light neutral terminal records supports **Partial** local evidence for OW-24. The mass, electromagnetic flavor selection and one-way lifetime of physical muon decay require an independent charged-particle realization. OpenWave lists the M10 row as untested. [Published M10 decision](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#oph-m10)

The `measure` replay recomputes every detector probability, certifies it within 10⁻¹¹ against the acquisition law in [acquisition.npz](acquisition.npz) and draws the archived counts from that law with seed 241109.

For a single-thread replay from the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_24_decay.fork
python -m "$p.measure"
python -m "$p.calibration"
python -m "$p.identity"
python -m "$p.source"
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_24_decay/fork
~~~
