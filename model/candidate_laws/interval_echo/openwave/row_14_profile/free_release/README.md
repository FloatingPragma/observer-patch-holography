# Which charge leaves the neutral source?

Open every contact of the federation after preparing the four-defect neutral source. At the first three registered appointments each of the three ports just outside its owner acquires **positive mean material charge**. The signed second moment about the launch owner is positive at all four appointments on the 20-, 80- and 320-owner boards. Reducing exchange by half gives a distinct, smaller response throughout the tested window. These are all-orders flight bounds from independent source measurements. The complete freely evolving four-body board was not sampled. [Finite-board proof](PROOF.md) · [Clock, pulse counts and bounds](results.json)

| Exchange-clock wait t | Native signed second moment / t² | Half exchange, same source and readers |
| ---: | ---: | ---: |
| 1/8192 | [0.230893, 0.267964] | [0.058003, 0.066711] |
| 1/4096 | [0.213479, 0.285377] | [0.053943, 0.070771] |
| 1/2048 | [0.177517, 0.321339] | [0.045573, 0.079141] |
| 1/1024 | [0.101044, 0.397812] | [0.027840, 0.096874] |

Distances are owner links. A separate local exchange dimer at the 1/4096 setting gave **1,042 flips in 2³⁶ simulated starts**. On its calibrated first half-period it reads a duration in [0.0002332836, 0.0002595380]. Propagating that entire clock interval gives [0.211256, 0.287600] for the native signed second moment divided by the actual squared duration. At 1/1024 the aggregate sign is resolved while the bounds for the three individual ports include zero. These time and population bounds have at least 99% joint sampling coverage on the simulated calibration sources.

## What each owner measures

The prepared cell contains two excesses and two vacancies. The rest of the board is filled. A single laboratory controller stops free flight at its chosen appointment. The terminal observations can occur sequentially on the held board. Each owner interrogates eleven internal seams that connect its twelve ports in a tree. A native echo writes a pair's signed difference in its two coin branches. Reading the signed value and undoing the echo realizes the corresponding diagonal sharp projection. A calibration checks all nine input columns and five signed outcomes of the 99-state seam.

From its eleven differences an owner reconstructs its port readings up to a common constant. At most four of its twelve ports are defects, so subtracting the most frequent level sets the filled host to zero and yields each signed material load. Both parties can share their local reading of an enrolled gateway. The reader recovers all 4,579 local words compatible with the total inventory. Its inputs are the actual differences recorded on its own seams, together with the calibrated majority-host convention. No reference loads need occupy the freely evolving board. [Local reader and complete census](reader.py)

Let Qᵢ be an owner's signed total and dᵢ its chart distance from the emitter. The coefficient of k² in Σᵢ E[Qᵢ] cos(kdᵢ) is −Σᵢ dᵢ² E[Qᵢ]/2. At t=1/4096 this spatial slope divided by t² lies in **[−0.142689, −0.106739]**. Its sign reflects the excess material launched past the first gateways. The source is net neutral at every appointment. The chart combines all twelve sites of one owner into one spatial bin. Changing the origin or using a body-centered ruler defines a different observable. This resolves an outward signed material response on release. A neutron electric radius calls for an intrinsic electric scattering response on a free baryon, and the neutron–proton mass ordering calls for a rest-energy comparison on identified partners. The full OW-14 criterion is Open.

## Checking the transport mechanism

The three outgoing swaps have distinguishable destinations. To leading order the signed accumulation at a particular receiving port is t²/4 times the initial signed charge at its gateway. The measured input ranges from +0.330414 to +0.334734 across all three conservative calibration intervals. To reach the jth neighboring shell a defect requires at least 2j−1 contacts. A bound on all such paths controls the later terms without truncating the board. Closing exchange preserves every sharp material word and yields exactly zero. With two planted excesses at gateways 1 and 3, the independent fourth-order word calculation on the 240 ports gives a signed moment per t² of +0.4999994 at t=1/4096. Exchanging these excesses and the two vacancies gives −0.4999994. Both Taylor truncation radii are below 0.0000138. [Independent word evolution](words.py)

A separate receiver measures the incoming load when only the three bridges are pulsed. The pulse at π/2 gives **6,752,920 excesses and 32,339 vacancies** at a uniformly chosen receiving port, plus 33,661,698 units and 68,679,029,779 source refusals. Half exchange at the same setting gives 1,979,381 excesses and 9,439 vacancies. A full 2π native pulse gives zero arrivals of either sign. These are draws from the exact three-bridge outcome law. The incident populations come from separate sign-sensitive and elastic probe responses at the individual gateways. An independent forecast uses the sharp source counts and calibrated dimer flips, then compares them with each receiving bank. The isolated pulses and the all-contact flight are separate experiments. The programme binds both source files by digest. [Simulation and bounds](measure.py) · [Frozen appointments](PROTOCOL.md)

From the repository root, reproduce the saved evidence with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_14_profile.free_release
nice -n 15 .venv/bin/python -B -O -W error -m "$p.measure"
nice -n 15 .venv/bin/python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_14_profile/free_release
~~~
