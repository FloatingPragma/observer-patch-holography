# An independently prepared neutral material pair

Set every port to one native load and every seam memory to blank. One actuator puts two loads at port $p$ and zero at a different port $q$. Relative to this material reference there is one excess, one vacancy and zero net charge. Preparing them changes the occupancy multiset. Free evolution conserves that multiset exactly. This experiment starts with an independently prepared material pair. The uniform unit-filled state supplies its formation-energy reference. [Full native conservation check](measure.py) · [Source, generator and receiver](neutral_pair.py)

## Finite formation-energy bound

Write $n_i=1+\delta_{ip}-\delta_{iq}$ and denote the physical degree of port $i$ by $d_i$. The undisturbed reference is an exact eigenstate of energy $E_{\rm ref}=\kappa\sum_i(6-d_i)$. Each equal-load seam has zero contrast and exchange energy. The 20-cell reference energy is 1,440 at $\kappa=8$. The squared seam differences in the pair configuration sum to $d_p+d_q+2\mathbf{1}_{p\sim q}$. The degree contribution above the reference is $\kappa[3(6-d_p)-(6-d_q)]$. Thus the diagonal agreement energy above the reference is

$$
\kappa[12+2(d_q-d_p)+2\mathbf{1}_{p\sim q}].
$$

Each exchange term $(I-\mathrm{SWAP})/2$ is positive semidefinite. Every sphere port has degree five or six. Consequently **every state of this blank-memory neutral-pair sector lies at least 80 native energy units above the uniform material reference**. The bound is size independent for boards retaining these degrees and the same coupling. It is a finite formation-energy bound with this choice of reference and native units. The isolated 12-port spectrum begins at 97.3786465 in those units. Holding the endpoint degrees fixed, adjacency raises the agreement energy by 16. This allows high-energy correlations even when low-energy attraction is absent. A comparison at different physical units or on a different background requires additional evidence. [Local spectrum and finite-board bound](pair_results.json)

## Crossing the boundary

Use the entire ordered basis of 57,360 positions $(p,q)$ on the 20-cell sphere. All 630 physical seams remain active. The generator exchanges each defect with a single-load neighbor and exchanges the two defects with each other when they touch. Independent checks project all 132 ordered columns on one cell from the full native local contact and compare 480 signed single-defect columns on the 20-cell board with the registered one-load generator after the unit-filling degree correction.

Each of the three ports in cell 0 that face another cell supplies a station. Use its lowest-index degree-five internal neighbor as the second source port. Prepare both sharp charge orientations and their equal even and odd superpositions. After an appointed free flight the facing cell reads the twelve possible positions of the excess or quiet and the same thirteen outcomes for the vacancy. There are 65,536 simulated trials per source, exposure and law, including every quiet or separated outcome. Sum the outcomes with both defects on a receiving seam to test transmission of an adjacent neutral pair. The comparison removes the agreement coupling everywhere while keeping the same carrier, preparations and ports. [Terminal counts and exact probabilities](pair_results.json)

At exposure 8 the largest calculated probability of an adjacent double arrival over the twelve original preparations was 0.0000185318. The corresponding launch put the vacancy on gateway 1 and the excess on port 0. The pair was still adjacent somewhere on the board with probability 0.994478, while the vacancy reached the receiving cell alone with probability 0.00135841. Reversing the launch orientation gave a global adjacency probability of 0.239300. In the matched vacancy-first exchange-only control the probability of adjacent joint arrival in cell 6 was 0.00829383. The native simulated early cohorts recorded eight adjacent joint arrivals in total across the tested source settings and exposures.

Follow the strongest early station to two later appointments, alongside its reversed preparation and changed-law controls. At exposure 32, adjacent joint arrival has probability 0.0000448707 with the vacancy at the gateway and 0.0000307684 with the excess there. Both exchange-only probabilities are 0.00756260. A direct Chebyshev propagation from launch covers the later appointments and reproduces the independent short-flight propagator at exposure 8. These calculations establish a small joint arrival channel for both constituents. Terminal adjacency alone allows separation and rejoining between the observations. The finite energy bound and the measured joint arrivals characterize a material dipole with a positive formation cost on the prepared medium. A freely propagating massless polarization remains a different target. [Earlier survey, held arrivals and complete records](pair_results.json)

Run from the repository root with a fresh output path.

```sh
.venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/photon_identity/native_emission/neutral_pair.py --output /tmp/native-neutral-pair.json
```
