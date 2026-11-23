# How much of the trapped orbital bit can an active boundary read?

An independently launched second load makes a trapped orbital population visible from the next cell under registered Interval Echo. A single receiving cell still misses a variation for which an internal observer measures a large matter-entropy slope. Combining two exterior gateways recovers that variation in ordinary joint occupation records with a coherent incoming probe. This obstruction holds for any terminal read of that cell at the two tested exposures, including a phase-sensitive read. The measured quantity so far is the entropy of a finite prepared matter subsystem. Identification of a boundary observable with proper area and a common fixed proper volume requires separate geometric evidence.

## Occupation and entropy measurements

The source is the same two-dimensional dark orbital sector as the [passive experiment](../hidden_entropy/README.md), with identical one-load energy 51 for both modes. An independent probe launches at one of the three receiving gateways. Each receiver reads the occupied subset of its twelve ports. All 79 outcomes enter the area-candidate tests. A three-axis regional tomography on independent preparations measures the source entropy. At central Bloch vector (0,0,1/2), the exact local entropy slope along z is −0.549306 nats per unit z. The signed finite difference at the actual displacement 1/8 is −0.554035. Reversing a real or imaginary coherence while holding the population fixed preserves regional entropy.

| Cells | Regional z slope | Binary capacity, nats per release | Largest held x error | Largest held z error | Regional y slope in phase challenge |
|---:|---:|---:|---:|---:|---:|
| 20 | −0.554008 | 0.000007156 | 0.00103 | 0.00067 | −0.554054 |
| 80 | −0.554046 | 0.000007258 | 0.00213 | 0.00162 | −0.554038 |
| 320 | −0.554034 | 0.000007254 | 0.00062 | 0.00231 | −0.554061 |

Each reported regional slope has standard error 0.0000147. The binary capacity optimizes the prior on the two pure z preparations at exposure 1/2, with gateway address retained. The changed law with exchange halved gives corresponding capacities between 0.000000676 and 0.000000681 nats. Capacity belongs to this specified channel, so varying the input mixture leaves the capacity fixed.

At 2³⁵ simulated releases per individual setting the measured z response of neither gateway-averaged empty-cell probability nor full terminal Shannon entropy reached the prespecified five-standard-error calibration threshold. The strongest reached 2.15 standard errors. The noiseless terminal Shannon derivatives along z are −0.00001281, −0.00000724 and −0.00000724, respectively. They have the same sign as the bulk entropy change. Interpreting this Shannon entropy itself as area with a positive constant cannot cancel the bulk derivative at this reference. On 80 and 320 cells the noiseless gateway-averaged empty-cell derivative is below 9×10⁻¹⁵ in magnitude. The full typed gateway records nevertheless discriminate the two pure populations.

The original three-axis exterior inverse has component errors as large as 0.33 on held inputs. The noiseless third weighted singular value is below 3.3×10⁻¹⁵ and its direction is the imaginary coherence y. The first two singular values are both above 0.0068. The x/z fit in the table is a diagnostic chosen after seeing the weak direction. It uses the original signed calibrations for its weights.

## Varying the invisible phase

A fresh three-to-one mixture of the two circular orbitals has Bloch vector (0,1/2,0). Changing y from 3/8 to 5/8 changes its entropy by approximately −0.13851 nats. Both mixtures are stationary under free one-load repair. Independent local tomography gives the signed slopes in the last column. For the original incoming gateway probe, the compressed model predicts total variation below 4.7×10⁻¹⁷ between their complete exterior occupation laws at every enrolled setting. The configuration-cutoff allowance for this comparison is 2.19×10⁻¹³ per release.

We also traced over every port outside each receiving cell. The resulting quantum trace distance bounds any one-shot terminal measurement of that cell. With a sharp gateway launch the largest distance is 4.61×10⁻¹⁶. Replacing the sharp load with the equally weighted coherent superposition of the gateway and an internal neighboring port, with relative phase +i, gives at most 3.15×10⁻¹⁶. Both values carry the same cutoff allowance and unquantified ordinary floating-point roundoff. For each launch and board, fresh signed circular calibrations selected one real or circular seam reader for a held test. The largest held difference was 1.42 standard errors. These statements refer to the tested preparations, receiving cell and exposures 1/4 and 1/2.

## Two exterior gateways recover the phase

We coherently prepared one incoming material load across two distinct receiving gateways with relative phase +i. Both receiving cells then reported their occupations, including quiet. The joint register has 301 exclusive outcomes. Across the three possible pairs, the predicted largest trace distance is about 0.00000526 under coherent incidence and 0.00000330 under an equal incoherent mixture of the two launches. Their direct joint occupation distances are about 0.00000386 and 0.000000000306 respectively. The primary pilot at 2³⁵ repetitions did not resolve its selected routed outcomes on held data.

For the replication we froze the first two gateways, exposure 1/2 and a routed circular reader. We used 2⁴³ simulated releases per setting, drawing aggregate counts directly. A fresh pure-source calibration selected an occupation event set before acquiring the three-eighths and five-eighths phase mixtures. The routed test required one load in the two receiving cells and zero in the auxiliary ports of their exterior connecting path. All other branches count as the complementary outcome. Its circular wave was compiled on the connected read graph using native splitters and degree-reference phases.

| Cells | Coherent incidence, joint occupation difference | Standard errors | Routed read, incoherent incidence | Routed read, coherent incidence |
|---:|---:|---:|---:|---:|
| 20 | +0.000003831 | 38.2 | 6.6 standard errors | 15.8 standard errors |
| 80 | +0.000003756 | 37.6 | 6.6 standard errors | 14.5 standard errors |
| 320 | +0.000003876 | 38.8 | 7.8 standard errors | 16.2 standard errors |

The incoherent occupation calibration selected an empty event set on all three boards. Joint phase analysis still resolved the changed interior preparation under incoherent incidence. The coherent control is needed for the strong ordinary-occupation read found here. Gate labels, preparation phases, appointments and the routing circuit belong to the instrument. No entropy stationarity or boundary-area normalization follows merely from this external distinguishability. [Pilot design](TWO_GATEWAYS.md) · [Replication design](JOINT_REPLICATION.md) · [Held joint counts](joint_replication.json)

## Native flight and truncation

The two-load generator is the registered unordered binary pair Hamiltonian up to a common shift. Each move has amplitude −1/2. The adjacency interaction and the exchange diagonal use the full board's degrees, including at the computational boundary. For each experiment we retain every pair whose sum of unconstrained site distances to its two input supports, minimized over the assignment of particles to supports, is less than 25. Any reachable omitted configuration requires at least 25 physical single-port moves. The primary runs retained 28,644 of 28,680 configurations on 20 cells, 141,023 of 460,320 on 80 cells and 182,345 of 7,370,880 on 320 cells.

In the interaction picture of the diagonal, each hopping operator has norm at most 6 by its maximum absolute row sum. The series on the retained and full boards agree on the source through 24 hops. The iterated integral remainder with its last factor the exact propagator is bounded by (6t)²⁵/25! for each evolution. Adding the two remainders bounds state-vector error by 1.093×10⁻¹³ at exposure t=1/2. Partial traces, mixtures and terminal measurements inherit this accuracy in trace distance. Two source hypotheses incur at most twice that allowance. Numerical linear algebra contributes separate floating-point error. On 20 cells a direct full-board replay of all original gateway effects differs from their restricted versions by at most 2.78×10⁻¹⁶.

From the repository root, run each module without arguments to regenerate and check its record. The accompanying tests independently check that the restricted operator equals the principal submatrix with its full diagonal, that the circular read has the claimed sign, and that the archived signed scores and selected settings follow from the integer counts.

```sh
python -m model.experiments.interval_echo_einstein_station.boundary_calibration.measure
python -m model.experiments.interval_echo_einstein_station.boundary_calibration.phase_receiver
python -m model.experiments.interval_echo_einstein_station.boundary_calibration.two_gateways
python -m model.experiments.interval_echo_einstein_station.boundary_calibration.joint_replication
python -m pytest -q model/experiments/interval_echo_einstein_station/boundary_calibration
```

The other records and tests in this station use distinct preparation and readout assumptions. This study probes one hardware cell on each board. It supplies no independently measured proper ball volume or geometric boundary-area variation across refinement. G4 and the Einstein consequence retain those premises. [Source and active-probe design](PROTOCOL.md) · [Phase challenge](PHASE_CHALLENGE.md) · [Quantum receiver comparison](PHASE_RECEIVER.md) · [Joint aperture protocol](TWO_GATEWAYS.md)
