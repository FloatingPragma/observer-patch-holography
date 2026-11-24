# Can two local probes distinguish a shared unit?

OpenWave lists charge quantization as untested for M10. On the registered Interval Echo law, two contact probes can distinguish a single vacancy shared between two distant cells from independent local sources with exactly the same mean response. Any measurement confined to the enclosing seams at a single free-flight appointment has identical outcome probabilities for the compact vacancy and the undisturbed host. This gives a specific obstruction to assigning the vacancy an electric charge through a Gauss flux measured on those seams. [OpenWave criterion and review](PROTOCOL.md) · [Native calculation](measure.py) · [All joint records](results.json)

## A shared vacancy and the enclosing charge

Start with one load per port and blank collars. Prepare a single vacancy in an equal coherent superposition of two compact spatial modes, one in cell zero and one in the lowest-address farthest cell. The two source cells are five, eleven and twenty-three owner steps apart on the 20-, 80- and 320-cell boards. Preparation of this state and the local probes are supplied services. All seams and degree references operate freely during the wait.

Let the input in each cell have amplitudes given by the first rational null vector of the local adjacency plus the identity, with zero amplitude on all three gateways. All supported ports have five neighbours. Applying the full graph Laplacian gives exactly six times the same vector. The free one-vacancy Hamiltonian has kinetic term equal to half this Laplacian. Relative to the unit-filled host the potential on its support is 32, giving energy 35 throughout the wait. The two localized alternatives and their coherent superposition are exact stationary modes; the largest numerical component residual for the vacancy was below 4 × 10⁻¹⁵. [Finite arrays and tests](test_split.py) · [Native seam law](../../../contrast_memory_spherical_v1/SPEC.md)

For each source cell measure first its three outgoing seams and on fresh copies the six seams leaving its one-step neighbourhood. Each signed echo reads the load contrast on that seam. Every cut gives the all-zero word with certainty for a vacancy on either side, the shared vacancy, and the undisturbed host. A vacancy placed at the first inner endpoint of a cut gives a minus-one there. There is also an exact limitation on *every* instantaneous measurement supported by these cuts. Every external seam joins two gateways. Both gateway loads equal their host values, and free evolution preserves blank collars. All of these registers factor from the compact state in their host state. Their reduced density matrix stays identical for the three vacancy preparations and the host at every free appointment. This holds for the boundary of any union of complete cells, even with coherent joint readout of the boundary. An intervention that couples into the interior over time is a different experiment.

## What the two internal probes can see

At each station, a fresh probe on two rails is split, touches one sample port, then recombines. An echo on the probe rails gives one of five signed contrast outcomes. All three ports belong to the native board. The touching and reference rails have zero vacancy amplitude in both source modes. Conditional on its cell holding the vacancy, the sampled port carries one quarter of its probability. Each contact lasts 1/32 of an exchange time. We retained the 25 joint outcomes for each of the four choices of the two probe signs. [Contact readout](../../../../experiments/material_inertia/pair_spectrum/PROBE_PROTOCOL.md)

The comparison uses an independent fair choice of compact vacancy or undisturbed host at each station. Each station then has exactly the same reduced state on its local detector as in the shared preparation. In the independent comparison, zero, one or two occupied stations are physically possible, with probabilities 1/4, 1/2 and 1/4. The shared vacancy occupies precisely one station. This tests whether single-station measurements of a half mean can identify the joint source.

For two positive probes score the terminal contrasts +1 and −1 by their signs and all other bins by zero. The table gives the mean product of the two scores. The initial bank used 4,194,304 launches per setting; only the 20-cell board resolved a difference at the specified simultaneous threshold. After inspecting this result we repeated the same measurement with a fresh seed and 268,435,456 launches per setting. The simultaneous 99% sampling half-width for any of the three differences is at most 0.000460. These are draws from the finite-model probabilities.

| Cells | Shared vacancy | Independent half-mean preparations | Difference |
| ---: | ---: | ---: | ---: |
| 20 | 0.096366 | 0.099116 | 0.002750 |
| 80 | 0.096400 | 0.099209 | 0.002809 |
| 320 | 0.096290 | 0.099260 | 0.002970 |

The exact total variation between the two full joint distributions for positive probes is 0.0013731. Removing the contrast energy from the measured contact, with the source, splitter, echo reader and ports fixed, reduces it to 0.000001967. The joint signals discriminate two sources that single-station histograms cannot separate. An incoherent mixture of the two exclusive vacancy locations gives the same contact histories as the coherent source. The independent half-mean control supplies whole vacancies intermittently. A persistent fractional excitation with its own dynamics would require a further preparation and test.

## Can the law select a unique elementary amount?

On an empty host the very same spatial modes are stationary for one load and for a doubled load. The respective energies are 51 and 195. At a port with weight 1/4, the same signed echo resolves a whole unit or a double against a reference neighbour. For example, the 20-cell single-mode sample contains 1,048,159 readings of +1 in 4,194,304 launches. The double-mode sample contains 1,048,535 readings of +2 and zero readings of +1. Both amounts persist indefinitely in a prepared localized stationary mode on every tested board. This excludes selection of single loads as the only stable local port amount on this blank material.

The reason extends to the free seam terms. Within the blank collar the generator only swaps loads or assigns energy. On the written subspace it acts as the compression of the memory generator to the complement of the native write range. In the projector onto that range, each input contributes a superposition of two balanced output orders. Both orders have the same two-port multiset. The projector and the memory generator therefore preserve the counts of 0, 1 and 2. The numerical check of both generating matrices uses all 99 local states. Addressed echo *writes* can balance unequal amounts and are separate interventions.

As a challenge to the detector, release a doubled load on an isolated seam at exposure 1/16. A terminal zero contrast at total load two means one unit on each side. Under the native term the probability is zero. A symmetric conversion term of strength 8 can couple a doubled load on either side to the separated configuration on the same two ports. Its predicted probability is 0.080864; 337,964 of 4,194,304 launches gave the split. This confirms the read can see conversion when the dynamics permits it.

The written-field extension previously measured a correlated boundary-head flux for a transported load. It introduces a different material hop and a different free field generator. On the compact preparation of the registered law, the interior contact correlations and the absence of any single-time boundary signature coexist. Recovering an electric topological unit requires a preparation whose enclosing records carry a conserved class with a separately measurable electrical response. [Written-field charge](../electron_field/README.md) · [Independent winding and flux comparison](../quark_winding/README.md)

## Replay

From the repository root, use the project NumPy, SciPy, SymPy and NetworkX environment.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_01_charge.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_01_charge
~~~

To save the full simulated record separately, append `--output NEW.json` to the measurement command.
