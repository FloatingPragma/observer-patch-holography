# A charged orbital loses the energy that an outgoing neutral orbital gains

A source owner prepares a negative vacancy in its **17.5** exchange-unit orbital and a neutral vacancy–excess pair in its **49** unit orbital. Every port of the receiving owner starts filled. After free evolution, the source can read the negative vacancy in its **16.5** unit orbital while the receiver reads the pair in its **50** unit orbital. At exposure 17,739, the two observers jointly record this transition in **648,729 of 4,194,304 simulated starts**. Preparing the lower charged orbital reduces the count to **15,102**. All local energies refer to the filled material on this bench. [Records](results.json) · [Executable apparatus](measure.py)

## What the observers do

Ten ports share a native blank-contact path, with six in owner 0 and four in owner 6. The unique bridge joins ports 1 and 82. Equal-load writes silence 33 off-path seams. A vacancy visits the six degree-five path ports to write their private reference at zero contrast. The emitter uses routed full exchanges to place two vacancies and one excess. A partial exchange and private phase pulses prepare the source energy orbitals. The two shutter openings close at unit load before flight. The candidate uses its registered contrast strength 8 throughout. [Pulse construction](../emission_channels/apparatus.py) · [Full-seam swap and silent-state checks](../conjugate_flight/measure.py)

Number the path positions from zero. The charged source uses positions 1–2. The initial pair occupies positions 4–5, and its receiving orbital is at 6–7. Each owner resolves the order parity in its own two-port window using an inverse native phase and exchange sequence and signed load reads. Each read checks a local unit-filled buffer and counts all outcomes, including a miss. On accepted starts the source has net charge −1 and the receiver has exactly one vacancy and one excess. The total inventory fixes all other ports to one. Independent copies supply a second detector, the load histogram over the whole receiving owner. Its two populations change only by current across the 1–82 bridge, verified by commuting each number observable with every other active term. At the stated exposure that detector reads one vacancy and one excess **1,419,965** times across all starts. This measurement integrates crossings between preparation and readout and does not distinguish individual vacancies.

| Preparation or control | Joint matching reads at 512 | Joint matching reads at 17,739 |
| --- | ---: | ---: |
| Upper charge and lower pair at source | 176,271 | 648,729 |
| Lower charge and lower pair at source | 13,424 | 15,102 |
| Close the only bridge | 0 | 0 |
| Prepare the outgoing states and reverse the experiment | 176,130 | 648,959 |

Each bank has 4,194,304 starts per exposure. The reversed measurement resolves the original source levels. The geometry and example at 17,739 were selected in an exploratory path-length and exposure scan. The complete distributions include nine appointments from zero to 32,768. These finite return and reverse probabilities give no exponential lifetime. [Full distributions](results.json)

## Where the unit of energy comes from

Local calibrations with the other excitation removed give a charged matrix with diagonal 17 and off-diagonal −½ and a neutral matrix with diagonal 49½ and off-diagonal −½. The lower orbital is even under order exchange and the upper orbital is odd. In the charge window the two squared load differences cost 16, the outer exchange contributes ½, and odd parity adds one. In the neutral window the squared differences cost 48, the outer exchanges contribute one, and odd parity adds one. The local parity meters read these levels independently. The full generator compresses to the sum of those same two local matrices on both the incoming and outgoing windows, including all nine bonds. The complete incoming and outgoing packets have mean **66.5** and variance **0.75** in squared exchange units. The matched local energy budget gives a one-unit drop and a one-unit gain. Because the windows leak into other positions, a regional read does not resolve a sharp global energy eigenstate on an individual run.

The coupling of these two specified windows is exact. Let the integer matrix be twice the free path Hamiltonian and let the two integer vectors be twice the normalized input and output. Their matrix-power moments vanish at orders zero through six and give **3185 at order seven** in the implemented pulse convention. The exponential transition amplitude therefore has a nonzero seventh derivative at zero. A floating-point recurrence or exposure selection cannot explain the existence of this channel. A bare lower neutral pair released on the same unit-filled, private-erased route has **exactly zero** amplitude in the upper receiving orbital at every exposure. Interchanging vacancies and excesses fixes that preparation, reverses the sign of the daughter, and commutes with its free flight. The surplus charge breaks this symmetry. [Exact integer certificate and independent exponential](test_transfer.py) · [Lean parity theorem](ObserverTransfer.lean)

The source supplies the pair’s 49 units of formation and orbital energy as part of the initial preparation. The charged-orbital change accounts for its additional unit in the resolved channel. The closed, real generator also supports reverse transfer. Identifying a heavy charged lepton requires a species preparation with a measured rest gap, a physical weak daughter tag and a decay law. The finite excited orbital here has no measured rest mass or weak tag. [Load inventory theorem](../neutral_ejecta/README.md) · [Comparison using a receiver with matter initially present](../resonant_delivery/README.md)

## Can neutral energy arrive in seam memory alone?

For each physical seam let F select unit loads at both endpoints. Compress its native term h to those loads, then extend the resulting head operator by identity on the loads to obtain B. The integer native matrices satisfy hF = BF and hB = Bh in all 99 local states. Every other term acts trivially on this head, hence B commutes with the full free generator even on a federation with shared material ports. Its spectral projectors commute with free flight. A final read with both endpoints filled measures precisely this conserved B for the seam energy. For any region whose counted seams finish with both endpoints filled, the joint support of their intrinsic energies at readout is confined to the support prepared initially. The path that matter takes between the two appointments can be arbitrary. Conditioning a record can reweight initially present spectral values. This proof concerns the energies of physical seam heads at filled endpoints. An excitation with transported matter, an energy-degenerate memory record, or a head still coupled to non-unit endpoints is outside the exclusion.

The filled-seam levels are 0 on blank and zero-contrast heads, 8 or 9 on unit-contrast heads and either 0 or 32/33 in the magnitude-two sectors. A native terminal echo distinguishes prepared energy 0 from 32 on the even-sign, magnitude-two subspace by whether the final physical loads are extreme. This confirms the distinction is readable on the specified analyzer states. [Complete integer matrices and energy analyzer](test_transfer.py) · [Lean terminal-energy identity](ObserverTransfer.lean)

From the repository root, use one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.level_transfer.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/level_transfer
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_11_spectrum/level_transfer/ObserverTransfer.lean)
```
