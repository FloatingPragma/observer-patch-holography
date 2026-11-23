# Two matched material rays at a stationary source

A written seam with 32 native units of energy changes both the source-cell entry and the relative spreading of freely released unit loads. Double loads with the same free propagation cross the same source as they cross a sham. A second source with identical stored energy gives sham transport to both tags. A separate local echo measurement predicts which flights respond.

## How the instruments agree in free flight

Each cell holds a stationary unit load on a nongateway port. Eight native memory preparations hold its five seams closed and close one seam for each of the other six non-neighbouring ports. Every mobile port loses exactly one exchange seam. All 30, 120 and 480 intercell links remain open on their respective hosts. The complete mobile orbitals have 220, 880 and 3,520 possible addresses.

On a blank seam adjoining an empty port, a load of magnitude q exchanges with amplitude −1/2. Its diagonal term is 1/2 + κq². Summing the open seams and the registered port term gives the same free operator L/2 for both tags plus 5κq² times the identity. Here L is the Laplacian of the active ports. The offset cancels from every terminal-position record. The two tags therefore match for every free input, release time and detector on this prepared background.

The selected energy insert replaces a guarded seam between a resident and an active gateway. Separate copies of this source undergo the inverse native echo, followed by native squared-contrast and reversal-parity reads. The two 32-unit inserts have the same sharp outcomes before a tag enters. With a unit tag at the gateway the even-coin insert measures zero energy. With a double tag it measures 32. The odd-coin insert reads 32 for either tag. Checking all 99 local dimensions gives a −32 diagonal shift at the gateway for the first encounter. All three other encounters have zero shift. The calibration uses local endpoint loads and source records and requires no flight data.

## Held terminal counts

One ray starts across the native bridge. Its companion starts on an adjacent port of that same outside cell in a separate preparation. At two native time units the central placement gives the following percentages. The separation read asks whether the two rays finish under different cell owners, using the first of two companion launches.

| Cells | Unit-tag entry, sham | Unit-tag entry, even | Double-tag entry, even | Unit-tag separation, sham | Unit-tag separation, even | Double-tag separation, even |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 27.586% | 0.068% | 27.584% | 50.076% | 27.675% | 50.068% |
| 80 | 22.490% | 0.051% | 22.492% | 31.662% | 7.465% | 31.660% |
| 320 | 22.494% | 0.051% | 22.493% | 31.670% | 7.463% | 31.665% |

Across all nine placements the even-coin insert reduced unit-tag source-cell entry by 22.44–27.52 percentage points at this appointment. Across the eighteen ray pairs it reduced different-owner outcomes by 18.62–27.56 points. All five appointments and both companion launches remain in the receipt, together with counts at each possible terminal-owner crossing distance. The largest observed departure from sham across the odd-coin control for both tags and the even-coin double-tag control was 0.00990 percentage points. The exact terminal laws in these controls coincide at every time. The 810 simulated primary frequencies have a simultaneous 99% Hoeffding half-width of 0.01057 percentage points each, or 0.02114 points for a difference. The largest observed error against the local calibration's flight prediction was 0.00646 points. The receipt uses 536,870,912 independent synthetic shots per entry count or pair histogram.

The result isolates a source-dependent material scattering channel. Intercell crossing counts measure coarse changes in relative transport. A curvature tensor would require a common proper-geometry read, and an Einstein response would require an independently reconstructed stress tensor on that same geometry. Equal seam energies constrain one scalar component of the preparation. Coin parity distinguishes the microscopic states and could change other components. This test leaves a possible common gravitational response in other preparations or at longer wavelengths open.

## Replay

From the repository root, run these commands with the scientific Python environment.

```sh
python -m model.experiments.interval_echo_einstein_station.native_tidal_response.measure
python -m unittest model.experiments.interval_echo_einstein_station.native_tidal_response.test_rays
```

The replay checks the saved [raw counts and exact outcome probabilities](results.json) against the [frozen protocol](PROTOCOL.md). The tensor-product test applies the unreduced native seam law with all three memories present on a four-port chain and checks the predicted flight and contact amplitudes.

The [stationary seam census and shared routing control](ROUTING.md) cover all sharp stationary seam-memory rays and demonstrate a common change in flight under a supplied re-routing.
