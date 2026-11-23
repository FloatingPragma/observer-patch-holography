# Three charged orbital bands on a physical ring

Three distinct coherent preparations move a signed load across cell boundaries on the same 48-cell guide. Each preparation occupies one of three dispersive bands associated with the three physical ports per cell. Their band-minimum energies have ratios **1, 4.561553 and 10.403882**. A collision with an independently prepared opposite charge can convert an upper orbital mode into a lower one. In isolated-end propagation the band index is exactly conserved. [Release and encounter records](results.json) · [Search design](design.md)

## Spectrum on the written guide

The fixture follows 144 distinct ports of the 320-cell federation with degrees 6, 5, 6 in each successive cell. It uses the proposed written-field extension at reference load one. Each allowed material contact has weight 4 and a charged degree-five port has local energy 8. The only cycle runs through all 48 cells, so the guide supports no elementary field flip. The neutral state has native energy zero. A sharp pair configuration contains one positive and one negative end and the consecutive arc of reversed heads between them. The executable checks local flux against the signed load for all 20,592 configurations. [Physical route](route.json) · [Model and carrier checks](measure.py)

Far from the compensating end, the charge Hamiltonian has diagonal entries 8, 16 and 8, internal hopping −4 and crossing hopping −4 times the cell translation phase. The only momentum-dependent term in its characteristic polynomial is proportional to cos(k). All three bands remain separated. Their extrema occur at zero cell momentum or the zone edge.

| Band | Energy at k = 0 | Band minimum | Momentum at minimum | Band maximum | Curvature at minimum |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.753789 | 1.753789 | 0 | 4 | 0.757464 |
| 2 | 12 | 8 | π | 12 | 2.666667 |
| 3 | 18.246211 | 18.246211 | 0 | 20 | 1.242536 |

Comparing all three at zero momentum gives ratios **1, 6.842329 and 10.403882**. The bottom and top zero-momentum energies equal $10-2\sqrt{17}$ and $10+2\sqrt{17}$. These are thresholds for a separated end in generator units. Their physical energy calibration and a common relativistic dispersion law remain open. The Koide statistic is 0.394540 at the band minima and 0.389849 at zero cell momentum. Removing the degree-five cost in a control gives precisely the three folded parts of the homogeneous chain dispersion, with crossings at the zone centre and edge. Its spectrum differs from the analytic folded chain by at most $9\times10^{-15}$. [Numerical bands and controls](results.json)

## Independently prepared signed pulses

Each source has a localized cell envelope with momentum width 0.23 and a fixed complex amplitude on each of the three ports. At a given central momentum the three source profiles are orthogonal to within $10^{-14}$. A low-band opposite charge is independently parked in a stationary packet 24 cells away. The joint wave is projected onto the allowed Gauss states and both ends then evolve freely. Prepare either charge sign, reverse the pulse momentum, and read the inclusive cell-load distributions. These full-envelope tests specify coherent initial wavefunctions. A second experiment below constructs finite sources from an already prepared sharp pair.

| Three-port source at k = +π/2 | Target band at time 4 | Signed projectile speed in cells per native time |
| --- | ---: | ---: |
| 1 | 99.0005% | +0.985746 |
| 2 | 99.0583% | −1.810491 |
| 3 | 99.9256% | +0.824745 |

Reversing momentum reverses all measured velocities. Exchanging the two charges gives the same speeds for negative projectiles under the verified exact symmetry. Across all fixed-profile tests the remote end moves by less than $2.5\times10^{-8}$ cells per unit time. Suppressing the intercell contacts gives speeds below $5\times10^{-16}$. Momentum-resolved pure-band packets separately reproduce the group velocities from the derivatives of the band energies. [All sign and momentum releases](results.json)

## An addressed circuit from one sharp signed pair

Start every run with the same excess and vacancy at route ports 18 and 90 and the intervening arc of reversed heads. Prepare each end in its own 39-port window, using two additional empty ports for phase control. Successive partial contacts distribute the initial amplitude. Each full positive-contact pulse swaps two adjacent addresses, including their seam write, with unit amplitude. A reversible train of swaps and two global degree-reference holds produces opposite phases at a selected port and an empty port of different degree. The two holds sum to the full degree period, restoring every spectator phase. Both circuit windows are disjoint throughout.

For each end the executable executes 38 partial contacts, 3,146 full contacts and 78 degree holds. This uses external isolation and timing of the existing terms during preparation. Free evolution then enables the full guide for four native time units. At positive pulse momentum the band-1, band-2 and band-3 projectiles respectively retain **98.8858%, 98.9453% and 99.9151%** weight in the selected band and move at **+0.982650, −1.805398 and +0.822750** cells per time. Negative input momentum again reverses all three. These circuits independently select the three orbitals of a pre-existing signed load. Preparing the initial written charged pair from the neutral vacuum is a distinct formation question. [Pulse synthesis and signed readouts](measure.py) · [All six finite-circuit releases](results.json)

## Conversion during encounters

Upper-band positive packets aimed toward a low-band negative packet initially ten cells away change their band populations. At native time 8 a band-2 projectile occupies band 1 with probability **70.5357%**, compared with **0.0566%** when the packets initially move apart. A band-3 projectile occupies the two lower bands with probability **34.6156%**, compared with **0.1149%** when moving apart. The record includes all times through 12. The faster apart-going pair begins to re-encounter around the far side of the finite ring late in that interval. Norm and total energy stay conserved to numerical precision.

Each isolated band projector commutes with the one-end Hamiltonian at all momenta. The full pair has a blocked contact at coincidence and reduced contact degree at adjacency. These boundary terms enable scattering between the two incoming modes. Charge number and active field content stay fixed. The experiment thus establishes several transportable charged orbitals and collision-induced transitions, while the identification of three charged lepton species and spontaneous decay with neutral daughters remains open. [Encounter series and conserved energies](results.json)

Run `python3 -B -W error measure.py --output /tmp/lepton-guide.json` from this directory with NumPy and SciPy installed. The executable uses the checked-in carrier and route only. The JSON includes SHA-256 digests of the measured code, input route and design.
