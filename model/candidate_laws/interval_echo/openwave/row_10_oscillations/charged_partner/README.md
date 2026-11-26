# What energy does a negatively charged daughter leave at the receiver?

A receiving observer transferred an arriving excess onto a unit-filled neighbour. On 510,459 simulated native launches the pickup changed from two units to one, the donor changed from one to two, and one vacancy remained on the local three-port material. The observer clocked this negative daughter without resolving its vacancy position. Changing the charged source arm changed the distribution of its measured energies at both receivers and both appointments. Each observer's cold and single-vacancy standards located the same three lines. [Frozen experiment](PROTOCOL.md) · [Native aperture and acquisition](measure.py) · [Full count banks](counts.json)

At owner 12 the negative charge occupies rails 152, 150 and 154. The excess on 212 moves to 207. At owner 17 the corresponding rails are 212, 207 and 206; the excess on 106 moves to 96. These are two independent copies of the 14-port guide. Every open contact is a native blank edge. A full pickup–donor pulse exchanges the sharp loads with unit probability. Both verifications and the total excess of the three daughter rails commute with the energy clock. During spectroscopy the outside contacts are held parked. The two sharp charges refer to their local unit-filled stocks. Charged inputs at the source, the original donor, both direction heads and the timing follow the [registered six preparations](../README.md).

The 4,096-path detector clock is referenced to a completely filled, stationary trimer. The registered sum of its two contacts on one vacancy has diagonal energies 8.5, 17 and 8.5, with exchange −0.5 on neighbouring vacancy positions. Reflection gives an odd line at 8.5. The even lines are (51 ± √297)/4. Independent simulated clock trials on the three sharp vacancy stocks resolve windows at 8.44158, 8.50000 and 17.05842 exchange units. The time aperture has a finite tail. Across all native captures, 95,175, 112,694 and 264,075 hits lay in these three windows, and 38,515 fell outside them. A clock on the fully filled trimer gave zero energy on all 4,194,304 commissioning starts. The initial and final material stocks, the Fourier delay clock and the addressed capture are supplied apparatus.

### Held charged-daughter records

At owner 12 and wait 768, the three addressed source arms gave 2,843 of 9,260, 17,609 of 24,553, and 5,337 of 9,358 captures in the highest window. These correspond to 30.70%, 71.72% and 57.03% conditional on the signed capture. Each bank began with 2⁵² simulated launches. Each of the 12 pairwise arm comparisons across the two sites and two waits resolved at least one of the three predeclared energy windows at simultaneous 99% confidence. None of the 12 within-winding comparisons was resolved. Opposite windings have the same remote density on this aperture by the retained-direction symmetry. The zero winding is a distinct preparation. Its unresolved differences have the precision of the enrolled data. [Count-only comparisons](readback.py) · [Remote symmetry](../theory/PROOF.md)

A sharp read of the daughter position before the clock changed the energy distribution in 17 of 24 native comparisons. Parking the second contact during the clock independently moves the lines to 0, 8 and 9; the enrolled native windows failed in all 24 comparisons with that apparatus. Halving only the neutral-flight contrast changed the measured distribution in all 24 comparisons with the original daughter clock. With pickup parked every capture refused. On the native banks 3,877,273 other transfers left the vacancy outside the three-port daughter. All categories, including undelivered precursors and refused sources, appear in the record. The total native exposure was 24 × 2⁵² original starts, about 4.72 jointly signed spectra per trillion launches.

The three addressed source arms therefore produce distinct occupation probabilities of the measured levels of the same negative charge on the same workstation. The shifted standards identify the levels as receiver orbitals. The observable establishes source-dependent charged-partner spectroscopy on a prepared material. Independently propagating charged daughter species, a weak-production flavour tag, and their masses remain open.

For a single-thread replay from the repository root:

~~~sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations.charged_partner
python -Bm "$p.measure" --acquire /tmp/oph-charged-partner.json
cmp /tmp/oph-charged-partner.json model/candidate_laws/interval_echo/openwave/row_10_oscillations/charged_partner/counts.json
python -Bm "$p.readback"
python -Bm pytest -q --pyargs "$p.test_partner"
~~~
