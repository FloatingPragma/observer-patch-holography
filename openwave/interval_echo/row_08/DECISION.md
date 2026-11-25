# Row 8 · Antimatter and annihilation

**Filing decision: HOLD.** OpenWave asks whether a particle and its conjugate annihilate into radiation. Positronium's two- and three-photon channels are listed as follow-up. In the simulated M10 ring, opposite order defects disappear on a measured fraction of launches and distant material receivers gain energy in one uninterrupted evolution. The same experiment still needs to establish the physical conjugacy of the incoming excitations and identify what propagates away as radiation. M10.5 nominates a one-load excitation for row 4. It has no established opposite electric charge partner. [OpenWave row 8 and the M5 precedent](https://github.com/openwave-labs/openwave/blob/bab6818b4a278ffef083c1e86c12134e5ff14ab5/MODELS.md#L110) · [Incoming-state comparison](../../../model/candidate_laws/interval_echo/openwave/row_04_conjugates/README.md)

## What was actually measured

At exposure 32 the opposite-wall run produced an ordered source ring in 2,066,989 of 67,108,864 original launches. Of those, 2,064,278 carried one excitation on each material wing. The one-per-wing event is 3.0760% of all launches and 99.8688% conditional on an ordered ring. It also produces measurable joint phase correlations. The passive receivers show energy uptake above sham. The resolved outgoing branch has a gap of 31.5 exchange units in the uniformly continued three-dimer-guide model. Before the ring is measured, 62.9782 of the original 65 excess energy units remain in the ring at exposure 32. These results support a material conversion and identify a concrete limit on the radiation interpretation. [Daughter spectrum and coherence](../../../model/candidate_laws/interval_echo/openwave/row_08_annihilation/daughter_modes/README.md) · [Energy budget](../../../model/candidate_laws/interval_echo/openwave/row_08_annihilation/README.md) · [Passive receivers](../../../model/candidate_laws/interval_echo/openwave/row_08_annihilation/passive_receiver/README.md)

The free blank-memory law has a stronger obstruction for literal excess and vacancy loads in a uniformly filled host. Its exact restriction is exchange plus diagonal contrast. Exchange permutes the two endpoint loads, and every diagonal term commutes with each count of zeros, ones and twos. The global free generator therefore preserves all three counts for every strength and exposure. A bare excess and vacancy cannot become the filled vacuum in this branch. Our symbolic check also detects a planted load-conserving conversion between (2,0) and (1,1). The obstruction applies to those literal loads. Collective wall numbers can change. [Exact check and instrument control](recheck.py) · [Reproduced numbers](receipt.json)

## Uninterrupted energy transfer

Two far receivers couple from launch until the terminal read. The ring evolves freely within the declared guards. Its sharp read and the two meter readings occur together at exposure 64. The source word has staggered order jumps +1 and −1. Both jumps vanish in either ordered terminal word. The full result includes all 16 ring words and all four receiver outcomes. [Measurement protocol](continuous_receivers/PROTOCOL.md) · [Population probabilities at all six exposures](continuous_receivers/results.json)

| Terminal event at exposure 64 | Opposite walls | Ordered sham |
| --- | ---: | ---: |
| Both far receivers excited | 0.00039741924 | 0.000000003836 |
| Ring ordered and both receivers excited | 0.00039735388 | 0.000000003807 |
| Mean energy stored in the two meters | 1.721103 | 0.000419 |

The meters take up 1.720684 excess units, or 2.6472% of the 65 supplied by reversing the two source dimers. The ring is ordered on 6.0358% of wall launches at this exposure. Cutting both source links gives the same unconditional meter distribution for both preparations. Removing the meter coupling gives zero uptake. No intermediate ring projection supplies energy or selects trials. These are model probabilities from unitary evolution. [Generator, energy check and independent propagation](continuous_receivers/measure.py)

## What would change the filing decision

Identify the opposite incoming excitations as a particle and its conjugate by a physical charge measurement or a demonstrated equivalent. The ring already supplies a local order observable for the signs. On the same collision, identify energy reaching the receivers in a freely propagating radiation channel with a diagnostic independent of the supplied material resonance. Positronium branching belongs to OpenWave’s follow-on test. The frozen bare excess–vacancy case and the collective-wall case have different dynamics. [OpenWave criterion and model comparison](https://github.com/openwave-labs/openwave/blob/bab6818b4a278ffef083c1e86c12134e5ff14ab5/MODELS.md#simplest-test-per-criterion)

From the RER root, run the following read-only checks.

```sh
python -m openwave.interval_echo.row_08.recheck
python -m openwave.interval_echo.row_08.continuous_receivers.measure
```
