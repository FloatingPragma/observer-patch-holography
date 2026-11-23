# Can observers measure cosmic acceleration?

Observers can recover accelerated, coasting and decelerated expansion from a calibrated counter of proper four-volume and an independently calibrated clock on a finite layered population. All 24 estimates classified the prepared history correctly. A change in the counter's production rate made decelerating dust yield exactly the same count distribution as accelerating de Sitter space. Comparison with the independent clock exposed the substitution on both of the later intervals in all four panels. [Measurements](results.json) · [Experiment and rejection thresholds](PROTOCOL.md) · [Finite record-density theorem](../../../paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex)

## Reading volume and elapsed time

The source supplies actual read and write histories from the two archived golden populations. Observers enroll the 39 or 105 local sites appearing in a receiving station's inbox. Every enrolled site wrote exactly once at each of four epochs. The comparison adds two records at those sites. A calibrated marker produces independent Poisson counts at fixed proper four-volume density. Repeating the exposures in inverse proportion to the comoving cell volume and conformal tick keeps their joint flat reference volume fixed across both populations and cadences. A separate comoving clock counts ticks between epochs. The apparatus supplies its conformal cadence, density convention and time unit. Only the producer receives the expansion profile. [Archived carrier writes](../../../evidence/source_net_causal_poset/README.md) · [Producer and observer inverse](measure.py)

The observer takes the fourth root of marker density relative to its flat reference to infer scale. Changes in that scale across successive proper-clock intervals give a first and second derivative. The dimensionless deceleration parameter is the negative of scale times its second derivative divided by its squared first derivative. Both interior epochs were measured at conformal spacings of 1/8 and 1/16 on each population.

| Prepared history | Reference deceleration parameter | Range of eight readings |
| --- | ---: | ---: |
| Dust | +0.5 | +0.478 to +0.521 |
| Coasting | 0 | −0.052 to +0.052 |
| de Sitter | −1 | −1.014 to −0.964 |

The histories share initial scale and expansion rate. An exact static control has zero first derivative and an undefined deceleration parameter. Varying the dust marker density in proportion to the fourth power of the ratio of the de Sitter and dust scales preserves the de Sitter marker law at every site. The clock still records dust proper time. The inconsistent records gave negative inferred deceleration parameters until the clock check rejected the interpretation.

## Measuring the expansion rate

In fresh runs with three de Sitter amplitudes, the observer used only the first interval's scale ratio and elapsed time to estimate a constant Hubble rate. The two later scale forecasts agreed within 0.042% across both populations. The same forecast missed the last dust scale by 1.70–1.73%. With light speed one and the pure de Sitter relation, the observer infers a cosmological constant of three times the squared rate.

| Prepared rate in the supplied clock unit | Inferred rate, 39 sites | Inferred rate, 105 sites | Inferred cosmological constant, 39 and 105 sites |
| ---: | ---: | ---: | ---: |
| 0.50 | 0.49701 | 0.50020 | 0.74105, 0.75059 |
| 0.75 | 0.75112 | 0.75096 | 1.69255, 1.69182 |
| 1.25 | 1.25099 | 1.24954 | 4.69495, 4.68403 |

The prepared expansion rates, proper-volume law, clock law and de Sitter relation are comparison inputs. The inferred magnitudes have units of inverse clock time squared. A source prediction would select the physical normalization and the amplitude independently.

## Which readings does Interval Echo supply?

On the separate registered 20-, 80- and 320-cell boards, home and traveller shared their positive and quiet elapsed digits. Home counted 13 positive digits in 33 services and the traveller counted 19, with the same seed and guide on each board. Counting those positive digits as proposed volume marks reproduced the elapsed-counter increments and the coherent hand's triggers on all 198 readings. Counting every publication instead measures the number of scheduled services. A freely evolving phase hand is another available clock. A physical four-volume rate or a capacity tied to an independently measured horizon for the same process remains necessary to interpret these records cosmologically. [Public clock measurements](../interval_echo_clock/PUBLICATION.md) · [Capacity and horizon relation](../../../cosmology/oph_cosmological_vacuum_and_structure_formation.tex)

From the repository root, replay the complete receipt and check the forecast and controls.

```sh
python -m model.experiments.observer_expansion_reconstruction.measure
python -m pytest -q model/experiments/observer_expansion_reconstruction
```
