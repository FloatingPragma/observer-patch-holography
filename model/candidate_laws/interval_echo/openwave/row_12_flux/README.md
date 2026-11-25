# Fractional material charge in a threaded belt

A period-three pattern on a ten-port belt can have a gap of four sites or a gap of one. When the gap moves three sites along a directly connected pair of wall configurations, one ordinary unit moves one site backwards or two sites forwards, respectively. An observer can test the signed fractions −⅓ and +⅔ by comparing the phases accumulated by each freely moving wall in opposite applied fluxes. A separate ordinary unit sets the sign and size of the probe. The flux is an external phase coupled to microscopic material transport. The belt, preparation, phase reference and detectors are supplied. [Wall preparation](../quark_winding/belt_walls/README.md) · [Protocol](PROTOCOL.md)

The two quadratures of the return phase give an energy splitting in exchange units. Ten zero-flux readings calibrate the dispersion. At four momenta and three held twist magnitudes, compare the maximum error of each forecast. Both fractional forecasts beat every integer competitor. Both exceed the frozen 0.01 tolerance for treating the wall as a rigid point charge. Removing contrast destroys the isolated, high-visibility wall response. [Calibration, counts, competing charges and changed-law measurements](results.json)

| Preparation | Fraction tested | Maximum observed splitting error | Best of charges −1, 0, +1 |
| --- | ---: | ---: | ---: |
| Single ordinary unit | +1 | 0.00006 | 0.00006 |
| Three-load wall | −⅓ | 0.01596 | 0.37831 |
| Four-load wall | +⅔ | 0.01373 | 0.33070 |

In the sharp-wall subspace the fractional flux phase holds exactly. The exchanges also take the wall temporarily through other configurations. Including just two-hop excursions predicts the observed splittings within 0.00769 and 0.01098 respectively. Including the entire complementary space gives an interval less than 0.000001 wide for each dominant splitting on 24 fresh combinations of momentum, twist and wall type. The simulated Ramsey frequencies lie within 0.00035 exchange units of those intervals. This prediction uses the specified microscopic Hamiltonian and detector. Assigning fractional *electric* charge to a quark requires an independently identified electromagnetic coupling to the same moving wall. [Virtual-motion calculation and held readings](dressing.json)

Over a full flux quantum the upper wall band stays separated from the other bands by more than 14.32 and 13.85 exchange units respectively. At the two tested momenta at least 97.75% of the wall preparations finish with a sharp wall, against less than one in a million for same-load compact preparations. Both give exactly the same final momentum shift. The change of gauge at the end multiplies each port word by a phase fixed by its total load and therefore carries both preparations from momentum index m to m−N modulo ten. A full-cycle momentum shift by itself cannot distinguish a quark from a compact cluster. [Flux-cycle estimate, joint momentum and wall counts](cycle.json)

Reproduce the count banks and compare one setting per wall with independent full-space propagation from the repository root.

```sh
python -m model.candidate_laws.interval_echo.openwave.row_12_flux.measure
python -m model.candidate_laws.interval_echo.openwave.row_12_flux.dressing
python -m model.candidate_laws.interval_echo.openwave.row_12_flux.cycle
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_12_flux
```

A calibrated local crossing counter resolves what the wall carries around a physical cut. The ideal four-load wall cycle transports two units per three gap revolutions. Free Interval Echo flight returns the identical wall word with odd current, and the isolated wall band has a certified half-flux spectral change of at least 0.00856873 exchange units. The rigid two-unit transport identity fails on this preparation. [Local current measurements and dressed-band spectra](winding_current/README.md)
