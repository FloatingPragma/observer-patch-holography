# Source-clock frontier

This package records the non-promoting, route-neutral producer frontier for
GitHub issue #633. It is the physical-unit interface required by issue #594
after a dimensionless pole packet exists.

The package does not derive a physical reference transition or an SI value of
the Newton constant. It does five narrower jobs:

1. records #634 as the sole hard issue dependency;
2. permits any target-clean physical clock route that satisfies the generic
   contract;
3. retains the cesium-133 construction as an optional incomplete candidate;
4. proves the exact interval inversion used after a positive clock-gap
   interval is supplied; and
5. refuses every physical-unit or gravity promotion while a source readout is
   absent.

The exact conversion is

```text
epsilon_clk = DeltaE_clk / E_star
DeltaE_clk = h * nu_clk
E_star = h * nu_clk / epsilon_clk
```

For `0 < epsilon_lo <= epsilon_clk <= epsilon_hi`, order reversal under
positive inversion gives

```text
h*nu_clk/epsilon_hi <= E_star <= h*nu_clk/epsilon_lo.
```

The committed arithmetic receipt uses a simple rational synthetic interval.
It is a theorem witness, not a physical clock payload.

## Dependency boundary

The machine-readable dependency contract is:

| Class | Issues |
| --- | --- |
| hard dependency | #634 |
| optional cesium-route owners | #32, #34, #317, #318, #425, #522, #545, #546, #569, #633 |
| downstream only | #334 |

Alternative physical clock routes are allowed. None of the optional
cesium-route owners is a hard dependency of #633. Issue #334 may consume a
completed source energy interval together with its separate gravity-scale
identification to form the Newton-\(G\) composition. It cannot be an ancestor
of the source energy interval.

The route-independent open gates are:

| Gate | Owner |
| --- | --- |
| route-neutral dimensionless clock observable | #634 |
| admissible physical clock route and process receipt | #633 |
| target-clean clock-to-SI attachment | #633 |

## Optional cesium candidate

The cesium-133 candidate has five incomplete component slots:

| Slot | Classification | Candidate-route owners |
| --- | --- | --- |
| `R_U` | conditional declared packet, not a physical source readout | #32, #545 |
| `R_alpha` | absent | #318, #425 |
| `R_e_abs` | absent | #34, #546, #569 |
| `R_QCD_nuc_133Cs` | absent | #317, #425 |
| `R_atom_133Cs` | synthetic theorem fixture only | #522, #633 |

The provenance graph sends these five slots only to the blocked
`cesium_candidate_clock_gap` node. There is no path from a cesium candidate
component to the generic `dimensionless_clock_gap` or to a physical source
energy interval.

The exact cesium frequency, Planck constant, speed of light, and joule value
of the electronvolt are SI unit definitions. They enter only the final unit
chart. They are forbidden as selectors or ancestors of the dimensionless
source gap.

The legacy gravity checksum skeleton and the synthetic Feshbach artifact are
hash-pinned as non-source diagnostics. The provenance graph gives them no
path to the source energy interval.

## Run

```bash
python3 build_source_clock_frontier.py --check-byte-exact
python3 check_source_clock_frontier.py
python3 -m pytest -q tests/test_source_clock_frontier.py
```
