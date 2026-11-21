# An excitation that stays inside one open cell

Four interior ports suffice to prepare a one-load state that stays inside its original echosahedron for all free times. On the 80- and 320-cell screens, some cells also support a packet whose observable shape oscillates while its entire probability stays in that cell. The source uses local gate exposures before release. During flight every seam on the sphere is active. [Protocol](PROTOCOL.md) · [All three boards, controls and terminal counts](results.json)

## Why the load stays

The single-cell icosahedron has six pairs of antipodal ports. Each port touches one member of each of the five other pairs. There are three external gateways per cell, so at least three whole antipodal pairs have no gateway. Put amplitudes +½ on the two ports of one such pair and −½ on another. At every unused port the arriving amplitudes cancel. At each occupied port the adjacency operator returns the negative of the original amplitude. No amplitude reaches an external contact. Thus the cell graph Laplacian has eigenvalue 6 on this vector, including when the cell is embedded in the full sphere.

On the one-load sector the registered generator is 6κI + L/2. This gives energy 51 in native units at κ = 8. An integer matrix check verifies the four-port eigenvectors on **all 420 cells** in the three fixtures. Keeping all external links and removing the uniform contrast offset would leave their shapes unchanged. Their confinement comes from interference on the fixed port geometry. The other row 4 experiment tests a pair whose binding *does* depend on the contrast interaction. [Native generator](../../../contrast_memory_spherical_v1/candidate.py) · [Two-load binding experiment](../row_04/README.md)

| Complete screen | Cells with two / five independent compact modes | Best certified home-weight floor for an unfiltered sharp port |
| --- | ---: | ---: |
| 20 cells | 20 / 0 | 1/3 |
| 80 cells | 56 / 24 | 5/8 |
| 320 cells | 216 / 104 | 5/8 |

The last column is a lower bound at every time for one localized port launch in a suitable cell. Its value is the squared projection of that sharp preparation into the complete local invariant subspace. In the oscillating cells even a uniformly random sharp launch among all twelve ports has mean protected weight 5/12. This gives an unfiltered way to populate the subspace. Unit retention applies to the coherent compact preparations.

## What the observer sees

The source starts from a single occupied rail. Partial exchanges and routed full exchanges distribute its amplitudes within the cell. Addressing the degree-reference term on nongateway ports sets the relative phases; this is supplied preparation control. A detector can read a single port by checking absolute written contrast one on two of its incident seams, undoing each query. On the one-unit sector both queries return yes precisely when the load is at their common port. The code checks the selective operation using the native 99-state echo. All reads occur on separate copies at the terminal appointment. [Source, detector and full-screen evolution](experiment.py)

* In cell zero the exact four-port state kept unit home weight at all five appointments through time 64. On 320 cells the same initial four-port **density** with all positive phases retained only 44.49% at time 64. A source with 10% of its weight moved to a gateway retained 92.73%, above its all-time 90% bound.
* In 24 cells of the 80-cell sphere and 104 cells of the 320-cell sphere, an antipodal gateway pair also permits two distinct compact energies, 49.381966 and 51.618034. Releasing their equal coherent superposition produces a local occupation beat with period **2π/√5 = 2.809926** native time units. In the selected 320-cell installation the patch saw **1,950/8,192** clicks at launch and **108/8,192** half a beat later, consistent with predicted probabilities 0.236803 and 0.013197. After 16½ periods the same local reader gave 115 clicks. For the incoherent mixture its probability stays at 0.125. The 20-cell screen has no cell with that gateway arrangement.

The excitation is a prepared, discrete compact mode. Its free persistence and internal motion have local readouts. The runs do not identify electron charge or spin, a spontaneously formed electron, or a relativistically propagating version of this mode.

From the repository root run `python -B -m model.candidate_laws.interval_echo.openwave.row_04_compact.experiment`.
