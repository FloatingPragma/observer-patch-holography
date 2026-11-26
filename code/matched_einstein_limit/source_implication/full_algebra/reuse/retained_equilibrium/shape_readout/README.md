# What the four observers measure when area falls

On a transfer of pairing at fixed reference survey volume, two of the six measured edges lengthen and four shorten. The rate of change of log length on edge 1–3 lies in [0.019763, 0.019929] and on edge 0–2 in [−0.020723, −0.020557] per unit increase in the twirl probability. These are finite differences across the two archived preparations. The very same observers measure a fall in total face area. [All six intervals](results.json) · [Original terminal counts and simultaneous uncertainty](../observer_banks/README.md)

Each endpoint contributes twelve directed return-time means, including timeouts. We keep the original clock-gain fit and finite settling, ticket and numerical allowances. The difference of the corrected log edge lengths divided by the separation of the two preparation settings gives each secant. All six intervals follow from one simultaneous 99% clock-and-head event for samples from the deposited numerical instrument. This is a retrospective reading of those synthetic banks. Its precision concerns that instrument. [Specified comparison](PROTOCOL.md) · [Independent interval arithmetic](measure.py)

**A shape change is resolved.** For a proposed common log dilation subtract its value from both extreme edges. Their separated intervals imply at least one absolute error of 0.020160, half the guaranteed spread. A static edge-dependent choice of length units cancels in each log difference. Choosing a dilation that varies with the source also leaves a common secant for all six edges.

**A flat realization exists.** Anchor station 0 at the origin. For stations 1, 2 and 3 form the Gram matrix with entries

$$G_{ij}=(d_{0i}^2+d_{0j}^2-d_{ij}^2)/2.$$

All three leading principal minors are strictly positive across each endpoint interval. The determinants exceed 274.92394 in the sixth power of the clock length unit. A Cholesky factor places three station vectors in Euclidean space with exactly these six lengths. Applying the same construction to any smooth positive-definite family gives a smooth flat realization. Even the two endpoints alone admit such an interpolation by convexity of positive-definite Gram matrices. Its tetrahedron has exactly the observed four triangle areas and Gram volume at each endpoint. This argument concerns the six enrolled lengths. Additional stations and other clocks provide independent tests of a common metric. [Five- and six-station tests](../../../../full_geometry/README.md)

On a smooth branch write the instantaneous affine displacement gradient as the sum of a skew matrix and a symmetric strain $E$. The volume derivative is $V\operatorname{tr} E$. Differentiating the cofactor that transports each oriented face gives $A_f'=A_f(\operatorname{tr}E-n_f^T E n_f)$ for its unit normal $n_f$. Thus at fixed volume, $A'=-\sum_f A_f n_f^T E n_f$. A shear can account for the full observed tetrahedron area change. At this resolution, identifying the response with curvature requires independent records of the radius, directions and boundary of a proper ball. This result does not assign any hardware coordinate to a station.

Run from the repository root with the project Python environment.

```sh
export PYTHONPATH=code:.
python -m matched_einstein_limit.source_implication.full_algebra.reuse.retained_equilibrium.shape_readout.measure
```

## A finite experimental control for shear

Twelve physically enrolled orientations of the same reference shape suffice to cancel any common first-order affine deformation. They can be specified from an independently calibrated orthonormal observer frame. Cycle its three axes and use the four sign choices with positive determinant. For any matrix $E$ their average of $R^T E R$ is $(\operatorname{tr}E)I/3$. Sign averaging removes the off-diagonal entries and cycling averages the three diagonal ones. The calculation holds for each face normal separately, so it holds for any shape and arbitrary face areas. [Exact check on every matrix unit](orientations.py) · [Tests](test_shape.py)

For each orientation measure $A_R,V_R$ under the same small source perturbation. Form the average of $\delta\log A_R-(2/3)\delta\log V_R$. The face and volume formulas above give **exactly zero** for any common affine deformation of all twelve independently registered figures. A rotation-independent additive response in that combination passes with unit gain. If each measured combination differs from the ideal one by at most $e_R$, the average error is at most $\sum_R e_R/12$. This supplies a control even for irregular tetrahedra. Performing the rotations on paper on a single set of measured ranges only predicts this zero. Testing it requires twelve independent physical enrollments on the same source with a common reference frame. A surviving signal can also come from nonaffine strain or instrument differences. Identifying curvature further requires proper-ball geometry and the scale and source controls of the small-ball theorem.
