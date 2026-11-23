# A register event that times a moving charge

**Multiplying the written guide Hamiltonian by the register's positive event effect puts the geometric hierarchy into three propagating charged thresholds.** Comparing the lowest orbital of each register character gives native energies **0.000105754442, 0.0218668794 and 0.367758199**. Their ratios are **1, 206.770316 and 3477.472837**, with Koide quotient **2/3**. The same rule multiplies each branch's propagation speed and its dispersion curvature by the corresponding register weight. The three branches fail the common-speed relativistic mass-shell test and cannot decay into one another under this generator. [Design](DESIGN.md) · [Source](measure.py) · [All 48 momenta, nine bands and signed launches](results.json)

## The coupling being tested

Use the single-load event operators and analyzer of the historical 432-state register. The three-dimensional compression of its tick has amplitude

$$
M=\frac{32I+16\sqrt2(zR+\bar zR^{-1})}{96\sqrt3},\qquad
z=e^{i\theta},\qquad \theta=\tfrac1{12}+\tfrac16-2\cdot\tfrac1{72}=\tfrac29.
$$

Here $R$ cycles the three register addresses. The coefficients are the historical tick's compressed integer matrices. Its exclusive count supplies the angle. The numerical experiment reconstructs this compressed amplitude, checks it on the Fourier characters, and forms the event effect $E=M^\dagger M$. It uses no measured particle masses. The historical single-load register supplies both balance and phase. A prepared measurement effect supplies a Hamiltonian rate only under the comparison postulate below. [Pinned register, literal compression and exclusive count](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_register_round/PROOF.md) · [Register reconstruction](measure.py)

Let $G$ be the existing signed-load guide Hamiltonian, with route contact weight 4, degree-five energy 8, and neutral energy zero. The comparison generator is $G\otimes E$ on a charged sector. This promotes the complete register effect to the clock for the material contact diagonals, seam reversals and degree potential. A character with effect eigenvalue $w_j$ evolves by $e^{-it w_jG}$. Its entire orbital dispersion is $w_j\epsilon_b(k)$, where $b$ is the original branch. Sharp register character preparation and readout are supplied apparatus. The neutral guide vacuum remains at its original zero in this comparison. [Guide construction and orbital dispersion](../bands/README.md)

| Register character | Effect eigenvalue | Orbital 1 energy range | Orbital 2 energy range | Orbital 3 energy range |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.0000603005590 | 0.000105754–0.000241202 | 0.000482404–0.000723607 | 0.00110026–0.00120601 |
| 2 | 0.0124683656 | 0.0218669–0.0498735 | 0.0997469–0.149620 | 0.227500–0.249367 |
| 0 | 0.209693556 | 0.367758–0.838774 | 1.67755–2.51632 | 3.82611–4.19387 |

The register has three characters and the route has three orbital bands, yielding nine branches. Comparing like orbitals is essential. Simply choosing the three smallest thresholds picks the three orbitals of character 1. The independent direct nine-dimensional eigenvalue calculation agrees at every allowed momentum, with error below $5\times10^{-15}$ on the initial run. The absolute thresholds are in native generator units. No calibrated rest energy in keV follows. [Per-branch energies and direct checks](results.json)

## Do the same loads move?

The verified 48-cell route contains 144 distinct ports. Prepare a lowest-band packet at momentum $\pi/2$ with its opposite signed end prepared 24 cells away. In the full 20,592-state allowed pair space, both ends and all attached written heads evolve. The same comparison register clocks this entire fixture. Swapping the two endpoints tests a negative projectile with the same spatial preparation. This experiment does not prepare a separate register character on each end.

| Character | Mean projectile speed at shared native time 4, cells/time | Group speed at $k=\pi/2$, cells/time | Lowest threshold × rest curvature |
| ---: | ---: | ---: | ---: |
| 1 | 0.0000610148 | 0.0000619896 | 0.00000000483039 |
| 2 | 0.0126161 | 0.0128176 | 0.000206518 |
| 0 | 0.212177 | 0.215567 | 0.0584130 |

The packet has finite momentum width, so its mean velocity averages the group velocity over the prepared distribution. Either charge sign gives the tabulated motion. Running each character for a native time $4/w_j$ also transports the prepared load across four cell spacings. These different elapsed times and the common-time measurements are both recorded. [Signed pair releases](results.json)

## Why the thresholds do not behave as a common relativistic mass family

For a shared clock and cell spacing, a relativistic low-momentum dispersion would satisfy $\mathcal E_j(k)^2=\mathcal E_j(0)^2+v^2 k^2+O(k^4)$. Its rest-energy–curvature product $\mathcal E_j(0)\mathcal E_j''(0)$ must equal the same $v^2$ across the family. In this clock model, both factors are proportional to $w_j$. Their product ratios are **1, 42,753.964 and 12,092,817.333**. One unit conversion for the whole apparatus preserves this discrepancy. Each character projector also commutes with the full free and pair generators. A charged orbital can scatter while the register character remains fixed. An autonomous weak transition, independent neutral daughter readout and an absolute energy standard remain open. The companion incidence experiment tests a different coupling designed to place the register in a mass term. [Incidence and mass-term comparison](../dirac_read/README.md) · [Original family ratio comparison](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_cycle/physical_comparison.json)

## Replay

Run `python -B -W error measure.py --output /tmp/clocked-charge.json` from this directory with NumPy and SciPy installed and a fresh output path. The receipt records digests of the experiment, its design, the shared dispersion implementation and the physical route.
