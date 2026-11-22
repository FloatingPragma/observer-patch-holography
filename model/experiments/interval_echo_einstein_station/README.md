# Can the boundary registers account for an observer’s entropy change?

A one-cell Interval Echo observer resolves an entropy change carried entirely by coherence between two of its ports. A double crosses a real intercell contact during the preparation. Every endpoint and memory on the three seams crossing the cell boundary has the same terminal joint state for all the tested preparations. Consequently any area reading derived from those instantaneous registers has zero variation, while the entropy of the cell changes. [Frozen protocol](PROTOCOL.md) · [Full probabilities, apparatus and block counts](results.json)

## What the patch measures

Cell zero on the 20-cell sphere uses internal ports A = 1, B = 0, spare port D = 5 and detector vacancy 7. A connects to port C = 82 in the neighbouring cell. A supplied sharp source puts a double at A, B or C with probabilities 4p/5, 4(1−p)/5 and 1/5. A native half contact splits A and B, then a π/3 contact transfers amplitude across A–C. The contact and detector use the registered 99-state contrast-memory generator and echo. The measurement uses four settings on fresh preparations. Two record A or B directly. Two combine an internal contact and a final occupancy read to reveal different components of the A–B coherence. At every setting the reader uses only ports of its own cell and retains quiet outcomes. The 80- and 320-cell embeddings exercise the same apparatus and local incidence. The finite samples below use one 20-cell realization of the design.

The conserved total number of doubles gives a local algebra with a vacuum block and a two-dimensional occupied block on the support of these preparations. If their weights are q₀ and q₁, its entropy is H(q₀,q₁) + q₁ S(ρ₁). Here the weights are 1/4 and 3/4. At p = 3/4 the two terms are 0.562335 and 0.434936 nats, giving 0.997272 nats in total. Both A and B have fixed single-port populations as p changes. Only the phase-sensitive readouts see the variation. The modular reference K = −log ρ is obtained from the four calibrated settings at the reference p on this fixed, faithful support.

A separate G2 comparison uses the bound pair that generates the two-gateway energy flux. The original source-mixture reference gives the wrong sign for central energy on a held pulse. A prepared cap reads independently metered weighted seam energy on that same pair. [Same-source energy comparison](PAIR.md)

| Source mixture p | Exact first-law slope tr(K dρ/dp) | Entropy slope at step 1/1024 | Expected finite slope at step 1/16 | Sampled slope at step 1/16 ± one SE |
| --- | ---: | ---: | ---: | ---: |
| 3/4 | −0.693669 | −0.693670 | −0.698376 | −0.689896 ± 0.007515 |
| 1/2, symmetric launch | 0 | 0 | 0 | −0.003230 ± 0.004832 |

Each sampled row uses 524,288 simulated preparations per setting and perturbation sign. The uncertainties use independent blocks of terminal readings. The expected finite slope allows for curvature across the sampled step. The finite first law is a matrix identity for faithful states. Here the additional evidence is its reconstruction from physical patch readouts, including the coherences that direct occupancy readings miss. [Executable instrument](__main__.py) · [Regression tests](test_station.py)

## Which area readings fail this test?

The two sharp internal launch alternatives acquire equal populations on A and B under the balanced splitter. Tracing over B therefore leaves the same boundary state before the crossing. The identical A–C unitary preserves that equality including coherence between its two endpoints. All four other boundary ports are empty and all three seam memories are blank. The trace distance of the complete final cut states for p = 0 and p = 1 is 1.48 × 10⁻¹⁵ in the calculation and zero by this argument. With the splitter bypassed the distance is 0.8. Thus the number of boundary links, their density matrix, every outcome distribution of an instrument on those registers, and the corresponding entropy all give a constant across the original launch family.

With the area defined by any of those terminal readings, the residual dS/dp + η dA/dp at p = 3/4 equals −0.693669 nats per unit p for every η. A boundary history, a larger collar including the interior port, or geometry that responds to the matter state can carry additional information. At p = 1/2 these same terminal readings are compatible with stationarity along this particular tangent. No physical vacuum is identified by this experiment.

G4 additionally requires the physical edge-centre area normalization and stationarity at **fixed proper volume** on the same geometric refinement tower used for the stress measurements. Here the hardware region is fixed, the tested centre is regional double number, and the gate windows and source mixtures are supplied. The conditional finite entropy and its derivative are measured. A physical area, proper volume and their response under the perturbation remain open. [Einstein realization premises and small-ball equation](../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex) · [Candidate specification](../../candidate_laws/interval_echo/SPEC.md)

From the RER repository root, replay the archived draws or run the exact controls with

```sh
python -m model.experiments.interval_echo_einstein_station
python -m pytest -q model/experiments/interval_echo_einstein_station
```

## Can a load probe tell how much energy the collar stores?

Prepare a single nonblank memory in the antisymmetric combination of its two coin states. Compare its zero-contrast mode with the even and odd combinations of the signed double contrast. A local interference measurement against the zero mode reads energies 32 and 33 in exchange units. With the agreement term removed the same apparatus reads 0 and 1. Preparation of the coin phase, coherent source interrogation and probe launch are supplied apparatus.

On these three memory states the complete seam operator is precisely the measured energy times the identity on all nine load inputs. To see why, write the repair isometry as $W$. For $d=0,2,-2$ its two coin branches have the same endpoint word, so their difference is orthogonal to the image of $W$ for every endpoint input. The echo acts as the identity on this subspace. Across the five occupied-memory energies $0,8,9,32,33$, the dimensions of the coin subspaces with scalar action for **every** endpoint word are $2,0,0,1,1$. Zero energy has both coin components because the written-image projection commutes with the record generator. These statements follow by integer rank calculations on the native write and an independent check of the full 99-state columns.

Other seams and degree references act trivially on the private coordinate of the marked seam. Changing its energy therefore multiplies the load evolution by a common phase. Every free load measurement, including a sequence of measurements on different patches that leaves this source intact, has identical statistics for the three preparations. On the 20-, 80- and 320-cell spheres, the observer across the first external contact compared its twelve port outcomes and the quiet outcome. A **blank** contact in the same position gave a total variation distance of 0.213–0.215 from the marked source at exposure 1.5. Thus the load probe senses the occupied contact while its stored energy can vary by 33 without changing the reading. A symmetric coin source can respond to two simultaneous loads on its marked endpoints. The measurement distinguishes a local material interaction from the common response of a calibrated geometry required by G6. [Operator proof and free-sphere rehearsal](backreaction.py) · [Patch probabilities and source energies](backreaction.json) · [Two-load comparison](../material_inertia/README.md)
