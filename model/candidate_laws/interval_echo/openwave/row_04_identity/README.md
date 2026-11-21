# Does the freely released pair behave as a colour singlet?

The same prepared pair that passes the stability test has a detectable response to the centre of the candidate's guarded SU(3) current. On the original 20-cell sphere, the two loads start on ports 0 and 1 and evolve with all 630 seams active. At native exposure 150, one particular seam in the home cell, (2,3), registered **12,051 times under sham and 1,576 times after a centre turn**, out of 262,144 simulated preparations in each arm. The verdict uses all 30 seam addresses and quiet. [Complete counts](results.json) · [Free-release experiment](../row_04/README.md)

| Free exposure | Centre versus sham, probability TV | Centre versus sham, sampled TV | Turn and inverse versus sham, sampled TV |
| ---: | ---: | ---: | ---: |
| 0 | 0.26802 | 0.26811 | 0.00000 |
| 32 | 0.21045 | 0.21019 | 0.00461 |
| 150 | 0.19636 | 0.20168 | 0.00655 |

TV is total variation distance across the 31 outcomes. The rejection threshold is 0.011996 for the six comparisons together. The unperturbed read matches the archived pair experiment at all three appointments. During each terminal test the apparatus protects the cell and applies the ideal guarded address services in 2,048 symmetric steps. The finite word differs from the intended centre operator by 0.00000104 in operator norm. The program retains the full unnormalized home-cell component. Branches with fewer than two loads in the home cell necessarily read quiet because every pulse conserves the cell's load count. The protected services act only after free flight. [Protocol and threshold](PROTOCOL.md) · [Executable](experiment.py)

### Which electron interpretation this tests

A singlet is invariant under every element of SU(3). A single registered element that changes the outcome distribution suffices to rule out that assignment for this preparation under these currents. The electron interpretation would have to explain why a different action represents physical colour, or identify a different excitation, such as a dressed state, that passes the same test. A spatially restricted current on a bare pair cannot exclude a dressing shared with other cells. Independently identifying the measured SU(3) with the colour interaction requires matter and interaction records. [Guarded pair currents](../../CURRENTS.md) · [Group action and its calibration](../../../../experiments/interval_echo_gauge/README.md)

There is a restriction even before identifying the interactions. The cell's two-load representation consists of a fundamental SU(3) triplet, a separate weak SO(3) triplet, and 60 spectator states. Half the sum of the twelve currents generates the primitive abelian action. Its unit-charge projector coincides with the projector onto the SU(3) triplet. Thus every SU(3) singlet has zero primitive abelian charge, independent of the choice of unit, and the weak triplet supplies no weak doublet. This rules out the usual charged weak-singlet lepton assignment using this abelian generator on the local pair. Physical electric charge can involve a weak generator as well; a response to this abelian pulse alone would not calibrate electromagnetism. [Representation and centre](../../../../experiments/interval_echo_gauge/README.md) · [Operator audit](results.json)

An electron assignment would require a reproducibly released excitation on which patches can identify unit negative electric response relative to a conjugate and a neutral control, measure a spin-half response against a phase reference, and establish colour neutrality with a reader sensitive to coloured controls. The low stable branch must also be distinguishable from the other charged leptons. Row 4 itself tests the narrower property of free persistence. The persistent pair supplies a common preparation on which to perform the identity tests. [Particle identity criteria](../../../../evaluation/particles.json) · [Public row definitions](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md#simplest-test-per-criterion)

From the RER repository root, reproduce the complete receipt with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_identity.experiment
```
