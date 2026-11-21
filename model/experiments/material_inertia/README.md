# Material inertia: binding, motion and rest energy

This experiment runs the reference seam law,
[A-contrast-memory-spherical-v1](../../candidate_laws/contrast_memory_spherical_v1/README.md), on rings and on the
federation, and asks what an observer who measures one binding energy can predict about the motion of clusters.
The law, its coupling, the apparatus and the ring are supplied; the observer reads only its own counters.

## What binding predicts about motion

On a 19-port ring an observer measures the energy released by binding two loads. From that one reading
and the droplet dispersion of the candidate law it predicts where clusters of two through five loads will
be detected, including every failure to find an intact cluster. Every complete position distribution
is within 0.00214 of its forecast in total variation. A competing law gives exactly the same binding
record. Its largest error across waits for each cluster size ranges from 0.82 to 0.94. [Derivation](../../candidate_laws/contrast_memory_spherical_v1/PROOF.md#9-binding-energy-and-the-motion-of-a-cluster) ·
[Detector records](records/flight/)

The same tests distinguish motion from rest energy. With the benchmark coupling of eight, clusters
of two, three and four loads have kinetic ratios of about 34, 1,155 and 39,236. Their rest-energy
ratios are about 1.686, 2.354 and 3.020. A common relativistic mass would make the two sets of ratios
equal. A common rescaling of energy or length leaves both sets unchanged. A separate change to the
memory of the law converts a doubled load into two single loads on the federation and permits the
reverse reaction. [Memory experiment](records/memory/)

## Run it

The material experiments use the [pinned model dependencies](../../requirements.txt). Set
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`. The short replay rebuilds the clock and pair-binding
inboxes and decodes the sealed forecast. The other commands regenerate the complete measurements.

```sh
python -B -O -m model.experiments.material_inertia.challenge --stage predict
python -B -O -m model.experiments.material_inertia.challenge --stage test
python -B -O -m model.experiments.material_inertia.flight
python -B -O -m model.experiments.material_inertia
python -B -O -m model.experiments.material_inertia.memory
```

Pass `--spectra-out PATH` to the base material-inertia command to export every binary momentum block's
eigenvalues and launch weights into a fresh directory. Rings have four fixed reference contrasts
per port. Translation-character preparation, coherent path switches, phased recombination and a
phase-stable counter are laboratory services. The baseline spectra use a parked phase reference.
The binding challenge compares two evolving arms of equal load, with separate registers for the cluster
and the uniform singles. Its clock counts local service intervals and compares two quadratures of a
one-load dimer. The longest direct flight uses 2,289,429,780 service ticks. This clock assumes a
uniform tick during the wait. Physical clock normalization and distances require further measurements.

The two outcomes of each interferometer setting use 4,096 evenly spaced counter thresholds. These
are deterministic simulated counts with bounded rounding error. The reader unwraps 30 doubling
intervals and checks a different wait. Its frequency error follows from the quantization bound and a
dominant spectral weight of at least 0.95, verified separately on the source. Broad or weak lines can
give an unresolved result. At a common exposure the reader subtracts whole phase turns exactly before
subtracting the measured residual angles. This resolves small kinetic shifts on large rest phases.
Position exposures use 65,536 thresholds per setting. All translated compact positions and their
complement are read in one ensemble. The choice of law, coupling, apparatus and ring is supplied.
For the sharp position read, a reversible query visits the occupied-run window, copies its conjunction
to a pointer and retraces the visits to clear its private scratch. The two-owner comparison of the law writes
the public flag. The [detector](detector.py) checks the complete partition on every
configuration before the Fourier calculation uses it. Readout assumes a switch that parks the source
for the queries. Across the five count archives the apparatus retains 191 arrays and 128,581,632
counter outcomes, including the competing laws.
