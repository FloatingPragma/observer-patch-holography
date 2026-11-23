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

## A mobile neutral excitation in Interval Echo

Fill every port of the spherical federation with one load. A local source on one seam prepares a signed contrast-two memory with symmetric coin and reverses the native echo. This emits one excess of load and one vacancy, leaves the collar blank and preserves the total load. The same operation on a zero-contrast reference restores the filled background. A phase coupled to the total load difference gives zero charge to the joint excitation and opposite responses to its separated controls. The phase source measures the native load number. Identifying it with an electromagnetic field requires an independent experiment. [Source and charge controls](neutral_defect.py) · [Finite records](dark_collar.json)

Let $x$ be the excess port and $y$ the vacancy. In this ordered sector the blank-memory swaps preserve both defects. Relative to the uniform background the free law is $12\kappa I+K+D$, where $K$ is half the graph Laplacian that exchanges neighboring loads, including $x$ with $y$. The diagonal is $D=2\kappa[d_y-d_x+\mathbf1_{x\sim y}]$, with $d_v$ the actual port degree. This follows by expanding the native contrast squares and the missing-contact references on the filled background. At $\kappa=8$, $0\le K\le12I$. The adjacent states with a degree-five excess and a degree-six vacancy have $D=32$. Everywhere else $D\le16$. Min–max therefore isolates an upper spectral band in $[32,44]$ above a complementary spectrum bounded by 28 in shifted units. Its ranks on the 20-, 80- and 320-cell screens are 284, 1,152 and 4,592. This is binding on a high-energy branch of the supplied medium.

For a sharp launch on one of these adjacent addresses, the shifted mean energy is 37 and its standard deviation is $\sqrt{10}/2$. Write $P$ for the projector on the high addresses and $Q=I-P$. On each of the three screens a high address has at most eight exchanges out and a low address has at most two exchanges in. The row and column norm bounds give $\lVert QHP\rVert\le2$. Conservation of the initial mean and variance, together with $QHQ\le28Q$, gives $9q\le\sqrt{10}/2+2\sqrt{1-q^2}$ for $q=\lVert Q\psi(t)\rVert$. Thus **at least 85.47% remains on these adjacent addresses at every free exposure**, for these sharp initial states. A separate two-quadrature comparison with the undisturbed filled state reads the unshifted launch energy from its short-time phase slope. The excitation mean is 133 exchange units above that medium. At exposures 1/4,096 and 1/2,048 the calculated slopes are 132.97662 and 132.90650 with finite-step bias bounded by 0.02726 and 0.10904. With 67,108,864 simulated trials per quadrature the readouts were 131.58093 ± 0.49974 and 132.63303 ± 0.24947 (one standard error). The first sampled estimate lies 2.79 standard errors below its outcome-law prediction.

On the 20-cell screen, a regular excess at port 0 and a gateway vacancy at port 1 were released with every seam active. At exposure 32 the exact probability of adjacency was 0.99440, versus 0.13872 with the exchange acting alone. The probability of an adjacent upper-address pair **with both constituents at different ports from launch** rose from 0.00465 at exposure 8 to 0.01641 at 16 and 0.06834 at 32. The larger reorientation probability of 0.81699 leaves the original vacancy fixed. A separate spectral projection assigns 99.516% of the launch weight to the isolated band and finds both-address transport there too. The raw joint detection at an upper address across the specified intercell seam is much rarer, with probability $5.74\times10^{-8}$ at 32. The projected preparation is a spectral diagnostic. The archived finite-trial local reads use sharp starts and include quiet outcomes. [Ordered-sector propagation and projection bound](neutral_defect.py) · [All appointments, controls and local counts](dark_collar.json)

The conserved signed number, positive excitation energy and transport make this a candidate for further dark-sector tests. A physical electrical charge meter, inertial dispersion for the same moving excitation, response of geometry and cosmic abundance remain open. The public OpenWave static row asks for a stable neutral massive state. [Public row definition](https://github.com/openwave-labs/openwave/blob/753232cac4976ab2ca208eb99d73fe7ac88ef333/MODELS.md#L121)

The lowest 24 lines of the same native bound band admit a free multi-cell flight on the 20-cell sphere. After one million exchange times, an adjacent regular excess and gateway vacancy are detected at least two cells from their filtered launch with probabilities of 6.8287%, 12.6035% and 18.8158%. All sharp-launch baselines lie below 1.6 × 10⁻¹² for the three admitted preparations after filtering. The cluster concentrates on twelve cells with the same local gateway motif. The 80-cell refinement has a connected 56-cell region of this type for a future dispersion measurement. [Energy lines, filter admissions and terminal counts](pair_spectrum/README.md)

## Reading a stationary collar by contact

At the same coupling an evacuated written seam has positive signed-contrast-two even and odd energy gaps of 32 and 33 in exchange units. An empty-port phase comparison reads both. The free Hamiltonian preserves the location of the occupied memory. Superpositions over different addresses of the same empty-port eigenmode all have the same energy. This gives a flat band for the bare collar, distinct from the mobile particle–vacancy branch above.

For a symmetric coin and binary loads the local generator on the contrast-two sign mode of energy $E$ reduces exactly to $E(1-n_an_b)$. The written zero mode has energy zero. A single load therefore sees the same motion in both backgrounds at every time. Two loads can interrogate the stored energy when they meet at its endpoints. At native exposure 2, a pair released on (0,1) was found there with probability 0.94057 on a written zero collar and 0.99321 on the two-unit even collar. On the held seam (84,89) the probabilities were 0.92937 and 0.99394. All other seams run freely and the finite patch retains other-home-seam and quiet results.

An active observer can read the stored bit perfectly in the ideal gate model while leaving its source intact. Park a load on one endpoint, split a second between the other endpoint and a spare port with a native half exchange, expose the marked seam for $\pi/32$ and recombine with the inverse exchange. A written zero collar returns the probe, an even two-unit collar gives the orthogonal exit and the independently read odd gap predicts a 0.002408 return. Removing the parked load gives a return for either source. An antisymmetric-coin control also returns the pair. These are local conditional phase measurements on prepared rails.

For a finite recovery test, a fair retained bit $Z$ selects the zero or two-unit source. A separately initialized native dimer flips that choice with probability $1/4$, giving the actual collar word $C$. The matter patch writes its terminal reading $B$. On these commuting records, the conditional recovery is $P_{ZB}P_{C|B}$. It minimizes relative entropy over classical reconstructions of $C$ from $B$. Define the recovery defect $\eta$ as the total variation from that reconstruction and the modular remainder as its log probability minus the actual log probability on their common support. The actual-to-recovered remainder variation has magnitude at most $2\log(3)\eta$, because the conditional probability of either retained bit is between $1/4$ and $3/4$. A blind reader gives $\eta=1/4$ and remainder variation 0.274653. The indirect two-load reading makes $C$ recoverable and both quantities vanish. The material energy marginal is preserved by every one of these recoveries. This realizes a finite information-recovery bound and its switch-off on the declared diagonal records. Assigning the remainder to gravitating stress requires the dark-sector paper's additional source identification. [Contact and recovery replay](dark_collar.py) · [Joint records and recovered distributions](dark_collar.json)

### A signed modular charge of the retained history

One cell uses seam (0,1) for the collar, (1,5) for the contact interferometer, (2,3) for its retained fair bit $Z$, and (4,9) for the independent flip. A sharp read of a one-load dimer after exposure $2\arcsin\sqrt p$ produces the flip with probability $p$. The controller uses the two bits to prepare the actual collar $C$, then parks the apparatus. The terminal probe leaves $C$ intact, so the patch can read its contrast after recording the probe outcome $B$. The rails, coherent source, sharp reads and addressing windows are supplied.

At fixed reference flip probability $p_0=1/4$, the measured operator is $E=-\log P_0+\log P_0^{\rm rec}$ on the common support of the three records. At contact exposure $\pi/64$ the bright return probabilities are $1$ and $1/2$ at the two collar words. Moving to $p=1/8$ or $3/8$ changes $\langle E\rangle$ by $-0.092480$ or $+0.092480$ nats. A fresh 131,072-trial reference reconstructed the operator. Independent 131,072-trial histories at each flip probability gave changes of $-0.09319\pm0.00156$ and $+0.09049\pm0.00174$ nats. These first-order standard errors include both the reference and the history samples. Selecting the $C=0$ or $C=1$ sector gives the respective exact magnitudes $0.095268$ and $0.089693$. Each sector carries half the probability in all three preparations. At contact exposure $\pi/32$ the outcome identifies $C$, giving $E=0$ on the full occupied support. [Native gates, signed comparisons and receipt](dark_collar.json)

The entire distribution of $(C,B)$ stays fixed as $p$ changes. Hence any fixed measurement of these two records gives the same distribution, including the stored energy reading with mean 16. On the partial reader the smallest possible operator-norm error in representing $E$ by a function of $(C,B)$ is 0.717542 nats. The change is accessible to a patch that reads its joint history. The paper's stress identification also requires this remainder to enter a physical boost split, a proper radius and a common response of the geometry. This gated experiment supplies no measurement of those quantities. [Contact experiment](dark_collar.py) · [Dark-sector paper](../../../cosmology/oph_dark_matter_paper.pdf)

Replay the experiments from the repository root with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.experiments.material_inertia.dark_collar`.
