# A conserved quantum field current and a source-sensitive transport curvature

Two independent simulated instruments follow local energy moving out of one triangle of the written quantum field. Changing the phase of this source preserves every history of its ring energies and currents and changes how independent material probes spread. The charged-probe curvature on this cell therefore requires information beyond the measured energetic history. [Primary protocol](PROTOCOL.md) · [Full-cohort and continuity protocol](REPLAY_PROTOCOL.md)

## Local energy leaving a face

The preparation acts on head (4,11) of the 458-word neutral component. It gives a coherent sign-dependent phase, then lets the field evolve for half an exposure unit under the existing positive written ring Hamiltonian. Each face has a positive local energy projection h_f and their sum is H. The transfer current from a to b is i[h_a,h_b]. For independent copies of the preparation, field energy can be read as a binary ring outcome and the current as an ideal spectral outcome. This latter reader is supplied by the experiment. The addressed-ring interferometer below replaces this spectral measurement with output bits.

At face (2,4,11), 149,351,740 of 2³² copies clicked at dwell 15/32 and 145,871,901 of 2³² clicked at dwell 17/32. An independent ensemble observed the spectrum of i[H,h_f] with a uniformly randomized appointment between those dwells. The integrated continuity equation predicts agreement between the current's window mean and the endpoint difference divided by the window length.

| Instrument | Measured transfer per unit exposure | Simultaneous 99% radius |
| --- | ---: | ---: |
| Endpoint ring-energy clicks | −0.0129634 | 0.0008733 |
| Window-averaged spectral current | −0.0128959 | 0.0000273 |

The signed rates resolve energy leaving this face under the finite Hamiltonian. The independent current read also resolves one nonzero bond contribution across faces (2,4,11) and (2,10,11). In the two opposite-phase preparations its mean readings are +0.0079014 and +0.0079021. [Energy and probe counts](results.json) · [Windowed current and complete cohort counts](replay.json)

## Counting current with an addressed ring interferometer

An independently simulated reader uses only exposure of existing rings and a supplied two-arm detector. On one arm, an addressed π pulse of ring b swaps each enabled pair of legal words. The other arm stays idle. Interference at the dark port consequently projects onto the antisymmetric ring amplitude. To read the current from ring a into b, first expose a to a π/2 pulse or a 3π/2 pulse on independent copies. Half the difference between the two dark-port rates equals the current i[h_a,h_b]. This follows by writing each unitary as I + (exp(−iα)−1)h_a. Only the three faces sharing a written head with b can contribute.

For the previously enrolled (2,4,11) region and time window, fresh output bits give net influx **−0.0129025 ± 0.0000893**. Its two resolved contributions come from the (2,3,4) and (2,10,11) neighbors. Independent endpoint clicks give an energy slope **−0.0128329 ± 0.000952**. The disjoint-ring control reads −0.00000465 ± 0.0000298. These conservative simultaneous 99% radii cover all ten settings. This implementation uses supplied coherent addressing and interference on separate copies after the autonomous field flight. [Interferometer protocol](INTERFERENCE_PROTOCOL.md) · [Legal-word gates and detector](interference.py) · [Detector-bit counts](interference.json)

## Which source data agree exactly?

Let g be the uniform ground state. Every face projector annihilates g. For a signed local rotation about the written-head orientation Z_e, the two freely evolved states are

$$
\psi_\pm = \cos\theta\, g \;\pm\; i\sin\theta\, v,\qquad
v=e^{-itH}Z_e g,\qquad \langle g,v\rangle=\langle g,Z_e g\rangle\in\mathbb R .
$$

Every ring operator and every free-flight unitary preserves the line spanned by g and its orthogonal complement. Each history made from their algebra, including conditional choices among these measurements, consequently has a positive history effect E with Eg = λg and real λ. The sign-dependent contribution to its probability is proportional to the imaginary part of ⟨g,Ev⟩, which vanishes. The same argument covers spectral projectors of i[h_a,h_b] and i[H,h_b] and arbitrary appointment spacings. A written-head orientation meter can distinguish the coherent phases. Physically the indistinguishability applies to the ring-energy/current history algebra. Pressure and shear have yet to be derived in an observer frame for this material.

## How independently launched material spreads

The calibrated price of transferring a unit between ports i and j is the total-variation separation of their charge-even vacuum endpoint profiles at dwell 1/8. At dwell 3/16, a pair of nearby launch ports gives two distributions of terminal records. The Wasserstein distance is the least cost of matching them using the vacuum prices. Discrete coarse Ricci curvature is one minus their matching cost divided by the calibrated launch separation. The experiment uses separate vacancy and excess charges and their equal mixture. The two calibrations, launches, and outcomes are independent draws in the simulated acquisition. [Ollivier’s definition](https://arxiv.org/abs/math/0701886)

For the matched opposite phases of head (4,11), the launch pair (10,11) was selected from the held heads using calculated laws and then measured on fresh simulated releases. It gave the following increase in discrete curvature when the source phase changed from negative to positive.

| Probe | Fresh change in curvature |
| --- | ---: |
| Vacancy | +0.092388 |
| Excess | +0.011183 |
| Equal-prior charge kernel | +0.042057 |

A common conservative simultaneous 99% radius is 0.000579, including cost-calibration uncertainty. Pooling the two signs means combining their endpoint distributions before taking the transport distance. A finite test using the primal transport problem reproduces the sampled pooled contrast independently of the dual solver used in acquisition. [Acquisition](measure.py) · [Regression checks](test_measure.py)

On incoherent mixtures of the source phases, one fitted coefficient predicts the change in mean curvature on the three edges of each face from the energy observed on that face. Twenty enrolled heads train the coefficient and ten are held out, with all twenty faces evaluated on each held head. Fresh acquisitions yield coefficient −0.116785. The relative held errors are 92.7% for vacancies, 88.7% for excess charges and 81.4% for the pooled kernel. The predeclared limit is 5% for each probe. This fit addresses only the specified pointwise scalar hypothesis. The exactly matched histories impose a stronger restriction on using any function of these energy and current records as the sole source for this specific discrete curvature. [Independent acquisition](replay.py)

The simulation operates on a single existing finite written-field extension, with no change to its Hamiltonian. The measured transport curvature depends on the material contact dynamics, source coherence, probe charge, appointments and vacuum distance calibration. Identifying spacetime curvature requires independent proper clocks and distances that agree across probes and scales. Testing gravitational sourcing then requires their curvature and an independently normalized field stress tensor on the same coarse geometry.

## Reproduction

Run from the repository root with NumPy and SciPy installed. Choose new output paths. The independent replay reads the archived primary result adjacent to the code.

```sh
python -m model.experiments.interval_echo_einstein_station.field_stress_curvature.measure --output /tmp/field-stress-primary.json
python -m model.experiments.interval_echo_einstein_station.field_stress_curvature.replay --output /tmp/field-stress-replay.json
python -m model.experiments.interval_echo_einstein_station.field_stress_curvature.interference --output /tmp/field-current-clicks.json
python -m pytest -q model/experiments/interval_echo_einstein_station/field_stress_curvature
```
