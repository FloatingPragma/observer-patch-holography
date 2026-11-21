# What the first pair reading changes

One observer on the 20-cell sphere can find an adjacent pair after the two loads start separated and evolve freely. Its probability of finding a pair anywhere in its own cell at duration five is 0.499425% under Interval Echo’s free κ = 8 interaction. Inserting a sharp reading of one nominated seam at duration one raises the late probability to 0.512724%, counting both outcomes of the inserted reading. The choice to read changes the observed population. [Frozen protocol](PROTOCOL.md) · [Complete patch probabilities](result.json)

| First service at duration one | Late pair inside the observer cell, κ = 8 | Exchange only, κ = 0 |
| --- | ---: | ---: |
| Bypass the reader | 0.00499425 | 0.23581465 |
| Choose read or bypass with a fair bit | 0.00506074 | 0.23591358 |
| Sharp read, keeping every trial | 0.00512724 | 0.23601251 |
| Condition on a sharp positive | 0.99322582 | 0.33867101 |

The first three rows use the same denominator of all enrolled starts. The fourth row conditions on the loading of one particular seam at the first reading. That seam gave a positive with probability 0.00018914 under κ = 8 and 0.01106690 under exchange. The freely evolving exchange-only control gives the larger terminal adjacency yield at this wait. After a sharp positive, the contrast interaction gives the larger probability of a repeated pair observation. The unconditional rows and the postselected row answer different experimental questions. [Prior two-read observation](../../evaluation/observer_encounter.json) · [Controls and quiet outcomes](result.json)

A separate experiment retains the native coherent repair on the same seam and tests the terminal yield away from its endpoints. [Retained-write experiment](MEMORY.md)

## What the patch can read

A patch of one echosahedron occupies cell 0. The 36 preparations put one load on each of two distinct, nonadjacent ports in that cell with equal weight. All 28,680 basis states of the 20-cell two-load sector enter the free propagation. The first binary reading tests seam (0, 1). The final local read resolves the 30 internal seams and a complementary quiet outcome. Quiet includes two separated loads and a pair beyond this observer’s aperture. Preparations, measurement windows, service decisions and waiting durations are supplied. The clock is idle. Each quoted probability is calculated from the complete finite record distribution. [Native sector](../../candidate_laws/contrast_memory_spherical_v1/sectors.py) · [Local records and compiler](../../candidate_laws/contrast_memory_spherical_v1/presence.py)

Let E denote the projector onto occupation of the nominated seam and Q = I − E. A sharp reading with both outcomes retained sends the input to EρE + QρQ. Deleting the record also loses the relative coherence between E and Q. Applying the local relative sign Q − E to half an ensemble and identity to the other half gives the same measured channel. With a sign on a quarter of the ensemble the channel matches the fair read-or-bypass service. Every final outcome in both interactions agrees to within 10⁻¹⁰. The sign comparison checks the channel identity mathematically. The random read-or-bypass service requires an ordinary sharp read and a fair service choice. [Readout comparison](__main__.py)

## What the free binding theorem explains

The dressed bound band is the isolated low-energy subspace of the contrast Hamiltonian. Its projector commutes with the free Hamiltonian. Its mean population therefore stays constant between interventions, even as the probability of adjacency on a particular seam changes. The simulator can evaluate the band weight. The patch detector directly measures local adjacency. [Binding theorem](../../candidate_laws/contrast_memory_spherical_v1/PROOF.md) · [Band diagnostic](result.json)

| κ = 8 ensemble at duration one | Dressed bound-band weight |
| --- | ---: |
| Separated input, then bypass | 0.00291210 |
| Sharp positive | 0.99556927 |
| Sharp quiet | 0.00290052 |
| Both sharp outcomes, weighted by their probabilities | 0.00308827 |

The all-outcome increase is 0.00017618 per enrolled start. A sharp positive prepares a state with 99.5569% weight in the bound band. The initial separated preparation also has a small bound component, which propagates freely. An unconditional disturbance and conditioning on a rare reading both enter the two-read result. An encounter-created durable record that predates the nominated read is an additional condition in OBS-04. [Observer criterion](../../evaluation/EVIDENCE.md#obs-04)

For the shifted κ = 8 Hamiltonian the binding theorem places the lower band at energies at most −10 and the upper band at energies at least 0. Positivity and the maximum hopping row sum put the full spectrum in [−16, 12]. With X = (H + 5I)/17, the spectrum obeys 5/17 ≤ |X| ≤ 1 and the projector is (I − sign X)/2. For a = 5/17, q = (1 − a)/(1 + a) and T = (1 + a² − 2X²)/(1 − a²), the Legendre generating function gives

$$
(X^2)^{-1/2}=\frac{2}{1+a}\sum_{k=0}^{\infty} q^k P_k(T).
$$

The spectrum of T is in [−1, 1], so truncation at k = 44 bounds the absolute operator error of the projector by q⁴⁵/(2a), below 2.43 × 10⁻¹² in exact arithmetic. A dense one-cell control independently diagonalizes the Hamiltonian and directly evolves the split density matrices. This checks both the numerical projection and the branch weights. [Dense check](test_encounter.py)

## Reproduce

From the RER repository root with the pinned model dependencies, run

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m model.experiments.interval_echo_encounter
python -m pytest -q model/experiments/interval_echo_encounter
```

The replay reconstructs the state amplitudes and compares all final seam probabilities, quiet outcomes and band diagnostics with the archived values. The protocol digest is checked before the run. [Implementation](__main__.py) · [Receipt](result.json)

## What a persistent pair identifies

The M10.3 release prepares two unit loads on seam (156, 157) of cell 13 on the 20-cell sphere. At a supplied exposure the observer reads one of its 30 internal seams or quiet. Keeping the launch and detector fixed gives the following probabilities of a home-seam reading. The κ = −8 arm reverses the sign of the contrast. The last column changes the missing-seam reference at κ = 8. [M10.3 source and trial records](https://github.com/FloatingPragma/observer-patch-holography/tree/4fe10f9ac3e81423155fb3ac7baaf0fdb4c0b11c/openwave/interval_echo/row_04) · [Frozen release comparison](confinement_protocol.json)

| Exposure | Native κ = 8 | Exchange alone | Reversed κ = −8 | Unit reference |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.990983 | 0.450002 | 0.988165 | 0.660699 |
| 7 | 0.989653 | 0.275769 | 0.989071 | 0.657885 |
| 37 | 0.980272 | 0.202906 | 0.983888 | 0.569245 |
| 74 | 0.983949 | 0.190458 | 0.975541 | 0.566346 |
| 148 | 0.982976 | 0.169847 | 0.979925 | 0.578064 |

Each exposure uses a fresh preparation. The local reader returns quiet for separated loads and for adjacent loads beyond its aperture. A complete terminal position census resolves those alternatives. For example, at 148 the native probability on an internal seam of another cell is 0.005346. Only the 30 home addresses and quiet belong to the original observer. The extended exposures are terminal distributions on the finite carrier. [Complete censuses](transport.json) · [Changed-reference censuses](centered_pair.json)

### Exact confinement and separated excitations

In the binary two-load sector the free Hamiltonian after its uniform shift is the exchange generator minus 2κ times the projector onto adjacent loads. Let P select either all 66 two-load positions in a given cell or only its 30 seams. Set A = 2PHP and B = 2(I−P)HP. A state initially in P stays there throughout an interval precisely when BAᵏ annihilates it for k from zero to dim(P)−1. This follows by differentiating the finite matrix exponential, then applying Cayley–Hamilton to A. Conversely, the common kernel is invariant under A and has zero outward amplitude. The calculation builds B from every allowed exit. Exact integer row reduction and an independent rank bound modulo 1009 give the following dimensions. [Exact calculation](confinement.py) · [Independent checks](test_confinement.py)

| Cell gateways up to icosahedral symmetry | Cells on the 20 / 80 / 320 boards | Confined single-load modes | Confined binary pairs in one cell |
| --- | --- | ---: | ---: |
| 0, 1, 2 | 0 / 24 / 104 | 5 | 0 |
| 0, 1, 4 | 8 / 0 / 0 | 2 | 0 |
| 0, 2, 8 | 12 / 56 / 216 | 2 | 0 |

The pair result holds in both apertures at κ = 8, 0 and −8. For the narrower aperture the contrast is scalar on P and diagonal across its boundary, so its common kernel is the same for every real κ. For the published sharp release the leading quiet probability is 5t²/4 for the seam read and t²/4 for leaving the cell. These coefficients follow from the squared norm of the outward column of H. A strictly positive escape probability is consistent with high retention. The result gives neither a disappearance time nor a decay width.

The one-position compact modes have twice-shifted kinetic energies 6 (multiplicity two) on each two-mode cell. A five-mode cell has 6 (multiplicity three) and 5 ± √5. A doubled load has the same position generator up to a scalar. Each compact mode vanishes at the gateways, because the external connections form a matching. Two such modes supported in distinct cells therefore have zero contact and evolve as a tensor product, even in neighboring cells. Their products span 760, 26,500 and 451,420 independent states on the three boards. These are lower bounds on the federation's confined two-load spaces with the loads in distinct cells. They require coherent modes; the sharp one-seam release has zero weight in the same-cell compact space. The energies use the model's exchange unit. [Spectra and multiplicities](confinement.json)

### Measuring the interaction sign

Sharp local persistence alone accepts both signs at the 80% population screen. The observer can distinguish them by splitting a pair coherently. Take x, y, z, w = 156, 157, 164, 159. Starting on xy, apply isolated π/2 contacts on yz and zw. On the three occupied paths xy, xz, xw the amplitudes become (1−i)/2, 1/2 and i/2. After the free exposure apply one more π/2 contact on yz and read the home seams. The binary columns of the native blank-collar operator give the same preparation and analyzer for the three registered strengths. [Sign protocol](sign_protocol.json) · [Native pulse and read](sign_read.py)

In those three input coordinates the exchange block has diagonal (4.5, 5, 5.5) and nearest off-diagonals −1/2. Only xy receives the contrast shift −2κ. The analyzer mixes the first two coordinates, so this block suffices for the initial slope even though the evolution uses all 28,680 positions. If F denotes the final xy effect, its slope is i⟨[H,F]⟩ = κ/2. The same incoherent mixture of the three starts has slope −1/8 at all three strengths. Thus the signed adjacent-pair energy contribution in exchange units can be read as minus four times the coherent slope for these settings.

| Probability of xy at exposure 1/16 | Native | Exchange alone | Reversed |
| --- | ---: | ---: | ---: |
| Coherent | 0.440067 | 0.123661 | 0.024072 |
| Dephased | 0.365418 | 0.364333 | 0.365925 |

The local phase control and terminal address distinguish the signs in a single specified ensemble comparison. Dephasing reduces the difference from 0.415995 to 0.000507 in magnitude. The full distributions include every other seam and quiet. [Terminal records](sign_read.json)

### A particle and a vacancy on the same background

Prepare one load at every port. Adding one load at a port gives an excess; removing one gives a vacancy. The free law conserves each one-defect subspace. The missing-seam term at port v of degree dᵥ is κ(6−dᵥ)(nᵥ−r)², with native reference r = 0. Relative to the uniform background, a defect of sign s = ±1 costs κdᵥ on the incident seams and κ(6−dᵥ)[1 + 2s(1−r)] on the references. Its position Hamiltonian is therefore

$$H_s = L/2 + 6κ I + 2sκ(1-r)\,\mathrm{diag}(6-d_v).$$

Both signs share the exchange operator. The three gateways of a cell have degree six and the other nine ports have degree five. This gives opposite detunings across a mixed-degree seam. A local observer prepares the defect on a gateway, splits it across two successive seams, and reads its address after a final half contact. The native excess and vacancy probabilities at the nominated gateway at 1/16 are 0.407586 and 0.029887 on each board. Their initial slopes differ by exactly 8. The detector resolves all twelve home ports and quiet. This refutes the simple identification of these two defects as conjugates with identical free position records under occupation reflection. Other encodings, including collective excitations, require their own records. [Common-background protocol](conjugation_protocol.json) · [All three boards](conjugation.json)

Under the sitewise reflection n ↦ 2−n, exchange and every seam contrast square are invariant. The change of the full reference energy is 4κ(1−r)Σᵥ(6−dᵥ)(1−nᵥ). Within this uniform-reference family, r = 1 is the unique value restoring this symmetry of the free law at nonzero κ on the federations. It gives identical excess and vacancy records and a probability of 0.101933 for the same gateway read. It also gives less than 80% home-pair retention at every tested sharp-release exposure in the table. The unit-reference control changes the law. A symmetry of its free generator establishes neither an electromagnetic charge nor conjugation of the complete clock and service program. Particle masses, spin, statistics and an identification with electrons or positrons require independent operational measurements.

### Replay and numerical accuracy

From the repository root, set the BLAS thread count to one. The short checks are

```sh
python -B -O -m model.experiments.interval_echo_encounter.confinement
python -B -O -m model.experiments.interval_echo_encounter.sign_read
python -B -O -m model.experiments.interval_echo_encounter.conjugation
python -m pytest -q model/experiments/interval_echo_encounter
```

The extended checks replay the published release, the changed reference and the independent error enclosure.

```sh
python -B -O -m model.experiments.interval_echo_encounter.transport
python -B -O -m model.experiments.interval_echo_encounter.conjugation --pair-regression
python -B -O -m model.experiments.interval_echo_encounter.check_transport
```

For the native, exchange-only and reversed releases, a degree-24 Taylor polynomial advances the centered generator by 1/8. The absolute row-sum bound controls the truncation and each sparse matrix multiplication. Accumulating the one-step error through the longest exposure and summing over disjoint terminal bins gives an L1 error below 1.2 × 10⁻¹⁰ for every archived distribution in those three arms. The centered-reference column uses a separate sparse exponential. All runs use prepared backgrounds, supplied isolation and exposure, and a terminal read. The simulation compares observer record probabilities; the exposure variable has no laboratory time calibration. [Numerical enclosure](check_transport.py)
