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

| Exposure | Native κ = 8, r = 0 | Exchange alone | Reversed κ = −8, r = 0 | Centered κ = 8, r = 1 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.990983 | 0.450002 | 0.988165 | 0.660699 |
| 7 | 0.989653 | 0.275769 | 0.989071 | 0.657885 |
| 37 | 0.980272 | 0.202906 | 0.983888 | 0.569245 |
| 74 | 0.983949 | 0.190458 | 0.975541 | 0.566346 |
| 148 | 0.982976 | 0.169847 | 0.979925 | 0.578064 |

Each exposure uses a fresh preparation. The local reader returns quiet for separated loads and for adjacent loads beyond its aperture. A complete terminal position census resolves those alternatives. For example, at 148 the native probability on an internal seam of another cell is 0.005346. Only the 30 home addresses and quiet belong to the original observer. The extended exposures are terminal distributions on the finite carrier. [Complete censuses](transport.json) · [Changed-reference censuses](centered_pair.json)

### Can a neighboring observer read the phase of a bound pair?

On the 20-cell federation, cell 6 distinguishes two coherent launches from cell 0. Both launches have the same local load probabilities. After 32 switches of the receiving clock, the minus launch gives nearly twice as many pair readings as the plus launch. A bound on the contribution of the complementary energy band lets the receiver attribute part of this difference to the moving bound excitation. These are finite measurement probabilities on a specified carrier. A physical species, dispersion in physical units and a rest mass require independent measurements. [Frozen experiment](gateway_protocol.json) · [Complete counts and sampling laws](gateway_counts.json)

Two unit loads start apart on ports 0 and 2. After one supplied π exposure the source reads seam (0,1). The probability of a positive is 0.0000563630. A negative terminates the run as a refusal. On a positive, an isolated native contact on (1,5) splits the second load between 1 and 5 while 0 remains occupied. Pulses of π/2 and 3π/2 prepare (e ± i e′)/√2 up to a common phase, with e = (0,1) and e′ = (0,5). Parking synchronizes the releases. An equipped auxiliary seam swaps a double load at each π pulse. At each pulse a probe reads the first endpoint against a blank, copies the yes or quiet outcome, and undoes its write. The receiver's seven-bit register advances when that recorded bit changes. Preparation switches, synchronization, read windows and fresh registers are supplied services. The free pair evolves on all 28,680 binary configurations. A compiled presence read gives one of cell 6's 30 seam addresses or quiet. Cell numbers label apparatus locations on the hardware. [Native pulse, clock and reader](gateway_source.py)

| Flight law at clock word 32 | Plus clicks / admitted | Minus clicks / admitted | Fair phase mixture |
| --- | ---: | ---: | ---: |
| Interval Echo, κ = 8 | 412,963 / 495,806,017 | 777,653 / 495,736,905 | 0.12015% |
| Exchange only, κ = 0 throughout | 0.69988% | 0.84994% | 0.77493% |
| Native enrollment, blocked exchange between adjacent and separated pairs during flight | 0 | 0 | 0 |

Each law, preparation and booking uses 2⁴³ simulated original starts. Seeded binomial draws aggregate independent trials. Each receipt includes all refusals and admitted quiet outcomes. The three menus share an enrollment law within each control. A Bernstein bound covers all 32 original-start bins and all 31 conditional apertures across 81 cohorts with simultaneous coverage of at least 99%. Conditional intervals include the rounding of the sampling law to multiples of 10⁻¹³. At word 32 the minus excess is between 0.07192 and 0.07523 percentage points. A count in any internal seam suffices for this comparison. [Analysis and replay](gateway.py) · [Results at every booking](gateway_result.json)

The two launches also have the same complete energy distribution. Up to a global phase, complex conjugation exchanges them. The free Hamiltonian is real symmetric, so every spectral projector is real and gives equal weight to the two launches, including at degeneracies. For any of the real receiving projectors the minus probability at t equals the plus probability at −t. The signed contrast is therefore odd in the flight parameter. Its detection distinguishes two preparations of equal bound-band weight.

**Which component travels.** In the energy convention shifted by 12κ the free Hamiltonian is H = K − 16P, where P selects adjacent loads. The binding inequalities and the min–max principle locate its two spectral bands in [−16,−10] and [0,12]. Let Π select the lower band, S = I − Π, U the free evolution, F a receiving aperture and ψ a launch. The complement contributes at most η = ‖FS‖ ‖Sψ‖ to the receiving amplitude at any time. Since Fψ = 0 at release, the lower bound √pₗₒ − 2η on ‖F(U − exp(iθ))Πψ‖ holds for every phase θ. At word 32 it gives positive travel bounds of 0.01371 and 0.02446 for the two signs. This time-uniform bound leaves the phase contrast within the band unresolved. [Binding hypotheses and proof](../../candidate_laws/contrast_memory_spherical_v1/PROOF.md)

A separate propagation bounds δ = ‖FUSψ‖ at each booking. The local low-band amplitude lies between √pₗₒ − δ and √pₕᵢ + δ, clipped to [0,1]. At word 32, δ is below 0.000644 for plus and 0.000616 for minus. The two **bound-band receiving norms differ by at least 0.00924** after counting and propagation error. For one individual receiving seam the lower bound is 0.00400. Seven of the nine bookings distinguish the two low-band norms on the whole receiving cell. A second calculation propagates the complete launches alongside their upper-band components, checks all native seam probabilities against the archive and subtracts amplitudes to check the inferred norm intervals. The patch observes ordinary clicks. Attributing their contrast to the lower band uses the declared source Hamiltonian. [Complement protocol](gateway_complement_protocol.json) · [All bounds](gateway_complement.json)

For the projector, X = (H + 5I)/16 has |X| ≥ 5/16 and ‖X‖ ≤ 17/16. The binomial series for sign(X) through degree 160 in I − X² leaves projector error at most (17/32)ρ¹⁶¹/(1 − ρ), ρ = 1 − (5/16)². The implementation adds 10⁻⁸ per unit column for arithmetic. Each normalized sparse multiplication has at most 13 terms per row and ‖|X|‖ ≤ 2. The error introduced by a recurrence step is below 10⁻¹³ in Euclidean norm for unit input. Summing the propagated errors of the contractive steps adds less than 10⁻⁸. Receiver compression uses a 30-column error allowance. For propagation, the Chebyshev expansion on (H + 2I)/14 through order 100 has a one-switch tail below 10⁻²¹. The Bessel series evaluated at 90-digit precision checks the weights. With absolute row sums at most one, each vector recurrence adds less than 2 × 10⁻¹⁴ error for unit input. The second-kind Chebyshev bound |Uⱼ(x)| ≤ j+1 on [−1,1] bounds its growth through order 100. The sum of absolute coefficients is below ten. The propagated amplitude allowance is the initial projection error plus 10⁻⁸ per switch. [Independent propagator](gateway_complement.py) · [Dense checks](test_gateway.py)

**How the interaction crosses the gateway.** Every move within PHP joins two sides of an internal triangle. The bridge (1,82) is isolated in PHP. Write D = (I−P)H(I−P) on separated pairs and C = −(I−P)HP. For E < 0, eliminating those pairs gives PHP − Cᵀ(D−E)⁻¹C. The inverse has a nonnegative path expansion since D is diagonally dominant with nonpositive off-diagonals and D−E is strictly diagonally dominant. A path from (s,u) via (s,b) to the bridge (u,b) gives a coupling of magnitude at least 1/[4(6−E)]. For each internal neighbor z of b, the path through (s,b), (s,z), (u,z) reaches (b,z), giving at least 1/[16(6−E)³]. Here u is the near gateway and s its nongateway neighbor. All these contributions have the same sign. At E = −13 the computed bridge entry is −0.01451027. The calculation retains the full column at three energy probes and checks it against a residual-bounded solve. The energy-dependent operator accounts for intermediate configurations in a stationary eigenproblem. The time experiment evolves the complete Hamiltonian. [Kernel protocol](gateway_kernel_protocol.json) · [Calculation](gateway_kernel.py)

The first phase-sensitive term also follows directly from the time generator. For each of the five internal seams touching b there are two four-exchange paths from e and four five-exchange paths from e′. Each exchange contributes −1/2 to H. The respective first amplitudes are t⁴/192 and it⁵/960. Their interference gives p₊ − p₋ = −t⁹/18432 + O(t¹¹) over the receiving cell. Diagonal waits cannot shorten a path. The enumeration includes every receiving seam on the 20-, 80- and 320-cell boards. This common leading coefficient in the supplied exposure t is a local statement across carriers. The finite-time rates above are for the 20-cell apparatus.

**Reading a message.** The source selects a phase after each successful herald. The receiver gets the live herald, synchronized booking and a sequence of local click or quiet words. It counts clicks across 131,072 admitted trials and reads minus for 155 or more. A frame that reaches 2³² original starts without enough admissions returns an erasure. The simulator samples frame counts in aggregate. The classical decoder can hold attempts in 33 bits and admissions and hits in 18 bits each, in addition to the seven-bit clock. Word 32 and the threshold use the gateway counts as design data. An independent seed generates 128 held-out frames for each sign and control with these choices fixed. The native receiver read 256/256 bits correctly with zero erasures. The exchange and blocked controls each gave 128 wrong bits, and erasing the sender's phase choice gave 119 wrong. The receipt preserves every frame including the number of original starts. With independent launches, binomial Chernoff bounds from the simultaneous gateway intervals give wrong-or-erased probabilities below 0.000295 for plus and 0.001115 for minus. This is a classical bit encoded in the relative phase of a costly, heralded preparation. [Message protocol](gateway_message_protocol.json) · [Frame records](gateway_message.json)

Run the frame replay and small independent checks with the pinned model environment. The two full-sector commands regenerate the detector probabilities, Schur columns and separate band bounds.

```sh
python -B -O -W error -m model.experiments.interval_echo_encounter.gateway_message
python -m pytest -q model/experiments/interval_echo_encounter/test_gateway.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.experiments.interval_echo_encounter.gateway
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.experiments.interval_echo_encounter.gateway_complement
```

### Exact confinement and separated excitations

In the binary two-load sector the free Hamiltonian after its uniform shift is the exchange generator minus 2κ times the projector onto adjacent loads. Let P select either all 66 two-load positions in a given cell or only its 30 seams. Set A = 2PHP and B = 2(I−P)HP. A state initially in P stays there throughout an interval precisely when BAᵏ annihilates it for k from zero to dim(P)−1. This follows by differentiating the finite matrix exponential, then applying Cayley–Hamilton to A. Conversely, the common kernel is invariant under A and has zero outward amplitude. The calculation builds B from every allowed exit. Exact integer row reduction and an independent rank bound modulo 1009 give the following dimensions. [Exact calculation](confinement.py) · [Independent checks](test_confinement.py)

| Cell gateways up to icosahedral symmetry | Cells on the 20 / 80 / 320 boards | Confined single-load modes | Confined binary pairs in one cell |
| --- | --- | ---: | ---: |
| 0, 1, 2 | 0 / 24 / 104 | 5 | 0 |
| 0, 1, 4 | 8 / 0 / 0 | 2 | 0 |
| 0, 2, 8 | 12 / 56 / 216 | 2 | 0 |

The pair result holds in both apertures at κ = 8, 0 and −8. For the narrower aperture the contrast is scalar on P and diagonal across its boundary, so its common kernel is the same for every real κ. For the published sharp release the leading quiet probability is 5t²/4 for the seam read and t²/4 for leaving the cell. These coefficients follow from the squared norm of the outward column of H. A strictly positive escape probability is consistent with high retention. The result gives neither a disappearance time nor a decay width.

For the one-position compact modes, twice the kinetic generator has eigenvalue 6 with multiplicity two on each two-mode cell. A five-mode cell also has 5 ± √5, and the multiplicity of 6 is three. A doubled load has the same position generator up to a scalar. Each compact mode vanishes at the gateways, because the external connections form a matching. Two such modes supported in distinct cells therefore have zero contact and evolve as a tensor product, even in neighboring cells. Their products span 760, 26,500 and 451,420 independent states on the three boards. These are lower bounds on the federation's confined two-load spaces with the loads in distinct cells. They require coherent modes; the sharp one-seam release has zero weight in the same-cell compact space. The energies use the model's exchange unit. [Spectra and multiplicities](confinement.json)

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

Both signs share the exchange operator. The three gateways of a cell have degree six and the other nine ports have degree five. This gives opposite detunings across a mixed-degree seam. A local observer prepares the defect on a gateway, splits it across two successive seams, and reads its address after a final half contact. The native excess and vacancy probabilities at the nominated gateway at 1/16 are 0.407586 and 0.029887 on each board. Their initial slopes differ by exactly 8. The calculation supplies an ideal terminal measurement of the defect position. Its outcomes are the twelve home ports and quiet. The runs compile the preparation and analyzer pulses on the native seam. Implementing the terminal position measurement as a patch service on the filled background is a separate instrument requirement. This refutes the simple identification of these two defects as conjugates with identical free position records under occupation reflection. Other encodings, including collective excitations, require their own records. [Common-background protocol](conjugation_protocol.json) · [All three boards](conjugation.json)

Under the sitewise reflection n ↦ 2−n, exchange and every seam contrast square are invariant. The change of the full reference energy is 4κ(1−r)Σᵥ(6−dᵥ)(1−nᵥ). Within this uniform-reference family, r = 1 is the unique value restoring this symmetry of the free law at nonzero κ on the federations. It gives identical excess and vacancy records and a probability of 0.101933 for the same gateway read. It also gives less than 80% home-pair retention at every tested sharp-release exposure in the table. The unit-reference control changes the law. A symmetry of its free generator establishes neither an electromagnetic charge nor conjugation of the complete clock and service program. Particle masses, spin, statistics and an identification with electrons or positrons require independent operational measurements.

### A symmetric reference that preserves the binary sector

There is a second way to realize occupation reflection. Keep the physical seam law, κ and the blank reference load. Write the missing-seam energy as κf(n). Preserving the native binary energies with the vacuum as zero fixes f(0) = 0 and f(1) = 1. Reflection fixes f(2) = f(0). Thus f(n) = n(2−n) is the unique normalized diagonal profile in this class. This is an experimental variant of the free reference service.

The native echo writes the signed difference d into its collar. On each missing seam, assign unit energy to the written labels with |d| = 1, leave every other label at zero energy, and undo the echo. Both coin branches of every blank input carry the same absolute contrast. The label projector therefore intertwines the native write with the diagonal test for |a−b| = 1 on all nine load inputs. During the middle phase operation the collar holds the written label. Applying the inverse echo returns it to blank for every exposure duration, including on coherent superpositions. A reference rail at b = 0 then implements the required f. The replay checks all 99 components of every blank input column. [Native collar construction](conjugation.py) · [Frozen tests](conjugation_protocol.json)

This gives occupation-reflection symmetry on the entire blank-collar free load space. Every binary configuration has exactly its native diagonal and exchange action, on any of the federations. Free propagation and the registered physical half contacts preserve the binary sector. Consequently **every probability in the original pair release** agrees at every exposure. A separate assembly sums the physical qutrit seam columns and missing-seam terms on all 28,680 binary basis states of the release board. It agrees entry by entry with the native pair generator and reproduces all 34 bins at the five registered exposures. On the common one-load background the excess and vacancy both have generator L/2 + κ diag(2dᵥ−6). Their terminal distributions equal the native vacancy distributions on all three boards, including the 0.029887 gateway reading. [Defect records](conjugation.json) · [Shared binary records](transport.json)

The doubled load above the empty background exposes the cost. Its free generator changes from L/2 + 24κI to L/2 + 4κ diag(dᵥ). The shifts relative to the native law are 0 at gateways and −32 on the other ports at κ = 8. Thus a single diagonal reference of this kind cannot preserve all three properties: the binary free energy, occupation reflection and native doubled-load propagation. The construction supplies the first two and gives a distinct prediction for the third. Equivalence of protocols using occupied collars, doubled clock fuel or interventions that split a double requires a separate calculation.

<a id="how-many-frequencies-can-a-pair-of-port-observers-see"></a>

## How many frequencies can one patch measure?

The observer occupies cell zero. Its two endpoint registers copy the same
predicate of one load at each port of seam (0, 1). A courier brings the
copies together to record agreement. On the 20-cell federation, the two
public readings bound the energy content of the state prepared by a
positive first reading. At simultaneous 99% confidence, **at least 36.72%
of its weight lies outside any two energy levels**. The full generator
supplies a separate exact result on each registered federation.

| Federation | Minimum native energies contributing to a sharp seam return | Minimum distinct positive frequencies in its continuous-exposure return probability |
| --- | ---: | ---: |
| 20 cells | 28 | 27 |
| 80 cells | 63 | 62 |
| 320 cells | 127 | 126 |

**The exact bound.** In the binary two-load sector every off-diagonal
transition exchanges one occupied port with an adjacent vacancy with matrix
element $-1/2$. The contrast term is diagonal. From a source configuration
$e$, a configuration at graph distance $j$ has zero amplitude in $H^i e$
for $i<j$. Its entry in $H^j e$ equals $(-1/2)^j$ times the positive number
of shortest paths. Along a shortest path of length $D$ the successive
shells therefore certify $D+1$ independent vectors. The spectral theorem
identifies the dimension of their full cyclic span with the number of
distinct energies having positive local weight. Any two-token
configuration is reachable on a connected port graph because edge swaps
generate all permutations of port occupancies. To go from $\{a,b\}$ to
$\{u,v\}$ costs at least
$\min(d(a,u)+d(b,v),d(a,v)+d(b,u))$ moves. Breadth-first searches on the
port table give witnesses of 27, 62 and 126 moves from seam $(0,1)$.
The proof uses adjacency to count powers of the operator. Physical distance
plays no part. [Source generator](../../candidate_laws/contrast_memory_spherical_v1/sectors.py)
· [Witness seams and carrier digests](return_lines.json)

After a positive read the state is $e$. Write $p_j=\|\Pi_j e\|^2$ for its
nonzero energy weights. The return probability is
$P(t)=|\sum_j p_j e^{-itE_j}|^2$. Each gap carries positive weight.
Subtracting the least visible energy from each of the others gives
the stated number of distinct positive gaps. The bound holds for every
real diagonal contrast potential with the same exchanges. It concerns
continuous exposures and places no lower bound on the weight of each gap.

**The local experiment.** Each cohort starts with an equal mixture of the
36 nonadjacent binary placements in cell zero. One reference interval later,
the two endpoints and their courier report whether both loads occupy
$(0,1)$. At a booked reading from zero through 16 intervals later they ask
the same question. All four ordered answer pairs are retained. A separate
reference seam carries a sharp double load. Each isolated native gate of
duration $\pi$ swaps the occupied endpoint. Addressing the endpoints
against blanks checks the alternating pointer with the native echo. A
five-bit register suffices for the 17 intervals including the first approach.
Matched gate durations, an equipped reference channel, sharp read services
with the material source held fixed, and the prepared cohorts are supplied.
The receiver uses the booking and both public words. It has no slot
amplitudes, eigenvalues, or calibrated physical time units.

For a real symmetric free generator one propagated column at the measured
seam suffices to calculate every launch-to-seam amplitude by transposition.
For launch $x$, let $a=\langle e|U_1|x\rangle$,
$c=\langle e|U_j|e\rangle$ and $b=\langle e|U_{j+1}|x\rangle$, with one
clock interval per subscript unit. The two joint terminal-yes probabilities
are $|ca|^2$ and $|b-ca|^2$. Subtraction from their respective first
branch probabilities gives both terminal-no outcomes. The code averages
over all 36 inputs and propagates the full 28,680-state sector. Switched
controls park the post-read state or leave only seam $(1,5)$ active. The
exchange-only arm uses $\kappa=0$ at both exposures.
[Protocol](return_protocol.json) · [Public counts](return_counts.npz)

| Simulated seam service | Moment rank from public counts, lower bound | Required positive gap frequencies, lower bound | Lower bound on weight outside any two levels |
| --- | ---: | ---: | ---: |
| Native $\kappa=8$ | 17 | 8 | 0.3672 |
| Exchange only | 17 | 8 | 0.5483 |
| Park after first read | 1 | 0 | 0 |
| Isolate one swap | 2 | 1 | 0 |

Condition on the positive first records within each independent cohort.
For $K_{ab}=P((a-b)\pi)$, $0\le a,b\le16$, each spectral difference adds
a positive rank-one matrix. Its sample uses the same conditional return
probability on each matching diagonal. Conditional binomial Hoeffding
bounds give radii $r_j$ simultaneous with the four unconditional bin
checks under a 1% family budget. $P(0)=1$ by repeatability.
The matrix error is at most $2\sum_{j=1}^{16} r_j$ by row sums.
If a negative eigenvalue exceeds that error in magnitude, the reader rejects
the free-return interpretation. A positive spectral measure requires a
positive semidefinite moment matrix.
In the native experiment this is 0.013095, versus the smallest measured
eigenvalue 0.157956. Thus at least 17 distinct sampled phases contribute.
A set of $s$ energy levels gives at most $1+s(s-1)$ phases, so **at
least five energies** contribute on the evidence of the public records
alone under the freely evolving sharp-state model. Sampling can fold
distinct frequencies onto the same phase. The count gives no individual
frequency value or rest-energy assignment.

There is also a bound on spectral weight. If the best $s$ levels have
total weight $1-\epsilon$, their mutual differences account for a positive
matrix of rank at most $r=1+s(s-1)$. The positive remainder has trace
$17[1-(1-\epsilon)^2]$. By the minimum trace on a complementary subspace,
the sum of the smallest $17-r$ eigenvalues of $K$ cannot exceed that
trace. Lowering each measured eigenvalue by the error bound and flooring
it at zero gives a conservative lower bound on $\epsilon$. For the native
preparation the bounds outside the best one, two, three and four levels
are 52.98%, 36.72%, 18.20% and 3.80%.

**Allowing free motion during the read.** The acquisition switches the
material source off during both sharp interrogations. For the finite-window
comparison, give each interrogation a fresh reader word in the same patch.
After the four forward stages, the private word is fixed by the public
yes/no result. Retaining it until both interrogations end gives the same
sharp measurement on the loads. Each reversible single-seam stage
has a Hermitian involution $G$. Running $\pi(I-G)/(2u)$ for duration
$u$ implements it. With the free source active throughout two interrogations
of duration $\tau=4u$ each, the total variation error of the two public
words is at most $\delta=2\|H\|\tau$ by Duhamel's formula and contraction
under readout. This compares equal free intervals *between* interrogations,
with the reference booked on those intervals. On the binary federation
the exchange is half the token-graph Laplacian, with degree at most 12,
so its spectrum lies in $[0,12]$. Subtracting $16$ times the adjacency
projector gives $-16I\le H\le12I$ in the registered energy convention.
For the native positive-herald probability $p>0.000060$ and
$\tau=2^{-32}$, every conditional return probability changes by less than
$2\delta/(p-\delta)<0.00025$. Thus the finite width can perturb the
17-by-17 return matrix by less than 0.00425 in operator norm. This is
smaller than the measured 0.14486 margin above the rank-error threshold.
The bound costs a supplied peak stage-generator norm $4\pi/\tau$.
It describes an instrument tolerance. The archived counts use the sharp
switched service.

Three returns alone give another test. Write $d_j=1-P(j\pi)$.
For a single gap of any amplitude or frequency,
$d_1d_3-(d_2-d_1)^2=0$ by the cosine triple-angle identity. The native
joint records give a lower bound of 0.48358 on that excess. The isolated
swap and parked controls satisfy the single-gap test. Across all arms the
acquisition represents $37{,}383{,}395{,}344{,}384$ independent simulated
starts from aggregate binomial draws. At the first native appointment
33,225,892 of $2^{39}$ starts gave a positive herald. [Detector, finite
sector and reader](return_lines.py) · [Independent checks](test_return_lines.py)

This is a spectrum of one locally formed **state**. A many-mode localized
preparation may belong to a single particle
species. Adding a constant to this entire number sector leaves every
reported probability invariant. Identifying a rest mass requires a
mobile dispersion and an operational energy or clock comparison linking
sectors. The [torus experiment](../../benchmarks/eight_slot_torus/formed_clocks/README.md#fringes-after-a-marked-birth)
uses the same principle of comparing local records with distinct dynamics
and distinct instruments.

### Replay and numerical accuracy

From the repository root, set the BLAS thread count to one. The short checks are

```sh
python -B -O -m model.experiments.interval_echo_encounter.confinement
python -B -O -m model.experiments.interval_echo_encounter.sign_read
python -B -O -m model.experiments.interval_echo_encounter.conjugation
python -B -O -W error -m model.experiments.interval_echo_encounter.return_lines
python -m pytest -q model/experiments/interval_echo_encounter
```

The extended checks replay the published release, the changed reference and the independent error enclosure.

```sh
python -B -O -m model.experiments.interval_echo_encounter.transport
python -B -O -m model.experiments.interval_echo_encounter.conjugation --pair-regression
python -B -O -m model.experiments.interval_echo_encounter.conjugation --singleton-pair
python -B -O -m model.experiments.interval_echo_encounter.check_transport
```

For the native, exchange-only and reversed releases, a degree-24 Taylor polynomial advances the centered generator by 1/8. The absolute row-sum bound controls the truncation and each sparse matrix multiplication. The polynomial replay assigns each distinct port pair to its terminal address independently of the sparse-exponential replay's detector masks. Accumulating the one-step error through the longest exposure and summing over disjoint terminal bins gives an L1 error below 1.2 × 10⁻¹⁰ for every archived distribution in those three arms. The centered-reference column uses a separate sparse exponential. All runs use prepared backgrounds, supplied isolation and exposure, and a terminal read. The simulation compares observer record probabilities; the exposure variable has no laboratory time calibration. [Numerical enclosure](check_transport.py)
