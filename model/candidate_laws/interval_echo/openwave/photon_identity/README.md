# Photon search with independent charge sources and receivers

**A parked charge changes what a distant material probe records. On the 80-cell sphere, the probes resolve one oscillating component at the lowest frequency. A second independently distinguishable component at that frequency is below 0.97% of the first at simultaneous 99% confidence.** This bound applies to the 24 tested source–receiver combinations in the positive ring-Laplacian written-field extension. Photon identification remains open. [Design and follow-up](DESIGN.md) · [Code](audit.py) · [Complete record](results.json)

| Finite spherical fixture | 80 cells | 320 cells |
| --- | ---: | ---: |
| Active spherical faces | 42 | 162 |
| Reachable field words from the prepared sector | 75 | 2,061 |
| Physical source seams × physical probe seams | 4 × 6 | 4 × 6 |
| Nonzero frequencies resolved from the generator | All 74 | Lowest 23 |
| Resolved frequencies with response rank one | 74 | 23 |
| Strongest absolute signed crossing contrast at the five survey times | 0.185596 | 0.025985 |
| Probe crossings after the two opposite source signs at that setting, out of 65,536 each | 35,743 / 23,514 | 19,111 / 17,145 |

At field time 128 the largest contrasts occurred between faces 19 and 21 on the smaller sphere and faces 89 and 101 on the larger one. The latter pair is five face adjacencies apart. Positive, neutral and negative sources were each read with separately prepared positive and negative probes. The signs of the two probe responses reverse. Switching off the intervening field evolution gives zero remote signed contrast throughout the survey. The retained four-outcome records include failed interventions and failed probe crossings. The figures in the last row are simulated counts. [Station maps, complete outcomes and control](results.json)

## What the material detector measures

A signed source on a private seam can change one head bit and move to its parking arm. Its sign determines which starting orientation permits that move. On the complete sphere the reference field is the uniform ground on all reachable words. A successful intervention blocks the source face during free flight. Every other originally admissible reversal persists. The simulation tests this statement on every reachable word and every selected source seam. Source and receiver share neither their material preparation nor their readout. The field carries the response between them. [Physical moves and apparatus checks](audit.py) · [Earlier written-field flight](../electron_field/README.md)

Let the head on the source seam be the zero-one column $x_s$ over the $N$ reachable words, and write $z_s=2x_s-1$. Define $z_d$ on the receiving seam in the same way. The successfully moved branch remains stationary because its selector is constant on every connected component after the cut. The other branch evolves under $U(t)=\exp(-itH)$. Including both source outcomes, the difference between positive-probe crossing probabilities for opposite source signs is

$$
C_{ds}(t)=\frac{z_d^T[I-\cos(tH)]z_s}{2N}.
$$

At angular frequency $\omega>0$ the coefficient of $1-\cos(\omega t)$ is $R_\omega=Z_D^T P_\omega Z_S/(2N)$, where $P_\omega$ projects onto that eigenspace. Its rank counts how many independent current-to-probe response patterns these electrodes distinguish at that frequency. All reported frequencies are simple at absolute resolution $10^{-8}$. Their residues have rank one. A calibration with two supplied, independently addressed copies gives rank two. Summing the two receiver tags returns rank one. A second tag therefore needs physically independent preparation and readout to appear in this test. [Eigenpair residuals, singular spectra and controls](results.json)

## What can be recovered from counts

On the 80-cell sphere, 256 irregular arrival times spanning a prescribed window of 4,096 field units give about 2.6 periods of the lowest oscillation. Each of the 24 station pairs has 65,536 simulated trials per source sign at each arrival time, for 805,306,368 enrolled trials. The fit receives the positive-probe crossing counts, both source signs and the frequencies of higher oscillations calibrated from the generator. Letting the lowest frequency vary recovers 0.004004446 per field-time unit. Its independent generator value is 0.004004700. The fit of the six-by-four response matrix imposes no rank restriction. [Outcome multiplicities, arrival times and frequency scan](results.json)

At the calibrated frequency the first two measured singular values are 0.646096 and 0.000734. For the fitted time weights $a$, each matrix entry is a weighted difference of two independent binomial proportions from $n=65{,}536$ trials. Hoeffding's inequality and a union bound over 24 entries bound the maximum entry error by $\lVert a\rVert_2\sqrt{\log(48/0.01)/n}=0.0011134$ with 99% simultaneous coverage conditional on the calibrated frequencies and ideal apparatus. Bounding the operator norm by the Frobenius norm yields an error of at most 0.0054544 on either singular value. The first component is at least 0.640641 and the second at most 0.006188. The finite counts support the one-component conclusion at the quoted sensitivity. The large sphere's 23 ranks come from its finite generator and the physical source/probe effects. Its lowest period is about 50,011 field units. The longest observation at 128 units covers about 0.26% of that period. [Unrestricted fit, exact instrument check and bounds](audit.py)

## Emission, flight and energy reception

The registered repair pulse writes a coherent neutral seam state from an excess and a vacancy. On uniformly filled matter the free Hamiltonian separates into independent memory operations. The original dipole and its phase can be recovered at the source. An independently calibrated remote receiver gets zero information from that emission. A mobile material defect can alter that same receiver. [Native emitter and receiver](native_emission/README.md)

The written-field extension has an additional neutral propagation law. Energy leaves the parked source and reaches two independent rings that both stay silent on the ground. On the 320-cell sphere, two disjoint receivers eight and ten face steps away gave 74 double clicks in 1,048,576 positive-source releases. Reversing the source leaves the energy probabilities exactly invariant while it reverses the signed material-probe response. These collective energy pointers read the local field Hamiltonian. [Remote energy and coincidence experiment](energy_receivers/README.md)

On the finite quantized Maxwell comparison, a charged circuit at one end of an open box heralds a packet with greater than 99.99999% one-quantum occupation. A distinct charged circuit at the other end distinguishes the two opposite phase windings by coherently comparing an early and a late interaction. These distributed records require the local emitter herald and retain both receiver misses and successes. The field-number counter can be left unread. The quantized field, charged gates and reference clocks are construction inputs. [Polarized model photons and occupation controls](quantum_transfer/README.md)

The globally neutral excess–vacancy pair is a different native preparation. On the blank seam background its entire spectrum costs at least 80 native energy units above the uniform unit-filled reference. Adjacent pairs can reach a neighboring cell with nonzero probability. [Neutral-pair formation and transport](native_emission/NEUTRAL_PAIR.md)

## Which signal could be a photon?

The registered Interval Echo ranging experiment sends a prepared material double along a reciprocally serviced route. Its departure, reflection, return and remote trigger are material records. They establish message delivery and optical timing under the supplied services. The two signs of the written-field probe instead report a local seam response following autonomous field evolution in the extension. Their opposite responses measure charge coupling. Circular polarization would require a calibrated rotation about a measured direction of travel acting on the same sourced mode. [Registered optical apparatus](../../../../experiments/interval_echo_optics/README.md) · [Written-field experiment](../electron_field/README.md)

The earlier guided experiment reads three oscillations from independent signed probes on two face chains. Replaying the original crossing vectors reproduces second-to-first frequency ratios 3.801938 and 3.941884. The measured branch has $\omega=1-\cos k$ in face and field-time units, hence quadratic small-$k$ dispersion. This branch has a different dispersion from a linear long-wavelength photon. The full-sphere rank result measures resolved modes in the selected finite sector. Curvature, finite-size splitting or inaccessible polarizations in another sector require additional apparatus and size scaling. [Guided replay and free-sphere survey](results.json)

The FZ-12 oscillator theorem gives a two-dimensional transverse space at each nonzero momentum and $\omega^2=\Lambda(k)$ for the declared Hamiltonian lift. Its independent replay passes. The observer test needs a physical map from these transverse coordinates to independently rotatable current preparations and polarization filters on the same law and clock. A one-quantum emission and absorption or a resolved nonclassical photon statistic would additionally test the particle interpretation. [Finite oscillator theorem](../../../../../Lean/Screen/SeamCurrentFreePhotonLift.lean) · [Independent replay](../../../../../code/a5_fingerprint/verify_fz12_free_photon_hamiltonian_independent.py)

The historical finite Wilson–Maxwell comparison supplies a quantized electromagnetic field coupled to matter. On a nonstationary free preparation at coupling $g=1/64$, the proved unconditional probability of one counted particle, one counted antiparticle and one counted field quantum exceeds $4.75\times10^{-8}$. Switching off the coupling gives probability zero. This is a probability bound for a supplied finite field and count instrument. It does not measure a photon detection in the registered Interval Echo trace. [Archived quantized interaction and complete count instrument](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/native_wilson_maxwell/PROOF.md)

## Replay

Run from the repository root with the project's NumPy, SciPy and NetworkX environment. Choose a fresh output path. The numerical run verifies every archived file against its SHA-256 manifest, replays the complete instrument and the guided dispersion records, rebuilds the full spherical sectors, and checks the recorded fit.

```sh
.venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/photon_identity/audit.py --output /tmp/photon-identity.json
.venv/bin/python -m pytest -q model/experiments/interval_echo_optics
.venv/bin/python -B code/a5_fingerprint/verify_fz12_free_photon_hamiltonian_independent.py
```
