# Does the arriving signed charge carry spin one-half?

A boundary seam head selects a charge arriving in cell 85 of the 320-cell board. The receiver compares two perpendicular half-turns in opposite orders. All three family preparations give the constructive output with order phase +1. A supplied spin-½ control gives the destructive output with order phase −1. A circulation encoded at the source can also send quantum information into the flux-positive receiver. Its calibrated address rotations give +1 on this preparation. [Protocol](DESIGN.md) · [Complete spin readouts](results.json) · [Field-loop records](circulation.json) · [Local tomographic counts](counts.json)

## What the observer actually measures

The source puts an excess or vacancy in a prepared written guide with a fixed compensator. Each local hop reverses its seam head. The receiver compares the entrance head with its prelaunch direction. A signed change is the arrival flag. Independently counting the material load and every interior incident head verifies this identity on the complete sharp basis, including the incidence grades. A single binary flux read keeps coherence between different locations inside the cell. A missed arrival stays in the trial record.

For the written-contact experiment, the entire eleven-edge path and eleven interior detector edges are open under a time-independent generator. For the family comparison, the first-order guide carries each character to the sole terminal port, and the prepared interior tree opens for readout. Conditional on the flux click, a specified controller isolates one contact at a time. Each 4-weight positive contact swaps its endpoints and flips its head at exposure π/8. Compiling the cell rotation from those pulses moves both the charge and its tree dressing. The apparatus covers all twelve sharp local charge addresses and all fifteen pairs of perpendicular axes. Its rotations are active terminal operations on the dressed charge. The anisotropic guide need not be invariant under them during flight.

| Source at exposure 8 | Probability of a signed flux click | Destructive probability given click |
|---|---:|---:|
| Written contact, sharp source, detector tree open | 21.7111% | 0 |
| Light character, positive-energy packet | 12.4820% | 0 |
| Middle character, positive-energy packet | 13.3130% | 0 |
| Heavy character, positive-energy packet | 18.6648% | 0 |

Both charge signs have these probabilities. In the equal-order sham the destructive probability is zero. In the spinor control it is one given a click. Cutting the entrance gives zero arrivals from the sharp launch, to numerical precision. The family register acts through the supplied mass term during flight and is unread at the detector. The measured +1 tests the implemented action on these dressed charges. A spin carried in a different material representation requires its own action and coupling to the apparatus. [Contact-level acquisition](measure.py) · [Guide generator and sources](../dirac_read/README.md)

## Can a field circulation travel with the charge?

Prepare the two sharp circulations of a triangular source face coherently. Open a triangular face at the receiving end as well. With both source tips of equal degree, changing the source from even to odd circulation raises flux detection from 15.7070% to 24.6265% at exposure 8. The source-tip reflection commutes with each part of the generator. After the herald it acts only on the exterior. The receiver cannot recover the even–odd coherence at any wait, even though it can detect their different arrival rates.

Open one existing source contact, 120–121, with the same weight. The receiver's sharp labels specify the charge port and its internal head word. Six source settings and 343 local receiver settings reconstruct the channel on a reference entangled with the two source circulations. The simulated counts give minimum partial-transpose eigenvalues −0.014302 for a vacancy and −0.014304 for an excess. A simultaneous 99% bound on the eigenvalue error is 0.004597 over both signs and both layouts. The mixed input reaches the receiver with probability 22.1126% at this wait. The eigenvalues use the unconditioned branch weights throughout. Thus the receiver can retain entanglement with a source reference even after the exterior records are discarded. The stopped-ring control also transmits coherence on the asymmetric layout. [Channel, changed geometry and rate control](circulation.py) · [Finite-shot reconstruction](tomography.py)

Dark counts from the ring-bearing cell require a separate calibration. At exposure 8 one pair of commands sends 56.2% of heralded events to the destructive output on the even, symmetric-source preparation. The 38 accessible sharp configurations distribute as either two or four at each charge address. None of the fifteen half-turns preserves these address ranks on the entire sector. For each axis pair the calibration keeps the largest support closed under both commands on which every sharp address turns correctly. Ten pairs have such supports of eight or sixteen configurations and five have empty support. The two orders agree on every retained complete configuration, including its heads. The dark port outside those supports measures path-dependent routing. It cannot establish physical spin-½.

## Replay

From the repository root, choose unused output paths.

```sh
python -B -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin.measure --output /tmp/arrival-spin.json
python -B -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin.circulation --output /tmp/arrival-circulation.json
python -B -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin.tomography --output /tmp/arrival-counts.json
```

The exposures and couplings are in model units. The coherent launch, controlled rotations and local tomographic analyzers are specified simulation apparatus.
