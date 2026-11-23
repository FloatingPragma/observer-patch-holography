# Does the neutral repair record leave its source?

**A registered Interval Echo repair converts an excess and a vacancy into a locally readable neutral record. Its sign and relative phase stay on the emitting seam during free evolution.** An independent observer in the neighboring cell detects no change for any choice of the source state, even when that receiver starts with a coherently oscillating memory. A separate mobile-defect control changes the receiving memory and leaves it changed after the defect has gone past. [Measurement design](DESIGN.md) · [Executable](measure.py) · [Complete settings and outcomes](results.json)

## A charged source writes a neutral state

Use one load per port as the material background on the 20-cell federation. Prepare an excess and a vacancy on seam 0–1 in either orientation. The source applies the registered 99-state echo. Input (2,0) ends at (1,1) with written contrast +2. Input (0,2) ends at (1,1) with written contrast −2. Both outputs have the same even coin state. Linearity preserves a superposition of the two orientations and its relative phase. The source intervention supplies the repair pulse. All material ports hold one after the pulse, and every physical seam is active during flight. [Native echo and free generator](../../../../contrast_memory_spherical_v1/law.py)

The two emitted states lie in the image of the write isometry. Their local seam energy is zero. Applying the inverse echo after any tested wait recovers the two input dipoles and their coherent superpositions exactly up to matrix roundoff. The same direct-write calibration on the receiver gives probability one for the matching sign and zero for the opposite sign. These signs label contrasting orientations of the original material dipole. The read calibrates their distinguishability locally. [Echo columns and sign effects](measure.py) · [Dipole readings](results.json)

## What every independent receiver can see

For each of the backgrounds with constant material occupancy $n=0,1,2$, apply four times the exact local Hamiltonian to all eleven input memory states with loads $(n,n)$. Every output entry outside those same two loads is zero. Consequently the whole uniform-material subspace contains arbitrary states of every seam memory, and its autonomous Hamiltonian is

$$
H\big|_{n\,\mathrm{uniform}}=E_{\rm degree}(n)I+\sum_e h_e^{(n)}.
$$

Each $h_e^{(n)}$ acts on just one memory. Its propagator factors into one unitary per seam and a global phase. This applies at all times, including densely prewritten backgrounds and quantum-correlated memories. A local trace-preserving source intervention followed by this factorized flight cannot alter an independent receiver's outcome law. For correlated initial memories, comparing conditionally selected source outcomes would use a different ensemble. Interventions that launch mobile matter can leave the uniform-material subspace. [Exact integer blocks for all 33 columns](results.json)

Prepare seam 72–82 in cell 6 as the receiver. A blank receiver stays blank with probability one. An active receiver starts in an equal superposition of two coin-odd modes whose energies differ by 32 native units. Its local Ramsey yes probability is $\cos^2(16t)$. This falls from one to zero at $t=\pi/32$ and returns to one at $t=\pi/16$ under all four dipole settings. The code checks the complete reduced receiver state, every memory outcome probability and the Ramsey yes/no counts. The measured station is sensitive to local evolution while the source choice has no effect on it. Each blank or absolute-contrast projector of any individual seam also commutes with the full free Hamiltonian on all material states. The counts of empty, single and doubled ports each commute with the full 99-state free term. A freely moving single excess or vacancy cannot create an additional material excess–vacancy pair on this background. The conserved-load current $i[h_e,n_b]$ annihilates all eleven memory columns on each uniform input. An independently prepared material dipole has a different initial state and a gapped formation spectrum on recovered seams. [Neutral material pair flight](NEUTRAL_PAIR.md) [Full receiver outcomes](results.json)

## What changes after matter passes through

Release one excess or one vacancy at port 1. The route crosses the cell boundary to port 82, traverses the prepared receiving seam to port 72 and continues to port 73. The boundary and outgoing seams start blank. The receiving seam starts in a sharp coin state and a symmetric superposition of contrasts +1 and −1. Every spectator seam starts in the balanced zero-contrast write. Its term annihilates all nine local load configurations on that word. The full native evolution therefore preserves the four defect locations and eleven receiving memory states used in the computation. The degree corrections stay on. The outcome reads the defect position, relative to the uniform background, together with the local memory state. All 44 outcomes are retained. [Invariant sector and exact reachable position–coin pairs](measure.py)

At exposure 4 the following joint event was observed. The defect had reached port 73, beyond the receiving seam, and the receiving coin had flipped. Both local ports at the receiver had returned to their baseline occupancy.

| Released defect | Initial receiving coin | Joint event probability | Simulated joint events in 8,192 trials |
| --- | ---: | ---: | ---: |
| Vacancy | 1 | 0.00666970 | 55 |
| Excess | 0 | 0.00020078 | 3 |

For the opposite receiving coin, the joint event is outside the exact connected subspace. Parking the initial defect using an idle written boundary or setting the whole law's contrast strength to zero also excludes it exactly for either sign and coin. This control shows why a changed receiver with baseline material occupancy at the final read is insufficient evidence of independent neutral flight. A mobile defect can pass, write the receiver and continue. [Complete count distributions at six exposures](results.json)

Reproduce the record from the repository root using a fresh output path.

```sh
.venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/photon_identity/native_emission/measure.py --output /tmp/native-emission-replay.json
cmp /tmp/native-emission-replay.json model/candidate_laws/interval_echo/openwave/photon_identity/native_emission/results.json
```
