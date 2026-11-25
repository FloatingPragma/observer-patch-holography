# A local impulse and two length controls on the three-arm charged doublet

The charged preparation that transmits a two-state record between the terminal triangles also responds to a longitudinal impulse. A load counter measures its displacement, and local seam contacts measure its energy change and signed current. Lengthening the arms moves every accessible pole away from the original formation energy. This identifies a positive-energy **guided resonance with measurable recoil**. An intrinsic rest energy would require a freely moving excitation and a separate spatial and time calibration. [Original formation and transfer experiment](../README.md)

## Holding the source fixed

The emitter, the three marked heads at each end, the logical input, and the published finite formation-clock setting are reused from the [original receipt](../results.json). On the original guide, the source and clock together accept about 0.000060365 of the starts. All measurements below use either this finite, prepared state or an explicitly labeled projection onto the ideal energy eigenspace. The two-dimensional space is the standard representation of simultaneous arm permutations.

Every ordinary backbone port has an ordinal along its arm. The two common hubs have ordinals zero and 19. These labels specify a detector for longitudinal load position, X. A terminal shortcut spans two ordinals. In particular, the coordinate has no independently established conversion to a physical length in an isotropic space.

## Applying and reading the impulse

A blank path carries any selected port's load to a shared work port and back while restoring every other load and all six marked heads. At the work port, an inverse echo reopens an exterior shutter for a filled port. A native swap with a staged ancillary vacancy, then a shutter contact, temporarily transfers the occupation to a degree-five private reference. Its private dwell supplies a relative phase. Reversing the contacts restores the reference, the staged load and the original shutter. The ancillary vacancy requires an extra balanced apparatus stock. The vacancy and its compensating excess each occupy a separately shuttered blank dimer. Each dimer uses two equal-degree ports and starts in a symmetric stationary state, present also in the clock host. A private quarter-turn and blank half contact localize the staged vacancy before use. The inverse preparation restores the parked stock after the pulse. The gate checks its full operator on the occupied and vacant inputs, including coherent superpositions. Applying the gate at each ordinal supplies the phase gradient exp(ikX). The observer supplies this impulse. [Native phase-gate check](impulse.py)

Separate banks measure three quantities on fresh starts. A relative quarter-turn phase and half contact followed by sharp load reads give the symmetric and antisymmetric populations, hence the seam energy. Omitting the quarter turn gives the current quadrature. On a marked contact a second full pulse and load read discriminate a transferring rail from a blocked rail. Simultaneously read contacts occupy disjoint ports. A third bank measures only charge position after free evolution. Every bank counts preparation refusals. The numerical observer banks use multinomial samples from the projected native probabilities with 2^50 virtual starts per setting. The simultaneous confidence level is at least 99% by a union bound. [Measurement implementation](measure.py) · [Simulated records](readings.json)

For impulses of +0.125 and −0.125 radians per guide step the measured currents are +0.02231 and −0.02226 steps per native time. Both intervals have half-width below 0.000566. Separate load counters measure displacements of +0.05925 and −0.05930 steps at native time one. The zero-impulse current and displacement intervals both include zero. The contact-energy banks resolve an increase for each sign. At +0.25 the energy increase is about 0.00568 native units. The matched frozen-load control suppresses every off-diagonal hop after preparation and gives zero current.

The ideal pole has a real spatial representative. Summing the native hopping terms at longitudinal separations one and two gives its exact impulse-energy dependence in the numerical representation,

$$
\Delta E(k)=c_1(1-\cos k)+c_2(1-\cos 2k),
\qquad c_1=0.56498534,\quad c_2=-0.09732012.
$$

The reciprocal curvature at zero impulse is 5.69136 in native guide units. This number depends on the defined impulse and coordinate. For the entire resolved doublet with projector P and energy E, the identity P i[H,X] P = i(E PXP − PXP E) gives zero. The kick commutes with all arm permutations, and its compression to this doublet is scalar. A current therefore requires coherence across different guide energies. A second use of the original finite clock accepts about 54.14% of the admitted +0.125 kicks and 2.238% of the +0.25 kicks. The original pole test cannot follow a continuously boosted state as a single unchanged finite energy.

## Changing the arm lengths

The detour rule chooses the earliest interior ordinary contact on all three arms with a spare adjacent actual port and replaces each of the three contacts by two blanks. Repeating it once more adds a second physical port to each arm. The old shortcut is shut, private references on newly included ports are erased, and every unused contact touching the active material is shuttered. The active arms remain connected by blanks. Both terminal triangles, their six marked seams, both hub ports, the logical encoder and the clock are identical. The calculation uses the complete standard row for the source on each material. [Detour construction](guides.py)

At arm lengths 19, 20 and 21, respectively, the frozen source-and-clock acceptance probabilities per original start are 0.000060365, 3.8822 × 10⁻¹³ and 5.4272 × 10⁻¹². The nearest source-row level to the original clock center lies 0.01899174 native energy units away for length 20 and 0.00123192 away for length 21. Numerical spectral exclusion allowances are below 0.000000027 in both cases. Small nonzero read rates follow from the finite clock response. Alternately embedding the same terminal controls thus rejects an arm-length-invariant interpretation of this particular resolved line. A mapping from guided resonances to a common infinite-volume mass, including independently calibrated momentum and confinement shifts, remains open.

Run from the repository root, with a single numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.recoil_search.measure
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/recoil_search
```

The measurement command prints the frozen-clock length comparison. Append `--write` to regenerate all counts. The tests check the stored and fresh observer classifications and compare native probabilities across numerical backends.
