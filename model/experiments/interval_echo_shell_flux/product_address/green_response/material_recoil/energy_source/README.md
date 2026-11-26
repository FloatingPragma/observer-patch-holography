# Source excitation measured by an independent guide-release counter

A prepared material excitation releases guide pulses at a rate proportional to its energy above the ground of its four-link circuit. The ordinary one-quantum source gave 0.250012 clicks per preparation. The held two-quantum source, using the other mobile material, gave 0.500139. Each click launched one source in the three-guide experiment. On its first outward face both mobile probes moved inward. Doubling the excitation changed the first probe's displacement per preparation from −0.01674 to −0.03331 receiver links and the second probe's from −0.01636 to −0.03305. The predictions were −1/60 and −1/30, with simultaneous 99% half-widths below 0.00098 links. Both ground preparations produced zero marks. Each bank has 134,217,728 starts. The terminal counts come from a seeded simulation of the native measurement probabilities. [Frozen protocol](PROTOCOL.md) · [All counts, independent spectral lines and release records](results.json)

## Two ways to read the same preparation

An observer prepares a fresh single load on four connected ports of equal degree. Write 00, 01, 11, 10 for their cyclic order. During unmonitored flight only the four ring exchanges are open. The full native diagonal and the degree-completing onsite term give the same energy at each of the four positions. The source operations precede stocking the guides. An isolated full exchange at the active stiffness interchanges the two sharp labels with their complex amplitudes preserved. One can temporarily park any component on an adjacent vacant degree-six port. A whole-board diagonal hold for one clock duration π gives this component a phase of −i relative to the other degree-five ports. Undoing the swap restores the source. The phase and mixer are identical for the two load species. These are native, clocked laboratory operations. The ring, the choice of four contacts and the conversion of one mark into one stocked guide release are supplied. [Native seam terms and both readers](instrument.py)

With the additive sector ground energy removed, the free generator is

$$H_{\rm exc}=I-\tfrac12(X\otimes I+I\otimes X).$$

The energy reader parks the second bit once and the first bit once, giving each component a quarter-turn phase for every occupied binary place. Both bits here label a single position. A single unmonitored ring flight for π/2 then sends each standing-wave sign pair to a sharp terminal address. Its Hamming weight reads an excitation of 0, 1 or 2. The signs are independent of the sector ground offset. Two additional interference preparations resolve the single and double gaps on a separate bank. Their even-return fractions at duration 0.6 give ordinary-material gaps of 1.00018 and 1.99965, and held-material gaps of 1.00011 and 1.99991. The first-branch intervals in the record contain the two independently assigned spacings. The end-point signs use inverse sorting of a sharp port. A half exchange, phase parking and full-swap routing prepare the two-site coherences.

The release reader draws just one of the four links. On its second endpoint it parks the phase, operates an isolated half swap, and records whether the load arrives at that endpoint. Other positions give a silent mark. For link $e=(a,b)$ this measures $P_e=(|a\rangle-|b\rangle)(\langle a|-\langle b|)/2$. Each site belongs to two links and each off-diagonal hopping appears once. Therefore

$$\sum_e P_e=H_{\rm exc},\qquad
p_{\rm release}=\tfrac14\operatorname{tr}(\rho H_{\rm exc}).$$

This identity holds for every density in either one-load sector, including mixtures and coherences outside the measured panel. The two readers use different measurement operations and independent copies of the same prepared stock. The measured single-quantum slope predicts the held spectra and source rates within their simultaneous intervals. In a control the first of the four draws is left silent. Antisymmetric dimers on that link and on the opposite link have equal excitation 3/2, while the control probabilities are 1/8 and 3/8. The full reader gives 3/8 on both. This distinguishes collection geometry from the energy of the material. [Checked four-link identity](EnergySource.lean) · [Native operator, spectra and held tests](test_energy.py)

## What the distal material measures

For each emitted pulse a stocked guide triple starts at a uniformly sampled central corner. Surviving appointments and the two orientations of the independently prepared receiver use the measured-step and ternary-contact protocols. A silent source starts no walk and contributes zero. For a source click probability $p$, the signed mean initial acceleration per source preparation is $-pJ_i$ in the matched receiver, where $J_i$ is the outward current per conditional release. The guide Green theorem supplies its inverse-square asymptotic in the interior address chart. At the registered finite shutter the conditional displacement on the tested face is −1/15. Both predictions count every attempted preparation. The observer reports a material terminal position and an independently read source excitation. [Common receiver](../common_receiver/README.md) · [Current and scaling domain](../../README.md)

Source and energy probability bounds use Hoeffding's inequality over at most 64 bank means with joint error 0.005. The mean sorted energy divided by two lies in the unit interval. The ordinary calibration interval propagates both the source-rate and independent energy error. The six receiver records use the complementary 0.005 error budget, conditioning the stopped-walk Bernstein bound on the measured release count, then adding the source-rate uncertainty. Integer cohorts sample the independent walks through drain, including idle and unmatched visits. The sampling grid in the source banks has spacing $2^{-40}$. The record contains every terminal bin and last-drain appointment. [Acquisition and bounds](measure.py)

An observer can change the source's excitation without changing its material tag and can compare two load species with a common excitation standard. A shift of the entire sector energy leaves these preparations, beat records and source marks invariant. A physical gravitating rest mass, the emission law of a freely formed body and proper source–probe distance therefore require additional records.

Replay from the repository root with a single numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.energy_source.measure
```
