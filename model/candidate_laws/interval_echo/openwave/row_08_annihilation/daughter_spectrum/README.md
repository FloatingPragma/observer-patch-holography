# Far-dimer spectroscopy of the phase-wall daughters

OpenWave records OW-08 as untested at [head 0df73dc6](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L443). The full local criterion is partial. On the existing phase-wall source, two far-dimer probes resolve downward transitions in the outgoing material after the central ring orders. At exposure 32 one pair of frequency settings gives **68,604 double clicks in 67,108,864 wall launches**, against **8** on the ordered sham. This adds an independently calibrated response of the same daughters to the contact and energy records. [Protocol](PROTOCOL.md) · [Complete spectra and detector counts](results.json) · [Measurement](measure.py)

## What energy can the daughters deliver?

The source, guide preparation, dimer guards and free flight are identical to the three-dimer-wing experiment. Both ordered central phases are retained. Once the ring is read, its observed orientation is protected during daughter spectroscopy. Each guide's observed boundary is zero or two. Its first three excitation gaps are 31.80075564, 32.51538098 and 33.21463444 exchange units. These differ from the 0.000945849 soft line of the separate converter's boundary-one guide. The native far-dimer swap couples the three hard modes to the ground with squared matrix elements 0.24976331, 0.49952663 and 0.24976331.

Each local detector coherently averages that swap over 32,768 time bins, then selects a downward-frequency clock output. The independent calibration records all eight wing levels and the ground. The quoted coincidence uses 32.51538098 on the left and 31.80075564 on the right. Its measured wall-minus-sham probability is 0.00102216, with simultaneous 99% Hoeffding half-width 0.00037201 across all 54 target/time/arm comparisons. Four of the nine pairs resolve a positive excess at exposure 32. At exposure 8 none clears that conservative threshold; all records are retained.

A separate complete wing-energy measurement gives **2,064,278 hard/hard records among 2,066,989 ordered-ring records (99.8688%)**. The native probability predicts 99.8668%. Its exact conditional mean excitation is 32.4411 per wing. This spectrum and the local detector use independent simulated launches. They measure first-band material transitions rather than counting every excited-wing flag as one quantum.

## What the controls distinguish

Both source-wing contacts are cut at release on a matched copy with the original guide input retained. This quench changes the guide energy convention. It leaves a small guide background, and the ring-order probability differs between sources. The detector response normalized to that probability is source-independent within 4.5×10⁻¹⁶, as the separated Hamiltonian requires. The two hard probes give zero double counts in both cut-source banks. Zero contrast on the same guarded routing also gives zero at these fixed hard-frequency settings. The complete control spectra are included; the absence of these clicks is a frequency-specific response.

The coherent clock, long time window, inverse wing evolution, protected source and energy-selective recombination are supplied apparatus. The pulse acts at the far dimer, but its Fourier analysis uses the whole guide's free evolution. A selected clock mode is an ideal transition readout; no autonomous absorber, detector work budget, photon number, massless dispersion or transverse polarization is derived. The result supports locally readable transport of the prepared material energy. Annihilation into a radiation field requires an independent identification.

## Replay and independent checks

From the RER root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_08_annihilation.daughter_spectrum.measure
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_08_annihilation/daughter_spectrum
~~~

The default replay reads the archive. An explicit --output PATH creates a separate first record and fails if it exists. The comparator uses explicit exceptions, including under Python optimization. Six tests reconstruct the wing operator from endpoint loads, sum a finite clock directly, decode both ring phases and all failures, check controls, reject invalid probabilities and optimized receipt mismatches, and verify spectral grouping under rotations inside an exactly degenerate test level. The raw probabilities are retained; their explicit sampling-complement corrections are below 3.4×10⁻¹⁶. Native control spectra retain unequal-load exchange diagonals and are not approximated by uncoupled dimer flips.
