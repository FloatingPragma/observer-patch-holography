# Where could a graviton appear in the observer records?

The pure Einstein branch has two classical, transverse, traceless polarizations around flat spacetime. An exact replay shows that the ten observers used to reconstruct the Einstein tensor give **identical source readings** in flat spacetime and in either vacuum polarization. The curvature differs and oriented tidal readings resolve both polarizations. On the registered Interval Echo board, a separate test finds two readable angular patterns of stored energy whose orientations make **exactly zero difference** to a remote material tracer. A quantum graviton identification remains open. [Test design](DESIGN.md) · [Complete outcome probabilities and simulated counts](results.json)

## Which result follows from the construction?

The gauge paper defines the quantum-particle conditions Q1–Q3. On the specified Einstein–Hilbert branch, the reduced quadratic action gives two massless classical metric modes. Particle identity requires positive-energy quantization, a physical two-point function with a positive-residue massless pole, and the asymptotic-state hypotheses of that phase. The string selector's conditional theorem assumes a positive-norm physical closed-string state and pole, agreement of its normalization and classical limit with the OPH metric, and a single visible tensor sector. The ten-observer Einstein theorem concerns regional entropy, area and stress in a common geometric limit. Its exact decoder recovers the symmetric tensor supplied to it. Establishing that input on Interval Echo is a separate physical requirement. [Classical carrier modes and the quantum-particle criterion](../../../../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex) · [Conditional tensor-carrier quantum completion](../../../../../paper/archive/observer_patch_holography_as_string_vacuum_selector.tex) · [Ten-observer theorem and hypotheses](../../../../../code/matched_einstein_limit/README.md)

## What the ten observers can see

Start with a linear metric perturbation on a flat background with zero cosmological constant. Calculate its Riemann tensor, contract to Ricci and the Einstein tensor, and apply the repository's ten independent timelike contractions. In parallel, calculate twice the relative acceleration per unit separation along five geodesic arms. Exact rational arithmetic gives the following responses at unit second derivative of the wave profile. The diagonal directions are the normalized spatial vectors along x+y and x−y. The contraction column lists the largest absolute reading. [Frame and inverse](../../../../../code/matched_einstein_limit/frame.py) · [All ten readings for each configuration](results.json)

| Supplied metric configuration | Largest absolute Einstein reading | Tidal x | Tidal y | Tidal diagonal + | Tidal diagonal − | Tidal z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Flat reference | 0 | 0 | 0 | 0 | 0 | 0 |
| Vacuum plus wave | 0 | 1 | −1 | 0 | 0 | 0 |
| Vacuum cross wave | 0 | 0 | 0 | 1 | −1 | 0 |
| Coordinate change of flat metric | 0 | 0 | 0 | 0 | 0 | 0 |
| Transverse breathing control | 9 | 1 | 1 | 1 | 1 | 0 |
| Plus perturbation off the null cone | 6 | 1 | −1 | 0 | 0 | 0 |

Both vacuum waves carry Weyl curvature while their Ricci tensor vanishes. The rank of the Einstein readings on this two-dimensional space is zero. The rank of the tidal readings is two. The coordinate-change and breathing controls distinguish a curvature response from a coordinate displacement or a scalar transverse expansion. The off-null control shows why measuring a plus-shaped tidal pattern alone leaves dispersion undetermined. The Lean calculation proves the vanishing Ricci tensor for every superposition of plus and cross. At any known nonzero wave scale s, the x arm and the x+y arm recover their amplitudes multiplied by s². These results apply to ideal linearized continuum inputs. Higher-order responses, separate regional observables and histories containing physical range or tidal information require their own treatment. [Formal curvature and tidal-read proof](../../../../../Lean/Geometry/VacuumTidalRead.lean)

## What the registered contacts do

Five seams around one pentagon in the home cell of the 20-cell federation store a stationary energy pattern. Each preparation puts an equal coherent superposition of energy 0 and 32 on one selected seam. The other four hold the zero-energy word. Changing the probability of selecting each seam produces independent cosine and sine components at twice the angular frequency around the pentagon. Both signs and the isotropic reference have the same mean total stored energy, 16 in native units. Sharp preparations, angular mixture weights, the ideal local Ramsey read and the exposure scale are supplied instruments. The energy and propagation use the registered contact generator at κ = 8. [Existing nine-load algebra](../electron_field/dark_contact.py) · [Native regional and tracer experiments](../../../../experiments/interval_echo_einstein_station/README.md)

| Source mixture | Expected angular coordinates | Coordinates inferred from local Ramsey counts |
| --- | --- | --- |
| Isotropic | (0, 0) | (−0.00047, +0.00073) |
| Positive cosine | (+0.5, 0) | (+0.50308, +0.00032) |
| Negative cosine | (−0.5, 0) | (−0.50023, −0.00113) |
| Positive sine | (0, +0.5) | (−0.00244, +0.49788) |
| Negative sine | (0, −0.5) | (+0.00049, −0.50129) |

The three receiving patches launch an ordinary material tracer at their respective gateway ports. Each patch resolves its twelve port outcomes and includes one outcome for a tracer anywhere outside the patch. All 13 probabilities agree for every energy preparation at each exposure 1, 2, 4 and 8. Each source calibration and remote setting has 262,144 simulated shots. The exact null follows from the contact matrix. Every local load word acquires the same scalar phase within a given frozen seam state. On the prepared sector the generator splits as H_material ⊗ I + I ⊗ H_memory. Tracing out the source therefore removes its phase for every material-only measurement history during free flight. An independent 2,640-dimensional propagation retains all eleven states on one seam and reproduces the reduced tracer probabilities to a maximum absolute difference of 7.94 × 10⁻¹⁵. [Registered-seam propagation and independent comparison](measure.py)

The sensitivity control replaces the selected occupied seam with a blank seam that transmits material. Under the same angular mixtures and observer settings, the maximum total-variation distance between opposite cosine preparations is 0.003126 at exposure 4, gateway 20. For opposite sine preparations it is 0.002997 at exposure 4, gateway 82. These are differences between the calculated outcome distributions. The frozen source is stationary throughout flight. This experiment tests whether the chosen local energy distribution influences a material path. A radiating, time-dependent stress source and a physical tensor detector are still needed. [Matched controls, settings and counts](results.json) · [Complementary material angular transport control](../mediator_controls/README.md)

A direct test of a classical carrier would prepare a varying source, then read differential free-fall separations or calibrated travel times in independent axial and diagonal arms in an initially quiet patch. Acquiring both orientations, the unchanged and scalar-source controls, and arrivals at multiple distances would test polarization and retarded propagation on the same trace. The earlier six-station return experiment already shows why its particular conditioned travel time needs further geometric calibration. A quantum assignment would additionally have to identify a physical spin-two operator and its positive spectral weight at a massless pole across refinement. [Existing clock and geometry controls](../../../../experiments/interval_echo_einstein_station/refinement.json)

## Replay

From the repository root, run the measurement with the project interpreter. It compares the saved numerical probabilities, counts, exact tensor components and design and program digests with a fresh calculation. The optional `--output PATH` writes a receipt to a new file.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/graviton_identity/measure.py
```
