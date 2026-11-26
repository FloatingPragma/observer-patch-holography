# Can a neutral pair carry three distinguishable source labels?

A distant receiver distinguished **three addressed sources of the same neutral pair** by the order of the loads arriving at three consecutive ports. Four independently calibrated direction settings on the emitting apparatus predicted all 1,296 held receiving probabilities. When the sender used the three coherent winding preparations, the two opposite windings gave identical remote-pair statistics. The freely moving pair retains information about how it was produced. Identifying neutrino flavours requires separately measured weak charged partners and oscillation angles. [Frozen flight design](PROTOCOL.md) · [Fresh mixture forecast](MIXTURE.md)

Two observers can also compare the distant arrival with a coherent read of the original source directions. In the held portable test, a decoder calibrated at one time and place recognized a winding preference at another receiver. Its native yield was 4.62–9.45 pairs per trillion original starts. The sharp source records must stay available during the flight. The physical weak tags, mass splittings and PMNS angles remain open. [Joint material read](joint_receivers/README.md) · [One-bypass theorem](theory/source_access/README.md) · [Fixed decoder at new appointments](portable_decoder/README.md) · [Test of a unitary flavour interpretation](kernel_check/README.md)

## Reading a whole receiving region

The source uses the established winding preparations and local conversion. A vacancy first travels to port 160. The local stock supplies an excess on 167, leaving one vacancy and one excess in the emitted pair. Both unknown writable heads are retained. The ring for the second freely evolving flight has 14 ports on five owners. The measurements cover all three ring sites at owner 12 (152, 150, 154) and owner 17 (212, 207, 206). The flight ends at one of three separately scheduled appointments. A receiver freezes the material, optionally gives either of its two internal blank contacts a half swap, and reads the three ternary loads. This distinguishes the six orders of one vacancy, one unit and one excess. Each original trial returns a word, a missed pair, a missed initial vacancy at the converter or a source refusal. Independent cohorts use each incompatible analyzer. [Existing emitter and geometry](../../../../experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge/neutral_ladder/release/README.md) · [Native contact and simulation](measure.py)

For the three addressed arms, the same receiver recorded the following native free-flight results at owner 17 after 2,048 exchange-time units. Each setting has 2⁵² original simulated launches. The displayed intervals use exact binomial tails with a simultaneous 99% allowance over 16,384 intervals. The conditioning is solely on both constituents reaching that receiver.

| Source setting | Complete pair arrivals | Word (2,1,0) | Fraction of arrivals, 99% interval |
| --- | ---: | ---: | --- |
| Arm 0 | 20,888 | 6,151 | 29.45% [27.89%, 31.04%] |
| Arm 1 | 65,426 | 51,123 | 78.14% [77.32%, 78.94%] |
| Arm 2 | 25,036 | 10,015 | 40.00% [38.46%, 41.56%] |

Eighteen separate receiver-local calibrations used the six sharp words in all three analyzers. A unit shutter, the native blank transposition, and a neighboring fifth-degree private reference allow an isolated local ternary read. In this source sector a three-port region always contains at least one unit through which the reading access can open. The signed stock, contact halt, calibration and reading time are supplied apparatus. On coherent received packets the phase-sensitive counts disagreed with an incoherent half-swap of the measured sharp-word mixture in 616 of 864 comparisons. Changing only the contrast from eight to four during the second flight resolved a change in all 108 tested whole-region arrival rates. Closing the two inter-owner contacts at release gave zero remote pairs. [Original integer banks](records.json) · [Count-based decoder](readback.py) · [Material and calibration checks](test_readout.py)

## What information reaches the receiver?

The two retained written heads initially admit four direction settings. On each written pendant, moving a nonunit onto or off the leaf flips its direction. The joint leaf-load and direction records keep the four initial settings distinguishable outside the remote aperture. When the whole pair reaches the receiver, both leaves hold units and the two source heads recover their original settings. Tracing the emitter therefore gives a mixture of four six-component receiving amplitudes. This reasoning permits interference between different pair words within a given setting.

At the preceding charged source, any three-bit receiving mask has two equal bits. Swapping the corresponding arms is an exact symmetry. It exchanges the two conjugate winding preparations up to a global phase while preserving that receiving mask. Summing over unread emitting masks gives the same probability of each receiving mask for both windings. Consequently they give equal weights to all four subsequent neutral-flight components. They give equal probabilities in **every measurement on the intercepted pair** at every time, including the arrival rate. A count comparison at the finite enrolled appointments found all 324 conjugate bin pairs compatible. Diagonalizing the full receiving difference, including off-diagonal amplitudes, gave conditional trace distances below 6.4 × 10⁻¹³ on both laws. This is the numerical residual of the symmetry and tracing argument. In the native 2,048-unit example the zero winding has trace distance 0.0390 from either conjugate, conditional on both loads arriving at owner 17. A joint measurement involving the old source heads accesses a different observation domain. [Independent symmetry and spectator check](audit.py)

The three addressed source arms and the three windings have different sharp source-direction populations in general. Independent ensembles measured the two head bits after conversion-source arrival. Independent sharp head stocks measured each of the four conditional word distributions. Another 216 held flight cohorts tested each predicted original-start probability. All **1,296/1,296** individual word comparisons agreed with their predicted intervals. Replacing the head proportions of arms 1 and 2 with the measured proportions of arm 0 failed **359/432**. The model supplies the sharp direction stocks for this calibration. The blind detector reads the same six words throughout. The receiver can distinguish all three addressed orbital preparations by their ensemble frequencies. At the same receiver the three winding recipes produce only two distinct distributions. PMNS mixing and a weak-flavour particle assignment remain open on Interval Echo. [Independent calibration and held banks](mixture.json) · [Forecast and interval propagation](mixture.py)

The exact bound also applies when a receiver observes all available active-guide sites outside the two source leaves and heads, including incomplete arrivals. In a single three-way decision with equal source priors the maximum success probability is at most two thirds. The charged preparation before conversion has a material analyzer that distinguishes both winding signs. A proposed joint read that revisits the source and also receives the distant pair yields a signed contrast forecast; it needs a coherent source instrument and control of the phase through parking. [Finite proof](theory/PROOF.md) · [Joint-read forecasts and charged calibration](ERASER.md)

On a clean ring the same sharp neutral pair can be prepared at three neighbouring source bonds. A remote receiver distinguishes all three positions. Moving that receiver together with the source tests the spatial interpretation directly. Independent measurements of one separation profile predict all 72 held outcomes on the second native setup. [Two translated-pair experiments](translated_pairs/README.md)

## Reproduce

Each file of synthetic integer records samples multinomials in aggregate. Run at the repository root with the scientific Python environment. Both count readers work without acquiring further data.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations
d=model/candidate_laws/interval_echo/openwave/row_10_oscillations
python -Bm "$p.measure" --output /tmp/oph-neutral-flight.json
python -Bm "$p.readback"
python -Bm "$p.mixture" --output /tmp/oph-neutral-mixture.json
python -Bm "$p.mixture" --decode "$d/mixture.json"
python -Bm pytest -q "$d/test_readout.py"
~~~
