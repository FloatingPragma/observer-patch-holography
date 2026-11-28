# Angular reads of the written four-face arrival

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

The fixed 14→24 receiver can distinguish positive and negative interference on a charge-and-doublet arrival. Addressed field operations give both signs without moving the receiving frame. The same fixture has a nonzero orbital response to an applied flux and measurable transfer between material and ring currents during flight. The local current probes read zero on every superposition of the captured field code. The magnetic, current and interference results use the original four-face preparations for both signs of charge. [Protocol](PROTOCOL.md) · [Records](results.json)

## Interference and the physical receiver

The first receiver outer-ring projector is the native generator when all other terms are shuttered. A duration-π pulse changes the sign of the first odd orbital and preserves the second. Independently supplied π phases on two outer rim head bits interchange the orbitals. The two complete actions have square $+I$. Their orders are opposite on all eight receiving code columns, including all four unobserved sender words. This verifies the phase operator before any ensemble is evolved. The balanced reference arm measures the sum and difference of the amplitudes. The analyzer is rerun with independently supplied $+I$ and $-I$ controls.

| Receiver at t = 27.0625, sender unobserved | Bright | Dark | Joint detector refusal |
| --- | ---: | ---: | ---: |
| Twice the native ring pulse versus hold | 20.14017% | 0 | 79.85983% |
| Native ring and supplied head pulses in opposite orders | 0 | 20.14017% | 79.85983% |

At t = 30.8125 the source-vacuum record can be required as well. The same comparisons give 13.15899% in the selected output and 86.84101% refusal for the equally mixed input. The receipt also retains the two input parities and both circular superpositions. Every launch from the ideal supplied source enters the denominator.

For a rigid physical turn of the receiving icosahedron, the observers hold the gateway aperture and its nine written seam addresses fixed. Of its 60 proper poses, 55 move the gateway and four preserve the gateway but move one of the nine seams to another address. Only the identity preserves the joint detector domain. The negative order sign uses field controls at a fixed pose. A spinorial spatial turn in this fixture would require a continuous physical action on the arriving charge and field whose intermediate poses and final read stay in a calibrated domain. [Geometry and phase operators](measure.py)

## Applied flux and native currents

On all eight physical triangles the observers can calibrate an applied phase against three enabled material hops and the closing ring flip. On each closed tour a phase applied to the material hop alone contributes $q\phi$ in the agreed sense, and the phase applied to the ring flip alone contributes $-q\phi$. Reversing the charge and complementing the written heads interchanges the positive and negative applied-flux arrival records. For example, for a positive load and equal mixture on the first receiving face (24, 31, 25), the local arrival at t = 27.0625 is 20.14017% with no phase, 5.87085% with a positive quarter-turn material phase and 5.59115% with a negative one. These phases are instrument settings and the face response is an output.

For a selected rim bond write $b_i\in\{0,1\}$ for the recorded head bit in configuration $i$. The two possible phase exponents on a transition sum to $b_j-b_i$. Applying them together is the exact diagonal conjugation $D H D^\dagger$, with $D_{ii}=e^{-i\phi b_i}$. The material-only Hamiltonian is likewise conjugate to the field-only Hamiltonian with the opposite applied phase. As a result their stationary energy records agree for all flux, although their fixed preparation and detector responses can differ. The replay verifies the conjugation using a rephased preparation, whose output rephases by the same diagonal action.

At zero applied flux both separate ground-energy slopes vanish. Complex conjugation reverses each flux and the ground is simple. The positive ground curvatures across all eight faces are **0.00600888–0.09788674** native energy per radian squared. The derivative from stationary perturbation theory agrees with independent finite probes at ±1/1024 radians. Positivity also follows from the connected graph. After dividing amplitudes by the strictly positive zero-flux ground vector, the ground-energy quadratic variation minimizes a sum of positive edge weights times squared differences between the applied phase and a gradient. The independently calibrated closed tour has nonzero period, which no gradient can cancel. A physical magnetic moment and gyromagnetic factor would require an independently established spin branch and conversion from applied flux to physical field.

At zero applied phase the time derivatives separate into material current and field ring current. Their sum is exactly $i[H,b]$. For the positive load with circular-plus input, the first receiving face gives integrated material current **+0.342827** and field current **−0.334731** up to t = 27.0625. Their sum **+0.00809620** agrees with the separately read net change in bond head occupation. The projections of both current operators vanish exactly on the prepared code for every source bond and on the captured code for every receiver bond, including arbitrary sender words. The pulse can therefore change phase records, and loops can contribute circulation during transfer, while these local probes give zero mean current after capture. A mechanical angular momentum measurement would require a physically coupled torque standard.

These outputs apply to the two-cell written-field fixture. The [separated-channel audit](../binding_audit/README.md) gives a separate one-level restriction on transporting the complete doublet below the parked threshold. The present experiments use simulated observers and ideal phase and shutter controls.

Extending the two-cell source to the symmetric receiver at port 116 gives an exact restriction on which field direction the receiver can record at any exposure. [Three- and four-fan terminal channel audit](TERMINAL_CHANNELS.md)

On the asymmetric routes of the 320-cell board, a finite coherent clock can preserve both admitted energy amplitudes and retain entanglement in the local receiver. [Finite coherent threshold measurement](COHERENT_FILTER.md)

The two individually resolved energies restrict a fixed-guide turn comparison. Both addressed phase controls change the energy records when applied to the entire admitted band. [Energy-preserving order implication](ENERGY_ORDER.md)

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/fan_rotation/measure.py
```
