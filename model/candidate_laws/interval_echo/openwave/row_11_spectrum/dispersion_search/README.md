# Momentum-first search for a light mobile charge

On registered Interval Echo, a vacancy and its conjugate carry the same resolved band through a prepared material guide. Both a cosine dispersion and a square root using the measured discrete translation momentum cover **17 of 17** held mode, length and sign comparisons after calibration on two vacancy modes. A quadratic dispersion and a relativistic square root cover **1 of 17** each. Their one success is the repeat of the calibration momentum at the opposite charge. The momentum here is the magnitude of the translation phase measured over one guide period.

## Selecting the charge

Use the 144-port route of the 320-cell federation. Terminate an active path after either 96 or 144 ports. In each even period $P$, alternate blank and even unit-contrast contacts and leave the final two contacts blank. The viable periods in this family that divide both lengths, retain at least six periods, and contain at least one live mark, are 4, 6, 8, 12 and 16. The source contains one missing unit or its conjugate excess against a host of ones; a remote stock compensates the charge. Private costs are erased and exterior contacts are closed. Each mark flips its direction as the charge crosses. Both conjugated native column tables give the identical finite generator after subtracting the host energy.

An addressed source compiles the lowest $N/P$ standing modes of each short guide. The position detector requires at least one quarter of the probability between the two end periods. Independent local energy counters give the lowest **interior** formation energies 0.484396, 0.234986, 0.135452, 0.061016 and 0.034425. Their simultaneous intervals separate all five. Period 16 wins this screen independently of physical lepton masses. On that guide two lower lines near 0.03419676 put only 0.08% of their charge in the interior. Counting them as freely propagating would lower the apparent threshold. Numerical modes supply the synthesis targets; the momentum, energy and transfer records are measured in independent simulated banks.

| Physical path | Resolved interior modes | Measured phase range, radians per motif | Lowest interior formation energy, negative-charge simultaneous 99% interval |
|---|---:|---:|---:|
| 96 ports | 4 | 0.635 to 2.520 | [0.03442265, 0.03443036] |
| 144 ports | 7 | 0.396 to 2.752 | [0.03442081, 0.03442850] |

A spatial interferometer reads position and the real coherence between ports one period apart. A separate energy reader counts gradients and material positions. Coherent free waiting against another mode of the same charge resolves recoil frequencies to a simultaneous radius of $4.44\times10^{-10}$ native energy units. The apparatus uses positive-time native nearest-neighbor controls, including the phase gates needed for emission and detection. Higher counts reduce statistical error in this simulated experiment. They require $2^{38}$ launches per absolute-energy setting, $2^{26}$ per phase or momentum setting and $2^{34}$ per current setting.

## Predictions on the held records

The training data are the lowest two interior modes of the negative 96-port source. Hold its fitted offset and speed fixed. The fitted relativistic speed is 0.01011 consecutive route-port spacings per native time. Each second length uses its independently read lowest interior line as the recoil reference. On the highest 144-port vacancy mode the observed recoil is $2.54596\times10^{-5}$. The fitted relativistic curve predicts $4.29551\times10^{-5}$, with simultaneous interval $[4.27094,4.32030]\times10^{-5}$. The cosine prediction is approximately $2.546\times10^{-5}$. Replacing the squared guide phase in the relativistic square root by the measured one-period difference $4\sin^2(\theta/2)=2(1-\cos\theta)$ also passes all 17 tests without a further fitted parameter. The largest residuals are $9.62\times10^{-9}$ for the cosine and $1.52\times10^{-8}$ for the discrete-momentum square root. The spatial-counter and calibration errors span their difference. Including both charge signs, eight different held momenta per separate sign fit follow the cosine curve. The fixed-vacancy training fit also predicts the opposite charge.

These data resolve a material transmission band. The shortest accessible phase is about 0.396 radians per motif. A long-wavelength continuum limit below that value, another excitation, a three-dimensional physical momentum calibration and a source of mass units require additional evidence.

## Two ways to emit and detect the same excitation

One emitter synthesizes a coherent pair of low-band modes and lets it fly into bare local load counters at two different owners. The reciprocal emitter prepares a sharp charge independently at either of those ports. The narrow analyzer on the return experiment inverts the pair preparation and reads the load at its output. The calibrated recoil frequency predicts both sets of held count fringes from their zero-time and incoherent-source records. The incoherent mixture fails the preselected second-wait comparison at each receiver for both signs and both lengths.

For example, on 144 ports after 524,288 native time units, the negative charge reaches port 1832 on 89,600 of 67,108,864 narrow-source launches. The sharp-source reciprocal analyzer records 90,271. In a fresh bank the material-position herald followed by an independent native reverse hop on the arriving direction head records 89,393. At port 3686 the same three readings are 89,889, 89,380 and 89,945. A direct cut-current splitter gives $+0.935\times10^{-6}$ in the positive coordinate direction for one phased vacancy source. Counting the increase of material probability beyond the same cut gives $(+0.933\pm0.041)\times10^{-6}$. The charge conjugate has the opposite signed current for the same phase preparation.

At any one sharp material port, the traversed marks uniquely determine the direction word in this prepared component. The registered reverse-hop reader checks that head on arrival. Its rank there is one, so this source supplies no two-state internal spin at a sharp charge address. The two-site or wider standing waves are spatial degrees of freedom. These measurements identify the charged band and its dispersion in the prepared material. They supply no electron assignment.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.dispersion_search.measure --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/dispersion_search
```

[Finite source and detector derivation](PROOF.md) · [Count banks and predictions](records.json) · [Native apparatus](apparatus.py) · [Measurement](measure.py)
