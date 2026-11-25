# Which interaction does a spatially limited observer measure?

A one-load probe on registered Interval Echo scatters from a single written seam. Two observers who keep different spatial wavelengths need different effective interaction matrices to predict the return of the same probe at the finer resolution. The change is detectable before either observer resolves the bound load around the seam. The effect establishes running of a material impurity interaction in a specified cutoff prescription. Its measured onset misses the core-scale test for OW-23. The OpenWave M10 cell is untested at the cited revision. [OpenWave criterion and decisions](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#simplest-test-per-criterion) · [Experiment protocol](PROTOCOL.md) · [Cutoff follow-up protocol](CUTOFF_PROTOCOL.md)

## A core and an independent probe

Write the native echo once at an empty emitter seam, then release a load with all federation seams active. An independent 256-position timing instrument tags the two bound lines below the blank one-load band with finite resolution. Locating the tagged load gives its root mean square distance from the written seam's midpoint. The tag accepts 97.96% of two endpoint launch banks on the 20-cell board and 97.97% on the 80-cell board. The radii are **0.53762 and 0.53768 links**. Conservative simultaneous 99% upper bounds are 0.53890 and 0.54043. The tagged state includes the finite filter's leakage. The scatterer in the second measurement is the same written seam. The separate bound-load measurement defines its operational core radius. The protocol supplies the write and the probing load. [Radius bins, all tag failures and measured spectra](results.json)

A second observer compares a freely evolving load with a stored coherent reference. Two choices of the recombination phase give the real and imaginary return amplitude from each of three ports and each equal two-port superposition. The third port is a neighbor of both written endpoints. The real symmetric Hamiltonian makes these six preparations sufficient to recover the three-port retarded response matrix from its complex diagonal readings in the six probe states. Controlled hold durations, preparation, blank uniform-mode phase calibration and spatial recombination are supplied metrology. With $z=E+i\eta$, the observer integrates each return against $-i\exp(izt)$ to obtain the resolvent; the same hold records supply every energy. The counts in [the first Ramsey bank](results.npz) are the numbers of plus outcomes out of $2^{36}$ independent simulated enrollments per quadrature and hold. Their axes are real/imaginary, hold index and source. The complementary count equals the enrollments minus the plus count. The eight sources are the three ports, their three lexicographically ordered pairs and the two filtered receiver impulses. [Source and readout](experiment.py)

## What changes with scale

The blank, centered one-load Hamiltonian is the graph Laplacian divided by two, $H_0=L/2$. On the two written endpoints, let $P$ restrict to those ports and $|{-}\rangle=(|a\rangle-|b\rangle)/\sqrt2$. The blank-to-written change of the registered generator is $P^T V P$, with $V=-8I-|{-}\rangle\langle{-}|$. This follows because the empty written contact annihilates its vacuum and one-load states. These formulas provide an audit. The interaction used to predict the held port comes from the independent endpoint observations. For the measured blank and written resolvents $G_0$ and $G$, put $b=PG_0P^T$ and $w=PGP^T$. The local interaction is measured by $b^{-1}-w^{-1}$; the scattering vertex, including repeat visits, is $b^{-1}(w-b)b^{-1}$. Removing the external blank propagation in the latter expression corrects for the source and detector geometry. Turning off the contrast while keeping the written contact leaves only the odd-channel cut. The inferred even scattering is consistent with zero in this control at every tested shell. [Endpoint and held-return counts](results.json)

If an observer retains only blank spatial modes with graph wave number $k=\sqrt{\lambda(L)}$ at or below a cutoff, omitted modes can carry the load away and back. Let $g_>$ and $g_<$ denote the endpoint blank resolvent from the omitted and retained modes. The interaction in the restricted space is $V_{\rm eff}=(I-Vg_>)^{-1}V$. Iterating it on the retained modes gives $(I-V_{\rm eff}g_<)^{-1}V_{\rm eff}$. The identity with $(I-V(g_<+g_>))^{-1}V$ explains the cutoff dependence. The physical test is a fresh return measurement. Prepare the third-port impulse projected into a spatial passband by the surveyed blank graph and then normalized. Infer $V$ once on the endpoints at the lowest frequency. At $z=-1+i/2$ the coarse and fine prescriptions use cutoffs 0.75 and 2.25 per link.

| Cells | Retained modes, coarse / fine | Error on fine return using coarse coefficient | Error using matched coefficient | 99% bound for matched error |
| --- | ---: | ---: | ---: | ---: |
| 20 | 20 / 80 | 0.097339 | 0.0000028 | 0.000286 |
| 80 | 80 / 320 | 0.097375 | 0.0000037 | 0.000286 |

The coarse prediction's own error allowance is below 0.000199 at the fine cutoff. The comparison uses the identical fine source, detector, linewidth and complex target for both forecasts. Thus repeating scattering through the eliminated modes changes a measurable effective coefficient, beyond a change of probe overlap. The eigenvalues, source weights and omitted propagator come from a supplied graph survey. The untruncated two-endpoint interaction also predicts the held third-port returns across all six probed spatial shells at both widths, including the cases in which a single encounter gives a large error. These are local finite-frequency response measurements. A scale dependence of the microscopic seam strength or a gauge interaction cannot be inferred from their variation. [Programs and counts](results.json)

## Does the change begin at the core?

At half-width 0.5 the magnitude of the even, externally corrected scattering vertex first changes by 10% at observed $kR=0.894$ on either board. Large blank spectral gaps precede those readings. At half-width 0.25 the first such sampled change moves to 0.299 and 0.330. This scan changes probe frequency and has limited resolving power for an onset.

A targeted follow-up holds the energy at $-1+i\eta$ and examines the **same cutoff-dependent coefficient** that predicts the held return. It enrolls all 103 and 399 distinct spatial shells on the two boards. The first resolved 10% changes of its even-channel magnitude relative to the uniform-mode cutoff have the brackets below. Both endpoints of each bracket are measured-mode cutoffs. Intermediate classifications, including unresolved shells, remain in the records. The post-crossing return was acquired with fresh third-port probes selected by the endpoint calibration alone.

| Cells | Half-width | Bracket in $kR$ | Post-crossing error with frozen coefficient | Error after matching |
| --- | ---: | ---: | ---: | ---: |
| 20 | 0.50 | 0.11348–0.18381 | 0.002802 | 0.0000022 |
| 20 | 0.25 | 0.11347–0.18381 | 0.010912 | 0.0000014 |
| 80 | 0.50 | 0.12462–0.13450 | 0.000618 | 0.0000015 |
| 80 | 0.25 | 0.12439–0.13450 | 0.001431 | 0.0000036 |

The frozen-coefficient errors exceed their individual simultaneous 99% allowances (0.0001021–0.0002261). All four matched errors lie below their own allowances. Even the upper confidence bound on the core radius puts each first certified crossing below $kR=0.185$. The declared core window starts at $kR=0.5$. A physical gauge charge, a dynamical charged scatterer, and an onset near its measured size require separate evidence. OpenWave's M5 example credits a classical core form factor in part; M10's Wilson-dressing rolloff used a selected coefficient. The present held tests resolve the effect of eliminating modes on a written material center, with a measured failure of the core-onset assignment. [All shell ratios and crossing decisions](cutoffs.json) · [Independent held counts](cutoffs.npz) · [Follow-up program](cutoff.py)

## Numerical checks and replay

Each complex read uses a six-node Gaussian rule on half-exposure intervals up to 80 with an explicit exponentially damped tail. The native one-load spectrum centered on the blank board lies in $[-9,6]$. Integration errors follow from the derivative bound of the exponential on that interval. Independent binomial banks and Hoeffding's inequality give simultaneous component errors, propagated through the two-by-two inverses and held predictions. Separate exact binomial upper bounds on the radial tails cover the core-radius check. The bank sizes quantify simulation sampling. [Independent full-board solve and record checks](test_coupling.py)

From the repository root with NumPy and SciPy installed, run the following commands.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_23_coupling.experiment
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_23_coupling.cutoff
python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_23_coupling
```

A two-load probe on the same seam provides an independent test of the effective interaction. One endpoint calibration on the 20-cell board predicts all 24 held native remote complex entries across two board sizes. Copying the single-load return fails on 20 of them. [Two-tag transfer, raw records and cutoff identity](tag_transfer/README.md)

An early/late contact on a separate freely moving pair isolates the effect of two-body flight on a guided scattering signal. The initial core radius and the flight delay can be measured separately. [Timed two-body scattering and all momentum bins](timed_scattering/README.md)
