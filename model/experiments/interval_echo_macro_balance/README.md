# Directional matter readings at coarse resolution

Twenty groups of observers compare twelve directional readings after a smooth load travels freely across the connected Interval Echo screen. Nine directions predict the other three if all twelve are contractions of one symmetric spacetime tensor in the specified cell chart. The withheld error drops in absolute units as the groups get larger. At 1,280 and 5,120 cells it is 42.5% and 73.3% of the directional variation. Allowing the best quadratic fit to *all* twelve directions leaves 8.0% and 14.1% unexplained on those two screens. [Source and measurement design](PROTOCOL.md) · [Measured results](results.json) · [Recorded clicks](counts.npz)

| Cells per group | Smooth-source kinetic energy | Directional RMS in nats | Withheld RMS in nats | Best twelve-direction residual / directional RMS | Largest sampled withheld difference in standard errors |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.05357 | 0.02224 | 0.10773 | 93.7% | 172.4 |
| 4 | 0.01150 | 0.00808 | 0.01144 | 29.8% | 19.6 |
| 16 | 0.00268 | 0.00641 | 0.00482 | 14.2% | 11.2 |
| 64 (holdout) | 0.000647 | 0.00292 | 0.00124 | 8.0% | 3.47 |
| 256 (holdout) | 0.000159 | 0.00105 | 0.000771 | 14.1% | 2.44 |

The RMS variation subtracts each group's angular mean. The table gives exact finite detector probabilities for the first four numerical columns and fresh multinomial samples for the last. Independent samples of antisymmetric seam outcomes give energies 0.05359, 0.01154 and 0.00268. The unreleased smooth source has zero angular variation at every resolution. On the 320-cell screen halving the kinetic exchange changes the twelve-direction residual to 20.5%. For a sharp load released from a single port it is 49.7%. At 5,120 cells the smooth half-exchange release gives 23.9%. All three residuals use the same chart, reference and detectors. There are six antipodal pairs. A quadratic on the sphere spans all six even readings but only three of the six odd readings. The irreducible angular discrepancy therefore has opposite signs at opposite ports.

The two holdouts preserve the instruments, grouping, 5% fractional target and acquisition budget. At 256 descendants per region, angular variation falls faster than the absolute prediction error. Both additional finite profiles miss the fractional target. [1,280-cell measurements](results_holdout_3.json) · [5,120-cell measurements](results_holdout_4.json)

An independent precision survey on 5,120 cells uses 64 times as many preparations per setting, selected after the original survey. Directly summing its fresh clicks gives an angular RMS of 0.0010463 nats and a withheld error of 0.00079695 nats, or 76.2% of the measured signal. The largest withheld difference is 15.00 propagated standard errors for native release, 2.39 for the unreleased input and 15.03 for half exchange. The RMS figures use the raw multinomial estimates and retain sampling noise. [Precision survey](results_holdout_4_precision.json) · [Independent count histories](counts_holdout_4_precision.npz)

The sources begin with one compiled coherent load, and all links act during flight. A detector measures the vacuum-centred negative logarithm of a six-port reference state, using sharp and balanced-pair settings. Each cell has the same supplied rotation-averaged reference. Its twelve labels and their alignment across cells belong to the apparatus. The cell observers sum their recorded outcomes by parent. The three withheld directions are tested on separately prepared loads without conditioning on arrival. The common angular chart would require physical transport between regional spacetime frames for a curvature comparison. [Local reference and angular calibration](../interval_echo_angular/README.md)

## A bound from the measured kinetic energy

In the one-load sector the kinetic energy above the uniform reference is $E=\langle L/2\rangle$. For the local projector $P$ onto the uniform mode and first triplet, write $F=I-P$. Every remaining local mode has internal kinetic energy at least 3. All intercell kinetic terms are positive, so the total weight $\sum_c\|F\psi_c\|^2$ is at most $E/3$ for the entire flight.

For each withheld direction let $R$ be its star operator minus the forecast from the nine calibration operators. Rotational invariance of the declared reference gives $PRP=0$. The largest $\|PRF\|$ across the three directions is 15.275 nats and the largest $\|FRF\|$ is 11.286 nats, rounded upwards. In each group choose the direction with the largest absolute discrepancy. Cauchy–Schwarz across the disjoint groups then bounds the sum of these discrepancies by $2(15.275)\sqrt{E/3}+11.286 E/3$. The bound uses the unconditioned, per-load normalization of the detectors.

For the specified envelope $f_c=\max(0,1-hd_c/2)$, adjacent cells differ by at most $h/2$ and the launch parent contains a twentieth of all cells with unit envelope. On any continuation with three contacts per cell and equal parent populations,
$E=\sum_{cd\ {\rm contact}}|f_c-f_d|^2/(24\sum_c f_c^2)\le5h^2/16$.
The absolute angular forecast error therefore tends to zero at least linearly in $h$ in these units. The same argument with the constant-mode projector and the first cell gap $(5-\sqrt 5)/2$ bounds the absolute angular signal itself by $O(h)$, since the constant mode has identical expectation at every pole. The measured fractional errors test whether coarse-graining preserves a predictive angular response above this common suppression. Stress in physical units, its covariant conservation and source-induced curvature require independent measurements on the same domains.

Run from the repository root to reproduce the release and detector counts.

The two large boards can be regenerated with the public OPH geometry builder. The deposit reader verified byte-for-byte SHA-256 agreement for both generated files. With archived counts the four tests below reconstruct every withheld discrepancy directly. [Public builder at the verified revision](https://github.com/muellerberndt/oph-physics-sim/blob/476aff7219fe05592563cb0a1281aba4cffc9d96/oph_defects/fedgeom.py) · [Carrier identities](../../executable_model/refinement.py)

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_macro_balance.measure --output /tmp/interval-echo-macro-balance
git clone --no-checkout https://github.com/muellerberndt/oph-physics-sim.git /tmp/oph-carrier-builder
git -C /tmp/oph-carrier-builder checkout --detach 476aff7219fe05592563cb0a1281aba4cffc9d96
export DEPOSITED_BOARDS=/tmp/interval-echo-carriers
PYTHONPATH=/tmp/oph-carrier-builder .venv/bin/python -m model.experiments.interval_echo_macro_balance.rebuild_carriers --output "$DEPOSITED_BOARDS"
for level in 3 4; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_macro_balance.measure --holdout-level "$level" --holdout-board "$DEPOSITED_BOARDS/FEDERATION_LEVEL$level.npz" --output /tmp/interval-echo-macro-balance
done
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_macro_balance.measure --holdout-level 4 --holdout-board "$DEPOSITED_BOARDS/FEDERATION_LEVEL4.npz" --precision --output /tmp/interval-echo-macro-balance
.venv/bin/python -m pytest -q model/experiments/interval_echo_macro_balance
~~~
