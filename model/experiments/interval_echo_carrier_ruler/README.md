# A material comparison length for the thirty-seam grating

A two-port standard defines the unit of length $L$. We declare the common chart edge of the twelve-port grating to have length $a=3L/2$. Two receivers on different owners infer $a/L$ from sign counts on the grating and the standard. The declaration supplies the numerical input required by PR-73. The counts check how the declared scale enters the preparation and the thirty-edge energy symbol. The unit's conversion to laboratory length requires an external calibration. [Prospective protocol](PROTOCOL.md) (frozen at commit 27c992335)

| Preparation | Expected scale in rod units | Result from both receivers and both gradients |
|---|---:|---|
| Declared grating | 1.5 | [1.4916694, 1.5081135] |
| Dilated grating | 2 | [1.9902556, 2.0090145] |
| Common phase gain reduced to 4/5 | 1.5 | [1.4856599, 1.5132132] |
| Folded grating | 8.5 | Principal-branch readings are incompatible. The separate wave hulls contain 1.5 and 5.5. |
| Receiver exchange shuttered | Gradient erased | The two reference phases are unresolved. Each apparent length ratio lies near 1. |

## What the receivers measure

A single double load occupies an otherwise empty board. Eleven native split pulses distribute its amplitude over all twelve ports of one owner. A common gradient sets its phase at port $u$ to $k a \,x_u\!\cdot n$, where the coordinates $x_u$ have unit edge length and $n$ is the selected seam direction. A separate rod uses ports 12 and 13 with phases 0 and $kL$. To set a phase, the controller transports the addressed amplitude to degree-five port 205, waits under the native private load-square term and reverses the transport word. A full blank exchange swaps the two double-load amplitudes and fixes vacuum. The routes consequently restore their intermediate ports, including occupied spectators. The chart, source calendar, phase gradients and shutters are inputs to the instrument. [Source and two-head layout](apparatus.py)

Fresh copies carry each pair to either ports (216,217) on owner 18 or (228,233) on owner 19. The receiver uses its private quarter-turn and a half exchange to separate symmetric and antisymmetric exits. Its blank reference neighbor serves the two sequential native write-and-inverse pointer reads. On all remaining addresses the load contributes to the off-pair outcome. For a unit-amplitude stable sign standard, this pointer reproduces the corresponding sharp exit. Each archived panel retains all three bins. The 1,100 independent cohorts of $2^{24}$ trials are seeded multinomial simulations of native-law probabilities. They represent 18,454,937,600 synthetic draws. [Count archive, routes and exposure words](records.json.gz)

For an equal-amplitude pair with phase difference $z$, the conditional antisymmetric probability is $p=(1-\cos z)/2$. On the enrolled branch $0\le z\le\pi$, the fitted length ratio is the recovered sample phase $\arccos(1-2p_s)$ divided by the recovered reference phase $\arccos(1-2p_r)$. This cancels a common positive phase gain. The decoder sees only the counts and order of the two registered gradients $kL=\pi/5,\pi/7$. At each gradient it evaluates every one of the thirty edges at both receivers. All sixty individual intervals contribute to the joint intersection. The reference must resolve the ordering of the two gradients at both receivers before reporting a joint interval. Reference failure and an empty scale intersection are recorded separately. We invert simultaneous Hoeffding intervals with total error at most 1%, budgeting 4,096 readings. Conditioning on the observed in-pair sample size gives a conditional binomial law. The in-pair probability is 1/6 or 1 on this enrolled cohort. The $64/2^{28}$ conditional allowance covers rounding of the two absolute exit probabilities. The unconditional allowance is $2/2^{28}$. [Inference and strict count replay](measure.py)

## From the comparison length to the coefficients

Write $w_e=x_v-x_u$ for the unit direction of edge $e=(u,v)$. The antisymmetric probability on a complete, uniformly illuminated grating is $(1-\cos(ka\,n\cdot w_e))/12$ per seam. The held experiment sums these unconditional readings at a separate fixed $kL=8/3$. In native exchange-gap units this gives $(5/12)\Lambda_1(ka n)=(5/12)a^2\Lambda_a(k n)$, with $\Lambda_1(z)=(1/5)\sum_e(1-\cos(z\cdot w_e))$. The exact equal-edge symbol fixes the entire Taylor packet. The finite directional counts test that symbol on held settings. [Underlying seam instrument](../../candidate_laws/interval_echo/openwave/row_29_maxwell/seam_symbol/README.md) · [Committed seam-moment identity](../../../Lean/Screen/SeamCurrentEdge30Moment.lean)

| Held direction | Energy interval at $a=3L/2$ | Energy interval at $a=2L$ |
|---|---:|---:|
| Port | [2.8818336, 2.9134637] | [2.5864604, 2.6170322] |
| Edge | [2.9817817, 3.0169372] | [2.9557532, 2.9905706] |
| Face | [3.0001583, 3.0342847] | [3.0128946, 3.0466894] |
| Oblique | [2.9263603, 2.9644896] | [2.7649820, 2.8025269] |

At the nominal scale all four exact seam-symbol values lie inside the recorded intervals. The face-minus-port interval is [0.0866946, 0.1524511]. A kernel on the twelve unit vertex directions with the same leading isotropic normalization gives −0.3882000. The dilated grating at the unchanged comparison wavenumber also resolves all four predicted means. [Independent route, pulse, count and kernel checks](test_ruler.py)

At $a=3L/2$ the edge coefficient slots are $C_4=-9L^2/80$, $B_0=27L^4/4480$, and $B_6=-9L^4/22400$. The quartic slot distinguishes the two positive scales. The three ratios $B_0/C_4^2=10/21$, $B_6/C_4^2=-2/63$, and $B_6/B_0=-1/15$ hold at every positive scale. Rescaling *both* lengths and inversely rescaling the phase gradient leaves both intensity distributions fixed. The declared reference length therefore remains essential. The checked theorem supplies the principal-branch inversion, common-gain cancellation, first-gradient alias, unit freedom, declared coefficients and dilation on the existing edge packet. [Scale.lean](Scale.lean) · [Scale and slot source](../../../Lean/Screen/DispersionArmingInterface.lean)

The supported PR-73 claim is the declared positive scale of this mounted finite scalar grating. The physical length of the reference is an independent input. The sector identification, common frame, cofinal attachment and positive exclusivity bound for a physical FZ-12 dispersion comparison remain independent. [PR-73 source and external-input qualification](../../../code/electromagnetism/fz12_arming/PREMISE_ROW_PROPOSALS.md)

Run the tests and replay from the repository root. Check the Lean file from the Lean project after building its DispersionArmingInterface import.

```sh
python -m pytest -q model/experiments/interval_echo_carrier_ruler/test_ruler.py
python -m model.experiments.interval_echo_carrier_ruler.measure
(cd Lean && lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_carrier_ruler/Scale.lean)
```
