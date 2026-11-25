# Keeping the material probe label changes the entropy maximum

On the freely formed 80-cell collar, toss a fair coin to send a vacancy or excess from port 1 to 38 and back, with one exposure unit per leg. Retain both the probe sign and its balanced comparison bit. The entropy of these four outcomes plus a positive coefficient times the independently surveyed area has a strict local maximum on a selected surface of constant survey volume, in all three admissible preparation coordinates. This is a population calculation on the finite apparatus. Its measured volume belongs to the four-station tetrahedral survey. [Raw source and clock records](results.json) · [Search and calibration](measure.py)

## The measurement and its ordinary calibration

Each repetition draws either the independently paired or the freely formed seam coins. A sharp read of the reversal sign with its result discarded populates the odd sector with probability one half. Let $u$ be the overall odd probability, $s$ the freely formed fraction and $v$ the probability of both. The four mixture weights are $1-u-s+v,s-v,u-v,v$. The conditional odd probabilities $v/s$ and $(u-v)/(1-s)$ are between zero and one half. These are the physical preparation bounds for the three-coordinate experiment.

For each probe sign $j$, compare completion with an independent even, independently paired reference trip and break a tie with a fair sign. Its mean is $m_j=p_j(u,s,v)-p_j(0,0,0)$, an affine observable of the four preparations. For the retained fair setting label and reply, the observed entropy is $S=\log2+[H_2((1+m_-)/2)+H_2((1+m_+)/2)]/2$. Here $H_2$ uses natural logarithms. Each individual binary entropy contributes a negative rank-one Hessian along its measured mean gradient. Their average supplies curvature in two directions when the gradients are independent. The area and the volume constraint have their own curvature.

A separate excess-probe four-station radar uses half-unit ticks and a common horizon of six. All missed trips are included. The directed ranges come from fresh copies of the same preparation. The survey fits four station clock gains, symmetrizes the corrected distances and evaluates the total triangular area $A$ and tetrahedral volume $V$. For each ordinary or pairing change, varying $v$ alone can hold the measured $V$ fixed locally. Ordinary calibration fixes $\eta=-(\partial_u S)_V/(\partial_u A)_V$. With $u=1/32$, selecting the independent conditional reference $v=us$ at $s=0.04129976557$ makes $(\partial_s S)_V+\eta(\partial_s A)_V$ vanish as well. Each parity-sign dephasing probability is then $1/16$. The paired condition participated in selecting the reference. [Geometry and source preparation](../intrinsic_area.py)

| Fixed-volume direction | Entropy derivative | Area derivative |
|---|---:|---:|
| Ordinary at fixed $s$ | $+1.29888\times10^{-6}$ | $-0.01520830$ |
| Pairing at fixed $u$ | $-1.31690\times10^{-7}$ | $+0.00154194$ |

The ordinary gain is $8.54058\times10^{-5}$ nats per survey-area unit. A separate 50-digit evaluation of the recorded corner populations gives constrained Hessian eigenvalues $-1.69022\times10^{-6}$ and $-1.39371\times10^{-9}$ for $S+\eta A$. It uses the Hessian of $S+\eta A-\lambda V$ projected onto the tangent plane, including the curvature of the constant-volume surface. Moving both ways by $10^{-4}$ along each principal direction and solving for the correlation probability to match the volume lowers the objective in all four cases. The changes are about $-8.45\times10^{-15}$ and $-6.97\times10^{-18}$. These are numerical differences on the stored population records. No finite-shot resolution is claimed. [Separate derivative and finite-displacement checks](verify.py)

## Search and held changes

The search uses 18 radar schedules (two tick intervals, three horizons and three probe-sign policies), five conditional odd rates in each pairing arm, and eleven initial pairing fractions from 0.001 to 0.999. It brackets gradient-determinant sign changes and refines each bracket. For single balanced bits, the two charges, four launch and four relay stations, and three leg times give 96 settings including same-station inspections. The 532 bracketed roots include 447 with positive ordinary gain and zero with two negative constrained eigenvalues. The count includes different preparations and duplicated reverse routes. A finite scan cannot exclude roots between samples or outside the scanned family.

For a retained fair label, each unordered pair of the six central-bridge charge/time choices defines a measured four-outcome record. Of 200 bracketed roots, 60 had a positive gain and two negative constrained eigenvalues. The displayed case has the most negative larger eigenvalue in this scan. Both the reference and the reported menu were selected on these data.

Holding $\eta$ fixed and erasing the sign label gives ordinary and pairing transfer fractions of 0.779 and 0.948, where cancellation requires one in both channels. Halving the formed fraction gives 1.816 and 1.462. Switching the independent radar to quarter-unit ticks gives 2.278 and −1.266. This establishes the strict maximum for the selected measured entropy and surveyed source family. Extending it to the native regional entropy, further admitted source directions, or the area of proper shrinking balls requires additional source identifications.

From the repository root:

```sh
PYTHONPATH=code:. python3 -m matched_einstein_limit.source_implication.stable_station.measure
python3 -m pytest -q code/matched_einstein_limit/source_implication/stable_station
```
