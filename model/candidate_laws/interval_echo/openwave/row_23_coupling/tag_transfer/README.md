# Predicting a second probe from one local measurement

An observer measures how a stored seam treats two different probes at its endpoints. That one measurement predicts how the second probe returns to a neighbor and to a shell two links away, on both a 20-cell board and an 80-cell board. All contacts are open during flight. The probes carry one or two units at one port. They have the same centered blank flight, which makes the comparison possible. The empty write, holds, uniform phase standards and coherent shell receiver are supplied. [Acquisition protocol](PROTOCOL.md)

At one frequency on the 20-cell board, the endpoint quadratures give a shift of **−24.00003 ± 0.00299** exchange-energy units. The non-scalar residual has operator norm 0.000095. At each held frequency, the remote forecast uses the single-tag record at that frequency and the one calibrated shift. The double-tag remote record is reserved for scoring. There are three complex entries per setting: the neighbor return, the coherent shell return and the transfer between them. The four frequencies share their time banks; each error bound covers those correlations.

| Condition | Held complex entries within forecast | Entries excluding an unchanged return |
| --- | ---: | ---: |
| Written seam, native contrast | 24 / 24 | 20 / 24 |
| Blank control | 24 / 24 | 0 / 24 |
| Written seam, half contrast | 24 / 24 | 20 / 24 |

Applying the full shift to the half-contrast control fails on 20 of its 24 entries. The distant shell and cross response at the lowest frequency on each board have insufficient resolution to distinguish the unchanged-return foil. The native fit uses no remote double-tag records. The comparison uses simulated independent cohorts of 2^36 per phase setting. [Raw counts](counts.npz) · [Calibration, bounds and every held outcome](results.json)

## Why one number suffices on this material

Let $P$ restrict the board to the two seam ports and $P_-$ project onto their antisymmetric combination. Removing the uniform phase $6\kappa q^2$ gives the same blank generator $L/2$ for either indivisible load $q$. On the written seam the local addition is $V_q=-\kappa q^2 I-P_-$. The simulation checks the corresponding columns against the registered local generator. Thus $H_2-H_1=P^T\delta I P$ with $\delta=-3\kappa$. For the full single-tag resolvent $G_1=(z-H_1)^{-1}$, write $B=PG_1P^T$, $X=RG_1P^T$ and $D=RG_1R^T$, where $R$ collects the two remote profiles. The finite resolvent identity gives

$$RG_2R^T=D+X(I-\delta B)^{-1}\delta X^T.$$

The two endpoint inverses measure $\delta I=(PG_1P^T)^{-1}-(PG_2P^T)^{-1}$. Real symmetric flight allows the entries of $G$ to be reconstructed from complex return measurements on individual and summed profiles. A separate sparse full-board solve checks every inferred matrix entry against its counting and integration allowance. The joint 99% budget uses Hoeffding bounds over independent shots at each quadrature, a six-point quadrature remainder on the centered spectral interval $[-33,6]$, and the omitted time tail. The implementation also allows for probability rounding before sampling. [Acquisition and inference](measure.py) · [Independent checks](test_transfer.py)

There is a second identity for a spatial cutoff. If $g_>$ is the omitted blank propagation at the seam, the effective interaction $U_q$ satisfies $U_q^{-1}=V_q^{-1}-g_>$. Hence $U_2^{-1}-U_1^{-1}=V_2^{-1}-V_1^{-1}$ at every common cutoff. For native contrast the even and odd entries of this difference are $3/32$ and $8/99$. The effective interactions can each change with resolution while this difference stays fixed. These two source tags and the material impurity have no independent gauge-charge identification. A running gauge coupling and a measured onset at its core radius remain open.

From the repository root, run:

```sh
python -m model.candidate_laws.interval_echo.openwave.row_23_coupling.tag_transfer.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_23_coupling/tag_transfer
```
