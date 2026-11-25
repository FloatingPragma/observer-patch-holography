# How light can a moving unit charge be?

A negative material charge crosses an owner boundary with a formation doublet at **0.484390 and 0.485307** above matched filling. That puts its lowest level a factor of 33 below the 16-unit cycle standard. A sharp vacancy launched at one end reaches the other end in 1,046,751 of 1,048,576 simulated starts. Closing the central contact gives zero opposite-tip arrivals. The splitting sets the appointment, and the full four-level dynamics predicts the terminal distribution. The exact energies come from the local generator. The transfer counts distinguish open and closed paths. Repeating the same construction with positive load gives a lowest energy of 16.484390. This double well supports tunneling between owners and has no bulk momentum measurement. [Protocol](PROTOCOL.md) · [Native experiment](measure.py) · [All counts](results.json)

A separate source prepares each of the four energy eigenstates of either sign. An adjacent partial exchange sets each magnitude, and routing that coordinate to the ordinary private reference and back sets its phase. After seven units of free evolution, a position reader estimates the diagonal cost. On three fresh banks, a quarter phase and half exchange give the odd component of each blank bond. Adding those four sample means measures the full cost above filling in the same head state. At simultaneous 99% coverage over eight lines, the vacancy doublet gives [0.484248, 0.484545] and [0.485129, 0.485427]. The lowest excess gives [16.484225, 16.484522]. The measured separation resolves the slow beat independently of the transfer appointment. [Native pulses and detector](meter.py) · [Eight energy readings](meter.json)

For a one-defect family, the earlier matched-host theorem bounds the heavy formation energy by 198 on every finite board. The observed tau/electron ratio would therefore require a positive light energy below 0.057 in these units. The following lower bounds apply to whole classes of preparations, including heterogeneous boards. The zero of energy always uses the same filled material and head states. [Upper bound and full local blocks](../marked_bands/README.md)

## Every passive transporting network

In a fixed head sector, write the vacancy generator on a connected component of blank, exchanging seams as
$$
H=\tfrac12 L+\mathrm{diag}(V),\qquad V_i=8d_i-8r_i+w_i.
$$
Here $d_i$ is the blank degree, $r_i$ is one for an ordinary private reference and zero otherwise, and $w_i$ sums incident written-head costs. On the single-defect sector magnitude zero costs zero. The magnitude-two relative term acts only on its own direction coin at either occupied endpoint. Its two coin eigenvalues are zero and 32 or 33 according to contrast parity. All these coin operators commute. Decomposing into their joint eigenspaces covers arbitrary superpositions and entangled physical heads of the stated magnitudes. The ordinary private reference for an excess contributes +24 in place of −8.

**Gap theorem.** On any blank connected component with at least two ports, every positive vacancy energy is at least $1/3$. Every excess energy is at least 8. A zero vacancy mode occurs only on an isolated two-port blank dimer with both private references ordinary and both written costs zero. Singletons have no transporting seam.

In a fixed vacancy block, $V_i=0$ is possible only at a degree-one tip with an ordinary private reference and $w_i=0$. Call these tips soft. Every other port has $V_i\ge8$. Outside the two-soft-tip dimer, each soft tip $i$ meets a hard neighbor $j$. Its edge gives
$$
\tfrac12|x_i-x_j|^2\ \ge\ |x_i|^2/3-|x_j|^2.
$$
The difference is $|x_i-3x_j|^2/6$. Each hard port pays for at most six tips and keeps at least two units of its own diagonal per unit probability. Sum the inequalities and drop the other nonnegative edge terms. On the exceptional dimer the two eigenvalues are 0 and 1. For an excess, every incident blank seam supplies eight units of diagonal, proving the second bound directly. This proof holds at every size. With a common positive electron/positron formation energy on passive transporting components and a one-defect heavy state, the largest possible ratio is $198/8=24.75$. Even using the vacancy alone gives at most 594. [Lean inequalities](../../../../../../Lean/SourceLaw/ChargedRegisterPole.lean)

## Live side heads on a ring

Keep the 144-port blank ring. At each of the 48 middle ports allow a magnitude-one pendant in either contrast parity. The outer port has no other live physical contacts. The backbone and pendant private references may independently be ordinary or zero-written. With one signed defect, the relative active two-state matrix of a pendant is
$$
-a\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad a=4\text{ or }4.5 .
$$
The two uncoupled coin/position columns have zero relative cost. Crossing an active side flips that coin. Returning restores it. Thus each component reachable from the backbone has a scalar backbone and one additional site for each active side. Other configurations include isolated pendant charges. The code checks both physical signs against all local columns. The argument also covers omitted, dark or passive magnitude-two sides, whose contributions on the backbone are nonnegative.

Partition the propagating component into the backbone, positive pendant diagonal sites and negative pendant diagonal sites. For a vacancy all active pendants have negative diagonal. The backbone costs at least $16-8-4.5=3.5$, using two blank neighbors and at most one live side per site. For an excess the backbone costs at least $16-4.5=11.5$. Any positive active pendant has diagonal at least $24-4.5=19.5$. Removing its coupling by the quadratic form inequality $2|uv|\le |u|^2+|v|^2$ costs at most 4.5 at either endpoint. The combined backbone and positive pendant block therefore costs at least 7 for an excess. For an eigenvalue $\lambda>0$, eliminating the negative diagonal block $D$ adds the positive operator $K(\lambda-D)^{-1}K^*$ to that combined block. A positive eigenvector must have nonzero projection there. This proves positive floors 3.5 and 7 respectively even for a nonperiodic arrangement. The same proof applies to a longer ring with the stated incidence.

As a numerical cross-check the homogeneous active rings give 1,536 exact-block eigenvalues per sign over 48 physical momenta, both head parities and four private-reference settings. The smallest positive energies are 5.574592 for a vacancy and 10.234637 for an excess. These are spectral calculations conditional on the prepared heads. The four-port experiment above uses equal-load writes and sharp launches with the native law. Active transporting seams on loops, multiple defects and comparisons that change the host sector are outside these bounds.

On the separately prepared fifteen-port all-live ring the signed formation operator is negative semidefinite for both charges. Replacing one contact by a blank or opening one side gives positive costs. Exact characteristic-polynomial counts exclude positive roots through 0.45 in every tested one-blank setting. [Active backbone and independent meter checks](ACTIVE_BACKBONE.md)

From the repository root run the following commands with one numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.light_charge_gap.measure
OPENBLAS_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.light_charge_gap.meter
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/light_charge_gap
```
