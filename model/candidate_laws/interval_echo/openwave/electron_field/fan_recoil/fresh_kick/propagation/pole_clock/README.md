# A rest pole on the sharp 1045 vacancy

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The original sharp 1045 arrival contains a measurable component of the lowest compensated formation line of the commissioned 48-cell route. The line is extended. It does not give an autonomous rest-frequency oscillation in the same one-vacancy sector. A separate long-wavelength comparison asks how well independently measured spatial and phase records obey the quadratic mass-shell approximation on this fixture. [Common propagation fixture](../README.md) · [Registered pole protocol](PROTOCOL.md) · [Long-wavelength protocol](LONG_WAVE.md)

## The measured premise and implication

Use the negative clock-bypassed arrival at native wait 30.8125. Keep both face-0 words 177 and 146 through the gate and trace the retained observer reference in the read. Its gate has probability approximately 0.0158436995 per original start. The first source enters the commissioned route sharp at 1045 and has no addressed standing-wave preparation. A fresh capture is moved by one full native exchange to 1052 for the held control. In the latter operation the grade-five offset is compensated during the pulse as in the route preparation. It acts at its native value during every flight and filter dwell.

Separate compensated Ramsey counts clock the extended standing formation line through the neutral-ring reference. Increasing waits resolve the phase branch to **[1.75378873136, 1.75378876373]** in neutral-ring units. The filter reads the measured centre in appointment units. In one stage prepare a uniform coherent register of $M=32768$ appointments, wait $j/8$ when the register records $j$, give that arm the opposite centre phase, and keep its return to the uniform state. The accepted amplitude on a line of energy $E$ is

\[
F(E)=\frac{1}{M}\sum_{j=0}^{M-1}e^{-ij(E-\widehat E)/8}.
\]

The six-stage positive effect has response $|F(E)|^{12}$. The read requires supplied coherent phase control and a matched compensated reference. The source measurement keeps the two-word gate refusal separately from later filter misses.

The exact route recurrence puts every other energy at least 0.0016 above the lowest line. For any off-pole line, the denominator of the geometric sum has sine magnitude above 0.000099, giving response at most $7.362\times10^{-7}$. At the pole, $1-|F|^2$ is the mean of $1-\cos((j-l)(E-\widehat E)/8)$ over appointments $j,l$. Using $1-\cos x\le x^2/2$ and the variance of uniform appointments bounds six-stage attenuation by $1.046\times10^{-9}$.

If $a$ is the calibrated pole response, $w$ the pole weight, and $q$ the conditional acceptance, then $q=aw+q_{\rm tail}$ with $q_{\rm tail}$ between zero and the certified off-pole bound. The existing finite theorem gives an interval for $aw$ from the count and tail intervals. Dividing its upper endpoint by a lower bound on $a$ gives the independent weight enclosure. A failure to observe nonzero weight, or an unresolved contrast with the swapped source, would fail the registered test. [Lean positive-filter weight theorem](../../../../../../../../../Lean/SourceLaw/ChargedRegisterPole.lean#L50)

| Input record | Pole probability given the two-word gate |
| --- | --- |
| Sharp negative arrival at 1045 | [0.0089890916, 0.0090035914] |
| Same gate, one full exchange to 1052 | [0.0028309651, 0.0028437797] |
| Filled compensated comparator | [0, 0.0000049596] |

The ideal line has $E_*=10-2\sqrt{17}$. Its mode weights at successive ports repeat $(1,1-E_*/4,1)$ after normalizing over the 48 cells. The wider and narrower acceptance on the two sharp sources agree with these local weights. Replay uses $2^{42}$ model-generated original starts per bank and includes source uncertainty, numerical evolution allowance, tail leakage and shot error. [Filter and calibration receipt](records.json)

## What the free clock cannot read

The active signed one-vacancy address layer has a real integer tridiagonal generator $H$. The characteristic recurrence $p_0=1$, $p_1=x-d_1$, $p_j=(x-d_j)p_{j-1}-16p_{j-2}$ computes $p_{144}(x)$ exactly in any coefficient ring. Degree costs and endpoints are those of the enrolled route. Use the homomorphism $\mathbb Z[\sqrt{17}]\to\mathbb F_{151}$ with $\sqrt{17}\mapsto 35$. The monic polynomial and its translate by $E_*$ retain degree 144. Their modular Euclidean algorithm ends at degree zero. Hence their resultant over the original ring is nonzero and there is no energy pair with difference $E_*$. Exact rational Sturm pivots also count one eigenvalue below 1.7554. [Modular and Sturm verifier](exact.py)

For any fixed read $M$ on the free charged sector, a time-dependent outcome probability is a sum of terms $\rho_{ij}M_{ji}\,e^{-it(E_i-E_j)}$. The same formula covers arbitrary charged inputs and fixed local or joint readouts with the parked spectators. The measured formation frequency lies at least 0.00014451 in appointment units from every Bohr line, using the finite spectral error bound and the independent phase interval. The external vacuum comparison may read the formation phase. A fixed free internal detector on this isolated branch lacks that frequency. This obstruction is scoped to the declared sector and shutters.

The typed free Dirac construction uses a separate CAR sector and family mass block. Its generator and receiver must both intertwine with a proposed charge-source transfer before its physical spin and mass conclusions apply. The compensated pole receipt closes the finite positive-weight premise. A localized state selected by the full dynamics, a laboratory energy conversion and a second independent calibrated quantity are required for row 2. The isolated autonomous clock fails the row-3 frequency premise. [Source and receiver interfaces](../../../../lepton_family/theorem_interfaces/README.md) · [OpenWave criteria](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#simplest-test-per-criterion)

## The long-wavelength shell

Reflect the open guide across each missing end contact. The even subspace of this mathematical double has the original open generator, because the boundary-neighbour difference vanishes on a reflection-even pattern. The lowest bulk eigenvectors at $k=n\pi/48$ and $-k$ combine to give exact standing modes of the original route. A one-cell shift on the double, averaged with its inverse, restricts to a symmetric operator on the original route. Its measured expectation is the real three-contact quadrature plus half the reflection quadrature at the two end cells. The prepared standing modes give $\cos k$ for that expectation. The signed spectral frequency comes from a separate compensated Ramsey cohort. [Registered long-wavelength comparison](LONG_WAVE.md)

For native rate four and degree-five offset eight the bulk determinant is

\[
P(E,k) = (8-E)^2(16-E)-16\bigl(2(8-E)+(16-E)\bigr)-128\cos k.
\]

It factors at $k=0$ with roots $E_*,12,20-E_*$. Set $D=(12-E_*)(20-2E_*)>160$ and $A=32-3E_*<28$. With $\delta=E-E_*$ the lowest-root equation reads

\[
D\delta=128(1-\cos k)+A\delta^2-\delta^3.
\]

For $|k|\le0.27$, the determinant is positive below $E_*$ and decreases on $[E_*,E_*+0.1]$. At the upper endpoint its value is negative. Thus $0\le\delta\le0.1$. The equation and $1-\cos k\le k^2/2$ improve this to $\delta\le0.41k^2$. The Taylor remainder gives $|128(1-\cos k)-64k^2|\le(64/12)k^4$. Inserting the bounds on $D,A,\delta$ into the root equation gives

\[
\big|E(k)^2-E_*^2-(128E_*/D)k^2\big|<0.4k^4
\qquad (0<|k|\le0.27).
\]

Using the independent bare transfer and degree calibrations and the formation anchor gives $c_{\rm eff}^2$ in **[1.3283578, 1.3285077]** in route and neutral-ring units. Three held spatial and temporal comparisons meet the registered $0.5k^4$ tolerance with ample count resolution. For two held modes the residual is resolved as positive. Removing only the intercell contacts preserves each owner projector and changes the compensated half-wait phase at fixed preparation. [Held count records](long_wave.json)

| Approximate wave number from fringe read | Squared-frequency residual enclosure | Registered absolute tolerance |
| --- | --- | --- |
| 0.1309 | [−0.00003290, 0.00009877] | 0.00014596 |
| 0.1964 | [0.00009456, 0.00023318] | 0.00074132 |
| 0.2618 | [0.00045122, 0.00059819] | 0.00234530 |

The approximation has a finite one-dimensional domain and an explicit quartic error. The parent route's larger-wave-number measurements reject a common exact square-root speed. A physical relativistic limit requires an independent length and speed calibration and a corresponding refinement of charged transport. The scalar transfer theorem in the Lean library uses a positive stiffness and scalar time transfer operator. This route's charged generator is a separate first-order address dynamics. A source map would have to preserve the tested clock and wavelength records before applying that theorem. [Finite scalar transfer](../../../../../../../../../Lean/SourceLaw/TransferStructure.lean) · [Larger-wave-number comparison](../README.md#the-square-root-test)

## Replay

Run from the OPH repository root.

```sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_recoil.fresh_kick.propagation.pole_clock.measure
python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_recoil.fresh_kick.propagation.pole_clock.long_wave
```
