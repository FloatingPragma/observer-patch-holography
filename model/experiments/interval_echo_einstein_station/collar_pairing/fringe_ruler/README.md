# A two-emitter phase ruler on the formed collar

Four observers recover the relative phase of waves sent from two bare ports on the 80-cell board. Strengthening the correlation in two parked head coins changes this phase, even though both head marginals and the parked energy remain fixed. The four observers read the same preparation with an excess and a vacancy. If the two phases measured a common change of path difference with one wavenumber per probe, the ratio of responses would be the same on every baseline. The first baseline requires a ratio between **0.104 and 0.127**. The second requires **0.274 to 0.319**. These simultaneous 99% intervals include the measured shutter uncertainty. [Counts and inference](results.json) · [Material-clock counts](timing.json)

## What is measured

A freely moving excess encounters the marked seams at 0–1 and 38–39 on a unit-filled federation. After four native exposure units discard its position and retain both heads. Their joint state is $\rho$. Independently prepared encounters provide the two factors of $\tau=\rho_L\otimes\rho_R$. A private source coin selects $C_s=(1-s)\tau+s\rho$ with $s=0,1/2,1$. Both end preparations have even reversal and parked energy 16. An independent ordinary standard dephases the left reversal of $\tau$ with equal even and odd weights, supplying one half unit of exchange energy. The arriving probe receives no indication of source choice. [Original source and ordinary standard](../ranging/README.md)

Mount the emitter at the first unmarked edge in address order of the first source-owner cell whose ends both have degree five. This selects ports 2 and 4. In the otherwise unit-filled board, prepare an excess or a vacancy on port 2. An isolated half exchange puts it on the two rails with amplitudes proportional to $(1,i)$. A private exposure on rail 4 shifts its energy relative to rail 2 by $8[(1+q)^2-1]$, where $q=+1$ or $-1$. The known sign of this contrast fixes the orientation of the phase for each probe. The program prepares $(|2\rangle+e^{i\phi}|4\rangle)/\sqrt2$. Both marked heads are parked during preparation. Release the entire board for four exposures, with no route, relay, intermediate measurement, or receiver-side phase standard. Each trial ends in one of the four original marked ports, or a miss. The clock, address selection, private exposures and terminal monitoring are supplied operations. [Native preparation and acquisition](measure.py)

Let $A_y,B_y$ be the four-by-four coin blocks propagating the first and second launch sites to receiver $y$. The count probability is

$$p_y(\phi)=D_y+\operatorname{Re}(e^{i\phi}F_y),\qquad
 F_y=\operatorname{Tr}(B_y C_s A_y^\dagger).$$

Four exposures at $0,\pi/2,\pi,3\pi/2$ read both quadratures. This determines the fringe directly from intensities. All 80 held probabilities at $\pm\pi/4$ agree with their count predictions. A mixture with equal source weights independently reproduces the midpoint of the complex endpoint fringes at every probe and receiver. Erasing the coherence between the two launch rails removes the phase modulation in every control. An arbitrary common launch phase cancels from the measured baseline phase $\Phi_{0y}=\arg F_y-\arg F_0$. A scalar phase of the free generator or an arbitrary phase convention at a receiving port also cancels. Since $F_y$ is affine in $s$, the measured endpoint disks exclude a zero on each interpolating segment: the phase difference along that segment has a unique continuous lift. [Count-only reconstruction](readout.py)

| Baseline | Vacancy phase change, radians | Excess phase change, radians |
|---|---:|---:|
| 0 to 1 | −0.009971 to −0.009860 | −0.001243 to −0.001043 |
| 0 to 38 | +0.1534 to +0.1546 | +0.0423 to +0.0489 |
| 0 to 39 | +0.747 to +0.760 | −0.0018 to +0.0034 |

The phases compare $s=1$ with $s=0$. The ordinary standard gives an independent finite calibration on each baseline. Dividing the correlation response by twice the observed ordinary phase response gives the equivalent energy in units of the ordinary exchange gap. On 0–1 the vacancy gives −0.0855 to −0.0843 and the excess gives −0.0184 to −0.0153. On 0–38 they give +2.29 to +2.42 and +0.087 to +0.104. This calibration is a finite ratio against the half-unit standard.

## Can the source change the wavelength?

Allow the source to change both the three path differences $L_b(S)$ and each probe's scalar wavenumber $k_q(S)$. Let $\beta_{qb}$ be an arbitrary static phase offset for each material and baseline. The model $\Phi_{qb}(S)=\beta_{qb}+k_q(S)L_b(S)$ requires differentiability only along the two source preparations used below. Any source-dependent phase common to all four receivers has already canceled.

At the **same fully formed stock**, interpolate toward the product stock or toward the ordinary half-dephased product stock. Both are convex preparations of the complete head state, including reversal parity. Label their three-baseline phase derivatives $g_q$ and $o_q$, respectively. The quadratures from the existing independent cohorts determine both derivatives because every receiver fringe is affine in the prepared density. In particular, its phase derivative toward stock $D$ is $\operatorname{Im}(F_D/F_{\rm formed})$. The two one-sided directions suffice, even if no opposite preparation is physically available.

For $c=k_+/k_-$, differentiating the proposed phase law cancels each unknown length derivative. Both $g_+-cg_-$ and $o_+-co_-$ must be multiples of the **same** vector of three reference lengths. Their cross product must vanish. Expanding in $c$ gives the three columns

$$M_0=g_+\times o_+,\quad
M_1=-g_-\times o_+-g_+\times o_-,\quad
M_2=g_-\times o_-.$$

The count-only 99% interval for $\det[M_0,M_1,M_2]$ is **$[3.9470,6.2744]\times10^{-5}$**, including the stable-clock uncertainty. The independent integer-propagator enclosure is **$[5.0549772193,5.0549772603]\times10^{-5}$**. Since the columns are independent, no $c$ can make $M_0+cM_1+c^2M_2=0$. If the vacancy wavenumber vanishes at the reference, its two response directions themselves must be collinear, contradicted by $M_2\ne0$. Thus these two source changes also exclude source-dependent scalar dispersion in this common-length phase law. A control with varying wavenumbers and genuinely shared lengths gives zero determinant. [Two-direction test](dispersion.py) · [Coefficient and determinant intervals](certificate.json)

## Local theorem and precision

The absence of a single calibration also holds at each endpoint as an infinitesimal statement. Wherever $F_y\ne0$, its phase has derivative $\operatorname{Im}(F'_y/F_y)$. The independent integer Chebyshev propagator encloses the complex four-tick amplitudes from both senders to all four ports, in both reversal sectors and for both probe charges. Its four-by-four operator error is below $10^{-19}$. The formed and product density errors in trace norm are below $10^{-18}$ and $2\times10^{-18}$ respectively. Quantizing each real and imaginary density entry to the nearest $2^{-80}$ changes its trace norm by less than $10^{-22}$. Contracting two blocks of norm at most one gives an error below $3\times10^{-18}$ per real component of $D$ and $F$. At the formed source the phase slope ratios on 0–1 and 0–38 lie within (0.115, 0.117) and (0.308, 0.310). The product reference also gives disjoint ratios. All 300 archived terminal-bin confidence intervals enclose the independently bounded native probabilities. [Integer replay](certificate.py) · [Derivative enclosures](certificate.json)

Each acquisition uses $2^{48}$ synthetic starts per setting. For each of the 300 terminal bins and the single-tick clock, the reader uses a two-tail Bernstein bound with $x=\log(80000)/N$. For observed frequency $h$ the radius $\sqrt{2hx}+4x$ follows by solving the Bernstein inequality for the unknown probability. A separate dimer with the same stable commands tests the flight tick and the $\pi/96$ preparation tick. Up to $10^{12}$ coherent firings resolve the phase alias in the final clock bank. These are simulated counts on identical fresh copies; the long trains assume the hand remains coherent and each repeated command has the same duration at the hand and apparatus. The largest additional error in any terminal probability is $2.55\times10^{-13}$. The per-clock-bank Hoeffding failure bound is $2e^{-16}$ for the twelve precision banks. Together with the 301 probability bounds the union failure probability is below 1%. [Timing and error propagation](timing.py)

In detail, the integer generator has norm below 32 after a common scalar energy shift. A four-tick timing error $d$ changes a flight effect by at most $64d$. It changes the formed state by at most $64d$ in trace norm, and its product of marginals by at most $128d$. Each beam splitter uses 48 preparation ticks with traceless norm $1/2$. Each private pulse uses at most 21 such ticks with relative gap at most 24. The total per-probability allowance is therefore $192d+2(48/2+21\times24)\delta$, for single preparation tick uncertainty $\delta$. The receiver choices and fixed source/flight commands define this apparatus.

This measures a transported phase and its source response. A proper path length also requires an independently characterized propagation law. Probe-dependent diffraction or a refractive response that varies by path can change these interferometric phases. The selection of a common proper metric and the OW-27 light-propagation criterion remain open.

From the repository root, use one thread and run

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
p=model.experiments.interval_echo_einstein_station.collar_pairing.fringe_ruler
nice -n 15 python -m "$p.measure"
nice -n 15 python -m "$p.timing"
nice -n 15 python -m "$p.certificate"
nice -n 15 python -m pytest -q model/experiments/interval_echo_einstein_station/collar_pairing/fringe_ruler
```
