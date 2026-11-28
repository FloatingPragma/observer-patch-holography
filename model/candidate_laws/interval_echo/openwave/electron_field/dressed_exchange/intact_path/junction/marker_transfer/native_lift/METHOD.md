# Pullback test for the central rotation phase

## Observable theorem premise

The [committed port cover][cover] maps the central element $-I\in SL(2,\mathbf F_5)$ to the identity spatial port action. A preparation can carry that lift if observers independently establish how the same physical rotations act on its coherent arrival readout. The [golden-sector representation][gold] factors through the port group, so its central action is $+I$. The cover itself allows either central response.

A more direct experiment compares the two orders of perpendicular spatial half-turns. Their spatial destinations agree. For the guarded vacancy orbital in one supplied flux quantum, [`MagneticVacancy.order_phase`][magnetic] checks a difference of ten twentieths of a phase cycle on every sharp port. Given faithful lifts on the same code and the same final exterior records, this forces opposite branch amplitudes and a dark outcome on every admitted input in a balanced phase comparator. Independently multiplying either turn by a constant phase leaves this order test unchanged. The orbital's supplied link phase is part of that premise. A sign under the exchange of two independently transported charges would require its own experiment. [Guarded orbital and controls][flux]

On the completed 1045 cell, an isolated term preserves $G_v=h_v-q\mathbf1_{\{v=x\}}$ at each port $v$. Here $h_v$ counts the internal heads pointing at $v$, $q$ is the signed load and $x$ is its location. A proper physical rotation of all charge and field words permutes the entire measured $G$ vector. The exact census leaves only the axial order-five subgroup. For instance, the proper half-turn that interchanges physical addresses 1044 and 1045 changes the six counted records at 1046, 1047 and 1052 through 1055. The paired negative order theorem thus has an unfulfilled physical premise on this source. A change of the external environment during a turn would need its own measured charge-and-field return.

## Coherent schedules

Each enrolled isolated term acts as $r_j(I-P_j)$, where $P_j$ transposes the two admissible sharp configurations and fixes inactive configurations. On the junction there are twelve material contacts and two native triangle terms. In the completed cell there are twenty triangle terms and, with mobile material, thirty material contacts. The proof uses the written-field Laplacian phase convention, including its zero on a blocked term. An additional address-dependent onsite energy or a change to the sign of an off-diagonal matrix element changes the premise.

Let $\Omega_N$ be the normalized vector which assigns the same amplitude to every one of the $N$ sharp configurations in a connected native orbit. The proof uses $\Omega_N$ to calculate a lower bound. It is not an outcome supplied to either observer. Each $P_j$ fixes $\Omega_N$. A sum of generators with any time-dependent exposure therefore annihilates it. The same is true when shutters disconnect the orbit during the schedule. The scheme includes adiabatic control of the unique ground line of a connected positive weighted native Laplacian. That line is $\mathbb C\Omega_N$ at every appointment, giving zero energy and zero Berry phase in the constant gauge. Excited sectors orthogonal to $\Omega_N$ require separate analysis. [Checked permutation and phase implications](NativeRecordPhase.lean)

For balanced arms $A,B$ acting on the same incoming records, the predicted bright amplitude is $w=(A+B)\psi/2$. Its squared norm, measured inclusively over charge-and-face accept and refuse, is bounded by
$$\Pr(\mathrm{bright})=\|w\|^2\ge |\langle\Omega_N,\psi\rangle|^2.$$
This follows by taking the $\Omega_N$ component and applying Cauchy–Schwarz. The bound is additive over native configuration orbits and over untouched references. A normalized symmetric sum of $m$ codewords in an orbit has weight $m/N$. One sharp input has weight $1/N$. For the archived coherently accepted arrival, the measured first-marker marginal gives $m_0+m_X$ per original launch. The independent second marker has marginal $I/2$. For both junction signs and both single-charge variants the inclusive floor therefore equals $(m_0+m_X)/N$. The source uncertainty propagates through this same linear combination. [Arrival preparation](../PROTOCOL.md) · [Orbit and floor receipt](results.json)

| Physical aperture | Configurations in one orbit | Bright floor on coherent code sum | Conservative inclusive floor per original launch |
| --- | ---: | ---: | ---: |
| 1045–1028 junction | 120 | $1/30$ | $0.00019480283$ |
| Whole 1045 cell, fixed hub charge | 6610 | $1/3305$ | $0.00000353542$ |
| Whole 1045 cell, mobile charge | 48820 | $1/24410$ | $0.00000047772$ |

Postselection can divert this bright weight into rejection. A separate bound applies when a unit-return reference is compared to an exact scalar negative amplitude $-a$, $a\ge0$, on the entire code. Let $P$ project onto the coherently accepted code. Its normalized equal sum $\psi$ satisfies $P\Omega_N=\sqrt p\,\psi$, with $p=m/N$. Swapping the two code labels fixes $\psi$. If $PA\psi=-a\psi$, the refused amplitude $d=(I-P)A\psi$ must carry $\langle\Omega_N,d\rangle=(1+a)\sqrt p$. The inequality $(1+a)^2p\le(1-p)(1-a^2)$ follows from $\|(I-P)\Omega_N\|^2=1-p$ and $\|d\|^2=1-a^2$. It implies $a\le1-2p$. On the junction an exact negative exchange limited to accepted events can therefore pass at most $196/225$ on this coherent input. The bound does not predict whether any schedule reaches it.

## Half-pulse computation

The path oscillator sets half the exposure needed for each isolated full exchange. The amplitude update is $[(1-i)I+(1+i)P_j]/2$. Integer real and imaginary numerators keep the acceptance/refusal effects exact across the two 37-step schedules. For each outcome the four-by-four matrix in the receipt acts on the original two-marker input, so any additional source-reference measurement can be computed on the same copies. The right-hand probability columns in the README exchange $A+B$ with $A-B$. That plant establishes phase sensitivity when refusals are counted. A single antisymmetric marker input or an additional phase service can produce a dark record as part of its calibration. Establishing intrinsic statistics requires the scalar minus comparison across the full spin code on an intact transport.

The positive lower bound excludes an exact full-code minus sign from this apparatus, even if the pulses often produce two separated physical charges and even if the two paths are programmed to reach the same final marker registration. It does not bound arbitrary selected bands, phase-dependent controls or a wider carrier with an independently measured route for the shared environment. The row and the physical electron assignment remain open on this domain.

[cover]: ../../../../../../../../../../Lean/Screen/SL2F5PortCover.lean#L76
[gold]: ../../../../../../../../../../Lean/Screen/GoldenSectorPSL2F5Representations.lean#L178
[magnetic]: ../../../../../../../../../../Lean/SourceLaw/MagneticVacancy.lean#L19
[flux]: ../../../../../../../../../experiments/interval_echo_angular/charged_orientation/flux_doublet/README.md
