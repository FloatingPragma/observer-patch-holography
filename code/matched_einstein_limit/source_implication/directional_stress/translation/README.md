# Local energy and translation balance across twelve cells

Both marked source heads lie in a connected twelve-cell region of the active 80-cell federation where a two-cell translation can be measured. Its signed generator has a flux through each of the fourteen boundary contacts. The marked heads supply a separately measurable force. At the four-tick stationary-source experiment the force and flux add to the measured rate of change of the translation observable. The native energy obeys its own continuity equation through exactly the same contacts. [Source and launch ensemble](../README.md) · [Complete chart, all cell rates and every directed edge](results.json)

## Registering the spatial comparison

The cells have three gateway ports apiece. Their distances inside the twelve-port cell distinguish inequivalent gateway configurations. A symmetry of an abstract honeycomb can therefore fail on the actual hardware. The search uses the planar cyclic order of the three incident macro-edges, the native internal icosahedral seams and the degree of every port. It maps each enrolled cell to a cell at graph distance two, without fixed cells, preserving the cyclic order of all three incident macro-edges. The shift across a pair fixes which three gateways of the first cell map to which gateways of the second. An internal port permutation must preserve all thirty internal seams and realize that correspondence. A consistent injection continues across a neighboring cell only when its gate assignments agree on the shared edge. The first compatible microscopic permutation in lexical order supplies the port comparison. The chart is fixed before reading a source response. [Chart construction](chart.py)

On the 80-cell carrier the shift sends the cells $(0,2,3,12,16,18,19,28,29,31,66,70)$ to $(19,1,18,14,17,31,29,30,71,69,64,68)$. Its center choice is cell 0 translated to 19. In this chart class the complete-ball search finds 24 oriented radius-one shifts on 80 cells and 96 radius-three shifts on 320 cells. The latter balls have nineteen cells. A deterministic connected extension witnesses up to twelve and 85 cells on those two boards, respectively. On the 20-cell board the specified equal-two-step shift does not extend over the incident edges of a single central cell. This census covers three finite boards. A displacement in native cell steps has no supplied proper length.

## The regional Ward identity

Let $S_u$ send the twelve ports of enrolled cell $u$ to their registered partners in $f(u)$ and vanish elsewhere. It acts as identity on the head factors. Set $S_R=\sum_{u\in R}S_u$ and $Q_R=(S_R-S_R^\dagger)/(2i)$. A plane wave on a homogeneous translation lattice would read a sine of its momentum in these units. On this finite apparatus $Q_R$ directly measures oriented coherences between the registered partners. It is additive under disjoint unions of enrolled cells, including when an image intersects the original domain.

Write $H=H_0+V$ up to an irrelevant scalar. $H_0$ is the actual mobile-sector Hamiltonian with the two marks replaced by blank contacts and the head factors left as spectators. $V$ is precisely the differing native action on the four marked ports. If $P$ selects the domain and $P'$ its image, the chart gives the exact intertwining
$$
P'H_0P'S_R=S_R P H_0 P.
$$
Subtracting the two restricted terms in the commutator therefore gives
$$
[H,S_R]=(I-P')H_0 S_R-S_R H_0(I-P)+[V,S_R].
$$
Taking the Hermitian part gives the real rate $i[H,Q_R]$. The first two terms record transport across the paired domain and image edges and the last gives $F_R=i[V,Q_R]$. Each edge read is local in the registered two-step comparison.

Specifically, on an oriented edge out of $u$, let $p,t$ be the near and far native ports, and $a,b$ the near and far ports of the corresponding image edge. The contribution to $\dot Q_u$ is
$$
J^Q_{v\to u}=\frac14\big(|a\rangle\langle t|+|t\rangle\langle a|
-|b\rangle\langle p|-|p\rangle\langle b|\big)\otimes I_{\rm heads}.
$$
Reversing an internal edge changes its sign. Its expectation is half the difference of two real pair coherences. Port-selective phases, blank swaps and interference read these quantities. The graph of blank routing seams with both marks removed is connected. Reversal parity and the two individual head Pauli measurements further resolve the force from the marked actions. Sampling a signed term with probability proportional to its coefficient and randomizing the sign on zero outcomes gives the recorded binary meter. The measured finite difference of $Q_R$ uses two independent uninterrupted flights bracketing the appointment. [Source operators and free-flight check](measure.py) · [Interferometer settings and sampling](counts.py)

For the energy meter each cell also owns its internal Hamiltonian, its degree-reference energy and half of each incident bridge. It uses the external bridge current defined by the independently calibrated contact pulses. On every enrolled cell the two simultaneous matrix identities are
$$
\dot E_u=\sum_vJ^E_{v\to u},\qquad
\dot Q_u=\sum_vJ^Q_{v\to u}+F_u .
$$
The proof is an operator calculation for all density matrices in each retained charge and parity sector. All thirty-six directed incidences, twenty-two internal incidences that cancel in pairs and fourteen boundary contacts are retained for each component. The unit coefficient of the force follows from the native change in the two marked actions. For the comparison Hamiltonian $H_0+\lambda V$ the same operator calculation gives the same boundary instrument and force $\lambda F_R$. A changed-law control takes the actual four-tick state, replaces only the two marked contact actions with the blank actions, and releases for two positive short durations. The observed shift-rate derivative then agrees with the boundary inflow alone within $10^{-8}$ for both species and both source arms.

At the registered four-tick appointment the two incoming species use exactly the tagged mixture of ten launches with all outside outcomes retained. The following outward bounds apply at the certified stationary source and use units of $10^{-6}$ per native time and per signed source variation. The ordinary and pairing sources each hold the registered radar volume fixed. [Fixed-point propagation and source enclosures](certificate.json)

| Species and source | Translation boundary inflow | Marked force | Free rate of $Q_R$ |
|---|---:|---:|---:|
| Vacancy, ordinary | [8.16, 8.30] | [−182.8, −180.1] | [−174.5, −171.9] |
| Vacancy, pairing | [2.42, 2.48] | [−53.46, −52.34] | [−50.99, −49.91] |
| Excess, ordinary | [−2.17, −2.11] | [204.9, 208.4] | [202.7, 206.3] |
| Excess, pairing | [−0.771, −0.754] | [61.07, 62.48] | [60.31, 61.71] |

At the nominal source the independently sampled ordinary boundary slopes are $8.243\pm0.106$ and $-2.149\pm0.106$ in the same units. The sampled marked-force slopes are $-181.436\pm0.528$ and $206.658\pm0.528$. The simultaneous 99% sampling bound covers both species, three meters and five preparations. Each setting draws $2^{76}$ simulated trials in binomial batches. On the same cut the ordinary energy-inflow derivatives are $+0.0294212$ and $-0.0199204$ in native energy per time per source variation. Independent differences of the free regional-energy readings agree to within $5\times10^{-7}$ on both ordinary variations. [Complete simulated counts](counts.json) · [Energy and shift readouts](results.json)

## What can absorb the force?

No redefinition of the directed momentum can absorb the force while retaining this measured boundary flux on all states. Let $h=2H+cI$ be the archived integral generator and $f=4F_R$. Multiplying the primitive sparse columns at the nonzero entries of $f$ gives, for each species and each reversal parity,
$$
\operatorname{Tr}(h^j f)=0\ (j=0,\ldots,6),\qquad \operatorname{Tr}(h^7 f)=64.
$$
The calculation uses integer arithmetic and gives $\operatorname{Tr}(H^7 F_R)=1/8$, independent of the scalar energy origin. For every operator $K$, cyclicity gives $\operatorname{Tr}(H^7 i[H,K])=0$. Hence $F_R$ cannot be the time derivative of any bounded observable of this fixed source sector, even a spatially nonlocal one. For example the strictly positive stationary state proportional to $I+(h/64)^7$ reads nonzero marked force exactly balanced by boundary transport. A balance for a larger material system could include support recoil, a change in the boundary flux or dynamics that translate the marks themselves. A coordinate momentum equation on curved spacetime can also have a connection term. Interpreting this measured force that way requires the connection from an independent metric on the same region. [Trace certificate, four sectors](bridge_limit.json) · [Direct column proof](bridge_limit.py)

There is a separate obstacle to directly interpreting the two finite conserved-current components as a symmetric spacetime stress. $Q_R$ connects pairs of cells at distance two. Every native bridge energy current acts only within a cell or between adjacent cells. Their Hilbert–Schmidt overlap therefore vanishes for **any** linear combination of bridge energy currents. Prepare the same stationary heads and, as a separate diagnostic probe, superpose port 0 and its translated partner with relative phase $+i$ or $-i$. The signed shift expectation changes by one. Every single native bridge energy-current expectation is identical. This excludes equality as operators on the full sector. It does not test an approximate relation limited to smooth states. [Exact support and two-mode witness](bridge_limit.json)

## A quantitative route to a coarse momentum

Choose the recorded blank-seam route from each port to its translated partner. The 144 routes use six to eight primitive hops. For wave amplitudes $z_0,\ldots,z_m$ along any one route, expand each amplitude as $z_0$ plus its consecutive differences. The imaginary coherence between the endpoints differs from the sum of the link coherences by
$$
-\sum_{i<j}\operatorname{Im}\langle z_{j+1}-z_j,z_{i+1}-z_i\rangle .
$$
Its absolute value is at most $(m-1)\sum_j\|z_{j+1}-z_j\|^2/2$. Purification gives the same inequality for mixed states. Summing over the routes yields a bound on the discrepancy between $Q_R$ and a registered sum of oriented local load currents. For each blank seam the squared-difference expectation is twice the antisymmetric-outcome probability, which gives an independent roughness test. At the nominal reference the actual discrepancies are −0.07950 and −0.01658 for a vacancy and excess, while the roughness bounds are 10.09 and 12.19. The prepared sharp-launch mixture does not provide a small roughness certificate.

Energy flux brings a second, independently measurable correction. Let $P$ denote mobile occupation in the twelve domain cells and let $H$ be the complete excitation Hamiltonian with energy referred to the parked background. The assigned regional energy obeys the exact formula
$$
E_R=\tfrac12\{H,P\}+\Delta_R,\qquad
\Delta_R=\frac{17}{4}\sum_{e\in\partial R}(P_{\rm outside(e)}-P_{\rm inside(e)}).
$$
Each crossing blank seam contributes diagonal contact energy $17/2$ and half is assigned to either side. With $J_N=i[H,P]$, differentiation and Cauchy–Schwarz give
$$
\left|\langle J_E-\dot\Delta_R\rangle-\langle H\rangle\langle J_N\rangle\right|
\le \sigma_H\,\sigma_{J_N}.
$$
At the reference the left side is 0.001880 and 0.007993 for the two species. The upper bounds are 0.06022 and 0.14520. Controlling wave roughness, energy spread and boundary storage in a matched regional limit is one concrete route toward identifying a common material momentum. These estimates have different current weightings. A symmetric tensor also requires compatible proper directions, lengths, time and a transfer of these estimates for the same physical regional weighting. [Exact operator identities and calculated error terms](bridge_limit.json)

## An exact weak-phase material transfer

There is an explicit preparation on which energy transport does select the same directed momentum to first order, without fitting a normalization. Retain the same stationary heads. In each native cell take the orbital amplitude equal at all twelve ports, and let $U$ denote the isometric embedding of these 80 cell amplitudes, including the four coin factors. At each unmarked **external** bridge $e$ the integrated energy and load currents obey
$$
U^\dagger J^E_e U=40\,U^\dagger J^N_e U .
$$
This is a matrix identity for arbitrary density matrices in the cell-uniform sector, including correlations between cell amplitudes and heads. An unmarked gateway has five internal blank seams and degree six. Acting on an internally uniform amplitude, each internal seam contributes exactly eight native energy units at the gateway and the degree reference vanishes there. Thus the incident internal action is $5\cdot8=40$ at both ends. The intra-cell part of the compressed bridge commutator vanishes and its intercell part is 40 times the load current. The integer calculation verifies all 119 unmarked bridges and rejects the same identity on the marked gateway $(1,38)$ in every charge/parity sector. All fourteen external edges of the enrolled region are unmarked, so the same exact relation holds for summed energy and load inflow of the whole region.

A second restriction recovers the longer translation read from these local currents. Put the same amplitude magnitude at every port and imprint $e^{i\alpha X_u}$ on each port in cell $u$, for any real profile $X$. Let $m_e$ count the net traversals of external edge $e$ in the 144 registered blank routes used above. None uses the marked gateway. Intra-cell phase differences vanish. The derivative of the end-to-end imaginary coherence is the sum of the consecutive phase differences, so at zero imprint
$$
\partial_\alpha\langle Q_R\rangle
=-\frac1{40}\sum_e m_e\,\partial_\alpha\langle J^E_e\rangle .
$$
The verifier checks the differential as a matrix on the heads for all eighty coordinate profiles. A constant profile induces zero response on both sides, leaving 79 independent directions. For $X_u=d(0,u)$ in the dual cell graph, eleven translated cells advance two distance levels and one stays on its level. The finite shift and signed routed energy divided by 40 are $-11\sin(2\alpha)/80$ and $-11\sin\alpha/40$, respectively. Both slopes are $-11/40$. Their finite ratio is $\cos\alpha$ away from its zero denominator. At $\alpha=1/64$ the relative discrepancy is below $1.23\times10^{-4}$. Separate simulated interference banks on the two mobile species reproduce the finite signals and distinguish their difference at $\alpha=1/16$. [Exact transfer, leakage check and count records](smooth.json) · [Operator proof and acquisition](smooth.py)

The phase-imprint direction differs from either registered matter direction. At fixed $\alpha$ both the shift and the routed energy measurements act as identity on the heads for this orbital preparation. Both signed matter derivatives vanish. More generally take any *common factorized* orbital density $\Omega$ with the same orbital state in both parity blocks and on both matter arms. Every unmarked bridge current acts as identity on both head coins and agrees between parities. On the marked bridge every position block has the exact form $L_b\otimes I+I\otimes R$, where $L_b$ can depend on left reversal parity $b$ and $R$ is the same in both parity blocks. The held unpairing preserves each parity-resolved left marginal and the parity-summed right marginal $\sigma_R$. The volume-matching rotation acts on the right about the Bloch axis of $\sigma_R$, so it fixes this same marginal and leaves all the left marginals invariant. Hence **every native bridge energy current has zero fixed-volume pairing derivative for a common factorized orbital density**. The marked current can respond to ordinary matter. At $\alpha=\pi/2$ its nominal ordinary slopes are $-4.79\times10^{-5}$ and $-3.03\times10^{-5}$, with pairing slope zero exactly on both charges.

The original source experiment evolves the probe together with the heads. At four ticks its nominal reference has only 0.08044 (vacancy) or 0.06337 (excess) probability in the cell-uniform orbital subspace. The uniform subspace itself immediately leaks under the active Hamiltonian. To see this, compare an unmarked gateway in cell $u$ to a distinct unmarked nongateway in that cell, with adjacent cell $v$. For an embedded vector $Uz$ the difference of the two $HUz$ outputs is $(\lambda z_u-z_v)/(2\sqrt{12})$, with $\lambda=33$ for vacancies and $-31$ for excesses. Choose distinct nongateways for all the unmarked incidences of a cell. The squared distance of two outputs from any common value is at least half the squared difference of the outputs. Summing over both ends of every unmarked edge, using $\|\lambda z_u-z_v\|^2+\|\lambda z_v-z_u\|^2\ge(|\lambda|-1)^2(\|z_u\|^2+\|z_v\|^2)$ and at least two such edges per cell, gives
$$
\|(I-UU^\dagger)HUz\|^2\;\ge\;\frac{(|\lambda|-1)^2}{48}\|z\|^2.
$$
The lower bounds are $64/3$ and $75/4$, respectively, on the initial squared leakage speed for unit input. No nonzero state remains in this flat-port subspace throughout an open time interval of native flight. These facts delimit the positive material transfer and explain why extending it to the nonzero collar response of the tagged four-tick source requires a different orbital sector or a dynamical error bound.

## The invariant transmitting sector

An exact cell-local sector that transmits along a bridge must retain its endpoint channel. In the orbital case take a reducing projection $P=\bigoplus_c P_c$. Between neighboring cells the Hamiltonian block is a nonzero multiple of $|g_c\rangle\langle g_d|$, with identity on the heads. The commutator vanishes only if the two gateway vectors are eigenvectors of the respective orbital projections with the same eigenvalue, zero or one. Admitting all three channels of a cell requires at least three independent orbital vectors. A single invariant orbital per cell cannot carry all three directions. For the full orbital and coin space, define the transmitting sector to be the smallest invariant block-local space containing all four coin channels at every gateway.

The exact construction takes the span of all powers of the native diagonal cell block $H_{cc}$ acting on its three gateway spaces. The orthogonal complement inside that cell avoids every gateway under all further internal evolution. Intercell hopping therefore acts only on the transmitting sum. The calculation gives 28 transmitting dimensions in each of 24 cells and 40 in each of 56 cells, hence 2912 of 3840 per charge and parity sector. The original tagged sharp-launch ensemble has exactly $4/5$ weight in this invariant sector for every prepared coin state. Each registered local shift intertwines its domain and image sector projections. All four marked port spaces are contained in the transmitting sector. The boundary shift flux uses gateway-to-gateway coherences, so both the boundary flux and the marked force are supported there. Conditioning the registered flux/rate measurements on this invariant sector multiplies them by $5/4$ for all five source settings.

Both obstructions survive removal of the dark modes. The two sign-controlled gateway preparations on ports 1 and 229 belong entirely to the transmitting sector. They give unit shift contrast and zero energy-current contrast on every bridge. The projected seventh force trace is $1/8$ as well. The dark Hamiltonian and all its powers are block diagonal in cell while the marked force has no same-cell block. The dark trace therefore vanishes. The positive weak-phase transfer has a precise preparation domain, while the nonzero collar force and these two obstructions persist on an exactly invariant sector that retains most of the native source. [Exact cyclic spaces, projections and weights](transmitting.py) · [Ranks and invariant source weight](transmitting.json)

From the repository root run `PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.directional_stress.translation.measure --write`. The modules `certificate`, `counts`, `bridge_limit`, `smooth` and `transmitting` accept the same invocation. The focused regression is `PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/directional_stress/translation`.
