# A flight-conserved collar charge with a signed local current

The same collar logarithm that responds to unpairing at the selected 80-cell equilibrium has an exactly conserved extension during native excess and vacancy flight. The extension measures spatial interference on six unmarked ports. It has the specified initial head charge at every one of those sharp launches. A local continuity equation gives a nonzero current inside the cell and zero total outflow from its six-port support. An apparatus using isolated native blank-seam pulses reads its expectation. [Integer proof and replay](measure.py) · [Pulse circuit](interferometer.py) · [Measured population records](results.json)

A second conserved extension transports the same initial charge through the source-cell boundary. The sign can be chosen independently of the initial calibration. [Intercell construction and limits](TRANSPORT.md)

## The exact extension

Freeze the faithful head reference $\chi_*$, including the left reversal block. Let $\pi=\chi_{*L}\otimes\chi_{*R}$ and $K=\log\pi-\log\chi_*$. The parked charge of the reference is $-D(\chi_*\Vert\pi)<0$. Erasing its pairing while keeping both marginals fixed has charge derivative $D(\chi_*\Vert\pi)+D(\pi\Vert\chi_*)>0$. The volume-adjusting right rotation contributes zero to this derivative. The stationary-source certificate encloses it in $(0.231079,0.231133)$ nats. [Reference and signed pairing certificate](../full_algebra/certificate.json) · [Parked charge and sharp-event obstruction](../remote_heads/README.md)

Use the first cell of the same 960-port board. Two integer position vectors supported entirely inside the cell are
$$
v_1=e_4-e_5-e_6+e_7,\qquad v_2=-e_5-e_6+e_8+e_{11}.
$$
Both have zero amplitude on the external gateways and marked ends $0,1,3,10$. To determine the largest cell subspace that stays zero there during internal propagation, let $A$ be the twelve-port internal adjacency and let $B$ select those four ports. The null space of the stack $BA^j$, $0\le j<12$, is exactly the span of these vectors, by rational elimination and Cayley–Hamilton. On this span $A=-I$. Independently, applying each **full** native integer flight generator to $v_i\otimes e_a$ for every coin column $a$ gives eigenvalue $-3$ for a vacancy and $29$ for an excess in both reversal sectors. These integers use twice the flight Hamiltonian with the sector's scalar energy zero removed. The invariant-subspace statement and commutators are independent of those scalar conventions. [Generator and scale](../full_algebra/native_precision.py)

Let $V$ have columns $v_1,v_2$. Their Gram matrix has diagonal 4 and off-diagonal 2. Thus $\Pi=V(V^TV)^{-1}V^T$ is an orthogonal projector with diagonal $1/3$ on $S=\{4,5,6,7,8,11\}$ and zero elsewhere. In each reversal sector put
$$
C_b=3\Pi\otimes K_b.
$$
Every full Hamiltonian commutes with $C_b$ because its range lies in the displayed eigenspace and its orthogonal complement is invariant. At a sharp port $p\in S$, $\langle p|C|p\rangle=K$ on the heads. For $p\notin S$ the compression is zero. These are operator statements for every head state. In particular the same signed change $\delta\langle K\rangle=-\eta\,\delta A|_V$ at the selected equilibrium survives free flight at every time for a launch in $S$.

For $P_p=|p\rangle\langle p|\otimes I$, $\|[P_p,C]\|=\sqrt{2}\|K\|$ on $S$. This noncommutation resolves the earlier obstruction for an extension commuting with the sharp event. A projective measurement of the dark subspace would succeed with probability $1/3$ on a sharp launch in $S$. Assigning $3K$ on success and zero otherwise gives second moment $3\langle K^2\rangle$. The bound $\|C\|=3\|K\|$ is optimal among dressings supported in this dark subspace because compression on a launch has norm at most $\Pi_{pp}\|C\|$. No optimality among *all* conserved extensions is assumed. The two-mode projector resolves a compact part of an energy eigenspace; it does not require measuring the entire degenerate energy band.

## Native currents and their initial sign

Let $H$ be the native flight Hamiltonian and $P_x$ the port projector. Choose the symmetrized site density $q_x=\{P_x,C\}/2$ and for each native hopping edge define
$$
j_{y\to x}=i(P_y HP_x-P_x HP_y),\qquad
J_{y\to x}=\{j_{y\to x},C\}/2 .
$$
Because $[H,C]=0$, the exact continuity identity $i[H,q_x]=\sum_y J_{y\to x}$ follows by inserting $\sum_yP_y=I$. The currents are antisymmetric. All $q_x$ outside $S$ vanish as operators. The net outflow through its boundary vanishes for every state, since $\sum_{x\in S}q_x=C$. An individual edge can have a nonzero response. This convention specifies an explicit finite current. Spatial density involves off-diagonal port coherences.

On a sharp launch from 5 the current from 4 into 5 is initially zero for any parity-block input $\rho$. Direct integer multiplication gives
$$
\left.\frac{d}{dt}\langle J_{4\to5}\rangle\right|_0
=-\tfrac38\operatorname{Tr}(\rho K).
$$
In detail, set $h=2H$, $D=P_4hP_5-P_5hP_4$, and $W=(6\Pi|5\rangle)\otimes I_4$. In each parity and for both carrier signs the four-by-four row block $-([h,D]W)_5$ is exactly $-3I_4$. The factor eight restores the two factors of $h$ and the doubled dressing. Taking the Hermitian part yields the identity for any $K_b$. The equal real scalar contact hopping and real spatial matrix give the zero initial current. At the true selected reference the current slope is therefore $+3D(\chi_*\Vert\pi)/8>0$. The fixed-volume unpairing changes this slope by between $-0.086675$ and $-0.086654$ nats per squared native exchange time per unit intervention.

At the archived numerical reference, the 4-to-5 signed current on an excess launch is $+0.01719641$ at time $0.5$ and $-0.01155533$ at time $2$. Per-site finite differences of the independent densities agree with the edge divergence within $1.3\times10^{-9}$. At time 2 the bare heads have changed their mean charge by $+0.0000914690$ nats. The diagonal six-port part of the dressing reads $-0.06394216$ and the spatial-coherence part reads $-0.04695108$, yielding the conserved $-0.11089324$ nats. These finite-time numbers use numerical propagation. The conservation and initial signed-coefficient results use integer identities.

## Reading the interference without a supplied phase tag

The six support sites and intermediate port 9 form a connected network of blank seams. Under an isolated edge pulse the symmetric and antisymmetric one-probe position modes have eigenvalues 8 and 9, for either sign. A probe elsewhere sees eigenvalue zero on this edge. The compiler checks these assertions against every row of the exact 99-state contact. A duration $\pi$ implements a swap with no spectator phase. Swaps route any two labels to adjacency. A duration $\pi/2$ followed by a signed difference of sharp populations reads $Y_{uv}=-i(|u\rangle\langle v|-|v\rangle\langle u|)$.

To measure the real coherence $X_{uv}=|u\rangle\langle v|+|v\rangle\langle u|$, pick a second neighbor $w$ of $v$ distinct from $u$. Under a pre-pulse on $v,w$ of duration $\pi/16$, the sum of the two subsequent $Y_{uv}$ and $Y_{uw}$ effects is $-(X_{uv}+X_{uw})$. A pre-pulse of duration $\pi/18$ sends their difference to $-(X_{uv}-X_{uw})$. This follows from the eigenvalues 8 and 9 on the relevant linear combinations. Adding the two settings with weight $-1/2$ recovers $X_{uv}$. Each $Y$ read uses a separate preparation. The endpoint not firing scores zero, including when the probe sits outside the read network.

One diagonal population read and four settings for each of fifteen pairs give 61 spatial settings. The head expectation of $K$ can be read after removing the probe conditional on the terminal port, using the independently established head tomography. Keep the outcome probabilities as weights, including misses. Summing these joint signed means reads $\langle C\rangle$ for arbitrary position–head correlations. The implemented seven-port matrix effects reproduce $3\Pi$ within $3\times10^{-14}$; the population-level combined read on the numerical flying source agrees with direct evaluation. This is a controlled terminal apparatus. The flight before it uses the unmodified native Hamiltonian. A finite joint shot bank for this new meter has not been sampled. [Head inference from sharp clicks](../remote_heads/README.md)

The same construction extends **every** parity-block head operator. Conservation of this compact extension cannot single out the modular logarithm as physical stress. Its bounded support and zero net outgoing charge restrict its role on macroscopic observer balls. The result supplies a finite charge, measurement and orientation-dependent graph current; identifying proper directional pressure, transported overlap stress, common null geometry and their large-scale coupling requires additional physical input.

From the repository root run

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/dressed_charge
```
