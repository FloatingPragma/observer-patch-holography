# Permanent clock pairs under refinement

The pressureless source admits a finite construction in which every previously enrolled physical pair, clock program and public outcome survives later refinements. On the matched regional sign records, a cofinal choice of cutoff gives an $O(R^5)$ error in the leading entropy of a proper ball of radius $R$ and $O(R^{12})$ relative entropy. The aperture takes one clock sample per selected endpoint. The radial response, dust and clock programs are inputs of the published-record candidate. [Dust geometry and independent stress](../PROOF.md) · [Two-clock channel](channel.py) · [Fixed-register reference and restriction](registered_reference.py)

## One source at every cutoff

Work in a fixed finite instrumented spacetime window in the regular homogeneous interior, with a buffer around tested apertures. Use separation bands $(D/2,D]$ with $D=D_0 2^{-j}$. Permanently assign each band a spatial grid of mesh $h_D\le D^2$, positive separation quadrature with interval widths $q_d\le D^2$, and positive normalized solid-angle weights $q_n$ from a partition of angular diameter at most $D$. Each service has a distinct address $(D,x,d,n)$ and endpoints $x,x+dn$, including when other services occupy the same location. One convention uses oriented services and total angular weight one.

With independently specified comoving mass density $\rho_0$, assign entropy exposure
$$
w=(8\pi\rho_0/15)\,h_D^3q_dq_n.
$$
For the unit initial sphere $8\pi\rho_0/15=2M/5$. Rounding each exposure to a rational with absolute error $O(h_D^3q_dq_nD)$ suffices. Freeze the finite clock grid, scale envelope, ticket precision, exposure and program upon enrollment. With $B_D\asymp1/D$ lag bins and maximum class $B_D-1$, every compact interval inside $(-1,1)$ is eventually covered with a margin. There are at most $2B_D$ lag classes. The largest class multiplier $1/[1-(k/B_D)^2]$ is at most $B_D$. The total allocation is $O(B_D^2 w)=O(D^6)$ for bounded clock scale because $w=O(D^8)$. Choose the finitely many coarser grids fine enough to pass the channel's quarter-ticket capacity check. This gives nonnegative normalized finite laws before the two readings.

The executable band uses $h_D=D^2/8$, separation spacing $D^2/8$, $8/D^2$ equal bins in the direction cosine and twice that number in azimuth. Even at a pole the angular diameter is bounded by a constant times the square root of the cosine-bin width, hence by $O(D)$. The exposure is an exact rational for the recorded stocks. The finite fixtures sample addresses from these much larger permanent populations. The bounds here apply to the complete population sum. [Band construction and sampled regional records](regions.py)

Preassign every independent sign and gate entry needed on the finite clock window. Deferred draws at previously unqueried entries realize the same law. Different services have independent seeds conditional on the stock. At the next cutoff take the product with the new band. Projection forgets only those added services. This preserves the previous hidden tape law and hence every previously available joint or adaptive public query. One product-space coupling reuses the literal old seeds and outputs. The blank and each admitted stock have this exact projective construction. Actual observers receive the selected outputs only.

## The same regional observable

For one service, $I\subseteq\{A,B\}$ denotes the observed ends and their actual clock settings. Keep every returned sign. Condition on the prepared stock, clock settings and geometric mask, selected independently of the sign tape. At the matched settings the blank law is $P_0(a,b)=(1+c_0ab)/4$ with $c_0=1/2$. Each individual sign is exactly fair at every setting under both the blank and $P_c(a,b)=(1+cab)/4$. Disabling shared lag on the same configured clock program realizes a reference rail at all settings. A query at two successive times of the same endpoint is available operationally. Its joint-history entropy is a different observable.

To compare physical apertures with different masks, every service retains a two-sign public register. It copies the observed signs in $I$ and draws the missing coordinates from $P_0(\cdot\mid I)$ with independent apparatus randomness. An empty service draws $P_0$, a singleton retains its actual bit and draws its partner, and a complete service retains both actual bits. For $N$ services the public entropy and relative entropy are exactly
$$
 S_{\rm pub}=N H(P_0)+\sum_{\rm complete}[H(P_c)-H(P_0)],\qquad
 D_{\rm pub}=\sum_{\rm complete}D(P_c\Vert P_0).
$$
The completed blank has entropy $NH(P_0)$ at every mask. Subtracting a separate raw blank at the held typed mask gives the same entropy difference. For an uncompleted mask with $n_1$ singletons and $n_2$ complete services the blank background is $n_1\log2+n_2H(P_0)$, which can change between held and flat physical apertures.

For nested observed sets $I\subseteq K$, keep the actual $I$ coordinates of a $K$ record, discard the remainder and complete conditionally from $P_0$. Reducing first to $K$ and then to $I$ gives exactly the direct reduction to $I$. The copied coordinates survive both steps and the final draws only use them. This proves composition for the product and for deleting newly added services. Overlaps of different tilted slices identify common event and clock labels. Finer rod surveys keep earlier appointments available verbatim. A newly resolved aperture can change the mask or time settings and is controlled by the estimates below. Private seeds are outside the published entropy. An unconditioned random mask would contribute its own entropy and correlations.

The mask identifies which signs were physically returned. Two observers reading a common physical tape entry agree on its sign. If an observer only receives $A$, its synthesized partner has $\mathbb E[A B_{\rm fill}]=c_0$. A larger held aperture records the actual partner with $\mathbb E[AB]=c$. Every coupling with the same physical $A$ therefore satisfies $\Pr(B_{\rm fill}\ne B)\ge |c-c_0|/2$. Compatibility for completion bits is the stochastic restriction law above. Treating such a bit as a prior physical reading of the distant end would contradict these correlations.

## Error bounds on proper balls

Fix a compact set of interior events and a bounded boost range. For small $R$ the exponential slices are smooth spacelike graphs $\eta(x)$ in comoving labels, with bounded second derivatives and $|\nabla\eta|\le v_*<1$. The projected ball $E_R$ has diameter $O(R)$, volume $O(R^3)$, boundary area $O(R^2)$ and boundary tube volume $O(R^2s)$ for $s\le cR$. These uniform estimates follow by rescaling the smooth exponential map about a regular point. Only $d\le C R$ can place both ends inside.

Put $z=[\eta(x+dn)-\eta(x)]/d$ and $f=\sqrt{a(x)a(x+dn)}/(1-z^2)$. The integral mean-value formula bounds the derivatives of $z$ in $x,d,n$ even at $d\downarrow0$. The same holds for $f$. On this lag interval the finite grid, tickets and separate scale gates approximate the linear entropy slope times the correlation shift by $wf+O(h_D^3q_dq_nD)$. The lag-grid error, including the two local tick roundings, is $O(B_D^{-1})$. The permanent ticket and gate errors can be made this small in relative units. For an additional observer-side conformal clock error it suffices to take absolute precision $O(Rd_{\min})$, where $d_{\min}$ is the smallest included separation. This gives lag error $O(R)$ for every included service.

* At fixed $d,n$, integrate $f$ on $E_R\cap(E_R-dn)$. Grid cubes meeting either boundary contribute error $O(R^2h_D)$. The interior Lipschitz error is $O(R^3h_D)$.
* Changing separation by $s$ changes this integral by $O((R^2+R^3)|s|)$. Rotating a direction through angle $\alpha$ moves an endpoint by at most $d\alpha$ and changes $f$ by $O(\alpha)$. The resulting error is $O((R^2d+R^3)\alpha)$.
* Each band has total separation weight $O(D)$. Summing $D\le C R$, spatial and separation errors are $O(R^2\sum D^3)=O(R^5)$. Angular and lag errors are $O(\sum(R^2D^3+R^3D^2))=O(R^5)$. Uniform geometric boundary displacement $O(R^2)$ costs $O(R^3 R^2)$ by the same boundary-cube count.
* If the first omitted interval is $(0,d_{\rm cut}]$, its contribution is $O(R^3d_{\rm cut})$. Choose $d_{\rm cut}$ comparable to $R^2$ and at most $R^2$. Every cutoff remains finite.

Write $a_0$ for the scale at the center and $\gamma=(1-|v|^2)^{-1/2}$. The tangent projected ball is $L B_R$ with $\det L=\gamma/a_0^3$ and $|L^{-1}n|=a_0\sqrt{1-(v\cdot n)^2}$. Its boundary differs from the true one by $O(R^2)$, while the true lag differs from $v\cdot n$ and the scale from $a_0$ by $O(R)$. These substitutions cost $O(R^5)$. The overlap of two Euclidean $R$-balls, integrated along their separation, equals $\pi R^4$ by integrating $4\pi R^3/3-\pi R^2s+\pi s^3/12$ from $0$ to $2R$. Set $s=d|L^{-1}n|$. The normalized angular mean of $(1-(v\cdot n)^2)^{-3/2}$ is $\gamma$. Combining the volume, line length, two-clock and scale factors gives
$$
 \Delta S_{\rm pub}={8\pi\rho_0\over15}\,{\pi R^4\gamma^2\over a_0^3}+O(R^5)
 =2\pi C_R T_{uu}+O(R^5),\qquad
 C_R=4\pi R^4/15,\quad T_{uu}=\rho_0\gamma^2/a_0^3 .
$$
The stress is also read from the independent shell counts and transported clocks. The equal-volume stop changes the leading fourth power by $O(R^6)$.

## Relative entropy and the two limits

For $\delta=c-c_0$ the pair identities and chi-squared bound give
$$
 H(P_c)-H(P_0)=-\operatorname{atanh}(c_0)\delta-D(P_c\Vert P_0),\qquad
 0\le D(P_c\Vert P_0)\le\sum_{ab}{(P_c-P_0)^2\over P_0}
 ={\delta^2\over1-c_0^2}.
$$
At fixed $d,n$ there are $O(R^3/h_D^3)$ full anchors. Since $|\delta|\le C h_D^3q_dq_n$, the sum of their absolute shifts across all bands is $O(R^3\sum_{D\le C R}D)=O(R^4)$. The largest contributing shift is $O(R^8)$. Thus $D_{\rm pub}=O(R^{12})$ uniformly in the number of additional finer bands. This includes all returned signs on this aperture and bounds the nonlinear entropy remainder.

At fixed physical aperture and fixed old clock appointments, additional bands retain every old complete law. The new relative entropy equals the old value plus nonnegative new contributions. As the ball shrinks, each previous band wider than its diameter has only single or empty services, whose held observables agree exactly with the blank. Thus $D_{\rm pub}/R^4\to0$ in the shrinking limit. Data processing explains the distinction. A scheme with total divergence tending to zero at a *fixed* aperture could not recover any old nonblank marginal of positive divergence by a stock-independent restriction. Permanent bands preserve the microscopic source and satisfy the shrinking regional bound together.

The argument covers the full recorded output of the specified single-time apertures on homogeneous pressureless stocks. Additional retained time histories, sources with pressure, and autonomous selection of this material by native Interval Echo require separate stress and entropy derivations.
