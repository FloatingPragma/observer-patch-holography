# Energy conservation for a moving region of the four-tick pump

On this 80-cell vacancy pump, the energies of all 960 ports obey an exact local continuity equation. It holds for every joint quantum state of the mobile, the heads and the strain register, including mixed states and preparations with an untouched external reference. The observer can change which fraction of each port belongs to a region while monitoring its energy. The energy carried by that change of membership enters the balance.

Let $P_p$ project onto mobile port $p$, leaving the internal degrees of freedom untouched. The source Hamiltonian $H$ assigns half of each boundary bond to either side. Port number and energy have observables $P_p$ and $e_p=(P_pH+HP_p)/2$, with $\sum_p e_p=H$. For distinct ports, orient the current from $p$ to $q$ and define

$$
j^N_{pq}=i(P_pHP_q-P_qHP_p),\qquad
j^E_{pq}=\tfrac i2(P_pH^2P_q-P_qH^2P_p).
$$

Both are Hermitian and antisymmetric under exchanging the ports. Summing directed currents with $p$ fixed gives $-i[H,P_p]$ and $-i[H,e_p]$. The energy identity follows by expanding $[H,P_pH+HP_p]=H^2P_p-P_pH^2$. It applies before any state is chosen. These currents supply one local representation of the regional flow. A choice of pair currents can differ by a divergence-free circulation.

For any differentiable memberships $0\leq w_p(t)\leq 1$, put $E_w=\sum_p w_p e_p$. Under native evolution the measured rate obeys

$$
-\frac{d}{dt}\langle E_w\rangle
=\sum_{p<q}(w_p-w_q)\langle j^E_{pq}\rangle-\sum_p\dot w_p\langle e_p\rangle .
$$

Replacing each energy operator with the corresponding number operator gives the number balance. A common energy-origin shift adds $c$ times the entire number balance. For nested regions, the membership of the annulus is the difference of the two memberships. The inner and annular accounts add to the outer account at every time. If the memberships partition unity, their corrected outward rates sum to zero. Partial memberships can be realized as averages of sharp enrollments by using a common uniform threshold on each fresh copy. This preserves nesting on each draw.

The code enumerates two-step paths in the actual integer matrix $40H$ and accumulates every coefficient of the 960 separate port conservation identities. The sparse operator comparison verifies all 276,704 nonzero entries of the resulting family. In this representation 2,520 native port pairs carry number current. The energy current also connects 3,600 pairs separated by two native hops. The regional energy includes contact terms, so this support follows from the product $H^2$.

## Independent finite-time check

The three enrolled regions are the moving inner region, the annulus and the exterior. Their five native port groups use graph distances fixed before reading the energy. At age four the outer region contains all left ports and half of the shared zone. As time increases, the two portions of that zone enter with velocities $+1/4$ and $-1/8$ per native tick. The inner region also changes the enrollment of the left neighborhood. The archived weights specify the complete rules for all three regions.

The endpoint cohorts measure port numbers and symmetrized energies on the two sides of age four, $4\pm 2^{-16}$. Separate midpoint cohorts measure the inter-group current operators and the regional inventory carried by the membership change. All envelopes include the source preparation error, error in the pump evolution and outward rounding. The integer Taylor polynomial for each independent endpoint has degree four. For a duration $d$, midpoint integration contributes at most $d^3\sup|F''|/24$. The current bound uses $\|[H,X]\|\leq 2\|H\|\|X\|$ twice and includes the first derivative of the moving weights. The sweep uses their constant velocities. These are certified population expectations from the executable law.

For the paired reference, the table gives signed readings on the outer region. Endpoint loss means initial inventory minus final inventory. A negative value records a gain. The corrected calculation subtracts the membership term from the outward bond current.

| Material reading | Endpoint loss | Integrated outward bond current | Integrated membership change | Absolute residual bound |
|---|---:|---:|---:|---:|
| Energy, native exchange units | $[-6.178477,-6.178319]\times10^{-5}$ | $[-4.651049,-4.650862]\times10^{-5}$ | $[1.527441,1.527444]\times10^{-5}$ | $1.730\times10^{-9}$ |
| Number | $[5.022342,5.022417]\times10^{-6}$ | $[6.998635,6.998814]\times10^{-6}$ | $[1.976344,1.976346]\times10^{-6}$ | $1.269\times10^{-10}$ |

Omitting the membership change misses the reference endpoint energy loss by more than $1.52727\times10^{-5}$. The corrected intervals contain zero on the inner region, annulus, exterior and outer region for all three actual preparations and both material observables. The three partition bins also sum to zero net loss.

This result supplies conservation of the native material energy through specified changing graph regions. A proper momentum frame, physical three-dimensional regional membership and the spatial components of covariant stress conservation require their own source measurements. The moving-cut comparison tests the proposed uniform sweep against the independently reconstructed face curvature. [Moving-cut comparison](../moving_cut.py) · [Finite records](results.json) · [Operator and evolution checks](measure.py) · [Regression tests](test_ward.py)
