# How finite repair can select stationary generalized entropy

A region can have one of three geometries, all with the same proper volume. Its observer reads the geometry sector and a material bit. Independent boundary observers read three, four or five binary tiles. Under a reversible repair law, entropy gained by redistributing the coarse states is exactly balanced to first order by the change in expected boundary area. Changing the repair probabilities breaks the cancellation even though the owners continue to agree on every shared register. [Protocol](PROTOCOL.md) · [Source and readers](measure.py) · [Observed counts](results.json)

The three sectors are centered geodesic balls in separately supplied ultrastatic, constant-curvature geometries. Each has volume $4\pi/3$ in the length standard of the source. Their curvatures were chosen to give three, four and five tiles of proper area $\pi$ each. An interior counter stops at a common volume threshold. A second interior counter verifies the stop, and a surface counter reads the area. These instruments use independent populations. [Proper-volume and area calibration](../fixed_volume/README.md)

A microstate is $(n,b,w)$ with tile count $n$, material bit $b$ and $n$ boundary bits $w$. Fair attempts flip the material bit or a visible tile bit, remain quiet, or move between adjacent geometry sectors. A newly exposed tile enters in state zero and an outgoing tile must be in state zero. Each attempted move is an involution. The 112-state transition kernel is symmetric, connected and has a nonzero holding probability, so it selects the unique uniform microscopic equilibrium from every launch. Its coarse weights are $p^*_{n,b}=2^n/112$. At fixed volume, with $A_n=n\pi$ and independently specified $\eta=\log 2/\pi$, every admitted coarse distribution satisfies

$$
S_{\rm bulk}(p)+\eta\sum_{n,b}p_{n,b} A_n
=\log 112-D(p\Vert p^*).
$$

Thus every mass-preserving probability variation has zero first derivative at $p^*$. Here $S_{\rm bulk}$ is the Shannon entropy of the observed sector and material bit. The edge entropy is independently readable through its conditional words. The equality between edge entropy and area comes from the supplied bit capacity per geometric tile. The stationarity follows from the repair kernel and this dictionary, without inserting an Einstein equation or an entropy-balance update.

After 512 attempts from either of two opposite sharp states, six coarse outcomes and an independent surface survey give a maximum residual below 0.0072 nats over all five independent first variations. A fresh binary-word survey checks the conditional entropy of each boundary to within 0.001 nats. Halving acceptance of an area-increasing move leaves the same metric sectors, bits, receivers and volume. The changed process selects a uniform coarse law and leaves first-variation residuals as large as 1.39 nats under the original area coefficient.

This gives a fixed-proper-volume stationarity example at one physical size and in one rest frame, on a declared comparison process. Deriving an Einstein response on the registered repair requires a measured edge capacity per proper area, the corresponding bulk stress, and agreement of those identifications along a common refinement. The scalar bit and the sector label here supply a classical coarse entropy.

Replay from the repository root with `PYTHONPATH=code python -m matched_einstein_limit.finite_entropy_region.measure`.

## One equilibrium across shrinking regions and moving observers

Five curvature sectors with two material-bit outcomes give nine independent probability transfers. At each flat reference radius, the source uses the first observer's equal-volume areas to determine which surface symbols are accessible. The same sector and bit counts are then read against ten separately surveyed rest spaces, including transverse boosts. Each survey stops an interior count at the flat volume target and checks the stop with another population. Surface counts have a separate random stream. The source prepares one ensemble per geometry family and radius. [Acquisition protocol](REFINEMENT_PROTOCOL.md) · [Source and reader](refinement.py) · [Counts and interpretations](refinement.json)

| Reference radius | Constant-curvature spacetime, reversible | Same spacetime, cleared contraction | Curved-plane product, reversible |
| ---: | ---: | ---: | ---: |
| 1/4 | 0.400 | 2.802 | 7.365 |
| 1/8 | 0.162 | 2.589 | 6.839 |
| 1/16 | 0.053 | 2.555 | 6.791 |
| 1/32 | 0.044 | 2.529 | 6.750 |

Entries are the largest observed absolute value of $(\delta S_{\rm bulk}+\delta A)/\ell^4$ across the transfers, observers and two launches. The finite bound specified at the smallest radius was 0.25. The largest area slope alone there was 2.515 on the maximally symmetric spacetime. All held volume counts on both preparations were within 3.69 standard deviations of their targets.

### Why the full refinement works in the symmetric example

Write $B=\ell^{-6}$. All sectors have the same number of registers and $B$ ordinary symbols per register. A full tile of area $1/B$ admits one extra symbol. Sector $i$ has $n_i=\lfloor B A_i(u_0)\rfloor$ full tiles. The inaccessible fraction of the last tile contributes zero symbols. A uniform refresh of the allowed symbols precedes a symmetric sector proposal. Contraction succeeds when every symbol on a lost tile is ordinary. On the coarse sector the acceptance is $(1+1/B)^{-\max(n_i-n_j,0)}$ for a move from $i$ to $j$. The bit is separately refreshed. Detailed balance and positivity select $p_{i,b}\propto(1+1/B)^{n_i}$ from either launch. Clearing obstructing marks on contraction makes the sector probabilities equal. The computation uses the exact coarse transition kernel and samples independent terminal experiments in batches. A test enumerates all microscopic words for a small alphabet.

For the transfer $j\to i$, let $\Delta A=A_i-A_j$. Its residual in the source frame is $\Delta A-(n_i-n_j)\log(1+1/B)$. The two floor errors differ by less than one tile. Integrating $1/(1+t)$ on $[0,1/B]$ gives $0\le 1-B\log(1+1/B)\le 1/(2B)$. Hence

$$
 |R(u_0)|\le \log(1+1/B)+|\Delta A|/(2B)=O(\ell^6).
$$

This bound holds uniformly on bounded areas for every normalized transfer. In a maximally symmetric four-spacetime, every observer's central exponential slice is a three-dimensional space form of the same sectional curvature. Hence $A_i(u)=A_i(u_0)$ at fixed proper volume, and the bound gives $R(u)=o(\ell^4)$ in every frame. The finite comparison uses a classical sector entropy, specified surface alphabets, calibrated geometry and dimensionless coefficient $\eta=1$. A separate Lean proof gives the sufficient bound $(1+|\Delta A|)/B$, the uniform-word transfer identity, and a lower bound on the errors of two frames sharing a scalar bulk readout. [Verified finite inequalities](../../../Lean/ObserverPatchHolography/EinsteinBranch/TileCapacity.lean)

### Why a transverse boost distinguishes the two geometries

On the product $\mathbb R^{1,1}_{t,z}\times M_K^2$, a rest-space vector projects onto the curved plane with principal stretches $\sqrt{1+u_x^2+u_y^2}$ and one. Write $b^2=1+u_x^2+u_y^2$, $f=b^2\cos^2\phi+\sin^2\phi$, $q^2=(1-\mu^2)f$ for spherical directions $(\mu,\phi)$. The angular exponential scales the component perpendicular to the projected ray by $s_K(rq)$, the ratio of the curved sine radius to the flat radius. The pullback metric has a rank-one correction $(s_K^2-1)ll^\top$, where $l$ is the pullback of that unit perpendicular, $|l|^2=b^2/f$, orthogonal to the radial direction. Thus the area density is $J=\sqrt{1+(s_K^2-1)b^2/f}$. Integrating $r^2 J$ over solid angle gives area and integrating over radius gives volume. Taylor expansion gives mean density $1-K b^2 r^2/9+O(r^4)$. Restoring the target volume then gives $\Delta A=-(4\pi/15)K b^2\ell^4+O(\ell^6)$ relative to flat space. The calculation uses the finite-radius integrals for all marker draws.

For a source referenced to $u_0$, the limiting residual divided by $\ell^4$ is $-(4\pi/15)K(u_x^2+u_y^2)$. For $K=1$ and $u=(3,2,2,0)$ this is $-6.702$. At the smallest radius the maximum absolute observed residual for this observer was 6.750. For the rest frame and both boosts along the flat line it was below 0.043. More generally two observers sharing a scalar bulk slope have residuals whose difference is exactly their area-slope difference. Changing that common bulk slope cannot cancel two different leading area slopes. An anisotropic realization would require observer-dependent bulk readings of the shared state, together with a physical calibration of the corresponding boundaries.

Replay the comparison with `PYTHONPATH=code python -m matched_einstein_limit.finite_entropy_region.refinement`.

## Calibrating the coefficient and withholding a curved sector

At proper volume $4\pi/3$, the independent surface counts give areas $9.42465$, $12.56648$ and $15.70791$ in the three sectors. The word detector observed every allowed word, giving capacities $\log 8$, $\log 16$ and $\log 32$. In a retrospective split, calibrate $\eta=\Delta\log\Omega/\Delta A$ on the three-to-four-tile transfer and withhold the fifth sector. The calibrated slope is $0.2206193$ nats per area unit. The four-to-five-tile transfer gives $0.2206472$ and a capacity-minus-prediction residual of $0.00008752$ nats. A simultaneous 99% interval puts that residual between $-0.000793$ and $0.000967$ nats. We use a gamma interval for each of the three stopping volumes, a Poisson intensity interval for each stopped surface count, and a union bound. The geometric derivative $dA/dV=2S'_K/S_K$ bounds the area correction from each stopping volume to the common target. The alphabets are exhaustive by construction, and the chart curvatures were chosen to assign area $\pi$ to a tile. This is an independent readout of that assigned relation. [Geometry and word counts](results.json) · [Coefficient calculation and intervals](held_calibration.json)

At $\ell=1/32$ a separate split calibrates on the $K=-1,+1$ transfer of the rest observer and withholds all nine other frames. Dividing the log ratio of the first launch's two terminal sector counts by the measured area change gives $\widehat\eta=0.9993684$ in maximally symmetric spacetime. The maximum held $|\Delta\log\Omega-\widehat\eta\Delta A|/\ell^4$ is $0.003749$. The corresponding entropy-plus-area residual is at most $0.044804$ across the two launches. On the curved-plane product, calibrating with the source's allowed symbols gives $\widehat\eta=1.0000246$, while the maximum held capacity residual reaches $6.70821$ in the same normalization. These held-frame values are point estimates. Both comparisons use supplied geometry and a surface alphabet assigned in the source frame. An Interval Echo determination requires that its independently measured proper areas and native completion capacities obey the same relation. [Refinement counts](refinement.json) · [Calibration replay](held_calibration.py) · [Regression](test_held_calibration.py)
