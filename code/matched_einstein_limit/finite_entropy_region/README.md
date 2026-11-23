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
