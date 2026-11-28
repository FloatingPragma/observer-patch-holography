# Does the parked source give all three signals a common delay?

**No deterministic change of elapsed time reproduces all three port curves near emission.** This conclusion applies to the registered exterior material tags, source preparations and marginal detector probabilities. It follows with exact rational arithmetic on the same 240-port candidate, with the exterior clock as the independent time standard. It tests a proposed interpretation of the earlier finite-appointment differences. The source comparisons and port choices were known when the test was designed. [Protocol and preparation](README.md) · [All four exact certificates](common_delay.json)

On every source word the local native column check gives \(H=L/2+\text{constant}\) for either load, where \(L\) is the integer Laplacian of all open contacts. The constant cancels in port occupation. For emitting port \(a\) and receiving port \(b\), let \(m_n=[(L-5I)^n]_{ba}\). Finite dimensionality and the real symmetric generator make
\[
 P_{ba}(t)=\left|\langle b|e^{-it(L-5I)/2}|a\rangle\right|^2
          =\sum_{k\ge0} c_k t^{2k}.
\]
Multiplying the two amplitude series gives
\[
 c_k=\frac{(-1)^k}{2^{2k}(2k)!}
       \sum_{j=0}^{2k}(-1)^j {2k\choose j}\,m_j m_{2k-j}.
\]
The executable weights these rational coefficients by the full pairs of sharp source-meter outcomes from the original archive. Its inverse ordering maps each meter outcome to its calibrated geometry word. It retains all 81 pairs at each of six placements. There is no fit to a time trace. Since \(\|L-5I\|\le5\), every integer moment through order 24 fits the checked integer range. The products and sums in the coefficient formula use arbitrary-precision integers. A two-port test checks every coefficient against the exact cosine solution. [Calculation](common_delay.py) · [Independent coefficient and retiming controls](test_signals.py)

Write \(x=t^2\) and let \(B_i(x)\) and \(S_i(x)\) be the fresh and held curves for near return, through transmission and far return. The common-delay hypothesis is that one continuous nonnegative effective elapsed time \(\tau(t)\), with \(\tau(0)=0\), makes all three held curves \(S_i(t^2)=B_i(\tau(t)^2)\) on a neighborhood of launch. The squared-time coordinate uses the time-reversal symmetry of each individual port probability. No prior differentiability or polynomial model for the proposed clock is needed. The reference near-return curve has derivative \(-5/4\) at \(x=0\), so it has a locally analytic inverse. Its observations force the unique candidate \(y(x)=B_{\rm near}^{-1}(S_{\rm near}(x))\), with \(y(t^2)=\tau(t)^2\).

Let \(k\) be the first order in \(x\) with any held-minus-fresh difference \(d_i=\operatorname{coeff}_{x^k}(S_i-B_i)\). Lower-order equality and the near slope imply \(y(x)=x+h x^k+O(x^{k+1})\) here, with \(k>1\) in every enrolled arm. Substitution in the other two reference curves then **requires** \(d_i=h B_i'(0)\) for all three channels. This provides a falsifier independent of any choice of delay size, gain, clock frequency or later-time fitting window. A synthetic control that applies \(y=x+x^2/7\) to all three curves gives a zero residual, as it must.

The three fresh slopes are \((-5/4,1/4,-5/4)\). The first discrepancy is at \(k=7\), or \(t^{14}\), for every stock. The exact leading vectors are \((d,0,0)\).

| Change from fresh | Exact near-return coefficient \(d\) at \(t^{14}\) |
| --- | ---: |
| Positive tilt | \(127/3417442234033766400\) |
| Negative tilt | \(-127/3417442234033766400\) |
| Fourth-order tilt | \(1/99166850541158400\) |
| Shared | \(1/231389317929369600\) |

Matching the near return forces the far return to change by \(d\), because their reference slopes are equal. Its actual coefficient is zero. The contradiction holds arbitrarily close to launch, independent of the sensitivity needed for a finite physical counting experiment. The coefficients characterize the exact ideal probabilities of these simulated observations. In particular, the finite-time throughput difference at appointment 16 cannot be promoted to a shared propagation delay for the complete pulse curves. The candidate's route-dependent exchange dynamics supplies all observed records.

A distribution of hidden elapsed times, a receiver-dependent change in the response to passage, or distinct travel times along distinct paths are different hypotheses. The argument also permits two curved-spacetime trajectories to have different gravitational delays. A physical Shapiro measurement would identify a light pulse, calibrate endpoint and geometrical path records, and compare travel times. The spherical experiment implements precisely that comparison on its specified light and clock action. [Matched-endpoint lens and radar records](../dust_exterior/README.md)
