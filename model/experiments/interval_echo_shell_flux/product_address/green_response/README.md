# What the current on three material guides measures

Three material doubles report their addresses on separately operated guides. After one release, a detector can count visits to its address, record signed crossings of its link, or report whether it was hit before the outer collector. An observer calibrates the single-step covariance on the guides and uses it to predict these records. An independent variation of the same eight-corner source changes its entropy and the angular distribution of visits. The experiment distinguishes a flux produced by repeated transport from acceleration of a second body. [Material apparatus](../README.md) · [Same-source entropy, energy and clock](../corner_source/README.md)

## Finite propagation from a measured step law

Let \(q_i\) be the probability per appointment of a step in either specified direction of guide \(i\). Both signs transfer on each of its eleven links under the native sharp exchange. The tag supply and shutter pattern are part of the preparation. The one-step mean is zero and the covariance is \(2D\), where \(D=\operatorname{diag}(q_i)\). With the outer collector at either end of any guide, let \(u(x)\) count expected visits per release, including the release address. For a corner distribution \(f\), conditioning on the last step gives
\[
 -\sum_i q_i\,[u(x+e_i)-2u(x)+u(x-e_i)]=f(x),\qquad u|_{\rm drain}=0.
\]
The finite Dirichlet solution is unique by the maximum principle. Expected signed crossings on the positive link are \(J_i(x)=q_i[u(x)-u(x+e_i)]\). Summing over a section counts final arrivals on its far side when all releases start on its near side. The probability of a first hit at \(y\) is \(u(x=y)/G(y,y)\), where \(G(y,y)\) is the expected number of visits under a release at \(y\). This follows by restarting at the first hit and includes time zero. The sine-transform solution is checked against all 1,000 vertex balances and the independent rational shell calculation. [Prediction and acquisition](measure.py) · [Cross-checks](test_green.py)

With two idle tags and the six single steps equally likely, the predicted fraction reaching either end of guide one is \(1/3\). Reassigning the two idle tags to the two signs of that guide changes its calibrated covariance and predicts 0.578927. The counts are 21,668 and 37,959 from 65,536 complete launches per setting. All fourteen probabilities for exit axis, final half-space and first visits agree within the simultaneous 99% Hoeffding radius 0.007782. Each setting has independent episodes. Multiple detectors in one episode can be correlated. Each kernel and all source probabilities are specified without detector counts. The same simulation draws the walk tags and executes the calibrated sharp transfers, so this checks the measurement protocol within Interval Echo. [Rates, counts and seeds](results.json)

## Large finite guides

The following existence construction uses the declared midpoint continuation of the triangulated screen. At face subdivision frequency \(m\), choose three different original faces. Just inside a side of each face, its \(m\) upward small triangles along the side alternate in a connected strip with \(m-1\) downward triangles. The dual strip is a simple path through \(2m-1\) distinct cells. A cell's three external matching edges use different gateways. Join the entrance to the exit inside each cell by a simple path in its connected twelve-port template. The resulting port paths are disjoint across the three faces. Reserve the source flag in a fourth face. Truncating the paths supplies three guides of any prescribed common even length for sufficiently large \(m\). Shutter unused contacts and use the measured sharp swap on each selected link. This constructs arbitrarily long addressed apparatus on that continuation, independent of the gateway labelling. The 20-owner acquisition uses three twelve-port guides. On the deposited 80-, 320- and 1,280-owner boards, independently checking each addressed seam gives three guides of 36, 84 and 180 sites respectively and a disjoint flag. All 132 ordered distinct endpoint pairs of the cell have an explicit twelve-port spanning path. Choosing these paths realizes exactly 12(2m−1) sites per guide on the stipulated continuation. [Constructed guides](guides.json) · [Transfer and topology checks](guides.py) [Carrier and deposited boards](../../../../executable_model/README.md) · [Midpoint construction](../../../../../code/support_refinement/README.md)

Write \(x\) relative to the centre of the eight source addresses, \(r=|x|\), \(r_D=(x^TD^{-1}x)^{1/2}\). At fixed positive \(q_i\) the centred, finite-range nearest-neighbour walk on the unbounded product lattice has
\[
 u_\infty(x)=\frac{1}{4\pi\sqrt{\det D}\,r_D}+O(r^{-3})
\]
for the symmetric source of unit mass. The random-walk Green theorem uses the independently measured covariance \(2D\). Averaging the eight bounded translations has no first-order term by source symmetry. On guides with the first drain at address distance of order \(L\), the difference \(u_\infty-u_L\) is harmonic in the interior with boundary values \(u_\infty\). It is uniformly \(O(L^{-1})\) there by the maximum principle. Differencing the leading term and bounding the two endpoint errors proves
\[
 J_i(x)=\frac{x_i}{4\pi\sqrt{\det D}\, r_D^3}
                +O(r^{-3}+L^{-1}).
\]
For \(r\to\infty\) and \(r^2/L\to0\) the vector error is \(o(r^{-2})\), uniformly in direction. With equal rates the leading vector is \(x/(4\pi r^3)\). This establishes the pointwise coefficient without uniform mixing of a finite shell. The norm is the diffusion norm of three independently chosen material addresses. The conditional mean increment of every freely walking guide marker is zero even where the current is nonzero. A gravitational interpretation would require source-dependent acceleration of a separately released test body, a common proper metric and calibration of the coupling. [Green asymptotic, Lawler and Limic, Theorem 4.3.1](https://www.math.uchicago.edu/~lawler/srwbook.pdf)

## Which part of the source reaches the field?

Write a source corner as \(a_i=\sigma_i/2\) with \(\sigma_i=\pm1\). Use
\[
 P_{b,s}(\sigma)=\tfrac18(1+b\sigma_1+s\sigma_1\sigma_2\sigma_3),
 \qquad |b|+|s|\le1 .
\]
The \(b\) derivative changes the first guide's marginal. The \(s\) derivative preserves every two-guide marginal. After tracing out the source tag, the entropy of the three orthogonal guide records at \(b=0,s=-1/2\) has derivative \(\partial_s S=\tfrac12\log 3\). The signed \(s\)-source is odd under each coordinate reflection. The Dirichlet solution inherits the three signs by uniqueness. Every reflection-invariant radial average has exactly zero \(s\)-response at every finite size. The total current through every enclosing surface also has zero derivative, since every setting releases one triple. Local detectors distinguish the settings. At \((6,6,6)\) the predicted first-hit probabilities at \(s=-1/2,+1/2\) are 0.291956 and 0.352713. Independent banks return 19,202 and 23,126 hits in 65,536 launches each. The reflected site exchanges the two predictions. [Angular detector records](multipoles.json) · [Acquisition and reflection check](multipoles.py)

For equal \(q_i=q\), set \(\Phi=q u\). In the infinite guide interior the Fourier denominator of \(\Phi\) is \(2\sum_i(1-\cos k_i)\). Relative to the source centre the two derivative numerators are \(-i\sin(k_1/2)\cos(k_2/2)\cos(k_3/2)\) and \(i\prod_i\sin(k_i/2)\). Their leading terms act as \(-\partial_1/2\) and \(-\partial_1\partial_2\partial_3/8\) on the continuum kernel \(1/(4\pi r)\). Consequently
\[
 \partial_b\Phi_\infty=\frac{x_1}{8\pi r^3}+O(r^{-4}),\qquad
 \partial_s\Phi_\infty=\frac{15 x_1x_2x_3}{32\pi r^7}+O(r^{-6}).
\]
For the error estimate, cut off a neighbourhood of the zero frequency. On an annulus of wave number \(h\) the remainders after subtracting the leading quotients and their order-\(j\) derivatives are \(O(h^{p-j})\), with \(p=1,3\) respectively. Integration by parts bounds the annulus contribution by \(C_Nh^{3+p}(1+hr)^{-N}\). Summing dyadically gives the displayed errors. Away from zero the original integrand, in integer address coordinates, is smooth and periodic. The centred phases are included together in this argument. Cutting off the homogeneous leading kernel changes only rapidly decaying terms away from its source. On the finite box the maximum principle transfers boundary errors of order \(L^{-2}\) and \(L^{-4}\) to the two interior fields. These are smaller than the respective natural field sizes at \(r=o(L)\). The correlation therefore has a measurable octupole, with zero monopole even when its entropy derivative is nonzero.

For an explicit shrinking family choose dyadic \(R\to\infty\), guide half-length \(L=R^3\), and \(\ell=R/L\). Use the clock ball of \(R\) address layers and the public source \(s=-1/2+\epsilon\ell^4,\ b=0\). All its corners lie inside every such ball. The surveyed boundary count is fixed, while \(\ell^{-4}\partial_\epsilon S=\log(3)/2\) exactly at \(\epsilon=0\). The total outward-current derivative also vanishes identically. Probabilities can be realized by fair bits at dyadic perturbations and the exact finite entropy gives the derivative. The experiment concerns a release lottery on a triple of guides. It exhibits entropy variation at the order used by a three-spatial-dimensional small-ball test even along a shrinking family. This rules out a uniform stationarity bound on this entire family of addressed preparations. Applying it to a fixed smooth physical perturbation requires identifying the same family with that perturbation in a common proper metric.

From the repository root, replay the two acquisitions with

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_shell_flux.product_address.green_response.measure
python -m model.experiments.interval_echo_shell_flux.product_address.green_response.multipoles
```
