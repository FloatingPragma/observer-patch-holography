# Entropy on overlapping balls of one homogeneous material

A comoving array of faithful endpoint pairs assigns a von Neumann entropy to every observed region. For a stock proper ball, observers record the enrolled worldline addresses \(E_g\). They read its material entropy change on copies by comparing excited and blank rotations on those **same** addresses, \(s_R=S(\rho_{E_g})-S(\sigma_{E_g})\). At one dust-clock epoch, all sufficiently small volumes in a compact, padded part of the homogeneous interior have
\[
 s_R=\frac{8\pi^2}{15}R^4\rho\,\gamma_u^2+o(R^4)
\]
for the ten specified timelike directions. Here \(R\) labels enclosed proper volume \(4\pi R^3/3\), \(\rho\) is the separately weighed rest density and \(\gamma_u=-U\cdot u\) at the center. Every center and radius uses the same material preparation at a given cutoff. The material's three angular responses were chosen to reproduce these ten directions. A further boost distinguishes it from pressureless matter read by an arbitrary observer. [Executable and rational enclosures](experiment.py) · [All finite readings](results.json)

## Angular response of common pairs

At the clock epoch with unit comoving scale, put \(v=\mathbf u/u^0\), \(\gamma=u^0\). In the tangent chart a moving rest ball projects onto the comoving endpoint labels as \(E=A B_R\), with \(A\) symmetric, eigenvalue \(\gamma\) along \(v\) and eigenvalue one across it. Thus \(\det A=\gamma\) and \(k(n)=|A^{-1}n|=\sqrt{1-(v\cdot n)^2}\). Its overlap with its translate by \(dn\) is \(\gamma I_R(dk(n))\). The cap-volume formula for two Euclidean balls gives
\[
 I_R(x)=\frac{4\pi R^3}{3}-\pi R^2x+\frac{\pi x^3}{12}
 \quad (0\le x\le 2R),\qquad I_R(x)=0\quad(x>2R).
\]
Integration gives \(\int_0^\infty I_R(x)\,dx=\pi R^4\). An array whose signed short-distance cut-entropy density per unit comoving volume and separation is \(-c\) in direction \(n\), with a compensating total \(+cL\) at distance greater than \(L\), therefore contributes \(2\pi cR^4\gamma/k(n)\) whenever the region diameter is smaller than \(L\). Indeed a cut is measured by the symmetric difference of the region and its translate. The two volume terms cancel between short and distant pairs. The remaining quantity is twice the integrated short-pair overlap.

Use the three unoriented cubic orbits given by the axis, face diagonal and body diagonal of the dust-frame triad. Within each orbit take equal weights and signed density \(c_i=(4\pi/15)\rho w_i\). Let \(M_{ui}\) be the orbit mean of \(\gamma/k(n)\). At rest, at \(u=(5/3,4/3,0,0)\) and at \(u=(3,2,2,0)\), solve \(Mw=(1,25/9,9)^T\). The remaining enrolled rows follow by cubic symmetry. Rational square-root and determinant enclosures in the executable give
\[
 0.06994924<\det M<0.06994925,\qquad
 w=(3.61011676,\ 13.76937831,\ -16.37949507)
\]
to the displayed precision for the weights. Signed increments are physically available about the same full-rank blank. The directions, negative weight, compensation and density conversion are design inputs. The archive material's dust action and the shell weighing supply the stress independently. They do not derive these angular conversion coefficients. [Timelike frame menu](../../../../../code/matched_einstein_limit/frame.py) · [Dust action and shell weighing](../PROOF.md)

## Finite source and errors

Prepare a padded comoving cubic mesh of spacing \(h\) before the surveyed epoch. All observers may restrict the material on independent copies made by the same preparation. At each lattice site and primitive unoriented vector \(b\), put separate endpoint registers on the sites \(x,x+jhb\) for \(1\le j\le m=\lfloor L/(h|b|)\rfloor\). Pair preparations and distinct direction labels are independent. In the code \(L=1/8\). Write \(\Delta=h|b|\), \(g_i\) for the number of directions in orbit \(i\). Assign each short pair the entropy increment
\[
 e_j=-h^3\Delta\,c_i/g_i\quad(j>1),\qquad e_1=-3h^3\Delta\,c_i/(2g_i).
\]
The first positive node also covers the half interval at the origin. One further pair at the lattice displacement \(\lceil 2L/\Delta\rceil hb\) receives \(e_{\rm far}=-\sum_j e_j\). Outside a sufficiently small ball its endpoints are never both enrolled. Keep potential partners through the entire outer displacement, including those with a base site just outside a measured ball. This requires a finite margin on a finite source. The finite array is covariant under lattice translations within this padded interior.

For each intended \(e\), apply the native rotation giving \(h_2(p(d))-h_2(3/8)=e\), where \(p(d)=1/4+\sin^2(\pi/6+d)/2\) and \(h_2\) is binary entropy in nats. The derivative at zero is \(\log(5/3)\sqrt3/4>0\). Both signs have a unique small exposure. All memories have the same full-rank blank with eigenvalues \((5,1,1,1)/8\), and each complete pair preserves its spectrum. The endpoint exposure is controlled by the weighed stock and pair label, before any spatial boundary is selected. The registers are held during the spacetime readout so that a crossing at a different dust time on a tilted ball still reads the same pair. The result concerns one fixed epoch of the homogeneous dust solution. Preparing and parking these correlations at all separations uses additional physical apparatus. [Native pair and entropy reader](../../../interval_echo_partition_repair/published_strings/memories.py)

For the set of registered lattice sites in the projected ball let \(N\) be its cardinality and \(N_j\) its overlap count after shifting by \(jhb\). There are \(2(N-N_j)\) short cut pairs and \(2N\) distant cut pairs for each direction. Thus the total change equals \(-2\sum_j e_jN_j\) exactly in real arithmetic, including centers off the grid. Complete pairs make zero entropy contribution. The executable obtains \(\sum_jN_j\) by grouping lattice points on parallel lines and counting their unordered pairs. If the extent in the primitive lattice parameter is at most \(m\), every such pair is counted. Adjacent pairs give the extra half-step weight at \(j=1\). It also evaluates the full product relative entropy of the complete and cut pairs by their independently calibrated native matrix elements. The compensating term uses the full physical distant-pair counts.

On the compact regular interior and a bounded set of velocities, the projected exponential map of actual observer rods has a Taylor expansion \(F(y)=x_0+Ay+O(|y|^2)\), uniformly in \(x_0\). The equal-volume rod stop is \(r=R+O(R^3)\). Its projection differs in volume from \(x_0+A B_R\) by \(O(R^4)\), since the boundary displacement is \(O(R^2)\) and its area \(O(R^2)\). On a sufficiently small neighborhood the projection is one-to-one. Both projected boundaries are smooth with tubular volume \(O(R^2h)\) for \(h\ll R\), by rescaling them to smooth perturbations of a unit ellipsoid. Changing the region changes any overlap by at most twice the symmetric-difference volume. Only separations \(d=O(R)\) contribute. Hence using tangent ellipsoids changes the integrated overlap by \(O(R^5)\). This estimate controls the true proper-ball projections uniformly over the interior centers. The numerical replay below uses tangent ellipsoids. [Proper clocks, projections and volume stop](../PROOF.md)

Digitizing any overlap costs \(O(R^2h)\) in volume since its two boundary layers lie in the layers of the original regions or their translates. Its continuous overlap as a function of separation is Lipschitz with constant \(O(R^2)\), by the same translation bound. The separation quadrature, weighted lattice counts and first half-step therefore give an aggregate \(O(R^3h)\) error, valid uniformly in center. For the ten enrolled frames this proves
\[
 s_R=(8\pi^2/15)R^4T_{uu}+O(R^5+R^3h),\qquad T_{uu}=\rho\gamma^2.
\]
With the same homogeneous rod and area law, \(a_R=-(4\pi/15)R^4G_{uu}+O(R^6)\), and the independent shell equation gives \(G_{uu}=8\pi G T_{uu}\). Thus the corresponding generalized entropy sum is \(O(R^5+R^3h)\) at the calibrated coupling. Choosing \(h\le R_{\min}^2\) in fixed body units suffices on a simultaneous menu of radii with minimum \(R_{\min}\). Centers and radii can be varied after preparing the grid. A different cutoff has more weakly driven pairs. This constructs common limiting regional readings, without an exact embedding that preserves every microscopic pair across cutoffs.

The faithfulness of the blank bounds both endpoint and whole-pair relative entropy by \(C e^2\) at small \(e\), using the inverse function estimate \(|d|\le C'|e|\). For complete pairs the exact value is \((\log5)\sin^2 d/2\). At most \(O(R^3/h^3)\) copies per separation touch the ball for \(h\ll R\). There are \(O(1/h)\) short separations with \(|e|=O(h^4)\), and one compensated separation per direction with \(|e|=O(h^3)\). Bisection of each gate to absolute angle precision \(h^6\) contributes at most \(O(R^3h^2)\) to entropy, since its derivative is bounded and \(O(R^3/h^4)\) registers can be cut. Rational calibration bounds and arithmetic precision can be tightened at each cutoff. The total relative entropy is
\[
 0\le D(\rho_B\Vert\sigma_B)=O(R^3h^4+R^3h^3)=o(R^4)
\]
on the same refinement. Each finite preparation is a product across pair registers, so these are entropies of actual overlapping regional marginals. The numerical replay uses their calculated probabilities. Finite sampling is possible at each cutoff by counting each of the finitely many distinct endpoint exposures on copies. The binary entropy derivative is bounded near \(3/8\). If a survey has at most \(Q\) split pairs, simultaneous endpoint probability half-width \(R_{\min}^5/(C Q)\) bounds the entropy error by \(O(R_{\min}^5)\) after blank subtraction. A union bound with summable failure budgets over cutoffs gives joint coverage along nested finite menus. Such counts are not deposited in this benchmark.

## What changes when the flat blank is enrolled separately

The raw difference of absolute entropies on independently measured proper balls contains the extra term \(S(\sigma_{E_g})-S(\sigma_{E_0})\). It can dominate the entire signal. Write \(H_0=\sqrt{2G(M_o+M_c)}\) for the initial magnitude of the dust expansion. A comoving central observer's spacelike ray has \(t(s)=H_0s^2/2+O(s^4)\), \(dx/ds=a^{-2}=1+H_0^2s^2+O(s^4)\). The fixed-volume radius is \(r=R+G_{tt}R^3/45+O(R^5)=R+H_0^2R^3/15+O(R^5)\). Its projected comoving radius is therefore \(R+2H_0^2R^3/5+O(R^5)\). The separately enrolled flat ball is a proper subset at nonzero density, with extra comoving volume of order \(R^5\).

Adding either endpoint to a blank region increases its entropy by at least \(\min\{h_2(3/8),H_{\rm pair}-h_2(3/8)\}>0\), where \(H_{\rm pair}\) is the Shannon entropy of \((5,1,1,1)/8\). On \(h=o(R^3)\) the lattice boundary-layer error \(O(R^2h)\) is smaller than the excess volume. The shell therefore contains \(\Theta(R^5/h^3)\) additional sites with \(\Theta(h^{-1})\) registers each. The raw blank increase is \(\Theta(R^5/h^4)\). This explicit divergence prevents using the uncorrected, separately enrolled absolute entropy in the three regional estimates. The common-address difference is an operational reference choice.

At \(R=1/64\) the independent flat and held proper-clock surveys recover comoving radii \(0.0156250000\) and \(0.0156254393\). The translated comoving observer at distance \(0.0156252\) from an endpoint address includes it only on the held ball. The finite receipt gives the number of blank endpoint registers there and a lower bound on the increased blank entropy. Spatial homogeneity allows the same rest-ray survey about this center. [Clock records and finite contrast](results.json)

## An independent velocity test

For a finite static signed angular measure \(w\) with \(\int dw=1\), the factor is \(\int\gamma[1-(v\cdot n)^2]^{-1/2}\,dw(n)\). Matching \(\gamma^2\) for all sufficiently small velocities would require
\[
 \int nn^T\,dw=I
\]
by expanding to second order in \(v\). Taking the trace gives \(1=3\). This excludes universal Lorentz response for this class of comoving independent pairs, even allowing signed angular weights. The specified cubic array at the untrained \(u=(2,1,1,1)\) gives \(-0.461805\) against dust's \(4\). Hence fitting ten directions can give an incorrect prediction on an eleventh, despite exact selection of the geometrically correct ten-component dust tensor on the enrolled menu.

For instance, on the held stock \(\rho=19/(44\pi)\), at the displaced center \((81,-43,31)/1021\) and radius \(1/64\), the \(h=1/1024\) survey gives entropy divided by \((8\pi^2/15)\rho R^4\) of \(1.03470,\ 2.83360,\ 9.07798\) at rest, axially and along a face diagonal. The unenrolled velocity yields \(-0.42444\). At this center and radius the largest \(D/R^4\) over the eleven readings is below \(2.64\times10^{-7}\). The archive contains two radii, three centers and three spacings. Two of the centers have overlapping rest balls with neither contained in the other. All data are model-generated. The extension uses a prepared comoving triad, a stock-dependent calibration, and remotely paired registers with prescribed exposure. Deriving those operations from unmodified Interval Echo or identifying this entropy response with isotropic physical gravity would require separate evidence.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m model.candidate_laws.interval_echo_published_records.spherical_dust.homogeneous_pairs.experiment
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo_published_records/spherical_dust/homogeneous_pairs
~~~
