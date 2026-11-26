# Seam correlation energy and its curvature under a supplied coupling

A four-word sector of the 99-state Interval Echo seam stores eight units of energy entirely in a change of correlation. Measuring each register separately gives exactly the same probabilities in both preparations. In a homogeneous gravitational model that couples to the measured seam Hamiltonian, the change in volume-weighted curvature is exactly eight times the coupling calibrated on ordinary loads. In the free Interval Echo law there are also correlations that change a modular reading while every permitted material probe gives the same records. The following calculations identify the physical premise that distinguishes the two cases. [Protocol](PROTOCOL.md)

## Which correlation measures energy?

Fix the physical Hamiltonian \(h\) on a finite tensor product and choose the faithful reference \(\sigma=e^{-\beta h}/Z\), with \(\beta>0\) measured in inverse energy units. Define \(K=-\log\sigma+\log(\sigma_A\otimes\sigma_B)\). For \(D=\rho_1-\rho_2\) with both partial traces zero,
\[
 \operatorname{Tr}DK=\beta\,\operatorname{Tr}Dh .
\]
Indeed the second logarithm is a sum of one-register operators, whose contractions with \(D\) vanish. The normalization of the joint Gibbs state adds a constant, whose contraction also vanishes. This proves the identity for coherent as well as classical preparations. A general correlated reference need not select the physical energy.

On an isolated native contact, use memory magnitudes 0 and 2 with even reversal parity and positive coin phase as rows, and loads \((2,0),(1,1)\) as columns. Their four tensor products are exact eigenwords. The table of native energies in row order is \(((0,0),(32,0))\). Reversing the coin phase gives \(((0,0),(32,32))\). On blank memory the vacuum and symmetric one-load word independently have energies 0 and 8. The integer test acts on the full native contact, including leakage out of every listed input. [Exact column checks](protected.py)

Take \(P_s=(1+sD_2)/4\) where \(D_2\) has diagonal \(+1\) and off-diagonal \(-1\). High uses \(s=-1/2\), low uses \(s=+1/2\). Both marginals are uniform. With these exact preparations the energy difference is the logarithm of the reference odds ratio, \(\log(\sigma_{00}\sigma_{11}/(\sigma_{01}\sigma_{10}))/(4\beta)\). It is 8 in the interacting sector and 0 with reversed coin. The source means are 12 and 4, or 16 and 16. The Gibbs reference is prepared on the specified four-dimensional invariant sector. Thermalization, a local vacuum boost and a physical regional collar normalization are additional physical identifications.

## How the geometry comparison works

The comparison places the seam energy \(E=\operatorname{Tr}(\rho h)\) uniformly in a comoving cell and supplies a semiclassical ADM geometry. One metric evolves for the expectation in a preparation. Averaging metrics over separately gravitating eigenword trials would implement a different prescription. Write \(q=\log a,\ W=V_0e^{3q}\); the expanding, flat branch has positive volume and unit lapse. Its Hamiltonian is
\[
 {\cal H}=-\frac{\kappa p^2}{12W}+\frac{\Lambda W}{\kappa}+E .
\]
The internal energy is independent of \(W\), hence its pressure is zero. The diagonal ensembles are stationary under the seam Hamiltonian. Hamilton's equations give \(\dot q=-\kappa p/(6W)\) and \(\dot p=-\kappa p^2/(4W)-3\Lambda W/\kappa\). Autonomy conserves the initial constraint \({\cal H}=0\). Substitution gives \(3\dot q^2=\Lambda+\kappa E/W\) and \(2\ddot q=\Lambda-3\dot q^2\). In the orthonormal comoving frame \(G_{00}=3\dot q^2,\ G_{ii}=-\Lambda,\ G_{0i}=0\). For any null direction \(k=(1,\mathbf n)\) with \(|\mathbf n|=1\), at every regular time,
\[
 \Delta[W(G_{00}-\Lambda)]=\Delta[W G_{kk}]
     =\kappa\,\Delta\langle K\rangle/\beta .
\]
Multiplication by each preparation's own measured volume permits comparison at common time even as their volumes diverge. At equal volume the local density change is \(\Delta\langle K\rangle/(\beta W)\). The correlation is part of the total Hamiltonian energy. Adding it to the total a second time would count the same energy twice. [Machine-checked scalar identities](../../../Lean/ObserverPatchHolography/EinsteinBranch/OverlapCoupling.lean)

A different observer reads two-way returns on a proper clock. Vacuum and an ordinary one-load source at the first appointment fix the background and coupling. The source meter uses a separate two-level thermometer and joint reference and preparation counts. At the finest spacing \(0.01\) the ordinary calibration gives \(\kappa=1.60000105\). Each count bank contains 4,194,304 simulated trials.

| Held source | Exact energy difference | Measured energy difference, one counting SE | Measured volume-weighted temporal curvature difference at three later times |
|---|---:|---:|---:|
| Interacting words | 8 | \(7.95234\pm0.02032\) | 12.799997, 12.799988, 12.799973 |
| Reversed coin | 0 | \(-0.00150\pm0.01203\) | 0, 0, 0 |
| Interacting words with gravitational charge held constant in this sector | 8 | Same counts as interacting words | 0, 0, 0 |

The finite-count forecast for the first row is \(12.72375\pm0.03251\) from source counting alone, about 2.35 standard errors below the radar result. The population forecast is 12.8. The three axial null estimates range from 12.79985 to 12.80041. The final control uses the same matter Hamiltonian, ordinary coupling, reference and counts, assigning the uniform mean energy 8 to gravity on each of the four words. Ordinary calibration cannot distinguish these two gravitational laws. The held test can. [Original counts and optical tapes](results.json) · [Forward source, meter, and radar comparison](measure.py)

## What an independent native geometry must measure

At a faithful state of a finite regional algebra, let \(Q\) represent the first variation of an independently read geometric statistic, in the units of \(W(G_{00}-\Lambda)\). Thus its derivative on a tangent \(X\) is \(\operatorname{Tr}(QX)\). For a differentiable statistic of record probabilities this is the pullback of its ordinary derivative through the measurement. At a fixed vacuum and calibration the exact criterion on **all** fixed-marginal perturbations is
\[
 Q-\kappa h=A\otimes I+I\otimes B .
\]
To prove this, decompose the Hermitian operators on each factor into the identity and traceless operators. The tangents with both partial traces zero are exactly the span of products of two traceless operators. At a faithful state each such tangent is realized by small positive and negative perturbations. Their orthogonal complement for the trace inner product consists of sums of local terms. This proves necessity and sufficiency. If the calibration or volume varies, its derivative belongs in the geometric statistic before applying this test. For a smaller admitted family, only orthogonality to its actual tangents is required.

The same criterion can test a candidate physical reference without assuming thermalization. For an independently prepared faithful \(\sigma\), its modular contrast measures \(h\) on every fixed-marginal tangent exactly when \(-\log\sigma-\beta h\) is a sum of local terms. A finite-precision version is available for either comparison. For the residual \(R\), the largest absolute error on tangents with trace norm at most one is
\[
 \sup_{\operatorname{Tr}_A X=\operatorname{Tr}_B X=0,\ \|X\|_1\leq1}
 |\operatorname{Tr}(RX)| =
 \inf_{A,B}\|R-A\otimes I-I\otimes B\|_{\rm op}.
\]
Each local sum annihilates all tangents, which proves the upper bound. For equality, extend the real linear functional on the Hermitian tangent subspace to all Hermitian operators with the same trace-norm dual norm using finite-dimensional Hahn–Banach. The representing Hermitian matrix differs from \(R\) by a local sum and has operator norm equal to that dual norm. A uniform bound on the right therefore controls every admitted correlation on that sector at once. For a raw cap statistic divided by \(\ell^4\), a bound \(o(\ell^4)\) before division suffices. On the diagonal four-word sector the tangent of trace norm one is the signed checkerboard divided by four, up to sign at the two endpoints. The absolute measured mixed residual is the exact maximum error.

On the diagonal four-word sector only one independent direction survives. In row order the contrast is \((Q_{01}+Q_{10}-Q_{00}-Q_{11})/4\); the target is \(8\kappa\). Ordinary vacuum and one-load calibration fixes the slope but puts no constraint on this mixed difference. A native proof needs a clock and proper-volume ruler that determine \(Q\) independently of the energy meter. For the proposed dark term the physical regional generator must also yield the claimed fourth-order coefficient on these same balls, with ordinary contributions accounted for and the required frame, pressure and conservation information. A global four-word thermal meter establishes only the energy identification on its declared sector.

## The direct test on a spatial collar

At a common proper volume \(4\pi\ell^3/3\), let \(\Delta A_\ell\) be the signed difference between two rest-frame rim areas, \(\Delta E_{C,\ell}\) their independently measured dimensionless regional collar charge, and \(\Delta T^{\rm ord}_{uu}\) the ordinary density difference. The Einstein paper's geodesic-ball estimator gives \(\Delta\widehat G_{uu}=-15\Delta A_\ell/(4\pi\ell^4)\); the dark-sector definition gives \(\Delta\widehat T^{\rm anom}_{uu}=15\Delta E_{C,\ell}/(8\pi^2\ell^4)\). Subtracting the calibrated sources from the area estimate gives exactly
\[
 \Delta\widehat G_{uu}-\kappa(\Delta T^{\rm ord}_{uu}+\Delta\widehat T^{\rm anom}_{uu})
 =-\frac{15}{4\pi\ell^4}\left[
   \Delta A_\ell+\frac{\kappa}{2\pi}\Delta E_{C,\ell}
   +\frac{4\pi\kappa}{15}\ell^4\Delta T^{\rm ord}_{uu}\right].
\]
Thus, when ordinary stress is matched, a positive regional charge must reduce the equal-volume rim area by \(\kappa\Delta E_{C,\ell}/(2\pi)=4G_N\Delta E_{C,\ell}\) to leading order in natural units. An error \(o(\ell^4)\) in that area comparison is necessary and sufficient for convergence of the *estimated temporal residual* to zero. Identifying it with physical curvature also requires the geodesic-ball expansion on the same frame and control of its error. Other frames and tensor conservation are needed for the full equation. The thermal meter on a homogeneous cell does not determine the radius-dependent regional charge in this formula. [Finite coefficient identity](../../../Lean/ObserverPatchHolography/EinsteinBranch/OverlapCoupling.lean)

In the homogeneous comparison the two preparations pass through the same cell volume \(V_0=(2\pi)^3\). There a calorimeter gives \(\Delta\rho=\Delta\langle K\rangle/(\beta V_0)\). Weight that density with the leading diamond kernel, whose integral is \(8\pi^2\ell^4/15\). This assigns a regional charge from a homogeneous energy density. A distinct instrument shoots spacelike geodesics from the observer, adjusts the aperture to match the interior proper volume and measures the rim with separate angular rods. Using the population energy 8, the couplings inferred from the rims alone are 1.62320, 1.60575 and 1.60143 at radii 0.4, 0.2 and 0.1. The ordinary radar calibration was 1.60000105. At the smallest radius the measured area change is \(-4.32691\times10^{-6}\). The finite-count leading forecast is \((-4.29729\pm0.01098)\times10^{-6}\) from source counting alone. The finite-radius geometric correction and source counting error are distinct. In the constant-gravitational-charge control both source preparations give identical rim records. The survey tests the sign and magnitude of an independently read area on the supplied metric. The regional weight is prescribed here; identifying it with a native spatial collar modular operator requires its own measurement or proof. [Rim and interior readings](rim.json) · [Geodesic aperture and regional energy comparison](rim.py)

## An exact invisible correlation on the native law

Take two private seams, each in the three negative-coin memory modes of energies \(0,32,33\). On every one of the nine endpoint-load inputs the local Hamiltonian restricted to mode \(i\) is \(E_i I\). The other seam generators and the port degree correction act independently of the marked labels. For fixed labels \(i,j\), any evolution generated by the free law or classically scheduled applications of these Hamiltonians differs on the remaining system only by a scalar phase. Label-independent material preparations and instruments therefore give identical outcome laws. Conditioning on the first material result and iterating proves the same statement for adaptive classical schedules, every later result, joint material probes and any function of their transcripts. Direct marked-memory measurements, a coherent control of service histories or interventions mixing the marked words lie outside this theorem.

Let \(T_{01}=T_{12}=T_{20}=1,\ T_{ji}=-T_{ij},\ T_{ii}=0\) and \(P_{ij}=(1+tT_{ij})/9,\ |t|<1\). All marginal distributions are uniform. Summing opposite ordered pairs also fixes the complete histogram of \(E_i+E_j\). A joint memory measurement distinguishes \(t\). If the reader pairs the two labels with probability \(f\) and otherwise takes independently drawn labels without exposing the pairing bit, it reads \(P_{ft}\). At the fixed reference \(t_0=1/2\), the derivative of the expectation of \(-\log(P_{ft_0}/P_0)\) is
\[
 -\frac{f}{3}\log\frac{1+ft_0}{1-ft_0}.
\]
It is \(-\log(3)/3\) at full inclusion. Even granting the regional law \(f=a\ell^2+o(\ell^2)\), multiplication by \(15/(8\pi^2\ell^4)\) gives the nonzero limit \(-5t_0a^2/(4\pi^2)\) for \(a>0\). The native material transcripts have exactly zero contrast at every size. This excludes interpreting this particular signed modular change as a temporal stress at nonzero coupling to a geometry obtained exclusively from those transcripts when ordinary stress and vacuum are fixed. It neither identifies a physical proper radius for the protected source nor rules out a different reference, stress, or metric instrument. The physical Gibbs state of the restricted additive two-seam Hamiltonian factorizes and gives zero interaction contrast. [Exact source and regional slopes](protected.json)

A genuine Gibbs-energy contrast can also be invisible. The protected memory contains two orthogonal zero-energy coin words at magnitude zero and both negative-coin reversal parities at magnitude two. Its four modes are a product of two logical registers, the magnitude bit and the within-pair label, with energy rows \((0,0),(33,32)\). The same high and low preparations have equal logical marginals and differ by exactly \(1/4\) unit of physical Hamiltonian energy. A thermal reference measures precisely that difference. The scalar-action proof on all nine loads still applies. On \(N\) such marked seams the total energy difference is \(N/4\) and all material transcripts coincide, including collective readings. If the same measured volume in both arms grows as \(Nv+o(N)\), the limiting energy density difference is \(1/(4v)\) for \(v>0\). No native proper-volume law is assumed in the finite indistinguishability proof. This encoding concerns two bits inside one memory; the two-seam example concerns a spatially separated pair with an explicitly different choice of reference. [Native columns, finite Gibbs calculation and energy histogram](protected.py)

Reachability from blank memory depends on the allowed interventions. Free evolution preserves the completely blank sector. A write on endpoint loads \((1,1)\) prepares the even-coin, zero-energy mode of the four-word example exactly. A write on \((2,0)\) gives equal coin amplitudes with contrast \(+2\) and balanced endpoint loads. Each protected odd-coin reversal parity has zero overlap with that pure write. A sharp coin read with its outcome discarded dephases the two terms and gives probability \(1/4\) in each of the protected energy-32 and energy-33 modes. This uses a direct instrument on the marked memory. On the odd-total input \((1,0)\) the two coin branches instead have orthogonal output loads; a coin-parity read after the write has odd probability \(1/2\). Thus a one-write image, a freely invariant sector and the support reachable under a specified instrument set are separate claims. Nonzero support alone does not prepare the two equal-marginal Gibbs ensembles used above. The finite-density formation experiment starts with two armed memories. [Exact integer checks](test_overlap.py) · [Formation preparation](../../../model/experiments/interval_echo_einstein_station/mobile_collar_source/README.md)

From the repository root the following commands check the saved counts, native column identities and both gravity arms.

~~~sh
PYTHONPATH=code python -m matched_einstein_limit.seam_overlap.measure
PYTHONPATH=code python -m matched_einstein_limit.seam_overlap.protected
PYTHONPATH=code python -m matched_einstein_limit.seam_overlap.rim
PYTHONPATH=code python -m pytest -q code/matched_einstein_limit/seam_overlap
~~~
