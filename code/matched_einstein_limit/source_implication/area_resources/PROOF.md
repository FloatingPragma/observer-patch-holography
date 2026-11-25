# Slope spread and the number of independent collar inputs

## Material comparison and admissible sources

Take the formed excess collar on the 80-cell federation. Its independent-head reference is \(\sigma\), the product of the two head marginals. Set \(D=\rho_{\rm paired}-\sigma\). An admitted fraction \(f=f(\ell)>0\) and an interior pairing strength \(s\in(0,1)\) give the unflagged density \(\sigma+s fD\). The source also supplies an independent comparison head and a fair tie bit in a fixed state. Keeping these auxiliary inputs independent of \(s\) preserves affinity. In the native two-leg return at time 2, the completed-return probability shift between the paired and independent preparations is \(z=-0.014370365536930434\). The registered zero-contrast probe gives a shift below \(10^{-17}\). The bit \(b=\pm1\) published after the comparison has probabilities \((1+b s f z)/2\). This uses the hit- and preparation-label-erased output. The supplied flat chart has \(f/\ell^2\to4\pi/3\) by its growing bank of placement marks. The theorem below permits any \(f/\ell^2\to a>0\); interpreting \(\ell\) as a common proper radius requires a separate measurement. [Return and admission](../../collar_identifiability/refinement/admission.py) · [Source and chart construction](../../collar_identifiability/refinement/README.md)

For fixed \(\sigma,D\), a symmetric chord \(\sigma\pm dD\succeq0\) exists whenever the reference is faithful on the support of both states. For the product of two marginals this support condition holds: a zero one-head probability implies zero support on that subspace in the positive joint state. On their support one may use any positive \(d<1/\|\sigma^{-1/2}D\sigma^{-1/2}\|_\infty\), with arbitrary positive \(d\) if \(D=0\). The numerical native input gives reciprocal norm \(0.155333476\), so \(d=1/8\) lies inside the measured chord. The endpoint states serve only as a mathematical decomposition; the hardware need not prepare them. The same positive \(d\) is required throughout the refinement for the full-input bound. On the published bit alone \(d=1/|z|\) works exactly when \(z\ne0\).

## One input and two endpoint slopes

Write \(k=fz\). The comparison-bit entropy in nats obeys
\[
 S'(s)=-k\,\operatorname{atanh}(sk),\qquad
 -S''(s)=\frac{k^2}{1-s^2 k^2}\ge k^2 .
\]
Thus, for fixed \(s_1<s_2\) with \(|sk|<1\), the slope spread is at least \((s_2-s_1) f^2 z^2\). The mean value theorem gives this estimate for either sign of the native contrast.

A fixed quantum procedure, including its internal state, adaptivity and classical postprocessing, gives an expected number linear in the input density. With one copy, the area mean \(A_\ell(s)\) is therefore affine in \(s\). Its derivative has the same value at the two endpoints. More generally let \(b\) bound \(|\eta A'_\ell(s_2)-\eta A'_\ell(s_1)|\) at a common fixed calibration \(\eta\). The endpoint residuals \(r_i=S'(s_i)+\eta A'_\ell(s_i)\) satisfy
\[
 \max_i |r_i|\ \ge\ \tfrac12\big[(s_2-s_1)f^2z^2-b\big].
\]
Subtract the two residuals and use the triangle inequality. For an exactly affine area \(b=0\). For \(A_\ell(s)=\ell^2F_\ell(sf)\) and \(|F_\ell''|\le H\) on the tested interval with \(H\) independent of \(\ell\), one has \(b\le|\eta|(s_2-s_1)\ell^2 f^2H=o(\ell^4)\). At \(z\ne0\) at least one normalized residual \(|r_i|/\ell^4\) is then bounded away from zero. This compares source strengths as well as orientations at a given strength. [Lean proof of the slope spread and residual inequality](../../../../Lean/ObserverPatchHolography/EinsteinBranch/AreaReadoutSpread.lean)

## Any bounded collective readout

Allow \(m=m_\ell\) independent copies of the unflagged input \(\sigma+xD\), \(x=sf\). Allow any joint operation, measurements and adaptive control, using \(s\)-independent apparatus states. Assign a number \(\ell^2Y\) to each complete outcome with \(|Y|\le M\) uniformly. A deterministic offset may first be subtracted from every outcome. Stops and failures are included; a variable number of copies up to \(m\) can be padded by discarding unused copies. The mean normalized read \(F(x)=\mathbb E_xY\) satisfies
\[
 |F''(x)|\le\frac{M\sqrt{2m(m-1)}}{d^2-x^2},\qquad |x|<d.
\]

**Proof.** Write each input as a mixture of \(\sigma+\epsilon_i dD\), \(\epsilon_i=\pm1\), with independent weights \(p_x(\epsilon_i)=(1+\epsilon_i x/d)/2\). Conditional on all signs, the whole procedure has a fixed mean \(y(\epsilon_1,\ldots,\epsilon_m)\in[-M,M]\). Only the product of the sign probabilities depends on \(x\). Its first logarithmic scores are \(u_i=\epsilon_i/(d+\epsilon_i x)\), with \(\mathbb E_xu_i=0,\ \mathbb E_xu_i^2=1/(d^2-x^2)\). Since each probability is affine, the second derivative of the product divided by the product is \(G=2\sum_{i<j}u_iu_j\). Distinct unordered pairs are orthogonal in \(L^2(p_x^{\otimes m})\), including pairs that share one index. Thus \(\mathbb E G^2=2m(m-1)/(d^2-x^2)^2\). Cauchy-Schwarz applied to \(F''=\mathbb E[yG]\) gives the bound. In fact \(\mathbb E G=0\), so one can replace \(M\) by the standard deviation of the normalized output at that \(x\), using the joint probability law of the signs and the measurement. A uniform variance bound is sufficient even for unbounded outcome labels when the conditional second moments are finite. This proof allows entangled measurements on the independent inputs.

With \(x_*=\max_i|s_i|f<d\), it follows that
\[
 b\le|\eta|(s_2-s_1)\ell^2 f^2
       \frac{M\sqrt{2m(m-1)}}{d^2-x_*^2}.
\]
For \(m=o(\ell^{-2})\), the right side is \(o(\ell^4)\). If both residuals are \(o(\ell^4)\), the same inequality even requires \(\liminf_{\ell\to0}\ell^2m\ge z^2d^2/(\sqrt{2}|\eta|M)\) for nonzero \(\eta\) and \(M\). For observations restricted to the balanced bit the factor \(z^2d^2\) can be set to 1. Uniform bounds on the chord and outcome scale are material assumptions. A recorded source flag, a state-dependent choice of apparatus, or conditioning by an \(s\)-dependent success denominator changes the input or observable and needs its own calculation.

## The copy rate can be attained for the scalar

Use two independent banks, each containing an odd number \(n\) of balanced material bits. Write \(J_1,J_2\in\{-1,1\}\) for their majority signs. At bit bias \(y=sfz\), let \(g_n(y)=\mathbb E_y J_1\). Differentiating the binomial tail and canceling adjacent terms gives
\[
 g_n'(y)=b_n(1-y^2)^{(n-1)/2},\qquad
 b_n=n{n-1\choose(n-1)/2}2^{-(n-1)}
       =\mathbb E_0\left|\sum_{i=1}^n b_i\right|.
\]
The last equality also follows by differentiating the product probability at zero. Stirling gives \(b_n^2=(2n/\pi)(1+O(n^{-1}))\). Since \(g_n(0)=0\), throughout \(n y^2=o(1)\) the estimates \(g_n(y)=b_n y(1+O(n y^2))\), \(g_n'(y)=b_n(1+O(n y^2))\) hold uniformly.

For fixed \(C>0,\eta>0\), choose an odd \(n=C/\ell^2+O(1)\) and the **designed scalar calibration** \(g=\pi/(4\eta C)>0\). Give each complete output the nonnegative area score \(g\ell^2(1+J_1J_2)\). Its mean and derivative are
\[
 A_\ell(s)=g\ell^2(1+g_n(y)^2),\qquad
 \eta A'_\ell(s)=2\eta g\ell^2 g_n(y)g_n'(y)fz .
\]
For bounded \(s\) one has \(ny^2=O(\ell^2)\). Therefore \(\eta A'_\ell(s)=s f^2z^2[1+O(\ell^2)]\) and \(S'(s)=-s f^2z^2[1+O(\ell^4)]\). These estimates give \(S'(s)+\eta A'_\ell(s)=O(\ell^6)\) uniformly on the tested source interval. Two banks of total size \(2n\) attain the order in the necessity statement. The positive offset keeps both possible area outcomes in \([0,2g\ell^2]\); its derivative is zero. [Exact binomial readout](response.py)

With \(C=\eta=1\) and endpoints \(s=1/4,3/4\), the native increasing-bank admission gives the following numbers. In the one-copy column the best possible common slope has at least half the displayed entropy spread. The limiting normalized entropy spread is \(0.00181168569\).

| \(\ell\) | Independent copies | Entropy slope spread divided by \(\ell^4\) | Larger endpoint residual divided by \(\ell^4\) for the majority product |
| ---: | ---: | ---: | ---: |
| 1/8 | 130 | 0.00165619 | \(5.8250\times10^{-5}\) |
| 1/16 | 514 | 0.00176760 | \(1.5527\times10^{-5}\) |
| 1/32 | 2,050 | 0.00179945 | \(3.9506\times10^{-6}\) |
| 1/64 | 8,194 | 0.00180937 | \(9.9304\times10^{-7}\) |

The construction concerns the mean of a repeated preparation at fixed \(s\). Its gain was chosen to counter the bit entropy, and no native proper-area reading was used to fix it. A single extended body whose global density is affine in \(s\) obeys the one-input theorem even if it contains many internal ports. A single body could supply independent changing subrecords if that product dependence and their physical aggregation are established. The collar data here do not supply a volume-matching variation, independent area calibration, null rays, the full ball entropy, or the matter coupling. The resource obstruction applies to the putative cancellation of this comparison-bit slope with all other varied contributions held fixed or explicitly included in the bounded area term. [Finite checks](test_area_resources.py)
