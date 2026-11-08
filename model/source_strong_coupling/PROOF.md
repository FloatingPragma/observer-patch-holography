# Color normalization and a complete finite coupling readout

## Source and scope

The source is the configured model in `native_quantum_gauss`, extended to
trees by `native_gauge_strings`. Its matter register is the declared
nonvacuum even exterior lift of the port-admitted current algebra. Its link
register contains the vacuum and all five spaces \(V_r^*\otimes V_r\), with
dimensions \(d_r=1,3,6,2,3\). The resulting link dimension is
\(1+\sum_r d_r^2=60\). The static dual reference, coherent Hamiltonian,
preparation and complete position detector belong to this configured model.
They are additional architecture beyond the unchanged repair rule.

The committed source specifications set the default hopping \(J=1\) and
electric coefficient \(\epsilon=1\). Their dimensionless ratio
\(\epsilon/|J|=1\) is exact, and observers allowed to read the simulator
can obtain it directly. The detector construction independently recovers
these specified values. Admitting other parameter choices in the same API
demonstrates freedom under gauge symmetry; it does not make the fixed
default source ambiguous. The physical identification problem concerns
the relation between this specified response and QCD.

No statement here identifies this finite model with the physical QCD
field. The finite algebra and its operational coefficient can be computed
without an empirical coupling input. Their relation to QCD requires the
source to select a spatial color action, state, physical clock and
renormalization map.

## Fundamental normalization on the same current

Let \(G_0,\ldots,G_7\) be the first eight anti-Hermitian matrices returned
by `native_matter_attachment.canonical_currents()`. They are admitted through
the actual twelve-port control map. Define

\[
 T_a=iG_a/2\quad(0\le a\le6),\qquad
 T_7=i(G_6+2G_7)/(2\sqrt3).
\]

On the color triplet these obey \(\operatorname{tr}_3 T_aT_b=\delta_{ab}/2\).
Their declared exterior lift satisfies
\(\sum_a T_a^2=C_3\), with the source-computed eigenvalues \(0,4/3\),
and \(\operatorname{Tr}_{15}T_aT_b=2\delta_{ab}\). Thus choosing the matter
trace instead of the fundamental trace multiplies the invariant metric by
four. The source's orthonormal six-dimensional trace frame has color
Casimir \(2C_3\). These are exact normalization identities, not values of
the strong coupling.

With the fundamental convention fixed, a curvature action can be written
\(S_E=(\kappa/2)\int\operatorname{tr}_3\mathcal F_{\mu\nu}
\mathcal F_{\mu\nu}\), where \(\mathcal A=g_s A\) and
\(\kappa=1/g_s^2\). Therefore \(\alpha_s=1/(4\pi\kappa)\) in that
declared convention. Ad-invariance permits every positive \(\kappa\).
It fixes the tensor shape and does not select its magnitude. Identifying a
finite coefficient with a renormalized coupling also requires a specified
scheme and energy scale. The definition and running convention here follow
the [PDG QCD review, equations 9.1--9.4](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-qcd.pdf);
no measured value is an input.

## Exact action on the full Gauss singlet blocks

Use the literal source states `reference_states(15, indices)` on a link.
The state \(A_r\) has its matter at the reference endpoint and the link in
vacuum. The state \(B_r\) has its matter at the other endpoint and an
occupied link carrying the compensating charge. Direct application of the
ambient source gives

\[
 \|A_r\|^2=d_r,\quad \|B_r\|^2=d_r^2,\quad
 H A_r=J B_r,\quad H B_r=d_rJ A_r+\epsilon B_r.                 \tag{1}
\]

All twelve actual Gauss generators at both endpoints annihilate both
states. These identities hold for every real \(J\) and every
\(\epsilon\ge0\). They follow from the literal raising link and its
adjoint, with electric energy \(\epsilon\) per occupied link. The exact
replay checks every state in all five blocks at three independent parameter
pairs; linearity in \(J,\epsilon\) proves the stated family.

In the normalized basis \(a_r=A_r/\sqrt{d_r}\), \(b_r=B_r/d_r\),

\[
 H_r=\begin{pmatrix}0&J\sqrt{d_r}\\J\sqrt{d_r}&\epsilon\end{pmatrix}.
                                                                    \tag{2}
\]

The same source position detector is \(P_0=|a_r\rangle\langle a_r|\),
\(P_1=|b_r\rangle\langle b_r|\). For \(U_t=e^{-itH_r}\) its complete
instrument is \(K_z=P_zU_t\), \(z=0,1\). Since
\(\sum_zK_z^*K_z=I\), the record isometry
\(\sum_z|z\rangle\otimes K_z\) preserves every input vector and every
untouched reference. Both branches and their backaction are retained.
Products of these Kraus operators define complete repeated histories.

Prepare the specified source state \(a_r\), keep the same clock parameter
\(t\), and put \(A=d_rJ^2\), \(\Omega^2=\epsilon^2+4A\). The exact
complete endpoint probability is

\[
 p(t)=\Pr(1|a_r,t)=\frac{4A}{\Omega^2}\sin^2(\Omega t/2).       \tag{3}
\]

Consequently \(p''(0)=2A\), \(p''''(0)=-2A\Omega^2\). For \(J\ne0\),
the records determine \(J^2=A/d_r\) and
\(\epsilon^2=-p''''(0)/p''(0)-2p''(0)\). The assumed
\(\epsilon\ge0\) fixes its sign. The sign of \(J\) is an endpoint phase
convention and cannot be identified from these position records.

## A finite-time error bound

The following construction needs two complete probabilities, rather than
unbounded numerical differentiation. Source specifications give
\(|J|\le J_{\max}\), \(0\le\epsilon\le e_{\max}\). Write
\(A_{\max}=d_rJ_{\max}^2\),
\(W_{\max}=e_{\max}^2+4A_{\max}\). Taylor's theorem and the bound on the
sixth derivative of \(\cos(\Omega t)\) give

\[
 |p(t)-At^2+Bt^4|\le A\Omega^4t^6/360,\qquad B=A\Omega^2/12.
                                                                    \tag{4}
\]

Thus the actual coefficients lie in outward intervals around

\[
 \widehat A=\frac{16p(t)-p(2t)}{12t^2},\quad
 \widehat B=\frac{4p(t)-p(2t)}{12t^4},
\]

with respective deterministic errors at most

\[
 E_A=\frac{A_{\max}W_{\max}^2t^4}{54},\qquad
 E_B=\frac{17A_{\max}W_{\max}^2t^2}{1080}.                    \tag{5}
\]

These constants include both times. Probability error intervals are
propagated before adding (5). If the resulting intervals certify
\(A_->0\) and \(B_-\ge0\), then

\[
 \Omega^2\in[12B_-/A_+,12B_+/A_-],\quad
 J^2\in[A_-/d_r,A_+/d_r],\quad
 \epsilon^2\in[12B_-/A_+-4A_+,12B_+/A_--4A_-]\cap[0,\infty).
                                                                    \tag{6}
\]

The program uses rational arithmetic throughout (4)--(6). Its witness
encloses (3) using an alternating cosine series with an explicit omitted
term bound. It does not label binary64 differences as error certificates.
Input probabilities can include a proved implementation error and a
statistical confidence interval. For independent identical preparations,
Hoeffding plus a union bound gives both frequency errors at most \(\eta\)
with probability at least \(1-\delta\) using
\(N\ge\log(4/\delta)/(2\eta^2)\) trials at each time. Such a statistical
claim requires those preparation assumptions. Taking \(t\to0\) and
probability errors \(o(t^4)\) contracts the coefficient intervals.

The source clock determines the units of \(J,\epsilon\). Replacing
\((J,\epsilon,t)\) by \((sJ,s\epsilon,t/s)\), \(s>0\), leaves every
history unchanged. The position-identifiable ratio \(\epsilon/|J|\)
survives this clock change.
Changing that ratio changes the oscillation contrast in (3), so its freedom
is not a common clock convention.

## The actual electric law excludes a Casimir-coupling substitution

At the configured \(\epsilon=1\), the five nonvacuum link sectors all have
electric energy one and the link vacuum has energy zero. Their exact
\((C_3,C_2,Y^2)\) rows, in source central-label order, are

\[
 (0,0,1),\quad(4/3,0,1/9),\quad(4/3,3/4,1/36),\quad
 (0,3/4,1/4),\quad(4/3,0,4/9).                               \tag{7}
\]

There are no real coefficients \(a,b,c\) whose quadratic Casimir energy
\(aC_3+bC_2+cY^2\) equals one on all these sectors. Equality on the second
and fifth rows forces \(c=0\); the first row then has zero energy,
contradicting one. A common identity term is also excluded by the
zero-energy vacuum. A separate nonvacuum-projector coefficient would
reintroduce the source's specified electric law and is not a Laplacian
coefficient. The exact replay also supplies a left-null witness for (7).

This excludes a direct identification of the configured electric law with
a pure color heat kernel or even a three-factor quadratic Casimir law.
It does not exclude a different emergent color action under a justified
coarse-graining. Such an action would require an additional derivation.
Likewise, the finite heat-kernel theorems explicitly assume the Laplacian
generator and its one-sided or two-sided trace convention; Gauss invariance
of this source does not provide those premises.

## Perturbative matching and the boundary condition

The inherited QCD producer supplies a decimal \(\alpha_3\), the D10
branch, a display clock and external heavy-quark thresholds. Its
perturbative evolution is conditional on these inputs. The beta function
cannot select its own integration constant: at one loop,
\(\alpha_s(\mu)^{-1}=\alpha_s(\mu_0)^{-1}+2b_0\log(\mu/\mu_0)\),
where \(b_0=(33-2n_f)/(12\pi)\) in the stated PDG convention.
Distinct positive input couplings remain distinct wherever both solutions
exist.

At \(\mu=m_h(m_h)\), the high-flavor to low-flavor matching is

\[
 \alpha_l=\alpha_h\left[1+\frac{11}{72}(\alpha_h/\pi)^2
 +\left(\frac{564731}{124416}-\frac{82043\zeta(3)}{27648}
       -\frac{2633n_l}{31104}\right)(\alpha_h/\pi)^3\right].    \tag{8}
\]

The inverse relation has the opposite displayed signs through this order.
The direct three-loop coefficients are given in
[Chetyrkin, Kniehl and Steinhauser, equation 23](https://arxiv.org/pdf/hep-ph/9708255).
Equation (8) and scale invariance of the integrated four-loop beta function
are checked independently. Perturbative truncation and threshold-bracket
variation are distinct uncertainties; a bracket over input thresholds
does not bound all orders of QCD or the missing source-to-field matching.
