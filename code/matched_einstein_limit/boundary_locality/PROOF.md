# Connected support of the regional logarithm

## Product reference and finite-order theorem

Let \(V\) be a finite graph with one finite-dimensional Hilbert space per vertex, and put \(H=\sum_X h_X\) on their tensor product. An interaction label includes its support \(X\). Multiple labels may share a support. Write \(R\subset V\), \(O=V\setminus R\), and let \(E=\operatorname{Tr}_O/\dim\mathcal H_O\) be the unital partial trace. Define
\[
 F(\boldsymbol z)=E\exp\left(-\sum_X z_Xh_X\right),\quad
 K=-\log F,\quad D=K-\sum_{X\subset R}z_Xh_X.
\]
The logarithm is analytic near zero because \(F(0)=I\). At physical couplings \(z_X=\beta\), the modular operator of the normalized regional density differs from \(K\) by a scalar. A trace-zero preparation derivative has the same response to both.

**Theorem.** For a multivariate Taylor coefficient of \(D\) of positive total degree \(m\), discard its scalar part. Any surviving monomial has a connected interaction-overlap graph and at least one interaction crossing the cut. Its operator is supported on the intersection with \(R\) of the union of its interactions. For graph edges and one-site terms, every vertex in this support has graph distance at most \(m\) to \(O\). If every crossing term has zero normalized partial trace on its exterior endpoint, the bound improves to \(m-1\). These bounds include repeated interaction labels. If \(O\) is empty, \(D=0\).

**Proof.** Set every unused coupling to zero. The coefficient is unchanged. If the remaining labels have disconnected overlap graph, partition them into two nonempty collections acting on disjoint tensor factors. Both the exponential and the normalized partial trace factor across these collections. Their logarithm is a sum of the two embedded logarithms. For real parameters this follows by diagonalizing the positive factors. Analyticity near zero gives the same Taylor coefficients. A monomial using a parameter from each collection thus has zero coefficient. This argument also proves support on the union, since all unused tensor factors carry identities. If every interaction is in \(R\), then \(K\) is exactly the linear internal Hamiltonian. If all are in \(O\), \(K\) is scalar. A connected collection reaching both sides must cross the cut. The algebraic cancellation of mixed finite differences with repeated directions is kernel checked in [MeanForceBoundary.lean](../../../Lean/ObserverPatchHolography/EinsteinBranch/MeanForceBoundary.lean). Tensor factorization and the passage to Taylor coefficients use the argument above.

For edges, a path through the cluster from an interior vertex to \(O\) has at most \(m\) edges. If the ambient distance is exactly \(m\), this path uses all \(m\) occurrences in the coefficient. Each edge occurs once, there are no one-site insertions, and only one crossing edge touches the exterior. Expand the logarithm through degree \(m\) in products of moments \(E(h_{i_1}\cdots h_{i_s})\). In the moment containing that crossing edge, every other factor acts only in \(R\). The bimodule identity \(E(ACB)=A E(C) B\) for internal \(A,B\) makes the moment zero. This excludes the distance-\(m\) case. Summing at common \(\beta\) preserves the support bound.

In a regular native binary carrier of degree \(d\), the term \(\kappa(6-d)\sum_i n_i\) is canceled by the declared chemical potential in the grand ensemble. Up to a scalar, the remaining contact is
\[
 h_{ij}=-\tfrac14[X_iX_j+Y_iY_j+(1+2\kappa)Z_iZ_j].
\]
Its one-end partial trace vanishes. On nonregular carriers any surviving one-site fields must be retained. The finite-order theorem also allows these fields. The native calculation here uses the binary blank-memory grand ensemble at the regular cell. Projecting first to a fixed total occupation changes the tensor-product partial-trace argument. [Native sector law](../../../model/candidate_laws/contrast_memory_spherical_v1/sectors.py) · [Grand ensemble and reduction](../../../model/experiments/seam_program/collective_rate/measure.py)

## An exact coefficient two contacts from the exterior

Take \(O=\{0\}\). The native cell contains edges \(0\,1\) and \(1\,8\), with no \(0\,8\) edge. Write \(C=h_{01}\), \(B=h_{18}\), and \(J_x=J_y=-1/4,\ J_z=-(1+2\kappa)/4\). In the subproblem with these two contacts, the coefficient of \(x^2y\) in \(-\log E\exp(-xC-yB)\) is
\[
 \frac16 E(CBC)-\frac1{12}\{E(C^2),B\}
 =-\frac1{12}E[C,[C,B]]
 =\sum_{j=x,y,z}\frac{J_j}{3}
   \left(J_j^2-\sum_i J_i^2\right)\sigma_1^j\sigma_8^j .
\]
For the first equality in this calculation, the cubic exponential contributes \(E(C^2B+CBC+BC^2)/6\), and the cross term in the logarithm contributes \(-\{E(C^2),B\}/4\). The last equality uses conjugation of one Pauli operator by another and the exterior Pauli trace. Terms containing only the cut contact have identity at port 8. Terms containing exactly one \(C\) give zero. The coefficient therefore gives the full cubic response on the measured pair for the two selected contacts.

The same pair projection can be calculated with **all 30 contacts** present. Write \(A\) for the internal Hamiltonian and \(C=\sum_j \sigma_0^j L_j\) for all cut contacts, with \(L_j=J_j\sum_{a\sim0}\sigma_a^j\). The cubic coefficient of \(K\) is
\[
  -\frac1{12}\sum_j[L_j,[L_j,A]]+\frac16 E(C^3).
\]
The last term is supported on the neighbours of 0, which omit port 8. Each commutator with a one-site sum acts on the support of its original term in \(A\). Hence only \(h_{18}\) can contribute on both ports 1 and 8. Among the sites in \(L_j\), only port 1 acts on that term. The two-contact answer survives at equal physical couplings on the complete cell.

| Coefficient in \(K\) on \(1\,8\), ordered \(XX,YY,ZZ\) | Native \(\kappa=8\) | Changed contact \(\kappa=4\) |
| --- | --- | --- |
| \(\beta^1\) | \((-1,-1,-17)/4\) | \((-1,-1,-9)/4\) |
| \(\beta^2\) | \((0,0,0)\) | \((0,0,0)\) |
| \(\beta^3\) | \((145,145,17)/96\) | \((41,41,9)/96\) |

The executable recomputation restricts the actual native contact to the binary block and checks its exact integer Pauli expansion after centering and multiplying by four. It multiplies that expansion on the graph, takes the unital partial trace, and expands the logarithm. A disagreement with the displayed coefficients on either graph, a surviving disconnected contribution, or a cubic change on a bond reaching distance three fails the tests. Related high-temperature coefficients for chains with nearest-neighbour interactions are studied by Burke, Nakerst and Haque. Their Hamiltonian of mean force divides the dimensionless logarithm by \(\beta\). [Exact cell calculation](coefficients.py) · [Tests](test_boundary.py) · [Phys. Rev. E 110, 014111](https://doi.org/10.1103/PhysRevE.110.014111)

**Fixed-temperature check on the two selected contacts.** At \(\beta=1/16\), retain the full exponential on the three sites \(0,1,8\). Subtract \((\kappa+1)I/2\) per contact and call the centered path Hamiltonian \(H_p\). Its norm is at most \(\kappa+1\). Global bit flip and occupation conservation make the symmetric and odd one-load words eigenvectors of \(F=E(e^{-\beta H_p})\), with eigenvalues \(f_+,f_-\). The coefficient of \(XX\) in \(K\) is \(-\tfrac14\log(f_+/f_-)\). The internal Hamiltonian gives \(-\beta/4\), so a strictly positive correction is equivalent to \(e^\beta f_--f_+>0\).

The executable evaluates the degree-12 Taylor polynomial of the centered matrix exponential in rational arithmetic. Each of the two expectation errors is at most \(u^{13}/[13!(1-u/14)]\), \(u=\beta(\kappa+1)\), by the scalar exponential series and norm contractivity. It multiplies the lower endpoint for \(f_-\) by the degree-four lower bound for \(e^\beta\) and subtracts the upper endpoint for \(f_+\). The resulting rational lower bounds exceed \(0.0012\) at \(\kappa=8\) and \(0.0003\) at \(\kappa=4\). Thus the correction on the gated two-contact system is positive at the declared temperature. The complete-cell statement above is the cubic Taylor coefficient.

## A convergent collar approximation at any fixed temperature on a finite graph

The formal support theorem alone gives no numerical remainder. The following reparametrization gives one for every real finite \(\beta\) on one finite graph. Scalar shifts of \(H\) and \(H_R=\sum_{X\subset R}h_X\) change the formulas only by scalar operators. Choose these representatives with \(\|H\|\le g\), \(g>0\), and keep \(\|H_R\|\) in the bound. If \(H\) is scalar, its regional logarithm is scalar directly. Set
\[
 z(w)=\frac1g\log\frac{1+w}{1-w},\qquad
 t=\tanh(\beta g/2),\qquad
 D(w)=-\log E(e^{-z(w)H})-z(w)H_R.
\]
For \(0\le |t|<r<1\) let \(P_m(t)\) be the degree-\(m\) Taylor polynomial of \(D(w)\) evaluated at \(t\). Define
\[
 L_r=\log\frac{1+r}{1-r},\quad a_r=2\arctan r,\quad
 b_r=\frac{(1-r)^2}{1+r^2},\quad
 M_r=(e^{L_r}+1)\frac{\log(1/b_r)}{1-b_r}
       +\frac{L_r+a_r}{g}\|H_R\|.
\]
Then, in operator norm,
\[
 \|D(t)-P_m(t)\|\le
 M_r\,\frac{(|t|/r)^{m+1}}{1-|t|/r}.
\]
Modulo scalars, \(P_m(t)\) is supported on the same finite collar proved above, of width \(m-1\) for centered native contacts and \(m\ge1\).

**Proof.** On \(|w|\le r\), the fractional linear map \((1+w)/(1-w)\) lies in the disk centered at \((1+r^2)/(1-r^2)\) with radius \(2r/(1-r^2)\). Its modulus lies between \(e^{-L_r}\) and \(e^{L_r}\) and its argument has magnitude at most \(a_r<\pi/2\). Consequently \(|g\operatorname{Re} z|\le L_r\) and \(|g\operatorname{Im} z|\le a_r\). The spectral theorem gives
\[
 \operatorname{Re} e^{-zH}
 = e^{-(\operatorname{Re}z)H}\cos((\operatorname{Im}z)H)
 \succeq e^{-L_r}\cos(a_r)I=b_r I .
\]
The unital partial trace preserves this inequality and contracts the operator norm. Thus \(F=E(e^{-zH})\) has \(\operatorname{Re}F\succeq b_r I\) and \(\|F\|\le e^{L_r}\). Its principal logarithm is holomorphic throughout the disk. The resolvent integral for \(\log F\) and the bound \(\|(F+sI)^{-1}\|\le1/(b_r+s)\) give
\[
 \|\log F\|\le(\|F\|+1)\int_0^\infty
        \frac{ds}{(1+s)(b_r+s)}
 \le(e^{L_r}+1)\frac{\log(1/b_r)}{1-b_r}.
\]
Adding the bound for \(zH_R\) gives \(M_r\). Cauchy's coefficient bound and the geometric tail prove the inequality. Since \(z(0)=0\), its degree-\(k\) composition coefficient only involves orders \(1,\dots,k\) of the original series. This preserves each finite collar. At \(w=t\), \(z=\beta\) exactly. The proof does not interchange a volume limit with differentiation or assume absolute convergence after summing individual cluster norms.

For example the two selected native contacts have spectrum in \([0,9]\) individually. Shifting each by \(4.5I\) gives \(g=9,\ \|H_R\|\le4.5\) for their sum and one internal contact. At \(\beta=1/16\), \(r=1/2\) gives \(M_r<10\), \(t/r<9/16\) and a degree-32 tail below \(1.30\times10^{-7}\). On the full 30-contact cell the sum of the half widths gives \(g=135\) and \(t=0.9995669\). In any regular binary grand ensemble in this normalization, the unshifted all-aligned word has energy zero and the mean energy in the normalized trace is \((1+2\kappa)|E_{\rm graph}|/4\). The spectral half width is therefore at least \((1+2\kappa)|E_{\rm graph}|/8\) for \(\kappa\ge0\). Even optimal global centering gives \(t\to1\) at fixed positive temperature as the graph grows. The bound does not give a uniform collar width on that refinement.

**Application to preparation derivatives.** Write \(C_m\) for the collar and let \(\delta\) be a trace-zero Hermitian regional derivative. If its marginal on \(C_m\) vanishes, then
\[
 |\operatorname{Tr}[\delta(K(\beta)-\beta H_R)]|
 \le \|\delta\|_1 M_r\,\frac{(|t|/r)^{m+1}}{1-|t|/r}.
\]
This follows by taking the trace against the collar polynomial and using trace duality on the remainder. For derivatives visible in the collar, their signed contraction with the polynomial remains in the formula. An observer small-ball proof needs these contractions and the tail to have the required order after stress normalization. A finite graph bound depending on its total spectral width gives no such scaling by itself. [Small-ball moment criterion](../source_implication/README.md)

An exactly commuting cut gives a useful endpoint. If ([H_R,H-H_R]=0), the exponential factors, (E(e^{-eta H})=e^{-eta H_R}E(e^{-eta(H-H_R)})), and the two factors commute. The second acts on the inner boundary alone for edge interactions, so the correction to the internal Hamiltonian is supported on that boundary at every real (eta). The native third-order coefficient on ports 1 and 8 excludes this exact boundary-only conclusion for the binary cell as an analytic identity near zero.

There are fixed-temperature results with different hypotheses. For finite-range one-dimensional spin chains, Kuwahara proves an operator-norm approximation to the actual regional logarithm in Eq. (17). In higher dimensions Eq. (19) controls an approximate density by exponentiating a local operator, a weaker conclusion for small eigenvalues. Kato and Kuwahara discuss connected coefficients and the separate convergence problem in Sec. VII. The native triangular cell and its federations require the hypotheses of any transferred estimate to hold in their own graph metric. [Kuwahara, Phys. Rev. X 15, 041010](https://arxiv.org/html/2407.05835v3#S4) · [Kato and Kuwahara, Sec. VII](https://arxiv.org/html/2504.02235v2#S7)


## A lower bound on conditional mutual information from a transport signal

Let the *whole* finite Hilbert space be \(A\otimes B\otimes C\), with faithful \(\rho=e^{-\beta G}/Z\) and real \(\beta\ne0\). Suppose two preparations have identical \(BC\) marginals. At native real time \(t\ne0\), let the difference of their readings on an effect initially supported in \(C\), \(0\le P_C\le I\), have magnitude \(p>0\). Define \(\lambda=\lambda_{\min}(\rho)\), \(\tau=\rho_{AB}\otimes I_C/d_C\), and choose
\[
 m=\lambda_{\min}(\tau)/\lambda_{\max}(\rho),\qquad
 M=\lambda_{\max}(\tau)/\lambda,\qquad
 Q=\frac{\pi p\sqrt{\lambda}}{4\sinh(\pi|t/\beta|)} .
\]
With natural logarithms,
\[
 I(A:C\mid B)_\rho\ \ge\
 \frac{Q^2}{4\log(1+4M/(mQ))}>0 .
\]
This applies to each finite buffer with a resolved propagation signal. The denominator, the Gibbs spectral floor, the modular time \(t/\beta\), and the signal itself determine the size of the certificate.

**Proof.** Let \(\mathcal N=\operatorname{Tr}_A\), \(\rho_0=\mathcal N\rho,\ \tau_0=\mathcal N\tau\). Write \(\widehat x=I_A\otimes x\) for lifting a \(BC\) operator. On the respective Hilbert spaces of matrices with Hilbert-Schmidt inner product, set
\[
 T(X)=\tau X\rho^{-1},\quad T_0(x)=\tau_0x\rho_0^{-1},\quad
 Vx=\widehat{x\rho_0^{-1/2}}\,\rho^{1/2},\qquad v=\rho_0^{1/2}.
\]
Both modular operators are positive self-adjoint. The trace identities \(V^\dagger V=I,\ V^\dagger TV=T_0\) follow by cycling under the full trace and then tracing out \(A\). For example,
\(\langle Vx,TVy\rangle=\operatorname{Tr}[\rho_0^{-1/2}x^\dagger\tau_0y\rho_0^{-1/2}]=\langle x,T_0y\rangle\).
Also \(Vv=\rho^{1/2}\), \(\|v\|_2=1\), and \(mI\le T\le MI\) because its spectrum consists of eigenvalue ratios of \(\tau,\rho\). Compression gives the same bounds for \(T_0\).

The relative-entropy difference \(D(\rho\|\tau)-D(\rho_0\|\tau_0)\) equals the conditional mutual information \(I\). Left and right multiplication commute, so \(D(\rho\|\tau)=\langle Vv,-\log(T)Vv\rangle\), and likewise for \(D(\rho_0\|\tau_0)\). Define
\[
 w_x=(T+xI)^{-1}Vv - V(T_0+xI)^{-1}v,\qquad x\ge0.
\]
Expanding the square and using \(V^\dagger(T+xI)V=T_0+xI\) gives
\[
 \|(T+xI)^{1/2}w_x\|_2^2=
 \langle v,[V^\dagger(T+xI)^{-1}V-(T_0+xI)^{-1}]v\rangle .
\]
Integration over \(x\ge0\) gives \(I\) by the resolvent integral for the logarithm. In particular \(\int_0^R\|w_x\|_2\,dx\le\sqrt{I\log(1+R/m)}\). For the tail, subtract \(x^{-1}I\) from both resolvents. Each difference from \(x^{-1}I\) has norm at most \(M/x^2\), so \(\int_R^\infty\|w_x\|_2 dx\le2M/R\).

For \(s=-t/\beta\) put \(S_s=\sinh(\pi|s|)/\pi\). The scalar identity
\[
 y^z-1=\frac{\sin(\pi z)}{\pi}\int_0^\infty
 x^z\left[(1+x)^{-1}-(y+x)^{-1}\right]dx,\quad y>0,
\]
follows from the beta integral for \(0<\operatorname{Re}z<1\) and extends to \(-1<\operatorname{Re}z<1\) by analyticity. Functional calculus with \(z=is\) and the two bounds above give
\[
 \|T^{is}Vv-VT_0^{is}v\|_2
 \le S_s\big(\sqrt{I\log(1+R/m)}+2M/R\big).
\]
The difference on the left is \((U-W)\rho^{1/2}\), with unitaries \(U=\tau^{is}\rho^{-is}\) and \(W=\widehat{\tau_0^{is}\rho_0^{-is}}\). Therefore \(\sqrt\lambda\,\|U-W\|_\infty\) is bounded by the left side. Because \(\tau\) commutes with \(P_C\),
\(U^\dagger P_C U=\rho^{is}P_C\rho^{-is}=e^{itG}P_Ce^{-itG}\).
The observable \(W^\dagger P_CW\) acts only on \(BC\). Its reading is identical in the two preparations. Their trace distance in the unhalved norm is at most 2, and the two conjugated effects differ in operator norm by at most \(2\|U-W\|_\infty\). Thus \(p\le4\|U-W\|_\infty\). Combining these inequalities and choosing \(R=4M/Q\) proves the bound.

**Native check.** For the two selected contacts \(0\,1,1\,8\), take \(A=\{0\}, B=\{1\}, C=\{8\}\). The initial preparations are the vacuum and one load at 0. At \(\kappa=8,\beta=t=1/16\), the occupation difference at port 8 is \(2.34630\times10^{-7}\). The lower bound is \(1.94754\times10^{-19}\) nats. The entropy difference gives \(8.08112\times10^{-6}\) nats. Direct numerical integration of the positive resolvent square agrees with the entropy calculation, also for \(\kappa=4\) and for inverted temperature. These are floating-point cross-checks of the finite-dimensional inequality. They do not estimate its sharpness on a federation.

On any connected native binary graph, the one-load sector has hopping matrix elements \(-1/2\) per edge. For occupied \(a\in A\), readout \(c\in C\), and distance \(d=\operatorname{dist}(a,c)\), the first nonzero amplitude coefficient of \(e^{-itG}\) is \((i/2)^d w_{ac}/d!\), where \(w_{ac}>0\) counts the shortest paths. All lower powers vanish and diagonal potentials cannot enter a length-\(d\) shortest walk. A sufficiently short nonzero time has \(p>0\), so the bound gives strict conditional mutual information at every nonzero finite inverse temperature on the full binary grand ensemble. If sites are traced out before forming \(ABC\), the generator of the marginal's modular evolution needs separate identification. Exact equality and its block structure are characterized by Hayden, Jozsa, Petz and Winter. [Transport and resolvent checks](transport_cmi.py) · [Hayden, Jozsa, Petz and Winter](https://arxiv.org/abs/quant-ph/0304007)
