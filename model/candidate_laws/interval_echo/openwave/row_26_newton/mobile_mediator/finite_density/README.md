# Two neutral sources in a finite-density vacancy medium

Two equal-energy neutral defects change the local density of a shared thermal medium on the unmodified Interval Echo board. This experiment sums all 512 vacancy configurations of a nine-site guide. At inverse temperature \(1/16\) the average source interaction is attractive at all five measured separations. Its curvature rules out every attractive shifted Newton potential in guide distance, even allowing an additive constant. On arbitrarily large boards there is also a high-temperature regime where this medium cannot transmit a Newton tail.

## Source and measurements

The guide follows nine consecutive ports on a geodesic of the native level-zero board. Its internal edges form a path. Other loads are clamped at one apart from two disjoint, mutually nonadjacent pairs of degree-five ports. Each pair is prepared in \((0,2)\) or \((2,0)\), always at the same static energy. Only its near end touches the guide, at a single site. These are independently switchable sources. Clamping enforces the stated compression. A free exterior exchange could otherwise carry a vacancy out of the guide. The Gibbs bath is a prepared state, and translating a clamped pair can change its local environment. On each guide site the vacancy bit \(n_i\) is one for load zero and zero for load one. Compression of the blank-seam native Hamiltonian, up to a constant, gives

\[
 H=\sum_i[(\kappa+\tfrac12)d_i-\kappa(6-d_i)]n_i
 -(2\kappa+1)\sum_{ij}n_i n_j
 -\tfrac12\sum_{ij}(S_i^+S_j^-+S_i^-S_j^+).
\]

Each sum over \(ij\) counts individual guide edges once. The contact with one near source end adds \(v n_i\), where \(v=-(2\kappa+1)\) or \(2\kappa-1/2\). A fair reversal thus has mean field \(-3/4\). The source energy drops out of the four-arm free-energy difference \(U=F_{AB}+F_0-F_A-F_B\). Each signed \(U\) is formed before averaging over the four orientations. The chemical and thermal settings are supplied controls, \(\mu=30,\ \beta=1/16\). The active stiffness \(\kappa=8\) gives 4.1101 mean vacancies. A changed-law control \(\kappa=4\) gives 6.5830 at the same settings.

| Guide distance | Active mean interaction in native units |
|---:|---:|
| 1 | −0.00730634 |
| 2 | −0.00453430 |
| 3 | −0.00083537 |
| 5 | −0.00004381 |
| 6 | −0.00001217 |

The interaction between particular orientations can have either sign. For the active law, the nearest target has vacancy probability 0.374526 with no source, 0.438617 with the near source end vacant and 0.334331 when it carries excess. The acquisition simulates \(2^{22}\) independent thermal preparations per setting, with nine empty-guide channels and three local target channels per bank. The simultaneous 99% Hoeffding half-width is below 0.000964 for each bank. Because the empty state has zero energy in every compressed arm, \(\beta U=\log[P_{AB}(0)P_0(0)/(P_A(0)P_B(0))]\). The distant, small averaged energies come from the exhaustive spectra. The finite count precision alone does not establish their sign.

The first three measured potentials have \(U(1)-2U(2)+U(3)>0.000926\). Every \(u-C/(d+c)\) with \(C>0,\ 1+c>0\) has second difference \(-2C/[(1+c)(2+c)(3+c)]<0\). The changed-law nearest mean is positive. Cutting the entire guide edge between its second and third ports gives \(U=0\) at distance three by exact grand-canonical factorization. Conditioning the same cut medium on exactly four vacancies gives a positive mean \(0.020817\). Fixing the global population couples the two halves even after their material interaction has been removed.

## First connected response at arbitrary filling

Fix fugacity \(z>0\) and reference filling \(p=z/(1+z)\). Take any finite blank, binary guide on the native graph. Two disjoint source fields \(v_A(a),v_B(b)\) first touch guide sites at separation \(d\ge1\). Denote the number of shortest paths between \(a,b\) by \(N_{ab}\). As \(\beta\) tends to zero at fixed \(z\),

\[
 U=-\beta^{d+1}(2\kappa+1)^d[p(1-p)]^{d+1}
 \sum_{d(a,b)=d}v_A(a)v_B(b)N_{ab}+O(\beta^{d+2}).
\]

To see this, expand the mixed logarithm of the normalized Gibbs trace about the product state with filling \(p\). Disconnected occurrences cancel because operators on disjoint sites commute and their expectations factorize. A surviving word needs one insertion from each source and at least \(d\) intervening edges. At the minimum degree its edges are precisely a shortest path. Write each interaction with the centered variable \(n_i-p\), absorbing its one-site and constant parts. A hopping edge cannot contribute at this order. At one end of the first consecutive run of hopping edges the trace contains a single off-diagonal operator among diagonal factors. For the remaining diagonal word, every nonempty proper subset of insertions has a vertex of degree one in centered variables, hence zero trace. The connected coefficient equals the full trace, with a variance \(p(1-p)\) at every vertex. All orderings cancel the exponential factorial. Multiplying by \(-1/\beta\) gives the stated result. Independent fair reversals yield a strictly attractive leading average when each nearest contact has the two native source orientations. With fixed finite \(\mu\), \(z=e^{\beta\mu}\) and the leading coefficient uses \(p=1/2\). Exact rational powers of the full hopping matrices check the first coefficient and all preceding cancellations on paths of two, three and four sites at three fugacities.

## A volume-independent range bound

In this same binary sector suppose the sources have disjoint support and are neutral pairs with clamped exterior. For the active stiffness one may use \(J=64\) as a common bound on every single interaction norm, after removing constant terms. Base one-site energies have magnitude at most 51, interior bond norms at most 17, and the source field from a pair has magnitude at most 34 on a site. Use one base term per site, one term per edge and one switched source term on each affected site. A label overlaps at most \(D=16\) labels including itself at degree at most six. The fugacity is in the product reference measure and imposes no restriction on this norm bound.

Set \(y=e^{\beta J}-1\). Each connected, nonempty set of occupied interaction labels \(S\) defines a polymer. Disjoint-support components factor in the Gibbs series: for fixed component words the sum over interleavings gives precisely the product of their separate exponential denominators. Summing words that use every label in \(S\) at least once gives an absolute activity bounded by \(y^{|S|}\), since the normalized product expectation has norm one. At most \((4D)^{m-1}\) connected sets of \(m\) labels contain a specified label, by a rooted ordered spanning tree with at most \(D\) choices per tree edge. For any polymer \(T\), the sum of \(|w(S)|e^{2|S|}\) over incompatible \(S\) is consequently bounded by \(D |T|e^2y/(1-4De^2y)\). This is at most \(|T|\) if \(5De^2y<1\).

Apply the Kotecký–Preiss cluster bound to activities weighted by \(e^{|S|}\), with size cost \(|S|\). The absolute weighted sum of clusters incompatible with a single label is at most one. In the four-log difference, only connected clusters that contain both switched fields survive. Each contains a label in the support of \(A\), and its total size is at least \(d+2\). Removing the exponential weight bounds the interaction by

\[
 |U|\ \le\ |\operatorname{supp}A|\,e^{-(d+2)}/\beta .
\]

The bound is uniform in board size, fugacity and orientation, and passes to the orientation average. The sufficient interval is \(0<\beta<\log(1+1/(5De^2))/J\), about \(2.64\times10^{-5}\) with these conservative constants. The nine-site measurements at \(1/16\) have their own finite spectral certificate. Neither experiment selects the physical temperature, a proper free-space radius or a freely moving gravitating source. [Kotecký and Preiss, cluster expansion for abstract polymer models](https://doi.org/10.1007/BF01211762)

## Spectral error and replay

The chemical term is included separately in each exact vacancy-number block of dimension at most 126. All matrix entries are exactly representable half-integers. For each eigenbasis the verifier checks the residual and orthogonality in Frobenius norm. For these symmetric blocks \(\|H\|_2,\|\Lambda\|_2<511\), \(\|V\|_F<11.9\); the computed residual is below \(5\times10^{-11}\) and the orthogonality defect below \(5\times10^{-13}\). With binary64 inner-product error \(\gamma_{128}<1.5\times10^{-14}\), the exact Frobenius errors are below \(2\times10^{-9}\) and \(3\times10^{-12}\). For the residual, use \(\||H||V|\|_F\le511\|V\|_F\) from the symmetric absolute row-sum bound. For the Gram matrix use \(\||V|^T|V|\|_F\le\|V\|_F^2\). These slack bounds include the elementwise multiplication, subtraction and norm rounding. The polar factor \(Q\) then differs from \(V\) by less than \(3\times10^{-12}\) in operator norm. Bounding \(HQ-Q\Lambda\) gives an eigenvalue error below \(10^{-8}\) by Weyl's inequality. Directed interval exponentials and logarithms then enclose each free energy and mixed difference. Saved endpoints in bounds.json are integer multiples of \(10^{-10}\).

Run from the repository root with one BLAS thread.

```sh
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.mobile_mediator.finite_density.measure
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.mobile_mediator.finite_density.certificate
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_26_newton/mobile_mediator/finite_density/test_density.py
```
