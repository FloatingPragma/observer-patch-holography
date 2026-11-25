# Intercell transport of the same parked charge

On the 80-cell board the full parked collar logarithm at port 5 admits an exactly conserved extension with nonzero current into a neighboring cell. Another extension with the same compression and conservation gives the opposite initial current slope. A third gives zero slope. The three use the same native flight, source and couplings. Their spatial operators differ. This proves that conserved extension and source calibration alone leave even the sign of intercell transport undetermined. [Exact integer matrices and verification](transport.py) · [Receipts](transport.json)

## Two independent matrix blocks on one flight

Fix an excess or vacancy and either left reversal sector. Let $H$ be its 3840-dimensional native flight Hamiltonian, $h=2H$ up to scalar, and $B=h-h_{55}I$, where the diagonal four-coin block at unmarked port 5 is scalar. Form the real polynomials $B^j$ for $j=0,2,3,\ldots,20$. For each power put the ten independent coordinates of its symmetric $55$ block in the first ten rows, and the ten coordinates of its $38,5$ block **plus its transpose** in the last ten rows. Port 38 belongs to the opposite owner at the marked bridge. Sparse arbitrary-precision propagation gives the following full $20\times20$ determinants modulo 1009.

| Flying load | Even left reversal | Odd left reversal |
|---|---:|---:|
| Vacancy | 756 | 886 |
| Excess | 637 | 254 |

Every determinant is nonzero over the integers. It follows that for **any** two real symmetric coin matrices $S,T$, there is a real polynomial $f$ of degree at most twenty satisfying
$$
 [f(B)]_{55}=S,\qquad
 \tfrac12([f(B)]_{38,5}+[f(B)]_{38,5}^T)=T .
$$
These polynomials depend only on the specified targets and the native Hamiltonian. The construction works in both probe-load sectors using the corresponding polynomials, equivalently as a direct sum over the conserved total load and reversal labels. A shared scalar polynomial independent of both labels is a stronger requirement.

Choose $S=\Re K_*$ for the exact frozen parked collar logarithm. Add the earlier compact charge $3\Pi\otimes i\Im K_*$ to $f(B)$ independently in each block. This addition restores **the entire original Hermitian charge** at launch and vanishes in the $38,5$ block. Both summands commute with the free Hamiltonian, which proves conservation for all flight times. Since $H$ only hops along edges, the polynomial couples position sites at most twenty native graph steps apart. This is a bound on the probe's position kernel on this finite graph; both stored heads enter the sector definition and collar calibration. [Compact projector](README.md)

## A signed current across the source boundary

Let $R$ be the receiving cell containing 38 and take $q_R=\{P_R,C\}/2$ for the conserved extension $C$. The expectation of its native boundary current equals $d\langle q_R\rangle/dt$ by the exact continuity identity. For any head state $\rho$ at sharp launch 5, the first path to a different owner is the unique two-hop path $5\to 1\to38$. Each hop has native matrix $-I/2$, so $(H^2)_{38,5}=I/4$. Expand the expectation as the real part of the inner product between $P_R U(t)|5\rangle$ and $C U(t)|5\rangle$, taking the head trace with $\rho$. The first factor starts at order $t^2$ and at that order has support only at 38. Hence
$$
 \langle q_R(0)\rangle=\langle J_{1\to38}(0)\rangle=0,
 \qquad
 \left.\frac{d^2}{dt^2}\langle q_R\rangle\right|_0
 =\left.\frac{d}{dt}\langle J_{1\to38}\rangle\right|_0
 =-\tfrac14\operatorname{Tr}(\rho T).
$$
The other receiving-cell boundary edges have zero current derivative at launch because neither endpoint is reached in one hop. This calculation transports a signed charge, for which $q_R$ need not be positive.

With $T=\Re K_*$ the *raw* unpairing tangent at the exact equilibrium changes the initial current slope by between $-0.057778$ and $-0.057764$ nats per squared native time per unit transfer (outward-rounded at the final decimal, see the interval in the receipt). To see this, write $D=\chi_{*L}\otimes\chi_{*R}-\chi_*$. Both partial traces of $\Re D$ vanish, so $\operatorname{Tr}(D\Re K_*)= -\operatorname{Tr}[(\Re D)\log\chi_*]$. The verified state radius $r=2\times10^{-7}$, the rational enclosure for the center logarithm, trace-norm contractivity of entrywise real part on Hermitian matrices, and the resolvent bound give $(0.23105604,0.23110963)$ for that response. The error allowance about the rational center is $[2/(0.02-r)+12\sqrt8]r$. At the numerical reference, the baseline current slope is $+0.02771495$, the raw unpairing changes it by $-0.05777071$, and its volume-matched version changes it by $-0.15431882$. The latter uses the numerical volume-adjustment angle. [State and fixed-volume transfer](../full_algebra/PROOF.md)

The earlier compact extension traps the charge within the source cell. The polynomials establish that a trapped charge and an intercell response both fit the same parked calibration. The interpolation also permits $T=-\Re K_*$ or $T=0$ while retaining that calibration. Consequently none of these chosen $T$ values predicts a physical flux sign. The degree-twenty solutions can be exceptionally large: at the nominal input a single matrix entry on another port already exceeds $9\times10^{12}$ in absolute value in all twelve reported sector and sign choices. The replay uses exact rational arithmetic for the fits to the binary-rational nominal matrices. These values give no uniform coarse-scale bound or new experimental meter. Independent transport selection, stable normalization on proper regions, and coupling to measured curvature remain physical requirements.

## Matching the ordinary calibration

The ambiguity survives a fixed reference reading and an independent ordinary-matter calibration of the intercell current itself. Let \(D_o,D_p\) denote the two tangents at fixed volume in the original six-unit survey. Use the real symmetric effect \(F\) for the reciprocal three-hop link from port 0 to 39, normalized so its length is \(\ell=(\operatorname{Tr}\chi F)^{-1/6}\). This link is the unique shortest route throughout the source enclosure. Freeze the length derivative \(L=-F/[6(\operatorname{Tr}\chi_*F)^{7/6}]\) at the reference. For the parked energy \(E\) one has \(E'_o>0.29756\) and \(E'_p=0\). The independent native link experiment gives \(\ell'_p\in[0.0357614,0.0362044]\) at this same fixed volume. [Interval-certified length and derivative](../common_front/certificate.json) · [Energy and preparation](../full_algebra/PROOF.md)

Define \(M=L-(\ell'_o/E'_o)E-cI\), choosing \(c\) so that \(\operatorname{Tr}\chi_* M=0\). Then \(\operatorname{Tr}D_o M=0\) and \(\operatorname{Tr}D_p M=\ell'_p>0\). Build a base cross block \(aI+bE\) to match any specified baseline current slope and ordinary response. Use the same launch block for three conserved extensions, with cross blocks \(aI+bE\), \(aI+bE+M\) and \(aI+bE-M\). All three match the complete parked observable, reference current slope and ordinary current calibration. Their fixed-volume collar responses are exactly zero, strictly negative and strictly positive. The response for the plus choice lies in \([-0.009052,-0.008940]\) in the native time and link normalization. The construction of \(M\) uses reference and ordinary data only. The held collar measurement certifies its sign. This also excludes choosing the flux by scalar ordinary calibration alone. A specified local source action can impose further constraints on its stress, including constraints on localization and scale.

Run from the repository root.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.dressed_charge.transport
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/dressed_charge
```
