# How far can one load mediate attraction between two written cores?

An observer comparing stationary energies on a 320-cell Interval Echo board sees a lowering when two nearby written cores share one mobile load. At survey distances seven and ten links, two different core shapes gave **zero dark counts in each of four banks of 1,048,576 launches**. The inverse-distance coefficient measured separately for each shape predicts a dark fraction above 99.88% in all four. A bound uniform in core shape, board size and orientation explains why an inverse-distance tail cannot return at larger separations in this sector. [Protocol](PROTOCOL.md) · [Code](measure.py) · [All simulated detector counts and energy certificates](results.json)

## What the clock measures

The two cores are sets of stationary empty native writes on disjoint port supports. One arm of a Ramsey interferometer contains the mobile load in the ground orbital of the two-core board. In the other arm it occupies the more strongly binding of the two separately characterized single-core boards. The companion board is empty and contributes zero energy. Free evolution, inverse spatial preparation and recombination convert the energy gain into a bright or dark endpoint. The source supplies the two stationary preparations, addresses, coherent clock and recombiner. The certificate selects the single-core reference. The measured quantity is the phase difference. The survey uses the shortest number of seam steps between the nearest touched ports, denoted by $d$. The parent experiment reports $r=d+1$ to account for the two half-links from seam midpoints. We retain the same convention for both shapes. Matching $r$ to a physical ruler is a separate hypothesis. [Native seam identities](../core_band_bound/README.md) · [Source and interferometer](../MEDIATOR_PROTOCOL.md)

Two near arrangements at $r=2$ gave 29,824 dark counts for the single-seam cores and 238 for bent two-seam cores, each in 1,048,576 starts at unit exposure. Separate clock quadratures at $r=3$ measured gains of **0.072197 ± 0.000664** and **0.030626 ± 0.000664** exchange units. The errors cover the eight calibration quadratures jointly at 99%. Adding either core can lower the ground energy by at most its operator norm, so $G\le\min_i\|W_i\|\le18$ for these amount-one preparations, including the bent cores. The short exposure of 1/8 confines the phase to $[0,2.25]$ radians, below $\pi$ even with the simultaneous clock allowance. This fixes the winding for the longer exposure. At each held distance the observer chose the exposure giving phase $\pi$ under the calibrated $1/r$ forecast. The numerical ground-energy enclosures put the actual gain between zero and $2\times10^{-9}$ exchange units on all four held boards. From the four zero counts alone the simultaneous one-sided 99% upper limit is 5.72 dark events per million per held setting. Removing the contrast cost gives zero lowering and a bright return at every measured placement. A $\pi$ shift at that receiver sends every launch to the dark port. [Count records and calibration](results.json)

## A distance bound for every pair of cores

Subtract the independently defined uniform blank-host energy. The free load of amount $a\in\{1,2\}$ then has generator $A=L/2$, with $L$ the intact graph Laplacian. For a written edge $e=(u,v)$, its removal from the native blank interaction subtracts the positive operator
$$
 W_e=k(P_u+P_v)+\tfrac12(|u\rangle-|v\rangle)(\langle u|-\langle v|),\qquad k=8a^2 .
$$
For core $i$, write $W_i=\sum_{e\in i}W_e$, let $P_i$ project onto its touched ports, and let $s_i$ be the largest number of its written edges touching any one port. The squared-difference inequality gives $0\le W_i\le (k+1)s_i P_i$. The graph has maximum degree six, so $0\le A\le6I$. At a touched port of core $i$, the single-port Rayleigh quotient of $A-W_i$ is at most $5/2-k$. Each single-core ground binding $t_i=-\lambda_{\min}(A-W_i)$ is therefore at least $g=k-5/2>0$. The argument allows arbitrary numbers of written edges within either core. [Exact local reduction](../core_band_bound/measure.py)

Put $t_0=\max(t_1,t_2)$, $t=-\lambda_{\min}(A-W_1-W_2)$ and $G=t-t_0$. Adding a positive attraction gives $G\ge0$. Let $R(t)=(A+tI)^{-1}$. The expansion around $3I$ in powers of $3I-A$ gives
$$
 \|P_1 R(t) P_2\|\le t^{-1}\left(\frac{3}{t+3}\right)^d .
$$
Indeed, $3I-A$ has norm at most three and propagates by at most one edge per multiplication. Every power of degree below $d$ vanishes between the two supports, and the geometric series bounds the rest.

At the joint ground energy the block matrix with entries $W_i^{1/2}R(t)W_j^{1/2}$ has largest eigenvalue one. To see this, factor the positive operator $A+tI$ out of the ground equation, then use the equal nonzero singular-value spectra of a map times its adjoint in either order. At $t_0$ each diagonal block is bounded above by the identity by the single-core ground inequality. Spectral comparison on $0\le A\le6I$ gives $R(t)\le(6+t_0)R(t_0)/(6+t)$, so the diagonal blocks at $t$ are each bounded by $(6+t_0)/(6+t)$. The norm of the off-diagonal part is the norm of its upper-right block. Combining the two estimates yields
$$
 0\le G\le(k+1)\sqrt{s_1s_2}\,(1+6/g)
                 \left(\frac{3}{g+3}\right)^d .
$$
For an amount-one mediator the decay factor is exactly $6/17$ per link. Even if the cores change with refinement, $s_1,s_2\le6$. Every fixed positive $C/(d+1)$ therefore exceeds the envelope for all sufficiently large $d$. For two single seams the entire 99% calibration interval from the 320-cell clock exceeds the envelope at every $d\ge7$. For two bent cores it does so at every $d\ge8$. These are uniform placement bounds. The measured gains at the four particular held placements are much smaller.

## Certifying the ground energies

The sparse eigensolver proposes intervals of width two billionths of an exchange unit. An independent rational sign test certifies that the bottom of the full material spectrum falls inside each interval. The matrix $J=6I-L$ has nonnegative integer entries with row sum six. We count all walks through length 24 between the touched ports by integer multiplication; $6^{24}<2^{63}$. At a positive rational endpoint $t$, the corresponding 24th-order resolvent polynomial has uniform operator error at most $(3/(3+t))^{25}/t$. On the support of $W$, the symmetric matrix $W-W R(t)W$ is positive definite precisely above the deepest binding, since $W$ is positive definite there. At the upper endpoint all leading principal minors are positive even when the truncation error is added to the resolvent. At the lower endpoint a principal minor is negative even when that error is subtracted. The recorded signs are exact rational decisions. This avoids mistaking a nearly degenerate excited eigenvector for the bottom of the spectrum. [Sign verifier](measure.py) · [Independent full-cell generator tests](test_range.py)

Each Ramsey count uses the same binomial quantile at the two ends of its probability allowance. The bracket includes the ground-energy interval and the phase error from the two compiled stationary sources, each with energy residual below $10^{-11}$. Every enrolled quantile selects the same count at both ends. The simulation keeps all outcomes. The clock and stationary compiler define the observable domain. The result bounds an indirect ground-energy interaction for one load on disjoint stationary writes. Moving sources, populated memories, gapless collective media, and standing-wave spectra in those sectors require their own measurements. The complete atomic criterion is Partial locally. OpenWave records the M10 orbital cell as untested. [Official criterion and M10 decision](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L468)

Run from the repository root.

```sh
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
.venv/bin/python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_31_orbitals.mediator_range.measure
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_31_orbitals/mediator_range
```
