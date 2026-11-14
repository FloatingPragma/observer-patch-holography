# Native pair record second moment

Let a finite connected simple graph have \(N\geq2\) vertices, \(M\) edges,
oriented incidence \(B\), Laplacian \(L=BB^T\), and cycle projector
\(P=I-B^TL^+B\). Each native event chooses an edge uniformly and a fair coin.
The two endpoint loads are replaced by the integer floor/ceiling split of
their sum. An initial excess \(+1\) at \(x\) and deficit \(-1\) at \(y\)
move or annihilate under this rule. No new pair is created.

The initial ordered pair is uniform over the \(2M\) oriented edges.
For a chosen edge with \(n_0=\pm B_e\), the initial integer current record is
\(F_0=\mp e\). Thus \(n_0+BF_0=0\). Each native transfer adds its signed
current to that edge's record, preserving \(n+BF=0\). Let \(\tau\) be the
annihilation time. Then

\[
\boxed{\mathbb E[F_\tau F_\tau^T]=\frac NM P},\qquad
\mathbb E\|F_\tau\|^2=\frac{N(M-N+1)}M.
\]

Neither regularity nor edge transitivity is required. The mean is zero for
the stated symmetric orientation preparation, so the second moment also is
the covariance. These are integer current records, with no modular wrap,
coherent extension, physical field energy, or radiation interpretation.

## Edge traversal counts

Before annihilation, the ordered load locations form the killed product walk
with generator \(A=L\otimes I+I\otimes L\), measured in \(2M\) native
attempts. Functions on the diagonal are set to zero. If \(g_e(x,y)\) is the
expected number of subsequent nonzero transfers on edge \(e=\{a,b\}\),
the finite Poisson equation is

\[
(Ag_e)(x,y)=1_{x\in e}+1_{y\in e},\qquad x\ne y.
\]

The right side includes the adjacent annihilation event: two endpoints on
the edge give transfer probability \(1/M\). Summing over all off-diagonal
ordered pairs gives \(4(N-1)\). Since every column sum of \(A\) is zero,

\[
\sum_x(Ag_e)(x,x)
=-2\sum_{x\sim y}^{\mathrm{ordered}}g_e(x,y)=-4(N-1).
\]

Consequently the stated initial average gives

\[
\nu_e=\mathbb E[\#\text{ nonzero transfers on }e]
=\frac{N-1}{M}
\]

for every edge, including edges in an irregular graph. Summing over edges
gives \(N-1\) expected transfers. Applying the same argument to
\(Ah=2M\), where \(h\) is expected annihilation time, gives
\(\mathbb E\tau=N(N-1)/2\). Connectedness and finiteness ensure absorption
with finite expected time; each nonabsorbed pair has a positive-probability
finite path to collision.

## Stopped cycle martingale

For the signed load vector \(n_t\), the fair endpoint rule gives

\[
\mathbb E[\Delta F_t\mid n_t,F_t]=\frac1{2M}B^Tn_t.
\]

Thus \(X_t=PF_{t\wedge\tau}\) is a martingale. Each nonzero transfer on
edge \(e\) has increment \(\pm Pe\). Summing the conditional second-moment
identity, with finite expected stopping time and bounded increments, gives

\[
\begin{aligned}
\mathbb E[X_\tau X_\tau^T]
&=\mathbb E[X_0X_0^T]+\sum_e\nu_e(Pe)(Pe)^T\\
&=\frac1M P+\frac{N-1}{M}P=\frac NM P.
\end{aligned}
\]

At absorption, conservation gives \(BF_\tau=0\), hence
\(X_\tau=F_\tau\). Taking the trace uses
\(\operatorname{rank}P=M-N+1\) and proves the norm formula.

## Scope, verification, and use

Fable's LA_L02 audit supplies the adjacent-lifetime argument and the cycle
martingale trace identity, with its closed form proved there for regular
edge-transitive graphs. The per-edge Poisson summation above removes that
restriction and supplies the full second-moment matrix. This attribution
distinguishes the extension from the prior result.

[pair_record_covariance.py](pair_record_covariance.py) checks every edge count,
the lifetime, and the cycle projectors in exact rational arithmetic on six
small graphs, including irregular examples. It also compares the complete
formula to the separately saved 30-by-30 Fable carrier moment matrix. Their
maximum entry difference is `1.571e-14`; both give trace \(38/5\).
[The receipt](pair_record_covariance_results.json) records the exact inputs.
The comparison reads archived geometry constants and matrix data; it does
not execute Fable's producer or audit code.

For a subsequent declared linear record-field map \(R\), the expected
quadratic output \(F_\tau^T R^T Q R F_\tau\) is
\((N/M)\operatorname{tr}(Q R P R^T)\). This replaces a product-state moment
solve by a graph projector when all the stated preparation and current-record
hypotheses hold. It does not derive the subsequent map, its energy scale,
field propagation, a Coulomb force, or radiation from the native process.

The uniform initial edge distribution matters. The code retains an irregular
triangle-with-tail case in which preparing only the first edge changes the
traversal vector. A different initial record changes the Gauss background;
a wrapped modular representative does not retain this integer quadratic
moment. Neither transfer is implied by the formula.
