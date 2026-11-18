# The binding theorem of the polynomial-sv model

This document proves that a pair of loads on a seam stays a pair under the contrast-memory law
(`A-contrast-memory-spherical-v1`, a candidate native repair law) on every federation of complete cells, with a
bound independent of the number of cells; gives the exact constants of one cell; and bounds the emitter and the
reader. The law, carrier and notation are those of [SPEC.md](SPEC.md). The constants are measured by
[theorem.py](theorem.py) and stored in [RESULTS.json](RESULTS.json).

## 1. Setting

$G = (V, E)$ is a federation: complete twelve-port cells, every internal seam kept, external seams forming a
partial matching on ports, so every port has degree $d_v \in \{5, 6\}$. On loads $n_v \in \{0, 1, 2\}$,

$$H = K + \kappa C,\qquad K = \tfrac12 \sum_{e} (I - T_e),\qquad C = 6 \sum_v n_v^2 - 2 \sum_{\{a,b\} \in E} n_a n_b .$$

Every summand of $K$ is a projector, so $K \ge 0$. $H$ conserves the multiset of loads. The sector of two loads
contains every pair of distinct ports and every doubled load; its dimension is $V(V+1)/2$. Its values of $C$ are

| Configuration | $C$ |
| --- | ---: |
| two loads on the two ends of a seam (a pair) | 10 |
| two loads on ports that share no seam | 12 |
| one doubled load | 24 |

since a binary configuration has $C = 6 \cdot 2 - 2\,[\text{adjacent}]$ and a doubled load has $C = 6 \cdot 4$. The
swaps conserve the multiset of loads, so the doubled loads form a reducing subspace; they stay part of the
complement below and of every outcome.

Let $P$ project onto the span of the pairs on seams and $Q = I - P$ in the sector of two loads. Write

$$H = \begin{pmatrix} A & B^* \\ B & D \end{pmatrix},\qquad A = PHP,\ B = QHP,\ D = QHQ .$$

## 2. The three inequalities

**Theorem 1.** On every federation of complete cells with a matching of external seams, and every real $\kappa$,

$$PHP \le (10\kappa + 6)\,P,\qquad QHQ \ge 12\kappa\,Q,\qquad \|QHP\| \le \sqrt{10}.$$

*Proof.* $C$ is diagonal, so $A = 10\kappa P + PKP$, $D = QKQ + \kappa\,QCQ$ and $B = QKP$.

*First inequality.* The diagonal of $PKP$ at the pair on $(a, b)$ is $(d_a + d_b - 2)/2 \le 5$. An off-diagonal
element of $PKP$ joins the pairs on $(a, b)$ and $(a, t)$ when one swap moves the load from $b$ to $t$; then $t$ is
adjacent to $a$ and to $b$, so $a, b, t$ span a face, and the element is $-1/2$. Two seams are so joined exactly when
they are two sides of one face. All faces lie inside cells, so an external seam has no such neighbour, and inside a
cell the seams with this relation form the line graph of the dodecahedron (faces of the cell are vertices of the
dual, seams are its edges). The adjacency matrix of a line graph is $M^{\top} M - 2I \ge -2I$ with $M$ the incidence
matrix, so the off-diagonal part of $PKP$, which is $-\tfrac12$ times that adjacency, is at most $1$. Hence
$PKP \le 6P$ and $A \le (10\kappa + 6)P$.

*Second inequality.* $QKQ \ge 0$ as a compression of $K \ge 0$, and $QCQ \ge 12\,Q$ by the table. Hence $D \ge 12\kappa Q$.

*Third inequality.* $B$ has entries $-1/2$. Its column at the pair on $(a, b)$ has one entry for every swap that moves
a load of the seam to a neighbour that leaves the two loads apart: at most $(d_a - 1 - m) + (d_b - 1 - m)$ entries
with $m$ the number of common neighbours of $a$ and $b$. For an internal seam $m = 2$ and the count is at most 6; for an
external seam $m = 0$ and the count is at most 10. The largest absolute column sum is therefore at most 5. The line
of $B$ at two separate loads on $x, y$ has one entry for every common neighbour $t$ of $x$ and $y$ and every choice of
which load came from $t$: at most $2 \cdot 2 = 4$ entries, because a nonadjacent pair of ports has at most two common
neighbours (inside a cell by the geometry of the icosahedron; across cells because the external seams form a
matching, so at most one common neighbour lies in either cell and none in a third). Doubled loads have no entry.
The largest absolute line sum is at most 2. By the Schur test $\|B\| \le \sqrt{5 \cdot 2} = \sqrt{10}$. $\square$

## 3. Persistence at every value of the evolution parameter

**Theorem 2.** Let $A \le a P$, $D \ge c\,Q$, $\|B\| \le \beta$ with $\gamma = c - a > 0$, $m = \operatorname{rank} P$, and let $\Pi$
be the spectral projector of $H$ onto its lowest $m$ levels. Then

$$\|Q\,\Pi\| \le \beta/\gamma,\qquad \|(I - \Pi)\,P\| \le \beta/\gamma,$$

and for every normalized state $\psi$ with $P\psi = \psi$ and every real $t$,

$$\|Q\, e^{-itH} \psi\| \le \frac{2\beta}{\gamma},\qquad
\Pr(\text{no pair on a seam};\, t) \le \min\Big(1, \frac{4\beta^2}{\gamma^2}\Big).$$

*Proof.* By min-max with the trial spaces $\operatorname{ran} P$ and $\operatorname{ran} Q$, $\lambda_m(H) \le a < c \le \lambda_{m+1}(H)$, so the
gap of the ordered blocks is not closed by $B$. Put $H_{\rm low} = \Pi H \Pi$ (levels at most $a$) and $X = Q\Pi$. From
$QH\Pi = X H_{\rm low}$ and $QH = QHP + QHQ$,

$$D X - X H_{\rm low} = -B\,P\,\Pi .$$

The spectra of $D$ on $\operatorname{ran} Q$ and of $H_{\rm low}$ on $\operatorname{ran}\Pi$ are separated by $\gamma$, with $D$ above. The
solution of this Sylvester equation is

$$X = -\int_0^\infty e^{-s(D - a)}\, B P \Pi\, e^{s(H_{\rm low} - a)}\, ds ,$$

and $\|e^{-s(D - a)}\| \le e^{-s(c - a)}$, $\|e^{s(H_{\rm low} - a)}\| \le 1$, so $\|X\| \le \beta/\gamma$. In the same way
$Y = (I - \Pi)P$ solves $H_{\rm high} Y - Y A = (I - \Pi)\,Q\,B$ with $H_{\rm high} \ge c$ and $A \le a$, so $\|Y\| \le \beta/\gamma$.

For $\psi = P\psi$ write $\psi = \Pi\psi + (I - \Pi)P\psi$. The first part stays in $\operatorname{ran}\Pi$ under $e^{-itH}$, so
$\|Q e^{-itH}\Pi\psi\| \le \|Q\Pi\| \le \beta/\gamma$; the second part has norm at most $\beta/\gamma$. The sum is at most
$2\beta/\gamma$, and the probability of the outcome "no pair", which is $\|Q e^{-itH}\psi\|^2$, is at most $4\beta^2/\gamma^2$.
Mixed states supported in $\operatorname{ran} P$ follow by convexity. $\square$

**Corollary.** On every federation of complete cells with a matching of external seams, Theorems 1 and 2 give
$a = 10\kappa + 6$, $c = 12\kappa$, $\beta = \sqrt{10}$, $\gamma = 2\kappa - 6$, and for $\kappa > 3$

$$\Pr(\text{no pair on a seam};\, t) \le \frac{40}{(2\kappa - 6)^2}$$

at every value of the evolution parameter, hence at every reading of every observer. The bound does not depend on
the number of cells or on the routing of the external seams. At $\kappa = 128$ it is $40/250^2 = 0.00064$.

The dressed band $\operatorname{ran}\Pi$ is invariant under the law. The bare outcome "no pair" also contains its virtual
dressing, so it need not vanish or decrease along the evolution: the pair may separate briefly and is read on
other seams.

**Measured constants.** The proof uses the bounds $a - 10\kappa \le 6$ and $\beta \le \sqrt{10}$. On the federations of the
model the largest level of the adjacent block of the exchange, $\max \operatorname{spec}(PKP)$, and the largest singular value of
the coupling, $\|QKP\|$, are

| Federation | $\max \operatorname{spec}(PKP)$ (bound 6) | $\|QKP\|$ (bound 3.1623) |
| --- | ---: | ---: |
| 20 cells | 5.7008 | 2.0821 |
| 80 cells | 5.7000 | 2.0811 |
| 320 cells | 5.7000 | 2.0822 |

With these constants the bound is $4\beta^2/(2\kappa - 5.70)^2 = 17.3/(2\kappa - 5.70)^2$: 0.163 at $\kappa = 8$ and 0.00028 at
$\kappa = 128$ on 20 cells.

## 4. One cell

On one cell without degree term every port has five seams, the pair costs 8, separate loads 10 and a doubled load 20
($C = \sum_e (n_a - n_b)^2$). There are 30 pairs on seams, 30 pairs at distance two, 6 antipodal pairs and 12 doubled
loads. The exact constants are

$$a = 8\kappa + 5,\qquad c = \min\big(10\kappa + 4 - \sqrt6,\ 20\kappa\big),\qquad \beta = 2.$$

- $\max\operatorname{spec}(PKP) = 5$: the diagonal is 4 and the line graph of the dodecahedron has smallest adjacency level $-2$.
- $\|QKP\| = 2$: the largest absolute line and column sums of $QKP$ are both 2, and the uniform vector on the 30 pairs
  attains the value.
- The bottom of the block of separate loads is $4 - \sqrt6$. The two orbits of nonadjacent pairs (distance two and
  antipodal) give the quotient $\begin{pmatrix} 3 & -\sqrt5 \\ -\sqrt5 & 5 \end{pmatrix}$, whose lower level $4 - \sqrt6$ has
  a positive lift; since the block has nonpositive off-diagonal elements, the Perron argument makes it the bottom.
- Doubled loads cost $20\kappa$ and the exchange on them is nonnegative.

The degree term of the law adds $\kappa \sum_v n_v^2$ on one cell: $2\kappa$ on every binary pair and $4\kappa$ on a doubled
load. Then $a = 10\kappa + 5$, $c = \min(12\kappa + 4 - \sqrt6,\ 24\kappa)$, and the gap $\gamma = 2\kappa - 1 - \sqrt6$ (for
$12\kappa + 4 - \sqrt6 \le 24\kappa$) is the same. The bound is $16/\gamma^2$ capped at
one: 0.1016 at $\kappa = 8$ and 0.00025 at $\kappa = 128$.

## 5. The emitter

Add one fuel state $|F_e\rangle$ at the seam $e = (a, b)$ with energy $E_f = \langle e|H|e\rangle = 10\kappa + (d_a + d_b - 2)/2$, and
during the pulse the term $g(|F_e\rangle\langle e| + |e\rangle\langle F_e|)$ for the duration $\tau = \pi/(2g)$, with $H$ acting throughout.

**Theorem 3.** The emitted state $\psi(\tau)$ satisfies

$$\big\|\psi(\tau) + i\, e^{-iE_f\tau}\, |e\rangle\big\| \le \frac{\sqrt{d_a + d_b - 2}}{2g}.$$

*Proof.* On $\operatorname{span}\{|F_e\rangle, |e\rangle\}$ the comparison dynamics with diagonal $E_f$ and coupling $g$ is an exact
resonant rotation, $\phi(s) = e^{-iE_f s}(\cos(gs)\,|F_e\rangle - i \sin(gs)\,|e\rangle)$, which reaches $-i e^{-iE_f\tau}|e\rangle$ at
$s = \tau$. The full generator differs from the comparison generator on $\phi(s)$ only through
$r = (H - E_f)\,|e\rangle$, the off-diagonal column of $H$ at $e$: $d_a + d_b - 2$ allowed swaps with amplitude $-1/2$ each, so
$\|r\|^2 = (d_a + d_b - 2)/4$. Duhamel's formula gives
$\|\psi(\tau) - \phi(\tau)\| \le \|r\| \int_0^\tau |\sin(gs)|\,ds = \|r\|/g$. $\square$

After the pulse the converter is off and the law preserves the error in norm. Residual fuel reads "no pair". With
Theorem 2,

$$\Pr(\text{no pair, doubled load or fuel};\, t) \le \min\Big(1, \Big[\frac{2\beta}{\gamma} + \frac{\sqrt{d_a + d_b - 2}}{2g}\Big]^2\Big).$$

On the federations $d_a + d_b - 2 \le 10$. At $\kappa = 128$, $g = 64$, $\beta = \sqrt{10}$, $\gamma = 250$ the bound is the exact
rational

$$10\Big(\frac{2}{250} + \frac{1}{128}\Big)^2 = \frac{64009}{25600000} = 0.0025003515625$$

on every seam of every federation, and the weight outside the dressed band, $(\beta/\gamma + \epsilon)^2$, is at most
$35721/25600000$.

## 6. The reader with the law acting

The reader of a seam (SPEC.md, section 7) runs as four stages $G_j$, each a reversible involution on the eleven
bits, with Hamiltonian $L_j = \frac{\pi}{2\tau_j}(I - G_j)$, so $e^{-i\tau_j L_j} = G_j$ on every word. Its only control
on the loads is the diagonal projector $\mathbf 1\{n_a = n_b = 1\}$, so every stage commutes with the diagonal part
$V = \kappa C$ (and the fuel energy) of the law, and stages on different seams commute. Write $H = V + K$. The comparison
evolution omits $K$ and is exactly $e^{-iV\tau_{\rm read}} \prod_j G_j$, the ideal reader at the start of the window.
Duhamel's formula gives $\|U_{\rm actual} - U_{\rm ref}\| \le \tau_{\rm read}\|K\|$. In the sector of two loads every line of $K$ has
nonnegative diagonal at most $\Delta = 6$ and off-diagonal absolute sum equal to it, so $\|K\| \le 12$, and the total
variation between the public laws is at most $12\,\tau_{\rm read}$, for any number of seams read in parallel. At
$\tau_{\rm read} = 1/16384$ this is $3/4096$, and the probability that the public records show no agreed pair after the
emitter at $\kappa = 128$ is at most $64009/25600000 + 3/4096 = 0.0032327734375$ at every reading. After one complete
read round the ideal reader acts as $\rho \mapsto \sum_e P_e \rho P_e + Q\rho Q$, each branch on a seam is a pair covered by
Theorem 2, and later "no pair" stays below $0.0025003515625 + 0.00064 + 0.000732421875 = 0.0038727734375$.

The finite readout of one cell propagates all 161,792 states (78 two-load states, the fuel and 2048 reader words)
with the law acting in every stage: the public total-variation error is $2.09 \times 10^{-7}$ against the bound
$\tau_{\rm read}\|K\| = 4.2 \times 10^{-4}$ on the cell.

## 7. Scope

The theorem concerns the sum of the seam terms with a continuous evolution parameter. As local events the same
binding holds when each port observer serves each of its seams once per round, and it needs small event angles
when seams act independently at random ([events.py](events.py)). The bound controls which seams the two loads are
read on; it states no speed of propagation, and no continuum limit.
