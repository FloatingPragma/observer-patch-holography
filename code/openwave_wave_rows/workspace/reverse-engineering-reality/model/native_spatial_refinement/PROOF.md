# A cubic family of native twelve-port cells

## Construction and complete records

Use the exact icosahedral ports `(0,±1,±φ)` and cyclic permutations, in the
ordering of `model.common_source.slow_geometry.raw_frame`, with
`φ=(1+sqrt(5))/2`. Squared distance four gives the thirty inherited internal
seams. In each coordinate-zero group, the two ports with positive product
of the other coordinates are the positive ports. Join each to the port in
the positive neighbouring cell obtained by flipping the cyclic last
nonzero coordinate. The three ordered pair lists are
`((5,7),(6,4))`, `((9,11),(10,8))`, `((1,3),(2,0))`.

Cells are indexed by `(Z/LZ)^3`, `L≥2`. Every port has five internal seams
and one external seam. There are `P=12L³` ports, `M=36L³` seams and
`D=2M` recoded seam/coin events. All twelve ports occur exactly once among
the external ends. The internal cells are connected and every adjacent
pair of cells has two connections. The global graph is connected.
These statements include `L=2`, where the two orientations join distinct
port pairs, so no physical edge is duplicated. This is declared cubic
routing. The geodesic federation uses a different macrograph.

The endpoint C/E instruments, their twenty-five raw slots, hold branch,
zeros, inverse recoder and retained quantum references are the literal
ones in `openwave_native_process` and `openwave_coherent_native_architecture`.
A one-load endpoint column transfers half its population across an occupied
seam and retains half. An unoccupied endpoint is an exact hold. Thus the
uniform seam law has the full population matrix

\[
T=I-\mathcal L/D,\qquad
D T_{yx}=\begin{cases}D-6&y=x,\\1&y\sim x,\\0&\text{otherwise}.\end{cases}
\]

The executable evaluates every incident raw branch for every input port,
retains all zero slots and factors nonincident events by their exact local
identity. It compares every entry of every full source column on the
96-, 324- and 768-port carriers. This population identity does not identify
the unread original quantum channel with a random swap channel.

For the coherent policy let `A|x>=D^(-1/2)sum_e|x,e>`. The controlled
involution `S` swaps the two endpoints on the recoded odd event and is the
identity on every other arc. The exact service `R=2AA*−I` reflects the
complete active event register. Set `W=SR`. Every `S` basis column is
checked on the three carriers. The reflection has diagonal `2/D−1` and
off-diagonal `2/D`; its column inner products are exactly one or zero.
Both operators are unitary on the full source/event carrier and tensor
with any reference. On power-of-two event padding, `S=I`, `R=−I`, and
`W=−I`, as in the common program. For padded event dimension `Dp`, this
retains `P(Dp−D)` additional minus-one modes. The terminal instrument
retains padding failure as well as every active source port. Uniform preparation and reflection are declared coherent
services; this package supplies no native pulse compilation for the full
cube. This distinction also applies to preparation of the cell packet.

The sequential coherent event word is a program ordering, not an inferred
global physical clock. Disjoint local endpoint operations commute. `R`
changes no source port and `S` moves at most one graph edge, so `W^n` has
exact source support within distance `n` of the initial support.

The inherited cyclic identity holds on any connected graph:

\[
\psi=Aa+SAb,\quad (a,b)\longmapsto(-b,a+2Tb),\quad
p_x=|a_x|^2+(T|b|^2)_x+2\Re(\bar a_x(Tb)_x).
\]

The second term uses the real stochastic `T` and includes every port.
Integer recurrences compute complete packet laws from equal amplitudes on
all twelve ports of one cell. The three periodic coordinate second
moments agree exactly. Their equality follows from the verified cyclic
coordinate symmetry, rather than an isotropic approximation.

For completeness, `S` has minus eigenspace dimension `M`. Since
`D>12`, `T` is positive, connected, and has a simple eigenvalue one.
The full `PD` dimensional walk has `W=+1` multiplicity `M−P+2`,
`W=−1` multiplicity `PD−M−P`, and `2(P−1)` dispersive eigenvalues.
The cyclic space has dimension `2P−1`; its orthogonal complement is
retained. Uniform event preparation has zero overlap with this
unprepared apparatus complement.

## Complete Bloch reduction

For `k_a=2πn_a/L`, the normalized cell characters tensor the twelve
internal vectors form a complete orthonormal basis. In this basis the
graph Laplacian is the direct sum of twelve-dimensional symbols

\[
H(k)=6I-A_{\rm ico}
 -\sum_{a=1}^3\sum_{(p,m)\in E_a}
   (e^{ik_a}|p\rangle\langle m|+e^{-ik_a}|m\rangle\langle p|).
\]

Each edge gives this identity directly by translation of its Fourier
character. All twelve branches are retained. The executable independently
applies the full graph matrix to every Fourier/internal column on
`L=2,3,4`; floating residuals below `2e−12` check the implementation,
while the translation identity establishes the exact statement.

At zero momentum, exact rational algebra gives spectrum
`0(1),4(3),6(2),8(6)`. Put `u=(1,...,1)/sqrt(12)`, `Q=I−uu*`, and
`G=(H(0)+uu*)^(-1)−uu*`. For axis `a`, let `j_a` equal minus one at
its two positive ports and plus one at its two negative ports. Then

\[
H(0)j_a=8j_a,\qquad j_a^Tj_b=4\delta_{ab},\qquad
j_a^TGj_b=\tfrac12\delta_{ab}.
\]

The first symbol derivative applied to `u` is `ij_a/sqrt(12)`.
The direct quadratic term is `|k|²/6` and its Schur correction is
`|k|²/24`, giving

\[
\lambda_0(k)=\tfrac18|k|^2+O(|k|^4),\qquad
\partial_a\partial_b\lambda_0(0)=\tfrac14\delta_{ab}.
\]

Two coordinate flips and cyclic coordinate permutation preserve all
internal and oriented external edges; their exact port permutations are
checked. Full rotational symmetry of the microscopic graph is not
assumed. A fixed alternative axis-zero partition selects the sign of
coordinate one and flips that coordinate across the external link. Its
ports retain degree six, but exact Schur algebra gives coefficient
`diag(4/33,1/8,1/8)`. This held-out architecture control resolves anisotropy
without changing any target or fitting a coefficient.

## Uniform fourth-order error

Let `r=max_a|k_a|≤1/2`, `V=H(k)−H(0)`. The external ends partition the
ports, so the varying part is a direct sum of two-dimensional edge
blocks. Therefore `||V||≤r`. Since `H(k)` is a sum of positive edge
squares,

\[
0\le\lambda_0\le a:=u^*Vu
 =\tfrac13\sum_a(1-\cos k_a)\le r^2/2,
\qquad \lambda_1\ge4-r.
\]

The lowest eigenvalue is simple. Write `b=QVu`, `C0=QH(0)Q` and
`E=QVQ−λ0 Q`, as operators on `Q`. The exact scalar Schur equation is

\[
\lambda_0=a-b^*(C0+E)^{-1}b.
\]

Split `V` into its real even part and purely imaginary antisymmetric odd
part. Their norms are at most `r²/2` and `r`; the odd part differs from
its linear Taylor term by at most `r³/6`. Thus `b=b_e+b_o`, with
`||b_e||≤r²/2`, `||b_o||≤r`, `||b_o−b_1||≤r³/6`, and `||b||≤r`.
Here `b_e` is real, while `b_o,b_1` are purely imaginary. Also
`||G||≤1/4`, `||E_e||≤r²`, `||E_o||≤r`. With `R0=(C0+E)^(-1)`,

\[
R0=G-GEG+GEGER0,\qquad
\|R0\|\le(4-r-r^2/2)^{-1}.
\]

The four errors relative to the quadratic Schur expression satisfy

\[
\begin{aligned}
|a-a_2|&\le r^4/24,\\
|b^*Gb-b_1^*Gb_1|&\le 7r^4/48,\\
|b^*GEGb|&\le r^4/8,\\
|b^*GEGER0b|&\le25r^4/864.
\end{aligned}
\]

For the second line the real/imaginary cross terms vanish. For the third,
the quadratic form of an imaginary antisymmetric matrix on a real
vector, or a purely imaginary vector, is zero; only the real/imaginary
cross terms of `E_o` contribute. The final bound uses
`(1+r/2)²/[16(4−r−r²/2)]≤25/864`. Their sum is `295/864<3/8`, proving

\[
\boxed{|\lambda_0(k)-|k|^2/8|\le(3/8)r^4.}
\]

Uniform internal preparation requires no exact spectral filter. Its
Rayleigh residual has norm at most `r`, and every other eigenvalue is
at least `4−r`, so its lowest-band weight is at least

\[
1-\frac{r^2}{(4-r-r^2/2)^2}.
\]

## Rational spectral intervals and refinement clock

The executable's outward interval certificates use only rational
endpoints. `π=16 arctan(1/5)−4 arctan(1/239)` is enclosed by alternating
series remainder bounds and rounded outward to sixty decimals. The
identity follows by applying the tangent addition formula to one quarter
of that expression: its tangent is one and its angle lies in `(0,π/2)`,
so the quarter expression is `π/4`. Sine and cosine use 24 Taylor terms,
the derivative-bound Taylor remainder, and the angle interval's
half-width; rounding to sixty decimals is outward. These operations
include all twelve symbol components.

A floating eigensolution selects a rational trial vector; it is not the
certificate. Exact interval arithmetic computes its normalized Rayleigh
interval `[q_l,q_h]` and a bound `σ²` for the residual at its Rayleigh
midpoint. This also bounds the variance at its true Rayleigh quotient.
With `β=4−r>q_h`, Temple's inequality gives the outward enclosure

\[
\max(0,q_l-\sigma^2/(\beta-q_h))\le\lambda_0\le q_h.
\]

No truncated eigenvalue is declared an exact pole. The receipts include
six held-out rational momenta and the first axis momentum on cubes of
sizes 16, 32 and 64. Those larger cubes use their complete Bloch symbol;
their global carriers are not enumerated.

On the cyclic space `cos θ=1−λ/D`, hence

\[
\frac{2\lambda}{D}\le\theta^2
 \le\frac{2\lambda}{D(1-\lambda/(2D))}.
\]

Declare cell spacing `ε=1/L` and time per coherent event `ε/sqrt(D)`.
For fixed Fourier momentum `p`, put `k=εp` and `Ω=θsqrt(D)/ε`. Then

\[
\Omega^2\longrightarrow |p|^2/4,
\]

with explicit error, writing `p_∞=max|p_a|`,

\[
\left|\Omega^2-|p|^2/4\right|
\le\epsilon^2p_\infty^4\left[
\frac34+\frac{1}{4D(1-\epsilon^2p_\infty^2/(4D))}\right].
\]

This follows by combining the fourth-order eigenvalue error with
`λ≤ε²p_∞²/2`. The optical branches have squared frequency at least
`2(4−r)/ε²`. The unprepared fixed event modes are retained separately;
no gap of the full event walk is inferred. For finitely many fixed
Fourier modes the displayed bound gives uniform dispersion convergence.
The constant internal preparation has vanishing optical weight in that
limit. A convergence theorem for all interacting fields or all apparatus
records is not asserted.

The dimensionless scalar speed is `1/2` under the declared clock and
spacing. Physical light speed, vector polarization, particle masses,
Einstein dynamics and a four-dimensional Yang-Mills construction do not
follow from this scalar dispersion theorem. The graph interface exposes
oriented seams, neighbours and cells for a joint charge/Gauss program;
that program's source, preparation and complete instruments require their
own transfer proof on this declared routing.
