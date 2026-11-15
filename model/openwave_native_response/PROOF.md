# The exact native response and its cell memory

Fix the supplied connected port-pair federation at level 0, 1 or 2, with
N=240, 960 or 3840 ports and M=630, 2520 or 10080 seams. Every attempted
event independently selects an existing seam uniformly, chooses C or E with
probability 1/2 and a fair endpoint coin, and executes the unchanged native
repair. The complete instrument has 25 local slots, including zero hold
and padding. The quantum carrier is the entire global one-load sector.
No per-cell occupation constraint or target Hamiltonian is imposed.

## Complete source operator

Let L be the unweighted graph Laplacian and A its adjacency matrix. The
exact source marginal after the raw event record is unread is

\[
 (\Phi\rho)_{ii}=[T\operatorname{diag}\rho]_i,\quad
 T=I-L/(2M),\qquad
 (\Phi\rho)_{ij}=\left[1-\frac{3(d_i+d_j)-2A_{ij}}{4M}\right]\rho_{ij}
 \quad(i\ne j).
\]

For one occupied endpoint, a fair repair stays or traverses the seam with
probability 1/2. On an offdiagonal matrix unit an edge meeting neither
address preserves it; an edge meeting exactly one address multiplies it
by 1/4 in the fair C/E mixture; an edge joining both addresses kills it.
There are M-d_i-d_j+A_ij, d_i+d_j-2A_ij and A_ij such edges. Summing them
proves the formula. The executable checks all 16 endpoint matrix units
for each of C, E and their mixture, including the empty and double-load
local complement, and all nine one-load endpoint/spectator columns. The
graph check accumulates every population column from actual branch maps
and checks the incidence partition of every offdiagonal column.

This is an identity of linear maps, hence remains exact on an arbitrarily
entangled idle reference. Kraus completeness is supplied by the actual
25-slot parent instrument. The source marginal does not replace or erase
the observer's retained event history. Independent seam selection makes
this marginal sufficient for the source response; history-dependent
feedback would require its joint source-memory channel instead.

The maximally mixed state I/N is stationary. The channel is already linear,
so its exact derivative at that state is Phi. The existing binary-sector
quantum mixing theorem gives its complete spectral decomposition without
a fresh diagonalization. Here the population eigenvalues are at least
1-d_max/M>0 by the graph Laplacian bound L<=2d_max I. The offdiagonal
multipliers are also strictly positive and smaller than one.

## Exact failure of a memoryless cell reduction

Let C sum the twelve populations in each cell and let S put a cell's
population uniformly on its twelve ports. Then CS=I, P=SC is the orthogonal
cell-constant projection and Q=I-P. Define A_c=CTS. An exact memoryless
intertwiner would require CT=A_c C on every source basis column.

Each cell has three boundary ports with one outgoing seam each. A boundary
port sends 1/(2M) into its particular neighboring cell; the uniform lift
sends 1/(24M) into each of the three neighboring cells. The exact maximum
column l1 discrepancy is

\[
 \|CT-A_cC\|_{1\to1}=\frac{11}{12M},\qquad
 2M\|CT-A_cC\|_{1\to1}=\frac{11}{6}.
\]

The code checks every column on each of the three carriers. A uniform
initial cell also produces a nonzero exact second-step memory term
CTQTS; its complete rational vector is included in the receipt. This
does not obstruct a source reduction with memory. Eliminating Q exactly
gives, for c_n=CT^n S c_0,

\[
 c_{n+1}=A_c c_n+
 \sum_{j=0}^{n-1} CTQ(QTQ)^{n-1-j}QTS\,c_j.
\]

Every complementary port mode appears in this finite expression. It is a
direct shortcut for subsequent coarse response calculations.

## Relation to the twelve-direction cosine object

The M10.1 arithmetic concerns the declared translation symbol
sum_u[1-cos(a u.q)]/(2a^2). Its exact icosahedral expansion is valid.
The native channel above supplies no translation representation, plane-wave
intertwiner or second-order wave-time law identifying that symbol with
its response. The failed natural cell intertwiner is a quantitative
counterexample to one proposed identification, rather than a theorem
excluding every observer readout or every OPH implementation.

No continuum, physical particle, mass, GeV or OpenWave criterion is
validated by this finite operator identity alone.
