# Electric energy and a native response geometry

The source is the twelve-port integer repair primitive and its compiled
proper-code quantum services. The charged sector, scalar cutoff, interaction
schedule, fresh preparations, pointer and scheduler table are specified
operations. The resulting response matrix and its spectral geometry are
derived from their complete outcomes. The construction uses one fixed
instrument for every matter input below.
The full 78-state integer instrument is an admitted native primitive;
its compilation into the six-level proper-code processor is not supplied.

## The fixed source instrument

Use the four-qubit register of `common_source.scalar_coupling` at scalar
dimension two. The charged Gauss-sector bits are 0 and 1; scalar bits are 2
and 3. Fix charged bit 1 and scalar bit 3 in their blank states. On this sector
the source electric energy is

\[
 E^2=(-n_0+n_1)^2=n_0,\qquad H_{\rm int}=\tfrac14 X_2 Z_0.
\]

Prepare scalar bit 2 in \((|0\rangle+i|1\rangle)/\sqrt2\). Execute the
existing scalar/charge interaction subword for time
\(4\alpha\), where \(\cos\alpha=4/5\), \(\sin\alpha=3/5\), and then the
existing native occupation pointer on bit 2. This is the emitted interaction
subword, with the other Hamiltonian subwords absent from this specified
protocol. It is not evolution under their sum. The exact unitary is
\(U=(4I-3iX_2Z_0)/5\). Its two terminal charge maps, with the scalar and
pointer output included, have charge factors

\[
 K_0=(4I+3Z_0)/\sqrt{50},\qquad
 K_1=i(4I-3Z_0)/\sqrt{50}.                         \tag{1}
\]

They obey \(\sum_bK_b^*K_b=I\). For any charge density matrix, including
coherences or an untouched reference,

\[
 p_1=\operatorname{Tr}\rho K_1^*K_1
     =\tfrac1{50}+\tfrac{24}{25}\operatorname{Tr}\rho E^2.       \tag{2}
\]

The output on branch b is the unnormalized \(K_b\rho K_b^*\), tensored
with its scalar and pointer output. The record and this work state are kept.
Both Kraus factors commute with \(E^2\), so the unconditional electric-energy
expectation is preserved. This is a driven measurement protocol; conservation
of an autonomous matter-plus-geometry Hamiltonian is not asserted.

These are complete outcomes for the exact specified proper-code controls,
whose complementary code-failure amplitude is zero. The Python control
contracts the raw six-level source pulses, compares their restricted work
columns with (1), and checks their numerical leakage through the compiler.
It uses no polar stabilization. Those binary64 comparisons are numerical
controls, not interval certificates for rounded analogue controls or a
physical service-failure instrument. Exact control identities establish the
ideal-source theorem; physical and rounded-control errors require separate
bounds and must retain their complement outcomes.

Eight fresh fair coins generate w in {0,...,255}; a ninth fair coin is the
native placement coin. Each coin has an emitted proper-code realization:
an X pulse of angle pi/4 on a blank bit, followed by its native CNOT pointer,
gives \((|00\rangle+i|11\rangle)/\sqrt2\). Keep all coin outputs, including
the placement output on idle attempts. Let the thirty native seams have
their source order, whose first two entries are E=(0,1), F=(0,5). Initially
assign eight consecutive words to each seam and words 240,...,255 to idle.
For b=0 reassign word 8 from F to E; for b=1 reassign word 0 from E to F.
The table depends only on the committed b and w, never on hidden native
loads. This gives conditional slot counts (9,7) or (7,9) on E,F, eight on
every other seam, and sixteen idle slots.

On a selected seam perform the unchanged native integer repair, retaining
its placement coin, signed current and final work state. On an idle keep
the work state and an explicit idle record. Every source work input has
1024 complete branches (two b values, 256 words, two placement values).
`source.py` executes all 79,872 such branches on the 78 total-two states.
Every branch has the same symbolic quantum-output reference (1) for its b;
none of the output charge/scalar states is substituted by its probability.

## Mean dynamics and its clock

For a native load column x, averaging the two placement coins on seam e
gives \((I-D_e/2)x\), where \(D_e=(e_a-e_b)(e_a-e_b)^T\). This follows
directly by averaging floor and ceiling of the endpoint mean, and applies
to every nonnegative integer load. Write \(n=\operatorname{Tr}\rho E^2\).
The source table and (2) therefore give

\[
 M_n=I-\frac{L_\theta}{64},\quad
 L_\theta=L+\theta(D_E-D_F),\quad
 \theta=\frac3{25}(1-2n),\quad 0\le n\le1.                 \tag{3}
\]

Here L is the actual unweighted thirty-seam Laplacian. Formula (3) is
checked on every total-two input from its complete native histories.
The mean law can be iterated when every attempt receives an independent
fresh probe with the same density matrix, independent of the native work
and all prior records. The native work can persist and become correlated
with its complete past. Fresh independent preparation is a specified
service. Reusing a probe, or providing probes correlated with native work,
does not satisfy this hypothesis; in that case the conditional joint law,
rather than a fixed marginal matrix M_n, governs the response.
For example, correlate charge 0 with two native loads at port 0 and charge 1
with two loads at port 1, with equal probabilities. The charge and native
marginals equal those of their product control, but the actual next load
mean differs by \((-3/1600,3/800,0,0,0,-3/1600,0,\ldots,0)\).
This is checked directly from the same native histories.

One attempt is one quantum readout, eight scheduler coins and one placement
coin, including idle attempts. Its interval is a specified model clock.
The mean load is \(M_n^t\mathbb E[x_0]\) for integer t under the fresh-probe
hypothesis. All thirty seam probabilities are positive. The active attempt
probability is 15/16, independent of n. Conditioned on activity, the local
pair-clock reconstruction in `relative_observer_geometry` recovers the
exact same seam weights from its 132 oriented count/coincidence laws.
The full attempt normalization additionally uses the retained idle record;
local endpoint clocks alone do not observe that record. Absolute loads
need a known preparation; the retained seam/current always determines
each relative load increment.

There is an exact response of the native kinetic quadratic form:

\[
 \frac{d}{dn}\langle f,(I-M_n)f\rangle
 =-\frac3{800}\big[(f_0-f_1)^2-(f_0-f_5)^2\big].          \tag{4}
\]

This is a source-derived repair response, not an identification of a
Lorentzian action or conserved spacetime stress tensor.

## A uniformly defined rank-three frame

The pinned native Laplacian has spectrum
\(0,5-\sqrt5,6,5+\sqrt5\), of multiplicities 1, 3, 5, 3. Put
\(D=D_E-D_F\). Direct multiplication gives \(D^3=3D\),
\(\operatorname{Tr}D^2=6\), and \(D\mathbf1=0\), hence
\(\|D\|=\sqrt3<2\). On the larger interval \(|\theta|\le1/8\),
\(\|\theta D\|\le1/4\). The three low positive eigenvalues lie in
\([5-\sqrt5-1/4,5-\sqrt5+1/4]\); the next lies above 23/4.
The separation is at least \(\sqrt5+1/2\). Thus one complete rank-three
band, rather than a choice of its potentially split eigenvectors, is isolated
for the entire source family.

Let \(P_\theta\) be this orthogonal projector. On \(\mathbf1^\perp\),
integrate its resolvent around the circle centered at \(5-\sqrt5\),
of radius 3/2. The unperturbed resolvent has norm at most 2/3 and the perturbed
one at most 4/5. The resolvent identity gives

\[
 \|P_\theta-P_0\|\le1/5,\qquad
 (P_\theta)_{ii}\ge1/4-1/5=1/20.                         \tag{5}
\]

Consequently
\[
 G_{ij}(\theta)=\frac{(P_\theta)_{ij}}
 {\sqrt{(P_\theta)_{ii}(P_\theta)_{jj}}}                 \tag{6}
\]
is an unambiguous positive rank-three Gram matrix with diagonal 1. At zero
it equals the pinned native icosahedral Gram. A coordinate factor of G is
unique up to an orthogonal transformation; no coordinate axes are supplied
to define its response. This carrier is a dimensionless response geometry.
The three eigenvalues can split, so a normalized infinite-time kernel that
keeps only a single slowest eigenvalue is a different object.

Let \(P_k\) denote the four unperturbed band projectors and
\(\lambda=5-\sqrt5\). The source-specific derivative is

\[
 P'_0=\sum_{\lambda_k\ne\lambda}
 \frac{P_kDP_0+P_0DP_k}{\lambda-\lambda_k},\qquad
 G'_{01}(0)=\frac{-33+29\sqrt5}{100}>\frac3{10}.         \tag{7}
\]

`geometry.py` evaluates all 144 entries over Q(sqrt5), rather than by
differencing eigenvectors. The matter-energy derivative multiplies (7)
by \(-6/25\), using the same fixed instrument (3).

## A rational finite comparison

For rational theta, define \(X=(2L_\theta-9I)/9\). All its eigenvalues
lie in [-1,1] and have absolute value at least 1/4. The negative subspace
consists of the constant mode and the three-dimensional band. Therefore

\[
 P_\theta=(I-\operatorname{sign}X)/2-\mathbf1\mathbf1^T/12.
\]

Apply \(x\mapsto(3x-x^3)/2\) nine times. If
\(e=1-x^2\), then \(e_{next}=e^2(3+e)/4\le e^2\), so the exact
polynomial's sign error is at most \((15/16)^{512}\). For matrix
calculation use integers scaled by \(Q=2^{80}\), rounding each entry
once at initialization and once after each cubic update. An entry error
at most 1/(2Q) has operator norm at most 6/Q. On the domain with total
rounding error at most 1/4, the cubic has operator Lipschitz constant at
most 4, including noncommuting matrix arguments. Induction bounds that
error after k updates by \(2(4^{k+1}-1)/Q\). The resulting projector
operator error is at most

\[
 \epsilon_P=\tfrac12\big[(15/16)^{512}
                    +2(4^{10}-1)/2^{80}\big].           \tag{8}
\]

All diagonal products and off-diagonal numerators are enclosed using
(8). Integer square-root bounds enclose the positive denominators in
(6). The complete 144-entry endpoint enclosures certify

\[
 G_{01}(n=0)-G_{01}(n=1)>3/40.                           \tag{9}
\]

No numerical eigensolver enters this certificate. An independent test
compares its enclosures with floating spectral projectors as a control.
All exact matrix inequalities, source execution and source hashes are
replayed by the package receipt.

## Finite observation and scope

The same contour bounds give \(\|P'\|\le48/25<2\). Positivity of P
and (5) then bound every \(|G'_{ij}|\) by 80. Estimating the matter
record probability p1 to delta and clipping to [1/50,49/50] changes
theta by at most delta/4, so the recovered Gram error is at most 20 delta.
For each of the two endpoint preparations take 12,000,000 independent
fresh-probe attempts. With delta=1/2000, Hoeffding and the union bound
give failure probability at most \(4e^{-6}<1/100\). Thus each recovered
G entry has error at most 1/100, and the endpoint contrast exceeds 11/200
with that confidence. A rational partial sum of the exponential series
certifies \(e^6>400\). This is a specified finite-cost experiment, not
a claim that those 24,000,000 attempts have been executed by the receipt.
Fresh probes, their preparation resources and all nine fair-coin services
per attempt are included in its interface; physical storage and laboratory
service errors need separate bounds.

This establishes an electric-energy-responsive native geometry, with a
common finite history, model clock and quantified readout. It supplies no
Lorentzian metric, spatial refinement, action-stress identification, regional
entropy/area stationarity or Einstein source instance. No state/algebra map
to the spatial fermion-core or defect packages is asserted. Those missing
maps and equations are necessary for a common interacting field and
responsive Einstein geometry.
