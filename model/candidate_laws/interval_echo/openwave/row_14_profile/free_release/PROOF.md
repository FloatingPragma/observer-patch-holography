# A signed radius bound for unrestricted free flight

The source has unit filling except for two ports at load two and two at load zero, initially in owner zero. All seam memories are blank. The terminal observer assigns q=n−1 to its ternary reading. Distances d are shortest-path distances on the graph of owners from the named emitter. Comparison of the observing stations' enrollment records fixes this chart. The signed second moment is the expectation of the diagonal operator M = Σ q d². The source has total signed charge zero, and M annihilates its initial sector. Any initial density matrix supported there is allowed, including one entangled with a reference. No support or sign measurement occurs during the full flight.

## The regenerated Hamiltonian

Let λ multiply exchange on every seam. Retain all native contrast terms at κ=8 and every private degree reference. The contrast contribution, relative to the uniformly filled host and in units of κ, is

$$
\sum_{ab}(q_a-q_b)^2+\sum_a (6-d_a)(q_a^2+2q_a)
=24-2\sum_{ab}q_aq_b+2\sum_a(6-d_a)q_a .
$$

Here d_a is the number of neighbors of a port. It is five or six. Four opposite-sign pairs and two same-sign pairs are possible and the last unsigned coefficient sum over the signed defects lies in [−2,2]. Thus contrast lies in [128,288]. Each word has at most 24 nontrivial swaps. Exchange is half the positive graph Laplacian on these words, bounded above by 24λ by row sums. Subtract the midpoint of the resulting energy interval. The shifted Hamiltonian H then has norm at most R = 80+12λ on the *whole* four-defect board, independently of its size. The shift is a common phase in this inventory.

Let P project onto the words with all defects in owner zero, Q=I−P, H₀=PHP+QHQ, and V=H−H₀. Only three outgoing swaps occur in QHP. A word reached via one of them has a unique outside defect and its gateway fixes the inverse swap. The images of distinct allowed moves are therefore orthogonal, even across initial words. It follows that ‖V‖=‖QHP‖≤b=λ√3/2. For a specified receiving port the retained block has norm at most λ/2. For any normalized initial vector, the leading expectation of its receiving charge divided by t² is λ² times the mean signed incident gateway load divided by four. The leading expectation for M is the sum of these three terms. The same conclusion holds for mixed sources by linearity.

Duhamel's formula for U=exp(−itH) against the block-diagonal U₀ gives

$$
\|QUP+itQHP\|
\le \int_0^t b[R(t-s)+Rs]\,ds=bRt^2 .
$$

The inequality uses ‖U(s)−I‖≤Rs for both generators. For any terminal diagonal observable O that vanishes on P, divide its expectation by t². Its deviation from the leading value is at most 2‖OQHP‖bRt + ‖O‖(bRt)². For a single receiving port these two norms are bounded by λ/2 and 1. For the total exterior charge they are bounded by b and 2. For M they are bounded by b and 2D², where D is the maximum distance on this board. The receipt retains this diameter-only bound.

## The outer shells

Let Qⱼ project onto words with at least one defect at owner distance j or greater. A path that has crossed to a new owner arrives on a port with only one exterior seam, which returns to the preceding owner. Before the next outward crossing it must make at least one internal swap. Consequently, to reach shell j requires at least 2j−1 generator applications. Specifically, every vector in VP has one outside defect on a first entry port and all other defects at home. On this range QⱼHⁿ vanishes for n < 2j−2. Diagonal terms and exchanges between defects cannot move a defect farther in fewer swaps. The statement follows by tracing any defect at the final occupied port back through its swaps.

Expand the left propagator in the Duhamel integral. The tail norm after removal of the zero powers is bounded by exp[R(t−s)] [R(t−s)]^{2j−2}/(2j−2)!. The block-diagonal right propagator preserves P and has norm one. With z=Rt<1, integration gives

$$
\frac{\|Q_j U P\|}{t}\le \frac{b\,z^{2j-2}}{(1-z)(2j-1)!} \qquad (j\ge2).
$$

The telescoping identity d²=Σ_{j≤d}(2j−1) writes M as total exterior charge plus weighted signed charges beyond each subsequent shell. A signed charge there is at most two in absolute value and vanishes outside Qⱼ. The sum of its absolute error bounds, divided by t², is at most

$$
\frac{6}{1-z^4}\left(\frac{b z^2}{6(1-z)}\right)^2 .
$$

Indeed the ratio of consecutive squared-tail terms including their shell weights is at most z⁴. This estimate is independent of D. The code intersects it with the diameter-only estimate. It substitutes the rational upper bound 0.866026λ for b and rounds the printed endpoints outwards. At zero exchange every word stays at its original sites at any time.

## What fixes the input and the slope

Independent sharp port banks from the guarded physical source give the three gateway means. Each has over four million accepted starts. For a score in [−1,1], the conditional two-sided Hoeffding error at radius 0.002 is at most 2 exp(−8) per bank, and the stored sampling roundoff is negligible at this budget. More precisely all three actual counts exceed ten million, so the same error allowance with radius 0.00199 is below 6×10⁻⁹ per bank. The extra 0.00001 covers the parent sampler's numerical normalization on these marginal bins. The archived sharp expectations also agree with these intervals directly. These are statistical input bounds for a numerical finite-source experiment. The propagation inequality itself is deterministic. The four first-half-period dimer clock counts use separate two-sided binomial intervals with one-sided tail 0.01/32, enclosing their individual probability-grid errors. The union for clocks and incident populations has coverage above 99%. No hidden evolution phase enters the clock decoder.

The owner-based cosine response is F(k)=E[Σ q cos(kd)]. Finite support and zero total charge give F(0)=0 and the derivative in k² at zero equals −M/2. This slope can be inferred from the spatially resolved ternary record. The graph ruler and launch center fix its convention. A three-dimensional electric form factor would require physical scattering momenta, an intrinsic center or rest frame and an electromagnetic coupling. Source mass is absent from the radius formula.

Two planted sharp four-defect preparations test the sign interpretation on the same free law. The independent sparse word expansion retains every path of length at most four on the actual 240 ports. Hermitian Taylor remainder gives a norm error bounded by (Rt)⁵/5!. Multiplying by the norm of the polynomial's exterior component bounds the error in the signed moment. This independent computation uses floating-point word sums and an analytic truncation bound. The main full-source enclosures above use rational arithmetic on the sharp count intervals.
