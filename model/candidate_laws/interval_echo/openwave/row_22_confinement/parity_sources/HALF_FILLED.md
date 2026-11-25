# A distance-uniform work bound with every native contact active

Fill exactly half the ports with double loads and the remainder with vacancies. Keep every seam memory blank and open **every native internal contact, external contact and degree-five private-reference service** of the registered Interval Echo material. This is an invariant sector because its blank generator only swaps two material loads or multiplies their amplitude. Coherent combinations of all words with that global inventory are admitted. No per-owner inventory constraint enters this argument.

Let H be this Hamiltonian, E₀ its ground energy and Eeven the minimum on the subspace where every surveyed face reads even. The latter subspace is nonempty, since each owner can park six doubles at private ports. For distinct faces f and g, their projectors are diagonal in the ordinary port load records. Define P to accept precisely when both f and g are odd, with all other face outcomes allowed. Define Epin as the minimum Rayleigh quotient of H on the range of P. This is the energy of the specified variational hard pin, equivalently the ground energy of its ideal projected generator. An unconstrained trajectory need not preserve P.

**Theorem.** On any such board with an even number N ≥ 240 of ports, maximum port degree six, and distinct 5- or 6-sided simple face meters sharing at most one exterior bridge, Epin − Eeven < 19,129 native energy units, independently of source separation. The argument also bounds Epin − E₀ by the same constant. Thus a static inter-parity ground energy growing without bound would require a more restrictive excitation sector or a different source read. The bound is deliberately coarse and allows an ordinary finite positive linear window.

### Pull the local source read back to a random same-number preparation

Measure bit 1 on load 2 and bit 0 on load 0. The parity of one face is the XOR of these bits over its gateway support F. Every boundary bridge contributes one bit at each end, all distinct within the face. The sizes of F and G are 10 or 12, their symmetric difference is nonempty and even, and s = |F ∪ G| ≤ 24.

Draw a uniformly random ordered sample of s distinct donor ports. Extend the bijection from the read support onto this sample to a permutation supported on their union, by closing the injection's open chains into cycles. The corresponding port permutation U (choosing its inverse if needed) pulls back the measured bits to the sampled donors and preserves global load exactly. Its support contains at most 2s sites regardless of how far the donors are apart.

For **every** initial state at half filling, averaging ⟨U†PU⟩ over the donor sample gives the probability of both faces reading odd in a uniformly sampled fixed-weight bit word. This holds first for each native basis word, because it has N/2 bits of each value. It extends to arbitrary coherent states since U†PU is diagonal for each permutation. No assumption of homogeneous occupation, finite correlation length or probability independence is needed.

For a uniformly sampled set of m ports the average sign product is

\[
 \chi_m=\frac{[z^m](1-z)^{N/2}(1+z)^{N/2}}{\binom{N}{m}} .
\]

It vanishes for odd m. For m = 2r it equals (−1)ʳ C(N/2,r)/C(N,2r). At m=2 its magnitude is 1/(N−1). The ratio of successive even magnitudes is (m+1)/(N−m−1), so |χₘ| ≤ 1/(N−1) for even 2 ≤ m ≤ N/2. The required three supports F, G and their symmetric difference satisfy this bound. Expanding the product of both odd projectors yields

\[
 \mathbb{E}_U\,\langle U^* P U\rangle
 = \tfrac14(1-\chi_{|F|}-\chi_{|G|}+\chi_{|F\triangle G|})
 \ge p_0=\tfrac14(1-3/(N-1)).
\]

Thus there is a finite-support permutation with success probability p ≥ p₀, even for a highly inhomogeneous coherent many-body ground ψ. Conditioning on P is a mathematical construction of an admissible trial vector. It does not require that the observer can address an unknown ground state or that the native unforced trajectory transports distant donors quickly.

### Cost of the projected trial

Put Y = PU. Its operator norm is at most one. The normalized Yψ has the correct total particles and both required measured parities. Since Hψ = E₀ψ,

\[
 E_{\rm pin}-E_0 \le
 \frac{\langle Y\psi,[H,Y]\psi\rangle}{\|Y\psi\|^2}
 \le \frac{\|[H,Y]\|}{\sqrt{p_0}}.
\]

The commutator vanishes for all terms disjoint from Y. Every port has six service slots. A degree-five port has five material contacts and one private-reference term κn²; a degree-six port has six material contacts. At most 2s × 6 terms of either kind meet the support. On the blank binary local material sector a native seam has spectrum 0, 0, 32, 33. Subtracting 16.5 times the identity from that term bounds the norm of its commutator with Y by 33. The private reference has binary eigenvalues 0 and 32 and contributes a commutator of norm at most 32 by the same centering argument. This step holds for a coherent, globally entangled ψ. Applying the variational inequality Eeven ≥ E₀ gives

\[
 E_{\rm pin}-E_{\rm even} \le E_{\rm pin}-E_0
 \le \frac{2s\cdot6\cdot33}{\sqrt{p_0}}\le 19{,}128.433.
\]

This is an upper bound. It neither estimates the actual work nor supplies a lower bound on Epin − Eeven when the two variational classes differ. The exact exterior-only minima and the much tighter dilute bounds are independent results. String breaking can therefore coexist with a finite apparent string tension. A bound on these binary hard-source measurements does not exclude physical colour confinement in another sector.

### Checks

The certificate surveys all 66, 861 and 13,041 distinct face pairs at the three registered board sizes. It stores each pair class, exact fixed-weight probabilities and the maximum local generator width. The per-board upper bounds are 15,940.37, 19,037.81 and 19,015.44. A small coherent many-body test independently builds every binary word of a six-port half-filled ring, propagates all source pulls over the complete permutation group, checks the mean against exact mask counts, projects the full native ground, and compares its Rayleigh quotient with direct constrained diagonalization. The uniform theorem only uses the geometry and inventory hypotheses stated above. Particle identification, a measured crossover and the other global inventory sectors remain open. [Frozen question](HALF_FILLED_PROTOCOL.md) · [Exact combinatorial certificate](half_filled.json) · [Certificate generator](half_filled.py) · [Tests](test_parity_sources.py)
