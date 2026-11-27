# Energy, connected support and two triangular cores

For binary loads on blank Interval Echo seams the native 99-state interaction reduces to half the identity minus half the exchange, plus \(\kappa(n_a-n_b)^2\). The reduction holds on the whole blank sector. Summing every edge and all \(6-d_a\) private zero references gives
\[
H=K+8(6N-2e),\qquad K=\tfrac12\sum_{ab}(I-S_{ab})\succeq0 .
\]
Here \(e\) counts occupied seams in a sharp word. All energy comparisons include the same absolute reference. A separated comparison has independently occupied disjoint copies of the same carrier and independent Hamiltonians added on their product. A disconnected word on one carrier concerns the components of its occupied induced graph. These two notions give the same lower bound by positivity, including superpositions of such words. [Blank interaction](../../../contrast_memory_spherical_v1/law.py) · [Graph certificate](certificates.py)

## Binding through unrestricted release

Counting the subsets of one icosahedron gives maximal occupied-edge counts \(0,1,3,5,7,10\) for \(1,\dots,6\) loads. If \(n\) occupied sites occupy several cells, with \(n_i\) in cell \(i\), intercell edges form a matching. Each consumes at least one site outside a largest cell. At most \(\min(\lfloor n/2\rfloor,n-\max_i n_i)\) intercell edges can therefore augment the within-cell sum. Enumerating the integer partitions through six preserves all six maxima. These bounds apply to every federation made from these cells with intercell edges forming a matching.

At total load six, every word has potential at least 128. For every proper partition of six the sum of the part maxima is at most seven, attained at \(5+1\). Thus the disconnected potential and the energy floor of every independent fragmented stock are at least 176. Specifically \(3+3\) gives at least 192. The polar cap on the single cell has ten occupied edges and ten crossing edges, so its diagonal trial energy is \(128+10/2=133\). On each deposited 20-, 80- and 320-cell federation the same cap has twelve crossing edges and mean energy 134 with all links active. The variational principle places an eigenvalue at or below those energies. At every free time a sharp cap has disconnected probability \(p\) obeying \(134\ge128(1-p)+176p\), hence its connected probability is at least \(7/8\). The single-cell bound is \(43/48\). This proof refers to connectivity of the occupied support at one terminal time and permits changes among connected shapes. Independence of the exterior graph check from any six-load propagation makes the free-time bound applicable to each certified deposited board.

With \(\kappa=0\) the uniform vector is a zero-energy ground state in each binary sector, hence every separation gap vanishes. For \(\kappa=-8\), \(H=-48N+K+16e\). The six singletons reach \(-288\). The connected configuration graph makes the kernel of \(K\) uniform. That vector has positive weight on words with \(e>0\), so six loads on one cell cost strictly more than the six separated singletons. Time-local persistence of the sharp cap under this law occurs high in energy. [Finite floors and controls](records.json.gz)

## Exact two-face bound for the entire lowest band

Let \(C\) project on the twelve caps, \(F\) on the hundred six-load words partitionable into two disjoint elementary faces, and \(R=I-C-F\). The configuration census gives \(CF=0\) and \(FHC=0\). The cap compression is \(133 I_C\). The complement \(D=I-C\) has potential at least 144, so \(DHD\) on its support is at least 144. By the min–max principle exactly twelve eigenvalues lie below 144 and the top of these is at most 133. Denote their spectral projector by \(P\) and the restricted energy by \(H_P\).

The matrix equation \(H_D DP-DP H_P=-DHC\,CP\) on the two spectral domains has the norm bound \(\|DP\|\le\|DHC\|/11\). Integrating \(e^{-tH_D}(-DHC\,CP)e^{tH_P}\) over positive \(t\) gives this estimate independent of the number of modes. The same argument on \(F\) gives \(\|FP\|\le\|FHR\|\|RP\|/11\) and \(\|RP\|\le\|DP\|\). A cap has ten outgoing edges, and a noncap meets at most one cap. From \(F\) to \(R\) the two maximum degrees are eighteen and six. Each nonzero matrix entry has magnitude \(1/2\). The geometric mean of the maximum absolute row and column sums bounds each operator norm. Consequently every normalized coherent state in the band obeys
\[
\langle F\rangle\ \le\
\frac{(10)(1)(18)(6)}{16\cdot 11^4}
=\frac{135}{29282}<0.004611 .
\]
The guarantee also holds for mixtures and when the preparation is entangled with a reference. The computed largest eigenvalue of \(PFP\) on its range is approximately \(0.000018394\). The rational bound is independent of the eigensolver. [Contact-by-contact test](test_binding.py)

## An oriented load quadrupole and a finite energy aperture

In the unit-radius icosahedral frame let \(r_j\) denote a port and \(u\) a unit polar axis. The cell-centred read is \(Q(u)=\sum_{j\ {\rm occupied}}[3(u\cdot r_j)^2-1]\). Inversion fixes each port weight and the sum of the twelve weights vanishes. Complementing the load bits reverses \(Q\). At half filling on a degree-five cell each word and its complement have the same number of occupied seams, since both have six ports and the same cut size. The Hamiltonian commutes with complement and with inversion. Their product fixes the sharp cap's word. It pairs equal-probability terminal records with opposite \(Q\) after every free delay and after any scalar energy aperture. Thus the *oriented* cap has mean cell-centred moment zero about every axis. The two-face source can break the symmetry. Its initial axis-1 reading is \(+2.4\).

At equal constituent masses the word centroid is \(\bar r=\sum_j r_j/6\). Recentring the same spatial moment gives
\[
Q_{\rm cm}(u)= Q(u)-6\{3(u\cdot\bar r)^2-|\bar r|^2\}.
\]
This added term can be anisotropic under an oriented cap preparation. Its mean has a negative sign on the registered polar axis both before and after the low-energy tag. The centroid convention and positive material charge enter this interpretation explicitly. A nuclear spectroscopic quadrupole would use a spin-polarized physical nuclear state and its electric charge about its mass centre.

For timing arms \(j=0,\ldots,8191\) the bright amplitude is their uniform coherent average at phase \(-j(E-f)/128\). Over the enclosed positive-contrast spectra the native Hamiltonian norm is at most 318 for the six-load cell and the combined norm on two separate cells is at most 348 (two copies with a 30-edge kinetic norm bound each). For either apparatus every eigenvalue from 160 to its respective maximum has bright probability at most \(1/[8192\sin((160-f)/256)]^2 <1.338\times 10^{-6}\) at \(f=132.924\), by the geometric sum. At 348 the sine is larger than at 160. The cutoff 160 is a conservative lower comparison threshold. It includes the proved fragment floors 176 and 192. A rare terminal success through this aperture alone has no guaranteed band assignment. The two-face run retains the full terminal word after every bright outcome and the number of dark refusals. [Timing and material readers](apparatus.py)
