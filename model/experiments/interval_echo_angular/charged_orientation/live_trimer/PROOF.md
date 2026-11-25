# Finite spectral and detector statements

## Full written-memory generator

Let $W_0$ be the integer contrast write pattern. Its normalized columns $W=W_0/\sqrt2$ are orthonormal. Swapping the two source loads intertwines with reversal $R$ of the signed written contrast, and the source squared contrast intertwines with the diagonal $D^2$ on written labels. Hence the blank and written regimes reduce the physical echo generator. At registered contrast $\kappa=8$, eight times its blank block is $4(I-T)+8\kappa C$, with $T$ the load swap and $C$ the diagonal squared difference. Eight times the full written block is the integer matrix
$$
(2I-W_0 W_0^T)[(I-R+2\kappa D^2)\otimes I_9](2I-W_0 W_0^T).
$$
For $m>0$, write a marked memory as direction coin $c$ and the sign combination $(|c,m\rangle+s|c,-m\rangle)/\sqrt2$, where $s=\pm1$. For $m=0$ there are only the two direction coins. Together with the blank block this gives $9+18+4(18)=99$ columns. The verifier uses the unnormalized integer columns $B$, evaluates the compression $A=(B^TB)^{-1}B^T(8H)B$, and tests $8HB=BA$ exactly on each block. It checks agreement with the full native generator and verifies orientation reversal $(a,b,c,d)\mapsto(b,a,1-c,-d)$ as an involutive symmetry. [Full matrix calculation](measure.py)

| Native sector | Every nonzero off-diagonal matrix element of $H$ |
| --- | ---: |
| Blank | $-1/2$ |
| Written $m=0$ | Empty |
| $m=1$, even or odd | $-4$ or $-4.5$, respectively |
| $m=2$, even or odd | $-16$ or $-16.5$, respectively |

For any fixed assignment of these sectors, tensor the local sign bases across all seams. Every term has nonpositive off-diagonal entries in the resulting global basis, including when two terms connect the same pair of configurations. All private references are diagonal. On one connected configuration component, choose $a$ larger than every diagonal entry. Then $aI-H$ is nonnegative and irreducible. Perron–Frobenius gives the simple positive ground vector, with no assumption about dimension or load composition. Any unitary preserving that line acts there by a phase. The commutator of two such unitaries is $+1$. This theorem concerns one connected component in a specified material sector. It gives no bound on degeneracies between components or in excited bands.

Contrast parity is a local conserved observable on the marked seam itself. Evolving a superposition of source parities produces two branches in orthogonal eigenspaces of that same head at all times. For an instrument supported outside the head, every off-diagonal matrix element between the branches vanishes. Consequently the destination channel from this parity encoding factors through measurement of that parity, even if motion transfers the entire net charge. The control in the test suite also evaluates this cancellation in the original signed-label basis for four input phases.

## Charge transfer and entanglement test

The six-port net-negative composition has 60 ternary words (two vacancies and one excess). With three live unit-contrast marks their direction coins give 480 configurations. The two ordinary orbital inputs occupy the same 168-word connected component. For each final direction-coin word $c$, let $K_c$ be the three-by-two matrix of amplitudes for placing the single excess at one of the three destination ports, both other destination ports vacant and the source filled. Tracing the three coins sums the channels $K_c\rho K_c^*$. Projecting to the two tip words defines the unnormalized receiver code. At each setting all complement events enter the counts. A detector can identify this subspace by the conserved inventory, one sharp gateway vacancy and the occupation of the two tips. The quadrature measurements on the tips resolve the needed coherent matrix elements after the corresponding readout pulses.

Let $B$ be the code Choi matrix, using a normalized maximally entangled two-rail input. If $B^{T_{\rm out}}$ has a negative eigenvalue, the accepted channel can transmit entanglement with the source. Each of the six prepared density matrices gives one acceptance statistic on $[0,1]$ and three Pauli statistics on $[-1,1]$. For $M$ total appointments and $N$ starts per bank, Hoeffding and the union bound give simultaneous errors $\delta$ and $2\delta$, respectively, with
$$
\delta=\sqrt{\log(2\cdot 24M/0.01)/(2N)}.
$$
The reconstructed two-by-two output density has operator-norm error at most $(1+2\sqrt3)\delta/2$. Expressing the input Choi matrix as $\sum_\mu\sigma_\mu^T\otimes\mathcal{E}(\sigma_\mu)/4$ bounds the norm error of the partial transpose by $(1+2\sqrt3)\delta$. Weyl's eigenvalue inequality gives the reported radius. The calculation uses original launches. Conditional normalization could conceal the unequal acceptance of the input rails. [Six-state readback](detectors.py) · [Long-time independent check](test_live_trimer.py)
