# A complete census for one written seam

Consider one mobile amount-one load, one arbitrary eleven-state core memory and blank memory at every other seam. Write $a,b$ for the core endpoints, $P$ for the projection onto those two material positions, $L_0$ for the graph Laplacian with that seam cut and $A=L_0/2-\kappa P$. The blank references bring every port to total contrast degree six. Subtract the independent uniform blank-load energy $6\kappa$ throughout. The energy of an isolated blank or zero-direction core is zero. This gives a common comparison across all memory states.

## Native block reduction

At all other seams the restriction to blank memory gives half the swap Laplacian plus the contrast of the local load word. Together with the degree references their sum, after the uniform subtraction, is $I\otimes A$. At the marked seam let $W$ be the integer pattern of the write, so $W/\sqrt 2$ is the actual isometry. Let $R$ reverse the direction of the memory and $D$ read its signed value. The local generator equals the native echo conjugate of the memory generator. Eight times its written-to-written block is

$$
(2I-WW^T)\{(I-R+2\kappa D^2)\otimes I_9\}(2I-WW^T).
$$

The blank subspace reduces the generator. In each coin state decompose each pair of signed directions $d,-d$ by reversal parity $s=+,-$. Set $\mu_{d,s}=\kappa d^2+\mathbf 1_{s=-}$. The complete Hamiltonian in this reference splits as follows.

| Memory sector | Hamiltonian | Number of copies |
| --- | --- | --- |
| Blank | Full graph Laplacian divided by two | 1 |
| Direction zero | $A$ | 2 (one per coin) |
| Absolute direction two | $A+\mu_{2,s}I$ | 2 for each parity |
| Absolute direction one | $I_2\otimes A+\mu_{1,s}(I-\lvert\chi\rangle\langle\chi\rvert)$ | 1 two-coin block for each parity |

Here $\chi=(\lvert0,b\rangle+\lvert1,a\rangle)/\sqrt2$. For an exterior material position the core sees word $00$, and a unit-direction memory evolves at $\mu_{1,s}$. At positions $a,b$ the local words are $10,01$. The integer multiplication on all six coin/word columns in this parity sector gives $8\mu_{1,s}I-4\mu_{1,s}vv^T$, where $v$ selects coin 0 with word 01 and coin 1 with word 10. All columns outside the embedding vanish. The zero-direction action vanishes, and each double-direction embedding gives $8\mu_{2,s}I$. These identities are checked separately for the native $\kappa=8$ and changed $\kappa=0$. The dimensions add to $11N$ for $N$ ports. An independent raw 132-state cell assembly tests the complete eigenvalue multiset and finite-time material distributions for both laws. [Exact integer action](core.py) · [Independent cell assembly](test_occupied.py)

## Negative index at the common separated reference

Take $\kappa=8$. The positive cut Laplacian gives $A\ge-8I$. The compression of $A$ to the two endpoints is diagonal with entries $(\deg p-1)/2-8\le-5.5$ because the endpoints no longer share a cut-graph edge and the maximum degree is six. Its compression to $I-P$ is positive semidefinite. The min–max principle therefore gives exactly two strictly negative eigenvalues of $A$, both in $[-8,-5.5]$. The blank sector is positive semidefinite. The double-direction sectors add at least 32 to $A$.

In each unit-direction block, $\mu_{1,s}=8$ or 9. On the subspace orthogonal to $\chi$ the compression is bounded below by $\mu_{1,s}-8\ge0$. The Rayleigh quotient of $\chi$ is at most $-5.5$. Adding the positive memory term to $I_2\otimes A$ leaves the lower bound $-8$. Each parity therefore contributes exactly one negative eigenvalue in $[-8,-5.5]$. Combined with both zero-direction coins, this makes **six negative eigenvalues counted with multiplicity**. No coherent combination can produce a spectral line in $(-5.5,0)$. This holds on each registered federation, including the 320-cell board. The bound counts all stationary lines, including lines dark to a particular source. An expectation value of a mixture may lie in a spectral gap.

## Binding relative to a conserved occupied sector

A physical occupied write and its endpoint phase reversal both maintain their signed-direction parity sector. A reference copy of the same written memory with the load parked away from the core contributes $\mu_{1,s}$ in addition to the independently measured uniform-load energy. Against this higher threshold the observed generator is $H_s=I_2\otimes A-\mu_{1,s}\lvert\chi\rangle\langle\chi\rvert$. Its restriction off the endpoints is positive semidefinite. The four-dimensional coin/endpoint compression has largest eigenvalue at most $-5.5$; its additional rank-one term is negative semidefinite. Hence there are exactly four negative eigenvalues per parity, all in $[-8-\mu_{1,s},-5.5]$. Even using the widest possible ground binding, the excited $n=2$ energy lies at or above $-17/4=-4.25$ and below zero. The $n=3$ energy lies still closer to zero. Both targets sit in the gap in their own conserved sector. The complete common-threshold result excludes them independently for every sharp line of the full single-core family.

Both statements use the operational separated references specified in the protocol. The theorem applies to the isolated unit load and arbitrary state of one stationary seam. Additional mobile species and multiple entangled occupied seams require their own Hamiltonians and separated references.
