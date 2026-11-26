# Which freely evolving ports can carry a rest line?

## Source and single-copy reader

The matched strand has one shared load on each physical dimer. Its two orientations give the eigenvalues of $Z_j$. Exchange contributes $-hX_j$, with $h=1/2$. A coherent connector has binary-endpoint cost $4J(1-n_a n_b)$ with $J=8$. Interior linear fields cancel. With virtual boundary orientations $Z_0=+1$ and $Z_{n+1}=-1$, the physical generator up to a scalar is

$$
 H=J\sum_{j=0}^n Z_j Z_{j+1}-h\sum_{j=1}^n X_j .
$$

The sign $w_j=Z_jZ_{j+1}$ equals $+1$ at an unsatisfied bond. All bond signs multiply to $-1$, so the number of unsatisfied bonds has parity $n\bmod 2$. For odd $n$ the source can contain one bare wall. A single read of $Z_1=w_0$ preserves all coherences within the selected eigenspace. At the terminal, $Z_7=\prod_{j=0}^6 w_j$. The simultaneous physical readings of $Z_1,Z_7$ retain both signs throughout the two-time record. The initial shape, route and commissioning are supplied. [Seam columns, end fields and dual generator](../../row_30_collective/experiment.py) · [Preparation rotations](../../row_30_collective/packets.py)

Each $X_j$ flips its two adjacent wall signs. Open-chain Jordan–Wigner ordering gives $w_j=iA_jB_j$ and a free Majorana generator with cross block $C=2(JI-hS)$, where $S$ is the one-step lower shift on the $n+1$ bonds. Its positive elementary energies are the singular values of $C$. The restriction to physical parity selects the orientation space. The physical state and fixed terminal effects determine the coefficients in a Bohr expansion of each word probability. Conditioning on a past record changes those coefficients. Its probability is fixed before the terminal delay is varied. Separate time references measure that delay. The small-chain verifier evolves complex one-wall preparations, both trigger outcomes and all four final words with the unreduced physical orientation matrix. [Independent original-basis tests](test_clock.py)

## A finite-resolution gap independent of length

Since $\|S\|\le1$, the triangle and reverse triangle inequalities give $15\|v\|\le\|Cv\|\le17\|v\|$ for every vector. Every constituent energy belongs to $[15,17]$. Within the exact one-quasiparticle space all energies differ by at most two. For a state supported there, any time-independent terminal effect sees this restriction. The experimentally prepared bare orbital can also populate higher odd sectors.

For an arbitrary physical state, any fixed function of $k$ simultaneously read bond bits is a linear combination of their parity products by Boolean Fourier expansion. Each product has even Majorana degree at most $2k$. In normal modes the time dependence of each Majorana factor is a sum of terms with time exponents $\pm i\epsilon_j t$. Each pair therefore contributes a frequency within two of $-32,0$ or $32$. Contractions can remove factors in pairs. Any frequency $\omega$ in the word probability satisfies

$$
 \operatorname{dist}(\omega,32\,\mathbb Z)\le2k .
$$

For $k\le7$ this excludes the entire open interval $(14,18)$. Marginalizing, mixing independent preparations, retaining a trigger sign and changing the initial state preserve this containment. It applies to the fixed physical terminal effects and free generator, including an initially entangled frozen spectator read by a time-independent effect. Tracking a different address or applying pulses as a function of delay changes the observable. [Lean arithmetic bounds](WallBand.lean)

## Exact rational frequencies with arbitrary terminal effects

A separate algebraic check motivated by the measured recoil compares all terminal observables with the exact bulk rest value 15. In site order $A_0,B_0,A_1,B_1,\ldots$, diagonal phases turn the Hermitian single-Majorana frequency matrix into a real symmetric tridiagonal matrix with zero diagonal and successive off-diagonal entries $16,1,16,1,\ldots,16$. It has $2m$ sites, $m=n+1$. Every eigenvalue is simple by the three-term recurrence and the nonzero couplings. At $m=14,30,50$ reflection gives two blocks of size $m$. The positive-reflection block has final diagonal $+1$ from the central bond, all other diagonals zero. Alternating site signs anticommute with the full frequency matrix and exchange the two reflection subspaces. The other block has the negated spectrum of the first.

Write $\lambda_1,\ldots,\lambda_m$ for the first block's roots. They sum to one, one root from each $\pm\epsilon_j$ pair. Its monic integer characteristic polynomial is irreducible over the rationals. For $m=14,30,50$ this follows by reduction at $p=3,11,3$, respectively. The continuant constructs the polynomial from the physical couplings. For each prime divisor $q$ of $m$, the code checks that its modular reduction $F$ has $\gcd(F,x^{p^{m/q}}-x)=1$. It also checks $x^{p^m}=x\bmod F$, the finite-field irreducibility criterion. [Exact certificates](rational_lines.py)

A many-body energy difference, including physical parity selection, is $\omega=\sum_i c_i\lambda_i$ with $c_i\in\{-1,0,1\}$ because each mode occupation changes by at most one. The Galois group of the splitting field acts transitively on the roots. Its average sends each $\lambda_i$ to $1/m$. If $\omega$ is rational this average fixes it, so $\omega=\sum_i c_i/m$ and $|\omega|\le1$. In particular, an exact line at 15 is absent for every fixed terminal effect and initial state on these three open strands. This arithmetic exclusion concerns exact equality. Frequency estimates from finite counts have finite resolution.

## A local reader at the boundary of the seven-bit theorem

The eighth physical dimer reads $Z_8=\prod_{j=0}^7 w_j$. It comprises the first 16 site Majoranas. Every eigenvector of the irreducible tridiagonal frequency matrix has a nonzero first component. Its $j$-th component is that first component times a polynomial of degree $j-1$ in the eigenvalue with nonzero leading coefficient. The determinant of the first 16 components of any 16 distinct eigenvectors is a nonzero Vandermonde determinant times nonzero factors. Site gauge phases and mode normalization preserve it.

Select the eight highest and eight lowest positive energies of the 49-dimer strand. Change the low group from occupied to empty and the high group from empty to occupied, leaving one identical occupied spectator in both states to enforce odd parity. A matrix element of $Z_8$ between these two states must use all 16 mode operations exactly once. Anticommutation makes its coefficient the determinant just identified, selecting the signed eigenvectors for the eight creations and eight annihilations. It is nonzero. A coherent superposition of the two states can therefore give this line at one physical dimer orientation. Preparing that superposition and resolving its amplitude would require a separate experiment.

On this strand $CC^T$ has first diagonal 256, remaining diagonals 257 and off-diagonal $-16$. Integer Sturm counts at squared rational thresholds bracket all 16 contributing positive singular frequencies. Subtracting the eight lower from the eight upper frequencies places their gap strictly between **15.23386 and 15.23402**. Thus the finite-width exclusion for seven wall bits fails at eight. This class-boundary calculation uses the native generator after acquiring the frozen records. [Exact Sturm brackets](eighth_port.py)
