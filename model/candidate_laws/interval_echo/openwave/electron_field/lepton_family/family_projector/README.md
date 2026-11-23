# Does the exact rank-three port band travel with a charged state?

The finite icosahedral theorem selects a particular rank-three port band. Its most direct coherent attachment to the fully written charged cell **leaks out in every direction for both signs of charge**. A separate whole-spectrum bound excludes the exact Koide value for any three native level energies on either tested cell sector. [Numerical record and exact certificates](results.json)

## Bind the theorem to the cell

The exact coefficient theorem resolves twelve port amplitudes into spaces of dimensions 1, 3, 3 and 5. The selected three-dimensional space has adjacency eigenvalue +√5. We read the theorem’s integer projector tables, match its incidence to the executable cell’s seams, and verify the projector polynomial identities with integer arithmetic. The two label sets are related by the physical port sequence 0, 1, 5, 7, 11, 9, 10, 8, 4, 6, 2, 3. At uniform degree the scalar port Hamiltonian preserves the +√5 and −√5 triplets exactly at native energies $4(5-\sqrt 5)$ and $4(5+\sqrt 5)$. The numerical residual is below $7\times10^{-15}$. [Lean coefficient theorem](../../../../../../../Lean/Screen/A5FamilyBand.lean) · [Executable port binding](measure.py)

Use the sending cell from the 3–48 matched seam of the 20-cell board. All 30 internal seams are written, and ring and material terms act freely. For every charge address $p$, the specified isometry $W$ maps the port basis vector to the equal positive superposition of all reachable sharp head words at that address. We retain their measured multiplicities when normalizing the columns of $W$. The charge and its distant compensator select a Gauss sector at preparation. All weights below sum over the unresolved head configurations. [Attachment design](DESIGN.md) · [Underlying charged-cell transport fixture](../transfer/)

For a port projector $P$, the three singular values of $(I-WPW^*)HWP$ measure how fast each independent input direction leaves the proposed charged triplet. The positive smallest singular value gives an escape probability proportional to $t^2$ for every initial triplet state. Values below use the native Hamiltonian, including the factor of one half in the archived integer matrix.

| Charge | Sharp configurations | +√5 leakage singular values | −√5 comparison | +√5 weight in three lowest cell eigenstates, respectively |
| --- | ---: | --- | --- | --- |
| Vacancy | 5,790 | 8.5098, 6.7461, 5.9673 | 8.6124, 7.4716, 6.1110 | 10.3432%, 0.7665%, 0.1549% |
| Excess | 1,954 | 7.5213, 6.5299, 5.0309 | 8.4280, 5.7722, 5.0931 | 9.0864%, 0.7949%, 1.4729% |

There is an exact check independent of these floating-point singular values. Start with amplitudes constant across head words at each individual port. Apply the integer Hamiltonian and subtract a reference output at the same port. Nonzero integer minors certify that these contrasts have ranks 11 and 10 for vacancy and excess. Direct multiplication verifies their entire kernels. The vacancy amplitudes must be equal on all twelve ports. The excess amplitudes must be equal on the eleven ports other than port 0. The remaining amplitude can be chosen freely.

Both theorem triplets are odd on every antipodal pair. To convert a constant sharp-word amplitude to a coefficient in the isometry $W$, multiply by the positive square root of the port’s word count. Two nonzero coefficients of the same sign cannot be opposite on a pair of antipodes. Applying this to an antipodal pair away from port 0, then to the pair containing port 0, leaves only the zero vector in either kernel. Thus every proposed triplet state immediately develops head-amplitude differences at a common charge address. Other coherent embeddings, including charge-dependent structured field patterns, require their own physical test. [Integer minors, kernels, singular values and level overlaps](results.json)

## Whole-spectrum Koide check

For three positive native energies the statistic is $Q=(E_1+E_2+E_3)/(\sqrt{E_1}+\sqrt{E_2}+\sqrt{E_3})^2$. The complete charged sectors are finite. For the lower bound, propose a positive integer vector from the numerical ground state, then check componentwise using exact integer multiplication that $H_{\rm doubled}x>2ax$. The matrix has nonpositive off-diagonal entries. Shifting its negative by a sufficiently large multiple of the identity and applying Perron comparison gives $E_{\min}>a$. The largest absolute row sum of the native Hamiltonian bounds $E_{\max}$.

| Sector | Certified energy interval | Upper bound on largest/smallest | Upper bound on Q for every triplet |
| --- | --- | ---: | ---: |
| Vacancy | 4.343 < E ≤ 55 | 12.66406 | 0.4745850 |
| Excess | 5.450 < E ≤ 53 | 9.72478 | 0.4475344 |

To maximize $Q$ in an energy interval $[a,b]$, vary the three square roots in turn. Each coordinate has a single interior minimum when the others are fixed, so a maximum occurs at a corner. The largest corner has one high and two low entries. With $r=\sqrt{b/a}$, this gives $Q\le(r^2+2)/(r+2)^2$. Hitting $2/3$ requires $b/a\ge(4+3\sqrt 2)^2\approx67.94$. The wider vacancy interval also contains the entire excess spectrum, so the exclusion covers every triplet formed from both sets of levels together. It applies to absolute energies above the specified zero-energy neutral vacuum, with no adjustable common offset. [Full bounds and positive integer margins](results.json)

Run `python3 -B -W error measure.py --output /tmp/family-projector.json` in this directory with NumPy, SciPy and SymPy. The run verifies the archived charged-cell code digests and reads the Lean tables directly.
