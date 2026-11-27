# Can excesses and vacancies form an energetically bound nucleus on blank material?

On any finite connected blank Interval Echo board containing both signs, the lowest mixed-stock energy exceeds the cost of putting all excesses on one copy and all vacancies on a second copy. The result holds for native contrast strength eight, and more generally for strength at least one half in the registered exchange units. The comparison uses identical contacts, stock counts and host-subtracted references on both boards. It rules out a mixed nucleus stabilized below all breakup thresholds by these contacts alone. It leaves recorded contacts, additional species, resonances and symmetry-restricted decay times available for separate experiments. OW-17's physical binding curve, nuclear levels and halo lifetimes remain open. [Frozen comparison](PROTOCOL.md) · [Contact proof](SignedFragments.lean) · [Preceding seven-defect measurement](../../row_16_binding/signed_pair/README.md)

## What the two cells measure

On the 12-port icosahedral cell, each observer measures one of 30 ordinary contact energies or nine private reference energies on a fresh stock. The contrast and exchange values on each blank pair come from the registered 99-state repair law. The 256 coherent delay settings of the local energy meter resolve every occupied quarter-unit eigenvalue at the enrolled strengths. The source makes a signed pair on ports 0 and 1, then either reads the sharp stock or puts these two loads in an even or odd exchange superposition. The corresponding sign on a second pair of cells is prepared as the appropriate independently sampled marginal stock. The preparation, shutters, coherent energy meter and inverse parity analyzer are addressed apparatus. Energy expectations refer to calibrated contact reads on independent cohorts. [Local acquisition and meter](measure.py) · [Complete terminal counts](counts.json)

For a single adjacent even pair at native strength, the joined stock cost **15.0073 ± 0.7575** exchange units more than its matched two-cell comparison. The independent contact census and parity measurement forecast **15.0000 ± 0.0022**. Four excesses with three vacancies, with the same even-pair preparation, gave **76.9717 ± 0.7575**, against **77.0000 ± 0.0022**. At strength one quarter the single even pair gave **−0.4985 ± 0.0421**, with forecast **−0.5000 ± 0.0022**. The exchange parity and contrast strength affect this matched-state energy difference. A negative difference against these two prepared fragment states leaves the comparison with the minimum fragment energies undecided at that coupling. All 36 inventory, parity and strength combinations agreed within the simultaneous 99% counting intervals. Each term reading in each bank retained 1,048,576 original simulated starts.

## Why separating the signs can lower the energy

For every sharp input word, place its excesses on board A and its vacancies on board B, filling the other ports with ones. This map J sends distinct orthonormal words to distinct two-board words. It therefore preserves norm and all relative phases. The image excludes simultaneous occupations at matching addresses. The two boards have independent generators with the common energy origin of a filled host. Relative private energy at port i is κ(nᵢ²−1), and it adds exactly under this embedding.

The blank two-port generator has exchange (I−swap)/2 and diagonal contrast κ(nᵢ−nⱼ)². Subtract the compression of the two independent generators from the joined generator. The difference vanishes on every input except the two unlike words 20 and 02. Its block on these two words is

$$
D=(2\kappa-\tfrac12)I-\tfrac12\,\mathrm{swap}.
$$

Indeed in the separate-board comparison, moving just one of the two signs between those occupied endpoints sends the product into a collision configuration outside the image of J. This removes the off-diagonal entry on compression. The joined contrast costs 4κ and the separated contrast costs 2κ. All 81 local matrix entries and the three private entries are checked exactly. The two eigenvalues of D are 2κ−1 and 2κ. Equivalently, on amplitudes x,y of two words related by the unlike swap, the quadratic form is (2κ−1)(|x|²+|y|²)+|x−y|²/2. Summing these terms, one pair per unchanged exterior word and per edge, gives

$$
H_{\rm joined}\ \geq\ J^\dagger(H_+\otimes I+I\otimes H_-)J
$$

for κ≥1/2. Taking the infimum on the left bounds it below by the unrestricted two-fragment floor. A graph rotation preserving the references commutes with J provided both copies rotate together. The same non-strict variational comparison holds within each common rotation sector. The finite-domain lower bound compares energies. A transition between the two preparations would need a communicating path; the comparison alone gives no breakup time.

## Why the unrestricted threshold is strict on a finite connected board

Each nonempty single-sign sector has a connected configuration graph: edge transpositions on a connected board generate all permutations of its sites. After a scalar shift, minus the fragment Hamiltonian is irreducible and entrywise nonnegative with positive diagonal. Perron–Frobenius gives a unique lowest eigenvector with a strictly positive amplitude on every allowed configuration. In the separated product ground g, the probability c that two nonempty fragments occupy any common address is strictly positive. In every normalized encoded joined input Jψ, that probability is exactly zero. Cauchy–Schwarz on the collision-free projector therefore gives |⟨g,Jψ⟩|²≤1−c. The finite spectral gap Δ of the independent pair of cells is the smaller of the two single-fragment gaps. Consequently the joined ground is at least **E₊+E₋+Δc**. This argument establishes strictness separately for each finite connected board. It supplies no size-independent gap.

In a separate experiment each observer sharply reads all twelve ports of an independently prepared fragment ground stock. Comparing the two words tests for collision at one or more matched addresses. On the native (1,1) fragments, 6,839 of 1,048,576 simulated pairs shared an address; the simultaneous 99% lower limit on c was 0.006287. The calculated smaller fragment gap was 0.067721 exchange units, giving a lower penalty estimate of 0.0004257 at that counting limit. For the (4,3) fragments the collision count was 158,773, the lower limit was 0.15036 and the calculated gap was 0.000002700. All eight comparisons on native and half-contrast material resolved positive collision probabilities. The gap values are numerical diagonalizations. These collision rates concern independently prepared model ground stocks; the general strict inequality follows from the finite operator argument. [Frozen cold-stock test](COLLISION_PROTOCOL.md) · [Collision counts](collisions.json) · [Decoder and fragment spectra](collisions.py)

From the repository root, run with the scientific Python environment:

~~~sh
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_17_nuclei.sign_fragments
python -Bm "$p.measure"
python -Bm "$p.collisions"
python -Bm pytest -q --pyargs "$p.test_fragments"
~~~

Both acquisition modules accept `--acquire PATH` to reproduce their complete count archives. In the Lean project run `lake env lean -j1 -DwarningAsError=true ../model/candidate_laws/interval_echo/openwave/row_17_nuclei/sign_fragments/SignedFragments.lean`. The Lean file verifies the local contact table, norm-preserving word recovery, private reference identity, real quadratic form, sum positivity and scalar collision-penalty implication. The configuration-graph, variational and Perron–Frobenius arguments are given above.
