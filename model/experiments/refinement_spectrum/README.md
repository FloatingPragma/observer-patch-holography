# Rest frequencies on a four-child register

The experimental member `four-child-incidence-CAR-v1` gives an observer three distinct light rest frequencies
at depth two. Their ratios in the observer's dimer clock are **1, 2.13582 and 5.15632**. Increasing the depth
produces additional light branches, and the ratio of two successive light frequencies is below 8/3 at every
depth. With any integer branching factor of at least two and any uniform parent-child hopping strength, the
ratio of two successive chain minima lies above the effective relative hop g and is at most g+1.

## Source and observation

The [construction manual, §§2.3 and 5](../../construction_manual/CONSTRUCTION_MANUAL.md) supplies the
20-, 80- and 320-cell carriers, the oriented seam-to-wedge incidence `C` and the hop magnitude
`a = 1/√2`. Each wedge joins its two seams with opposite signs. The external one-quantum operator
`D` has off-diagonal blocks `a C` and its transpose and zero diagonal. It is the hop of the scalar channel
without the level term `m` of the manual's section 5.3. The fourfold subdivision of the carrier
motivates a separate internal register. Its parent-child links, two sublattices, grading, and controlled
interferometer are design choices. These internal links represent additional coherent hardware.
The [official executable](../../candidate_laws/contrast_memory_spherical_v1/SPEC.md) has its own repair and measurement laws.

Attach two fermionic modes A and B per node of a rooted four-child tree of depth L to every external mode.
Write S for the map sending a parent to the sum of its children and let T = a(I − S).
The internal Hamiltonian h has off-diagonal blocks T and its transpose. The external hops on sublattice B
have the opposite sign to those on A. With Γ = diag(I, −I), the one-quantum generator is

    H = h ⊗ I + Γ ⊗ D,        H² = h² ⊗ I + I ⊗ D².

The number-conserving CAR Hamiltonian is the second quantization of H. All nonzero microscopic hops have
magnitude a. The same-node A–B hop supplies a common rest scale even with the children disconnected.
The anticommutation of h and Γ proves the squared identity on every state in the one-quantum sector.
Thus a pair of internal and external rates μ and p has squared frequency μ²+p². Here p is an external
oscillation rate measured by a reference experiment. Assigning a physical wavelength requires a further observable.

A detector couples its fuel to a declared profile v by a balanced splitter. It parks the fuel for the
source evolution, applies one of four phases, reverses the splitter and reads the fuel occupation.
For return amplitude A = ⟨v|exp(−itH)|v⟩ and phase θ, the fuel-to-fuel probability is
s = |exp(iθ)+A|²/4. On a half-occupied mixed source, the final full-fuel probability with initial fill z
is (1−s)/2 + zs. Subtracting the empty-fuel trial from the full-fuel trial gives s. Opposite phases give
Re A and Im A. All operations preserve number during a trial. For a fixed total parity one can prepare an
equal mixture of even source occupations (at least two modes), with one-body covariance I/2, and pair
the fuel with an oppositely filled spectator. Quadratic evolution and a single occupation read depend
only on that covariance.

An independent half-occupied A–B dimer supplies the clock. Its splitter selects the negative eigenmode
as the reference orbital. Filling that orbital in every trial would change the raw fuel counts.
The four phase records give
an unwrapped angle in units of a. The spectroscopy schedule supplies a positive increment of 0.375 radians
per sample. The reader checks contrast and accepts increments
between zero and π/2. A complete turn between samples would alias. The schedule, phase rotations, stopped
hops during switches and repeated initial preparations are supplied services.

## Records and controls

The root and a difference of two siblings at each remaining layer give one preparation for every distinct
internal chain length. Spectroscopy combines each with a uniform seam profile, a port profile or a four-seam ring.
Their measured external rates are 0, √6 and √10 within the frequency gate. The federation experiment
tests the uniform and ring profiles on all three carriers. For each transferred profile the code checks the full return of C times
the partner to the original vector, including its orthogonal complement. A transfer certifies the same
record law on that preparation. Federation copies provide no additional independent count evidence.

The reader estimates lines and weights from the first 704 clock readings and predicts the remaining 320.
It sees fuel counts and clock counts. A separate check compares every visible line and its weight with the
source projection, including the bulk bands. The maximum held-out RMS error is 0.002161 across all
30 spectra and references. Control runs use the same apparatus. Disconnecting all children gives one
rest frequency. Giving both sublattices the same external sign produces |μ±p|. The maximum measured
residual |ω²−p²−μ²| is 0.00011345 for the declared law and 17.13 for the equal-sign control (nearest rest line).

The transfer experiment prepares one excitation on an otherwise empty source and uses eight rounds of a
balanced evolving and parked path at fixed dimer interval 1.
The passing Kraus operator is (I + exp(−ih/a))/2. Each first failure has its own record. After preparation,
a local seed at the sender evolves on all 3,330 external modes. The receiver groups the final detector
outcomes into sender, neighbor, intercell seam and elsewhere. The preparation columns give the source
probability and its projected light weight. The peak neighbor fraction is an observer count divided by
**all** trials, including failed preparations, over five specified waits.

| Chain length | Lowest frequency in dimer units | Preparation success | Light weight among successes | Peak neighbor fraction |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.999999 | 0.12377 | 1.000000 | 0.006026 |
| 2 | 0.414213 | 0.60412 | 0.99999998 | 0.025970 |
| 3 | 0.193937 | 0.73077 | 0.99996424 | 0.031340 |

The stopped-exterior control has zero receiver probability. The receiver's dimer reads 0, 1.00008,
2.00048, 4.00009 and 7.99999. Transport concerns the simulated detector records. Cell labels identify
the readout groups and supply no physical distance or speed calibration.

## Complete spectrum and the depth limit

The normalized child-sum S/2 propagates a radial unit vector one level. At every parent, the three
orthogonal child differences seed independent copies of the same shift. This decomposes the whole tree
into one radial chain of length L+1 and 3·4^(L−n) chains of each length n = 1,…,L. In each chain
T/a = I−2J, where J is the truncated forward shift. The positive frequencies are the singular values
of T, paired with their negatives in h. At depth two, lengths 1, 2, 3 have multiplicities 12, 3, 1.
The six distinct positive frequencies include the three light lines and the bulk lines 1.90321, 2.41421
and 2.70928 in dimer units. Multiplicity-weighted counts give all 42 internal states.

The squared singular values are the eigenvalues of the tridiagonal matrix with diagonal (1,5,…,5)
and off-diagonal −2. Deleting the first row and column leaves the constant tridiagonal matrix, with
eigenvalues strictly between 1 and 9. Strict interlacing (all off-diagonals are nonzero) puts every
line except the lowest above 1. The norm of I−2J is below 3 at finite length. Leading principal interlacing
strictly decreases the lowest rate with each added layer, starting at 1 for length one.

Let uₙ be the lowest rate in dimer units. The inverse of I−2J has entries 2^(i−j) below and on its diagonal.
Its largest entry is 2^(n−1) and its squared Frobenius norm is (4^(n+1)−3n−4)/9. Bounding its spectral
norm between these gives

    3 / √(4^(n+1)−3n−4) ≤ uₙ ≤ 2^(1−n),        uₙ₋₁/uₙ < 8/3   (n ≥ 2).

Even the three lightest branches at arbitrary depth have largest-to-smallest ratio below 64/9.
Picking widely separated layers would leave intermediate light branches in the same source.

### Allowing arbitrary uniform hopping

Take any branching factor b ≥ 2, depth L ≥ 0, nonzero common dimer hop a and uniform parent-child hop c.
Put g = √b |c/a|. The orthogonal child sums and differences give the same chain matrices I−gJ with
b−1 sibling differences per parent. Measure rates in units of |a|. If uₙ is the least singular value of
the length-n matrix, write Kₙ = (I−gJ)⁻¹. For g > 0,

    g < uₙ₋₁/uₙ = ‖Kₙ‖/‖Kₙ₋₁‖ ≤ g + 1       (n ≥ 2).

Indeed Kₙ = I + g J Kₙ and J Kₙ is an isometric embedding of Kₙ₋₁ with one zero column.
The triangle inequality gives the upper bound since ‖Kₙ₋₁‖ ≥ 1. For the lower bound, delete
the top row and last column of Kₙ. The remaining entries dominate gKₙ₋₁ entrywise. Apply this
submatrix to a positive maximizing right singular vector of Kₙ₋₁ (its Gram matrix is strictly positive).
Restoring the top row adds a strictly positive component. All successive ratios therefore lie in one
interval of width one. Principal interlacing of the squared tridiagonal gives uₙ < uₙ₋₁.

For L ≥ 2 and 0 < g < 2, the second distinct positive frequency is at most uₗ/uₗ₊₁ times the first,
hence the first adjacent ratio is below 3. For g ≥ 2, deleting the first row and column of
(I−gJ)(I−gJ)† leaves a tridiagonal with spectrum strictly above (g−1)² ≥ 1. Interlacing puts
every nonminimal chain frequency above 1. The three smallest distinct frequencies are therefore
uₗ₊₁, uₗ and uₗ₋₁. Their two adjacent ratios differ by less than one. At L = 1 the three lines are
u₂, 1 and 1/u₂ because the length-two matrix has determinant one. At L = 0 or g = 0 there is one rate.
These cases exhaust the uniform family.

With the uniform relative hop and the branching free, the three lowest distinct lines of the family
therefore have a first adjacent ratio below 3 (0 < g < 2 and L ≥ 2) or two adjacent ratios that differ by
less than one (g ≥ 2, or L = 1). Selection rules that hide intervening states, nonuniform couplings and
interactions fall outside the proof.

### Fixed root in the infinite tree

At the fixed root of the infinite chain, T†T/a² = 5I−2(J+J†) has the sine-transform spectrum [1,9].
TT†/a² also has a zero mode with normalized components (√3/2)·2^(−j). Since T†e₀/a = e₀, the positive
squared-rate spectral measure at the root equals the sine-transform measure divided by its spectral variable.
The normalized root return in dimer time t is

    A(t) = 3/4 + (2/π) ∫₀^π sin²θ cos(t √(5−4cosθ)) / (5−4cosθ) dθ.

The fixed detector sees a stationary weight 3/4 and a band from 1 to 3. It has no isolated positive subgap
line in this limit. The alternating hopping path first reaches the omitted link and returns in 4n hops.
The power series for the unitary gives a uniform amplitude error at most
min(2, 2 Σₖ≥₄ₙ (3|t|)^k/k!). Records at lengths 1, 3, 5, 9, 17 and 33 use the same root detector.
At length 33 the bound is 3.37 × 10⁻⁴² on [0,8]. The largest recorded difference for dimer clock
readings in that window is 0.00177, including clock and counting noise. Over all 257 queries it reaches
0.04895. The last query is scheduled at time 64 and the dimer reads 64.0012. The theoretical
length-33 rate is 1.74623 × 10⁻¹⁰. Those finite clock records do not resolve it.

These statements concern the register coupled to the scalar incidence. The squared dispersion identity is
read at the three profiles of the menu, and the records carry no calibrated propagation speed.

## Resources and replay

For N tree nodes, E seams and W wedges the source has 2N(E+W) logical modes and
(2N−1)(E+W)+4NW undirected hops. The depth-two cell uses 6,300 modes and 16,230 hops. The 20-cell
flight uses 139,860 modes and 363,330 hops. The largest profile transfer checks 2,237,760 logical modes
by invariant-subspace reduction. The length-33 root calculation propagates 66 radial modes. Preparing the
literal depth-32 full register would require (4³³−1)/3 nodes per external mode.

The spectroscopy budget is 2³⁷ simulated binary trials aggregated by setting (2²⁰ shots per setting).
Flight uses 30 multinomial settings, 33 four-phase reference queries and eight possible first-failure
records per flight. A balanced parked path needs a companion for each active mode it carries. An arbitrary
profile on M connected orbitals can be brought to the root by at most M−1 adjacent two-mode rotations
with controlled phases on a spanning tree. Conjugating a root-fuel splitter and reversing the sequence
costs at most 4(M−1)+2 such gates for the two splitters. The calculation assumes these ideal controls.
It supplies neither their native pulse errors nor thermodynamic work and archival costs.

Binary acquisition digitizes each probability to the nearest multiple of 2⁻³⁰. For a multinomial partition
the cumulative edges use the same grid, giving per-bin error at most 2⁻³⁰. Each binomial draw uses an
independent 52-bit uniform midpoint and the inverse binomial distribution. At each integer threshold the
resulting CDF differs from the exact digitized binomial CDF by at most 2⁻⁵³ in exact arithmetic.
A fuel quadrature combines
four sample proportions. For K simultaneous binary proportions with N shots each, Hoeffding's inequality
and a union bound give a quadrature error bound 4√(log(2K/δ)/(2N)) at confidence 1−δ, plus 2⁻²⁹
digitization error. Discretizing the inverse draws adds at most 2K·2⁻⁵³ to the tail probability.
Frequency and dispersion gates are checks on the retained sample and its reserved tail.
The spectral bounds and the bounds on successive ratios follow from the source.

From the repository root, with the [model dependencies](../../requirements.txt) and one BLAS thread:

```sh
python -B -O -m model.experiments.refinement_spectrum
python -B -O -m model.experiments.refinement_spectrum.limit
python -B -m pytest -q model/experiments/refinement_spectrum
```

Each replay verifies source hashes and regenerates the complete count receipts. The first also scores
the held-out spectra, control laws and all transfer outcomes. The second compares fixed-root records to the
infinite measure. The tests propagate all 139,860 modes on the labeled federation, compare every preparation
failure and detector region, read the fuel on Fock space and check the smallest chain root against a
high-precision characteristic polynomial. The measurements establish the three rest-frequency ratios,
the squared dispersion at the profiles of the menu and finite transfer to a neighboring cell.
