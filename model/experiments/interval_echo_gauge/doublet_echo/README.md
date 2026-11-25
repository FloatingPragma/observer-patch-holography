# A material doublet returns through two interfering apertures

An observer can return any state of the registered [pair doublet](../bright_line/README.md) to its home cell with one calibrated rotation and one success probability for the entire plane. Both flight legs use the free Interval Echo law on all 20 cells. A neighboring cell turns the pair around with guarded local contacts. The neighboring observer preserves two distinct paths through different port pairs at turnaround. Their amplitudes interfere on return. At exposure 32 per leg the terminal counts give their relative phase as **+0.458 radians**, in a simultaneous 99% interval from **+0.360 to +0.563**. A separate midpoint measurement records the aperture and removes that interference. [Protocol](PROTOCOL.md) · [Integer records](records.json)

## What the observers do

Two unit loads inhabit an otherwise empty, blank 240-port screen. The home gateway is port 1. The observer launches a second load at either port 0 or 5 while the other occupies port 1. The six source settings are the two individual words and equal superpositions with relative phases 0, π, π/2 and −π/2. A guarded exchange and a private port phase prepare these single-edge superpositions. The home cell's independently calibrated current service gives a common phase to the plane orthogonal to its bright pair line. Comparing this tag coherently against idle and collecting the difference arm selects that same plane. The current service, coherent switch and its phase readout are supplied. Both arms of this probabilistic preparation count.

At the remote gateway the observer checks for a unit at port 82 and emptiness of every port outside 82, 74, 83, 78 and 79. The remaining unit may occupy either edge (74,83) or (78,79). The detectors leave these four alternatives coherent. A private phase of π on the second port of each edge followed by a full guarded exchange implements the signed swap J on that edge. A further private phase controls the relative phase of the two edges. All four active ports have the required native degree-five reference. The intact board opens again for the return. The same difference-arm filter at home selects the pair plane, and an inverse local source rotation reads the target and its orthogonal logical word. The archive counts both accepted words, all other returned words, failed home collection or filter, missed midpoint, and failed launch filter. The experiment with a recorded aperture also retains that label on its successful midpoint histories.

For each phase in the sharp-load experiment the six inputs contribute a total of **6 × 2⁷⁸ synthetic launches**. These probabilities are very small. The acquisition samples pooled independent integer batches from the finite law; a wall-clock duration or an empirical particle yield cannot be assigned to this number.

| Remote intervention at exposure 32 per leg | Target returns, native contrast 8 | Target returns, contrast 4 |
| --- | ---: | ---: |
| First aperture alone | 1,486 | 65,948 |
| Second aperture alone | 3,479 | 397,354 |
| Both, with the aperture recorded | 4,888 | 462,976 |
| Coherent, relative phase 0 | 8,768 | 785,146 |
| Coherent, relative phase π | 824 | 140,513 |
| Held phase π/3 | 5,107 | 648,286 |
| Held phase 4π/5 | 482 | 219,213 |

The two held probabilities lie inside the intervals forecast using only four independent phase-calibration banks. The inferred phase on the contrast-four flight is −0.084 radians, with interval [−0.096,−0.072]. The reference gates were recalibrated for this control, and the home plane instrument kept its native calibration. A plain swap without the relative sign produced **54,109,609 orthogonal logical outcomes** on the six sharp native preparations. The signed swaps produced zero orthogonal outcomes. The native normalized-plane bank at exposure 8 produced zero signed target returns at every phase, leaving the phase unresolved there.

## Why every input undergoes the same rotation

Let V be the real orthonormal two-column frame of the calibrated plane. Let Wᵢ collect the two sharp configurations of remote edge i. During each flight U = exp(−itH), with H real symmetric in the native pair basis. The one-way map is Tᵢ = WᵢᵀUV. Reciprocity gives the return map VᵀUWᵢ = Tᵢᵀ. Both legs run forward in time. For J = [[0,−1],[1,0]], direct multiplication gives TᵀJT = (det T)J for every complex two-by-two matrix. The full collected map is therefore

**[det(T₁) + exp(iθ) det(T₂)] J.**

This identity holds for every complex input state and for a doublet entangled with an untouched reference. The detection probability is the squared modulus of the scalar. At exposure 32 the two individual path weights are 1.084 × 10⁻¹⁷ and 2.550 × 10⁻¹⁷ for a normalized pair-plane input. Their phases differ by +0.45752 radians. Recording which edge was occupied at the turnaround replaces the square of the amplitude sum by the sum of its two squares. The one-way map to the adjacent bare doublet has strongly unequal transmission probabilities. The present equality concerns a **collected round trip**, with a different remote instrument. [One-way comparison](../doublet_flight/README.md)

For sharp inputs let E embed the two occupied-port words. The difference port of the local home tag has amplitude αP, where P = VVᵀ and α = (exp(iφ)−1)/2 is determined by the independent home service. Set B = VᵀE. Collecting both filters and returning to the two sharp words multiplies the map above by α² det(B), since BᵀJB = (det B)J. The resulting success probability is 0.00007365485 times that on normalized doublet inputs. This accounts for the preparation, analyzer and terminal word overlap without dividing away their failures.

There is also a constraint on any attempt to extract the gauge group from this family of loops. Allow a receiver aperture of any finite dimension and any antisymmetric turnaround K on its coherently collected modes. On the two-dimensional home plane, TᵀKT is antisymmetric and hence a scalar multiple of J. Sums of coherent alternatives obey the same identity. Consequently all successful operations of this class, with the same home frame, commute. After removing their scalar amplitudes their projective products lie in the two-element set generated by J, since J² = −I. The measured route phases affect return counts, but this family cannot distinguish a noncommuting doublet connection or select an order-six global identification.

An integer computation certifies that both det(Tᵢ) are nonzero analytic functions of the exposure on the native board. It uses the real, unnormalized basis obtained by spectral projection of differences of pair stars. In the integer spectral coordinate u = −it/2 both determinants have nonzero coefficient 85 modulo 101 at degree 17. The antisymmetric combination d₁(−u)d₂(u) − d₂(−u)d₁(u), where dᵢ are those determinants, gives coefficient 87 at degree 37. The coefficients use the irreducible native band quartic and its simple root 69 modulo 101; a nonzero residue certifies a nonzero algebraic coefficient at the physical root. Thus both routes contribute and their amplitudes have a non-real relative phase except at discrete exceptional exposures. Changing the real home frame multiplies both determinants by the same nonzero real factor. The trace of each collected operation is fixed by the signed turnaround. A physical gauge connection or a particular central quotient would need independent loop and matter-representation measurements. [Exact certificate](certificate.py)

For confidence intervals, every pooled count k in n launches has binomial mean μ = np. The binomial moment-generating function is bounded by exp[μ(exp(t)−1)] for every real t. Optimizing the exponential tail bound gives exp[−μ+k−k log(k/μ)], with the continuous extension at k = 0. Invert it separately above and below k with tail allocation 0.01/4096. Linear interval propagation from the four quadratures bounds each held probability, including the sign of its sine coefficient. The tests independently check every interval by exact binomial zero-count tails or Poisson tail inversion with a binomial approximation error at most np². For vanishingly rare outcomes the simulator computes the first binomial mass with log1p and inverts successive masses; rounding 1−p to one would suppress arrivals at these exposures. [Acquisition](measure.py) · [Independent returns and statistical checks](test_echo.py)

From the repository root, use an unused destination for the replay.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_gauge.doublet_echo.measure --output /tmp/doublet-echo.json
.venv/bin/python -B -m model.experiments.interval_echo_gauge.doublet_echo.certificate
.venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/doublet_echo
~~~
