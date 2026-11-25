# A formed two-head source driving a reciprocal bridge

Ordinary rail energy and a change of head correlations both move the bridge qubit on the 20- and 80-cell federations. The two changes point in different directions in its measured three-component state. Each of the 32 size, time and charge comparisons excludes an ordinary-calibrated scalar transfer to the entire qubit. This distinguishes a mechanism that senses correlations from one that couples universally to modular charge.

## Preparation and source measurements

Use the active Interval Echo native single-mobile sector during formation. The region consists of a marked seam in each of two adjacent cells and their whole blank connecting bridge. A vacancy or excess starts on the outer left port and evolves for four native ticks. Tracing out that mobile leaves the two native even heads in a state ρ. Choosing the two heads from independent repetitions gives σ = ρL ⊗ ρR and preserves both one-head marginal states. The head pair has four basis words. The receipt gives the four outcome probabilities for each of its nine local Pauli settings. These are finite-model population probabilities. All charges and arrival sites are retained in their respective preparations.

A fresh mobile vacancy or excess loads the symmetric and antisymmetric rail eigenstates on the two ends of this same bridge. In the native energy units the bridge action is

    b = (17 I − Xrail)/2  on the two ports, tensored with both head identities,
    b = 0                  elsewhere.

Keep the rail weights independent of the head state. The reference uses weights (2/3, 1/3) and the paired heads. The ordinary arm uses (1/3, 2/3) and the same heads. The decorrelated arm uses (2/3, 1/3) and independent heads. The reference defines both source contrasts. The initial strain is |+X⟩ in all arms.

The ordinary excitation increases the initial expectation of both b and the bare native action by exactly 1/3. The only inter-rail matrix element between these two ports is the blank hop −I/2. The two rail mixtures share their diagonal. All other one-head action terms cancel in the rail contrast. Conversely, decorrelation changes neither energy. Each term of the initial native action meets at most one head, whose marginal is identical in the two preparations. Bare energy values use the integer-centered native generator divided by two. Its constant offset cancels from all contrasts.

On the incoming eight-dimensional support the reference modular generator is K = −log(τ ⊗ ρ), where τ has the 2:1 rail weights. Its two contrasts in natural units are

    ordinary       ΔK = log(2)/3
    decorrelation  ΔK = Tr[(σ−ρ)(−log ρ)] = D(ρ‖σ) + D(σ‖ρ).

The last identity uses Tr[(σ−ρ) log σ] = 0, since log σ is a sum of one-head operators. All four formed head states have positive least eigenvalue. For 20 cells the decorrelation charges are about 0.001924 and 0.004173 nats for vacancy and excess formation. For 80 cells they are about 0.200825 and 0.293848 nats. The tomography settings and rail counts determine these charges independently of strain measurements. The factor log(2) between ordinary energy and its modular charge is fixed by the chosen rail weights.

## Bridge response

The extension adds one bridge qubit and makes the **whole** connecting action b reciprocal.

    H = I ⊗ H_IE − Ω X ⊗ I + (exp(−ε Z)−I) ⊗ b,
    Ω = 1/2,     exp(ε) = 5/4.

The initial expectation of the full coupled action includes −Ω + (cosh ε−1)⟨b⟩. The ordinary change is 41/120 and the decorrelation change is zero. Forming and pump charges independently take each sign. A declared pump time of 1/2, 1, 2 or 4 ticks precedes tracing out the entire mobile and both heads. Each fresh reading then receives the resulting two-dimensional state. Its density matrix is (I+xX+yY+zZ)/2, with y = −2 Im ρ₀₁. The Python interface exposes the complete input isometry. Cross terms between different source coordinates are retained until the partial trace.

For a short-time check, the Heisenberg equations on the initial product state give Y′(0) = −2 sinh(ε)⟨b⟩ and Z′ = −2ΩY. The ordinary contrast has ΔY′(0) = −3/20 and ΔZ″(0) = 3/20. Both quantities vanish for decorrelation despite its positive ΔK. At ε = 0 the two factors evolve independently and |+X⟩ is stationary as a density operator. The uncoupled simulation also checks the projection onto |−X⟩ at every enrolled time and pump charge.

At 80 cells, vacancy formation followed by a vacancy pump for four ticks gives

| Contrast | ΔX | ΔY | ΔZ |
| --- | ---: | ---: | ---: |
| Ordinary | −0.12224 | 0.094731 | −0.009987 |
| Decorrelation | 0.005815 | −0.001222 | −0.002107 |

Multiplying the ordinary line by the independent modular charge ratio 0.86919 fails each component. More generally the 2×2 minors of the two response columns certify noncollinearity even if their scalar ratio is fitted after measurement. Across all enrolled cases, the distance of the collar response from the ordinary line exceeds 1.99 × 10⁻⁶ in Bloch-vector norm. The largest component of the modular-calibrated residual exceeds 1.01 × 10⁻⁴ in every case. The result concerns the entire local qubit at each fixed preparation. A selected scalar reader can discard components. Local geometric null balance would require independent identification of proper null focusing, stress and their common coupling on macroscopic sources. The extension and its finite population tests supply no such identification.

## Enclosures and replay

The exact source-to-qubit map uses an eight-column launch isometry V. For each qubit index a define Uₐ = (⟨a| ⊗ I)exp(−itH)V. Then the output entry (a,b) is Tr[Uᵦ† Uₐ(τarm ⊗ ρarm)]. The nominal receipt stores all three Bloch moments, both real and imaginary density matrices, head tomography, and the independent source readings. The separate certificate bounds every stored component by replaying the native evolution with integer arithmetic. Half of the Bloch component error bounds the corresponding real or imaginary off-diagonal part and the normalized diagonal populations.

For formation, a two-limb integer Chebyshev replay bounds the distance of the normalized rational head center from the native reduced state. The reused log certificate encloses the eigenbasis and matrix logarithm. During the pump, 40H has integer entries. The Taylor replay takes 32 steps per tick, keeps degree 20, and rounds real and imaginary parts to the nearest multiple of 2⁻⁵¹. The maximum absolute row sum bounds the spectral norm since H is real symmetric. Each substep has L = ‖H‖/32 ≤ 1. A rounded Taylor term costs at most √(2d)/(2·2⁵¹) in column norm. Its accumulated contributions to higher terms sum to at most exp(L) times this error. The degree-20 Hermitian exponential remainder costs at most L²¹/21! times the input norm. Unitary evolution preserves the accumulated error, so the running bound below 10⁻³ also keeps every numerical column norm below 2. Overflow checks precede each sparse product and integer addition. The final Pauli effects are integer Gram matrices, evaluated in short limbs. For each rail the isometry error is bounded by twice the largest of its four column errors. Polarization and tensor-product trace-norm bounds propagate both pump and formation errors to each response.

These commands replay the source from the repository root.

```sh
PYTHONPATH=.:code python -m matched_einstein_limit.rate_geometry.seam_strain.source.measure --write
PYTHONPATH=.:code python -m matched_einstein_limit.rate_geometry.seam_strain.source.certificate --write
python -m pytest -q -W error code/matched_einstein_limit/rate_geometry/seam_strain/source
```
