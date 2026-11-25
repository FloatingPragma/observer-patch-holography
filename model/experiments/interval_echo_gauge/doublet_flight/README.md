# Can the pair doublet cross an owner boundary?

Yes, on the all-open 20-cell Interval Echo screen. Give the [registered pair doublet](../bright_line/README.md) a phase at launch, let both loads travel freely, and read the seams of the adjacent owner. With the receiver untouched, 29 of its 30 seam frequencies respond at exposure 32. Both independent doublet modes have a nonzero transmission amplitude into the receiving doublet for all but isolated exposures. Their efficiencies differ substantially. If the two cells communicate through their single bridge in isolation, symmetry blocks the doublet exactly, even though other pair preparations cross. [Frozen protocol](PROTOCOL.md) · [All simulated detector banks](records.json)

## What the receiver measures

Use blank memory and two unit loads on otherwise empty ports. On the first inter-owner edge of the 20-cell carrier, the gateway is port 1 of owner 0 and port 10 of owner 6. An isolated cell's first retained odd pair triplet contains the bright projection of the pair star at its gateway. The doublet is the orthogonal plane *in this triplet*. The cross-axis preparation is the ordinary registered odd pair star at source port 0. The entire preparation goes into the experiment, including its components outside the triplet. The on-axis control uses the ordinary star at port 1.

For the tag, apply the two calibrated current pulses at the gateway and its antipode before free release. Their generators commute, and their sum has eigenvalues β on the bright line, α on the plane and zero on the remaining modes. Each pulse lasts 2π/β = 8.727524 in the service units. This gives the bright line and every spectator phase one; only the two plane amplitudes acquire the phase exp(−2πiα/β). Use a sham or the same tag at the receiving endpoint after free flight. The receiver then reads its thirty internal seams and one other outcome. During endpoint service the apparatus uses the isolated-cell current action. The flight uses the ordinary blank-sector Hamiltonian. The exposure is a scheduled appointment in exchange units.

| Active law, exposure 32 | Receiver tag | All receiving-seam clicks, launch sham | All receiving-seam clicks, launch tag | Seams resolving the difference |
| --- | --- | ---: | ---: | ---: |
| Two cells, one open bridge | None | 62,427,028 | 62,431,745 | 0/30 |
| Full federation | None | 426,559,088 | 400,799,285 | 29/30 |
| Full federation | Applied | 426,312,368 | 400,640,818 | 30/30 |

With the all-open receiver untouched, the source tag lowers its total seam-click rate by 6.04% relative to sham, or 0.00002343 per launch. Each setting draws one multinomial bank of 2⁴⁰ simulated launches from the calculated distribution. All missing pairs are counted. The per-seam two-sided Clopper–Pearson intervals have joint coverage of at least 99% across 2,880 prespecified frequencies. The on-axis and zero-exposure controls give no resolved differences. At exposure 8 the untreated all-open receiver already resolves 12 seams. With contrast set to zero and both meters fixed to the original calibration, the untreated all-open receiver resolves all 30 at both nonzero appointments. Geometry and exchange alone can therefore convey this information. The counts measure sensitivity to the launch phase on an ordinary pair source. Resolving the weakest singular direction separately would require preparing and discriminating the corresponding coherent modes.

## Why the protected bridge blocks the doublet

On the two-cell one-bridge fixture, a fivefold turn of the source cell about its gateway fixes the bridge and the entire receiving cell. It commutes with the free Hamiltonian because it permutes the exchange terms and preserves both adjacency and degree. Let R be this turn and A = (I + R + ⋯ + R⁴)/5. In the first real triplet R fixes the gateway bright line. It acts nontrivially on that triplet, as can be checked from the orbit of a projected neighboring star. A real order-five action in three dimensions with a nontrivial rotation has precisely one invariant axis. Thus A kills the whole doublet. Any state with both loads outside the source cell is fixed by R. For the projector Q onto these states and the free propagator U, QU = QAU = QUA. Hence QU vanishes on the source doublet at *every* exposure. The corresponding independent receiver turn also forbids populating the receiver doublet from a pair wholly outside it. This argument continues to hold for a source with only its axial bridge open regardless of the geometry beyond that bridge, provided the turn remains a symmetry. A single load may leave while the other still carries the source's rotation character.

On the full screen the other gateways of each owner break its independent fivefold symmetry. The following squared singular values are the maximum and minimum transmitted probabilities for a normalized input in the source doublet and an ideal receiving-doublet analysis. They come from the full transfer amplitude with no renormalization.

| Screen and flight | Larger probability | Smaller probability |
| --- | ---: | ---: |
| Guarded, any exposure | 0 | 0 |
| Full screen, active contrast, 8 | 1.79083 × 10⁻⁹ | 1.82031 × 10⁻¹³ |
| Full screen, active contrast, 32 | 5.70793 × 10⁻⁷ | 4.05908 × 10⁻¹² |
| Full screen, exchange only, 32 | 1.84776 × 10⁻³ | 4.45245 × 10⁻⁴ |

## Exact rank and unequal transmission

A receiver could try to undo this unequal transmission by filtering the two collected modes and keeping only successful trials. There is a sharp limit for recovery of an arbitrary input, including its phase relative to a reference. Write T for the measured-mode transfer and s for its smaller singular value. For the weakest input the total collection probability is s², so every subsequent physical filter succeeds with probability at most s² on that input. When T has rank two, the ideal filter s T⁻¹ is a contraction and its product with T is s times the identity. Thus the largest uniform success probability for exact recovery from these two modes is s². At exposure 32 the calculated ceiling is **4.06 × 10⁻¹² per launch**, about one success per 246 billion trials even with an ideal receiver. This is a bound on a proposed decoding instrument. Collecting additional output modes changes the transfer being bounded.

An integer calculation certifies the two effects independently of numerical time evolution. In an isolated cell put K = 2H and use the [quartic of the native odd stars](../bright_line/certificate.py),

$$F(t)=t^4+34t^3-128t^2-6904t+41936.$$

Its lowest root λ selects the first triplet on those sources. Apply Q(K;λ) = [F(K)−F(λ)]/(K−λ), understood as polynomial division, to two successive differences of neighboring odd stars going around the gateway. The cyclic source obeys F(K)s = 0, so Q equals F′(λ) times its first spectral projector there. The differences have zero fivefold average and therefore lie in the doublet after projection. Call the two-column results V at launch and W at receipt. Their Gram matrices are nonsingular, as the certificate verifies. On the full board use the integer K = 2H of the free two-load law. Put u = −it/2 and T(u) = Wᵀ exp(uK)V.

The Taylor coefficients are computed without diagonalizing the full board. Reduce the integer entries modulo 101 and substitute λ = 69, a simple root of F modulo 101. The polynomial F is irreducible over the rationals. A nonzero reduced coefficient is thus a nonzero algebraic coefficient at the physical root; the factorial denominators through the tested order 64 are invertible modulo 101. At degree 18 the coefficient of det T is 69 modulo 101 for both active and exchange-only flight. The determinant is analytic and not identically zero, so both independent input amplitudes transmit for all real exposures except a discrete set.

A unitary doublet link with one common efficiency would preserve all norm ratios. In these unnormalized coordinates put G = VᵀV, D = WᵀW and L(u) = T(−u)ᵀ adj(D) T(u). Common efficiency would require L proportional to G. At degree 18 the coefficient of L₁₁G₂₂ − L₂₂G₁₁ is 98 modulo 101 for the active law and 85 for exchange alone. Thus unequal efficiencies also persist outside a discrete set. This concerns the free transfer between these calibrated endpoint modes on a single board. Establishing a weak gauge connection requires a common matter action and compatible transports on further paths and scales. [Integer certificate](certificate.py) · [Native flight and detector](measure.py)

From the repository root, use a fresh path for the generated records.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -m model.experiments.interval_echo_gauge.doublet_flight.measure --output /tmp/doublet-flight.json
.venv/bin/python -B -m model.experiments.interval_echo_gauge.doublet_flight.certificate
.venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/doublet_flight
~~~
