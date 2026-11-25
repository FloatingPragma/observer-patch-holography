# Which energies transmit the charged orbital qubit?

**The qubit crosses the live trimer boundary through a superposition of different energies.** On the two dynamical components accessed by the prepared launch heads, all 336 energies are distinct. At one exactly resolved energy the receiver can retain only one source amplitude. The intrinsic transfer of an unknown qubit, verified on the unfiltered preparation, requires coherence across energies. This statement concerns the two written triangles, their three live coins, and the actual pair of source orbitals. [Original orbital transport and tomography](../README.md) · [Integer certificate](certificate.py)

## Exact single-energy test

Take the original 20-cell preparation, including the complete sector reachable from the two tip rails. The two choices of launch-head contrast parity each give a 168-dimensional connected component containing both prepared rails. In each, form the columns $x,Ax,\ldots,A^{167}x$ from the first sharp source rail, with $A=8H$ the integer native generator. Modulo the prime 65521, their determinants are 59318 (even head) and 19952 (odd head). The two column systems therefore span their components over the reals. Repeating the certificate from the second sharp source rail gives determinants 36521 and 25337. Each of the two sharp inputs thus has strictly positive weight on every eigenline of its component. A superposition can have zero weight on a particular line. A cyclic vector for a real symmetric matrix contributes at most one independent column direction per distinct eigenvalue, so each component has 168 simple eigenvalues. Expressing $A^{168}x$ in the same cyclic basis gives its monic characteristic polynomial modulo 65521. The two polynomials are coprime there. By the monic integer characteristic polynomials and Gauss's lemma, the two components also have disjoint spectra over the complex numbers. Arbitrarily small splittings are kept distinct by this proof. [Characteristic polynomials and certificates](results.json)

At a simple energy $E$ let $v_E$ be its normalized eigenvector, $C$ the source isometry and $R_c$ the two-tip receiver for final coin word $c$. The incident coincidence map is
$$
 K_{E,c}=R_c^*v_E\,v_E^*C.
$$
Every $c$ has the same input bra. For a reference initially maximally entangled with the input, a successful selection at this energy gives a product between the reference and the receiver plus its coins. Discarding the coin word preserves separability. Forgetting the outcome of a perfect energy measurement likewise yields a mixture of separable outputs. The finite source weights $C^*P_EC$ and all eight coin-resolved real $2\times2$ matrices are deposited for each of the 336 lines. They are numerical amplitudes in the displayed energy order. Very close individual numerical eigenvectors require more accuracy than a bound on the ordered eigenvalues alone. The exact simplicity and rank statements are independent of their numerical orientation. [Line arrays](lines.npz) · [Reconstruction of both original unfiltered channels](test_energy.py)

At the even launch parity the exchange reflecting the two cells and reversing the direction coin of the bridge commutes with $H$. Resolving this additional symmetry stabilizes the closely spaced opposite-reflection branches. A single native energy projection anywhere on the joined triangles would itself be a nonlocal postselection; subsequent free evolution multiplies that line by a scalar. Thus an arrival recorded immediately after an ideal projection cannot be counted as a flight caused by that projection.

## Which frequencies reach the receiver?

The first source rail concentrates near 66.8725 units; the second concentrates near 101.8273. The two mixture components still have overlapping spectral weights. The leading transfer from the first rail to the two receiving tips involves opposite cell-reflection lines at approximately 66.872509274 and 66.872528962. Their computed separation is 0.00001968813 native units. The exact ordered-eigenvalue certificate places the separation between **0.000019687735 and 0.000019688532**, rounded outward. Energies use the native trimer's displayed formation convention. The difference frequency is invariant under a common shift of energy.

For this certificate round the orthogonal diagonalizer to a dyadic matrix $X$ and the energies to a dyadic diagonal $D$. Exact integer arithmetic bounds $g=\lVert X^*X-I\rVert_2$ and $r=\lVert HX-XD\rVert_2$ using row and column sums. Since $g<1/2$, the orthogonal factor $U$ of $X$ obeys $\lVert U-X\rVert_2\le g$. Thus $\lVert H-UDU^*\rVert_2\le r+(\lVert H\rVert_2+\lVert D\rVert_2)g$. Weyl's eigenvalue inequality gives maximum ordered errors of $1.99\times10^{-10}$ (even) and $2.60\times10^{-10}$ (odd). The comparison with the measured beat does not require resolving all lines to the smallest pair splitting. [Dyadic verifier](certificate.py)

On fresh repetitions prepare the first sharp source and measure only the two tip arrival outcomes, with every other outcome retained in a third bin. Average two independent uniform timing jitters from 0 through 4096 over each of 1024 appointments, spaced by 4096. The two jitters multiply the contribution of frequency difference $g$ by $\mathrm{sinc}^2(4096g/2)$, where $\mathrm{sinc}(u)=\sin u/u$. They suppress fast beats while preserving the large sublevel transfer. Each appointment has 262,144 simulated starts. A discrete Fourier peak followed by a one-tone least-squares fit in that frequency bin gives **0.00001968742** in the first 512 appointments and **0.00001968765** in the other 512. This read uses unfiltered local load records. The remaining physical tones account for deviations of up to 0.00987 in the bright probability in the numerical spectrum. The largest observed deviation from the dominant cosine is 0.01095 with a simultaneous 99% shot allowance of 0.00483. The fitted values also include interference from the remaining tones and the finite duration of each window. An independent native dimer clock measures 1.000000014 per scheduled unit. [Complete count banks, settings and dimer record](results.json)

A separate receiver measures total signed load in the destination triangle, retaining all four possible values from −2 through +1. Its mean is **0** initially, **−0.36398** at 65,536 and **−0.89102** at 131,072. The only active seam crossing the boundary is the bridge. The change in total signed load measures its integrated net current, including launches where the three defects arrive separately. The tip coincidence interrogates the full outgoing pattern. [Virtual observer instruments](measure.py)

## How much bandwidth does the example use?

At the previously tested exposure 131,072, insert the mathematical projector onto the indicated half-open interval before free evolution, coherently retaining all its lines. All probabilities below include initial weight outside the band as a failure. The eigenvalues of $C^*P_{\rm band}C$ bound the weight of every normalized source input, and the partial transpose uses the unnormalized successful channel.

| Even-parity window | Lines | Minimum and maximum input weight | Mixed arrival | Smallest partial-transpose eigenvalue |
| --- | ---: | ---: | ---: | ---: |
| [66, 67) | 3 | 4.18 × 10⁻¹¹ to 0.79967 | 0.29523 | −1.29 × 10⁻¹² |
| [101, 102) | 8 | 1.07 × 10⁻¹¹ to 0.95890 | 7.68 × 10⁻¹⁰ | within 10⁻²² of zero |
| [81, 84) | 14 | 0.03345 to 0.06573 | 0.000863 | −0.000423 |
| [66, 84) | 55 | 0.04054 to 0.99363 | 0.42188 | −0.012332 |

The last window gives a rank-two source effect and reference-preserving flight across distinct energies. These are calculations for an ideal coherent band projector. No such gate is assumed by the local beat or charge measurements. Implementing the preparation and measuring a freely propagating two-component dispersion relation are separate requirements for a particle mass assignment. [All per-parity windows](results.json)

From the repository root run the following commands.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.live_trimer.energy.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/live_trimer/energy
```

An unrecorded uniform free wait before the same flight provides an operational finite-resolution test. Entanglement remains detectable for a timing width of two units. At width 64 the accepted reference–orbital negativity is at most 0.00111 while the dominant charge-transfer beat retains more than 99.99999% visibility. The synchronized comparison at the same mean wait retains negativity of at least 0.02741. [Random-wait experiment and full outcome banks](WAITING.md)
