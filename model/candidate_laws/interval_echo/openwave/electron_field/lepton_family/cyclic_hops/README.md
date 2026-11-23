# Carrying a three-state register with a signed load

**A three-state register can travel with the material contact on the verified guide. Its three characters give at most two distinct lowest energies, and that splitting vanishes on a long guide.** At 48 cells, the largest threshold ratio is 1.00041112 and the Koide quotient is 0.333333336462. Each character is conserved, so this coupling supplies no decay between them. This is a comparison extension of the proposed written-field model. The three register states and their coupling to a material hop are additional physical assumptions. [Design](DESIGN.md) · [Executable and controls](measure.py) · [Numerical record](results.json)

## A register on every allowed hop

Keep the 144 route ports, the native hopping weight 4 and the degree-five cost 8 of the separated-end calculation. The neutral reference stays at zero. A forward hop on seam $e$ increments the mobile three-state label by $a_e\in\{0,1,2\}$ modulo three, and a reverse hop undoes it. Denote the one-step label shift by $R$. The off-diagonal forward contact becomes $-4R^{a_e}$ and the reverse contact is its adjoint. The positive contact diagonal remains 4 per incident seam. This construction adds no species mass, new hopping scale or selected lepton phase. [Verified route and original guide](../bands/README.md) · [Physical route ports](../bands/route.json)

Relabeling the fibre at each successive port by its accumulated shift removes all $a_e$ up to the final seam. That seam carries $h=\sum_e a_e\bmod 3$, the net register change in one circuit. All $3^{144}$ edge assignments consequently fall into three spectral classes. A reader fixed in a specific local register basis must be transformed under the same relabeling to preserve its measured response. In character $j=0,1,2$, the boundary contributes the phase $e^{-2\pi ijh/3}$. The original three-port Bloch Hamiltonian is sampled at $k=2\pi(\ell+jh/3)/N$ for $\ell=0,\dots,N-1$. The full Hamiltonian commutes with the global shift $R$. With a sharp register label, crossing a seam may change the label. With a sharp character, the transition probability into either other character is zero at every time. [Full matrix and symmetry checks](measure.py)

| Net shift around the physical ring | Character 0 threshold | Character 1 threshold | Character 2 threshold | Koide quotient |
| --- | ---: | ---: | ---: | ---: |
| 0 | 1.753788749 | 1.753788749 | 1.753788749 | 1/3 |
| 1 | 1.753788749 | 1.754509769 | 1.754509769 | 0.333333336462 |
| 2 | 1.753788749 | 1.754509769 | 1.754509769 | 0.333333336462 |

All energies are native generator units. The lowest original band has its minimum at $k=0$. For either nonzero net shift, the two nontrivial characters see opposite phases. The guide is real in the sharp-label basis, so complex conjugation equates their spectra, including the minimum. The program verifies the complete 432-dimensional real-space spectrum character by character, and verifies relabeling and global-shift invariance by exact integer comparisons. Sizes other than 48 below repeat the same motif on mathematical rings. Only the 48-cell route was checked on the physical carrier. [All nine sector spectra checks and size series](results.json)

| Cells | Nontrivial minus trivial threshold |
| ---: | ---: |
| 12 | 0.011528607 |
| 24 | 0.002883696 |
| 48 | 0.000721021 |
| 96 | 0.000180261 |
| 192 | 0.000045066 |
| 384 | 0.000011266 |

The nearest allowed nontrivial momentum has magnitude $2\pi/(3N)$. With the measured lower-band curvature 0.757464, the splitting is $1.66130527/N^2+O(N^{-4})$. A sharply localized observer can distinguish preparations in the added register given suitable apparatus. These finite-size threshold differences do not provide three bulk rest energies. [Underlying dispersion](../bands/results.json)

## Adding a balanced register energy

A separate comparison uses $H_{\rm field}\otimes I+s I\otimes M^\dagger M$ for any $s>0$. Suppose the circulant amplitude $M$ obeys the geometric balance and lies in the positive chamber, so that its three squared eigenvalues $p_j$ alone have $Q(p)=2/3$. The three absolute charged thresholds above the same neutral reference are then $E_0+s p_j$, where $E_0=10-2\sqrt{17}>0$ is the guide threshold. Since the terms act on independent factors, the propagator factorizes. Each character survives, with zero amplitude to convert into another. The squared eigenvalues are supplied by the chosen register operator. [Geometric identity and positive chamber](../../../../../../../extra/koide_identity_from_positive_c3_face_circulants.tex) · [Additive check](measure.py)

The quotient is homogeneous. Write $c=E_0/s$, $x_j=p_j+c$, $S=\sum_jx_j$, $A=\sum_j\sqrt{x_j}$ and $B=\sum_j1/\sqrt{x_j}$. Differentiating gives $dQ/dc=(3A-SB)/A^3$. For these unequal thresholds,

$$
SB-3A=\sum_{i<j}(x_i-x_j)(x_j^{-1/2}-x_i^{-1/2})>0.
$$

Thus $1/3<Q(E_0+s p)<2/3$ for every finite $s>0$. Increasing $s$ raises the quotient continuously toward the geometric value. Subtracting $E_0$ restores the isolated register's identity but changes the measured quantity from absolute charged energy to energy above its common guide floor. For illustration, $a=1$ and $\delta=\pi/24$ give total $Q=0.375533$ at $s=1$ and $Q=0.639965$ at $s=100$. The chosen phase illustrates the proof and uses no mass inputs. A physical family with distinct persistent rest energies and decay would need dynamics beyond either tested coupling. [Numerical scales and isolated register](results.json)

## Replay

From this directory run `python -B -W error measure.py --output /tmp/cyclic-hop-family.json` with NumPy and SciPy and a fresh output path. The run checks the canonical carrier, samples all three bands and records SHA-256 digests of its code, design, route and band implementation.

The counted phase of the 432-state single-load register is 2/9. With its own compression normalization, its positive effects sum to 2/9 and give ratios 1 : 206.770316 : 3477.472837. Adding these effects to the native guide floor with illustrative unit coupling gives ratios 1 : 1.007075 : 1.119528 and Q = 0.333559530. The all-coupling inequality also covers this phase. The unit coupling is a comparison choice. [Historical compression](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_register_round/PROOF.md)
