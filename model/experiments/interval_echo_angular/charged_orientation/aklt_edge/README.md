# Fractional spin at the end of a blank qutrit chain

A spin-one chain can have a two-state boundary degree of freedom on which perpendicular internal half-turns anticommute. This experiment asks whether the three load levels of registered Interval Echo generate such an endpoint. It compares the registered bond with the spin-one AKLT projector on identical physical port addresses. Both interactions preserve the two tested global half-turns and the total load. The AKLT interaction and its spin interpretation are supplied comparison inputs. [AKLT construction](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.59.799)

## What the registered bond permits

On a blank seam the matched generator is
\[
h_{ij}=(I-S_{ij})/2+8(n_i-n_j)^2 .
\]
Both terms are positive. A zero-energy vector has zero amplitude on every load string that differs across that seam. Conversely all three equal-load pairs cost zero. On **any finite connected graph** of such bonds, with the degree-correction rails erased, the entire ground space therefore spans the three constant strings. This includes coherent superpositions with entanglement across an arbitrary cut.

Let an operator act on any proper subset of the ports. At least one port outside that subset retains its original load. The operator's matrix element between different constant strings vanishes at that port. Its compression to the whole ground space is diagonal. Two such compressed operations commute on every input, even when their supports overlap. This rules out a locally controllable projective endpoint in this sector. It applies also to endpoint-supported operators which temporarily leave the ground space and return within the same local intervention. A sequence that uses the rest of the chain as a nonlocal control has a different domain. [Lean proof of the connected-support and compression steps](EndpointAlgebra.lean)

There is a second test at fixed total charge \(Q=\sum_i(n_i-1)=\pm1\). With one defect in a filled host the restriction of the registered bond sum is \(8.5 D-0.5 A\), where \(A\) and \(D\) are the host adjacency and degree matrices. Each sign has a unique sector ground vector on a connected host by irreducibility and Perron–Frobenius. On an open chain the uniqueness extends to **all** states of that charge. A histogram with all three loads pays at least two unit contrasts. A histogram with only zeroes and twos pays at least one contrast of magnitude two. An isolated defect at the end has expectation 8.5, strictly below the other histograms' lower bound of 16. Histogram preservation under exchanges completes the comparison.

The two lowest open-chain one-defect orbitals can be close in energy. On six physical sites their splitting is 0.0000118903 native energy units for either sign. A sharp arrival of the specified charge at the left end projects onto exactly one material word, including all coherences between those low orbitals. The finite energy filter selects the lowest mode with probability exceeding 99.994% conditional on admission. About half the selected defects are seen at the left end. These records distinguish a pair of end locations from an internal two-state record at a single charge arrival.

## An internal spinor that the detector can resolve

The alternative bond is the projector onto total spin two of a pair of spin-one records, with \(S_z=n-1\). The source contracts a triplet at each site with oriented singlets on adjacent auxiliary indices. Full diagonalization independent of that construction gives four ground states on each tested open chain. The tensor amplitudes span these ground states. Keeping the right boundary labels as two levels of an ordinary reference qutrit gives an isometry from the left boundary input into the chain and reference.

Measuring the terminating qutrit along either Cartesian basis element \(x\) or \(y\) shortens the chain by one site. Contracting the **physical** measured amplitudes with the shorter isometry reconstructs two boundary operations. Each outcome occurs with probability one third on every logical input, and its normalized operation is unitary. Acting with the corresponding spin-one half-turn on every physical qutrit factorizes in the measured boundary frame as \(U\otimes\overline U\). The commutator of the two normalized left operations is \(-I\). This is an internal spin action of the alternative interaction. A physical rotation of observers' spatial frames, a charged mobile excitation, a fermionic exchange and a rest mass require independent observations.

## The two-port read

The same detector checks the relative sign between the two endpoint outcomes \(xy\) and \(yx\). In the load basis the negative combination is the antisymmetric \(02,20\) state. The positive combination lies on \(00,22\). For the charge-invariant source mixtures the charge-dephased positive detector is exact and accepts either equal extreme pair with probability one half. The negative detector uses a balanced reference arm and one native blank π-contact, which is the full SWAP on all nine input columns. An odd swap result and final unordered load pair \(\{0,2\}\) give the negative event. Other results are retained. Both comparisons use this same supplied contact and interferometric read service. The two-path read avoids conditioning on the final sign.

The source is a supplied uniform mixture of load strings, with compensating charge in the external stock. The chains are actual shortest paths of three through six ports on the level-zero federation. All unused touching seams are dark zero-mark fixtures; the same native columns erase degree-five private costs. A coherent clock preparing either material uses the same zero center, tick 1/32, and \(2^{28}\) delays. Exact ground-projector controls are recorded separately. Each Monte Carlo bank has \(2^{22}\) original simulated starts and retains preparation refusals.

On six sites, the finite clock and registered bond gave **5,762 positive-order and zero negative-order records**. The AKLT bond gave **zero positive-order and 5,189 negative-order records**. The admitted rank-three registered mixture has positive probability one third in the ideal ground. The ideal AKLT mixture has negative probability about 0.22223. The finite clock leaves excited fractions below \(2\times10^{-12}\) in both six-site mixtures. The sign comparison applies to those source mixtures; the registered local-algebra obstruction applies to every ground state. [Source, readout and physical paths](measure.py) · [Recorded probabilities and counts](readings.json)

Run from the repository root.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_angular.charged_orientation.aklt_edge.measure
python -m pytest -q -W error model/experiments/interval_echo_angular/charged_orientation/aklt_edge
(cd Lean && lake env lean -j 1 ../model/experiments/interval_echo_angular/charged_orientation/aklt_edge/EndpointAlgebra.lean)
```

Append `--write` to the measurement command to regenerate the probability and count receipt.
