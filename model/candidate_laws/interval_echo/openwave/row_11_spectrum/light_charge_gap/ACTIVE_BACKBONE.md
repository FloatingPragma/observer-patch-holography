# Formation energies on the fifteen-port winding route

A negative charge winds a physical ring while changing its seam memories. Its sixteen reachable words have one zero formation level and fifteen negative levels relative to the same memories at unit filling. The nearest level below zero is **−0.153717757** in native exchange units. The conjugate excess has the same spectrum. The motion is observable in the local charge and head records. The sign of the formation energy specifies that inserting this charge lowers the energy relative to the prepared material. Counting its magnitude as a positive rest energy would require a different physical host and a preparation and removal measurement for that host. [Winding apparatus](../../../../../experiments/interval_echo_angular/charged_orientation/winding_memory/measure.py) · [Native sector and spectra](active_backbone.py)

## Why the negative sign persists

In the positive contrast-parity sector each flippable directed bond between two reachable words contributes minus four times the projector outer product of their unnormalized sum. The remaining words have zero relative contribution from that bond. On the winding source these bonds form a sixteen-word open path. Alternating the sign of every other word changes its signless Laplacian to the ordinary path Laplacian, giving all sixteen formation energies as
[
-8[1-\cos(k\pi/16)],\qquad k=0,\ldots,15 .
]
The negative-semidefinite conclusion extends to any single-defect material graph made exclusively of magnitude-one transport links and zero-written shutters with erased private references. Each active native two-state block is minus 4 or 4.5 times the same positive outer product, according to its contrast parity. This statement also covers superpositions across the invariant contrast parities. Both physical load signs obey the local identity. [Matched native blocks](measure.py)

## One blank physical contact

The following minimum positive levels range over all fifteen positions of the changed link, starting at the winding source with all fifteen direction coins zero, or the corresponding remaining coins when a mark is omitted. The charge-conjugate measurement uses the same head-direction stock. The extra side port, when used, is an actual degree-five neighbor outside the ring. All its other physical contacts are sealed, and the five private contacts on the ring stay erased.

| Material change | Vacancy | Excess |
| --- | ---: | ---: |
| Replace one marked ring link by a blank link | 5.333333994 | 5.333333994 |
| Open a blank side link, erase the outer private reference | 3.635255691 | 3.635255691 |
| Open a blank side link, retain the ordinary outer reference | 0.456783148 | 3.666312076 |

All 90 perturbed matrices are assembled from the native local integer columns on the component reached by the preparation. Each has at most eighteen states. An exact characteristic-polynomial sign count excludes every positive eigenvalue through 5 for a replacement, through 3.5 for an erased side, and through 0.45 or 3.5 for the respective signs on an ordinary side. The calculation removes zero roots first. All roots are real by symmetry, so the coefficient sign changes count all roots above zero. Comparing that count with the polynomial shifted by the rational bound proves the exclusion. Each variant therefore misses the one-defect light threshold of approximately 0.057. The table gives finite spectra conditional on the prepared heads. The four-port blank double-well is the apparatus with independent absolute energy readings. [Four-port meter](meter.py)

To replay from the repository root, run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.light_charge_gap.active_backbone`.
