# Field energy reaches two independent receivers

A signed material intervention injects energy into one face of the written-field sphere. The charge is parked before the free flight. Two rings, owned by different cells from each other and from the source, can both register energy. Their pointers stay silent on the uniform ground and when the source face is disabled. Reversing the charge gives exactly the same energy probabilities. The separate signed material probe reverses its response. [Design](DESIGN.md) · [Complete survey and counts](results.json) · [Signed material-probe response](../README.md)

| Sphere | Source and receiver faces | Field time | Face distances from source | Joint energy probability per release | Simulated double clicks, positive / negative source |
| --- | --- | ---: | --- | ---: | --- |
| 80 cells | 19 → 16 and 21 | 128 | 5 and 2 | 0.000225622 | 220 / 231 |
| 320 cells | 89 → 101 and 120 | 32 | 8 and 10 | 0.0000843550 | 74 / 101 |

The table gives the largest calculated coincidence for each sphere in the declared survey. At these settings the source supplied 0.0400000 and 0.0194081 mean field energy units, respectively. Faces sharing no source cells held 0.0304831 and 0.0185042 units. Each count uses 1,048,576 fresh simulated releases per sign. The same receivers gave zero clicks on the sham and source-face-cut controls, and double clicks on every trial of the supplied two-ring calibration. The stored record includes the source stay or move outcome, both receiver outcomes, the tested face pairs and every time. [Acquisition and controls](measure.py)

Each detector is the local positive face term already used inside one donor cell. It projects a pair of legal circulating configurations onto their antisymmetric combination. A coherent pointer across the named face is an ideal supplied instrument. All free field faces of the spherical fixture remain active. The full sector has 75 or 2,061 words. For two disjoint faces the common positive subspace has one dimension per legal four-word square. We checked this representation against both orders of the two original projectors. Between 503 and 535 disjoint pairs per source were tested on 80 cells, and between 11,484 and 11,635 on 320 cells. Many pairs have zero joint support for every preparation. The two displayed pairs have ranks one and eighteen and admit a directly verified double-click calibration. [Donor-cell detector definition](../../electron_field/lepton_family/decay/README.md)

The two opposite charge preparations have unnormalized stay amplitudes that add to the field ground state. Free flight preserves this sum. Any of the energy projectors annihilates the ground, so it maps the two amplitudes to opposite vectors of identical squared norm. On the moved branches the source face is blocked and the remote pointers are silent. The same argument works for the product of two commuting energy projectors. Charge sign can therefore be read by the signed probe while the energy pointers measure the sign-even energy deposition. [Formal sign identity](../../../../../../Lean/SourceLaw/WrittenChargeSupport.lean)

These pointers project local ring energy and leave a post-measurement field. Inferring photon number from their coincidences would require a calibrated map to absorption of a field quantum. In laboratory single-photon anticorrelation measurements, a triggered source and two output detectors after a beam splitter supply that different apparatus. [Grangier, Roger and Aspect](https://doi.org/10.1209/0295-5075/1/4/004)

A further comparison tests whether one material traversal implements the field operation used by a quantum Wilson circuit. For either charge and traversal direction, one of the two written orientations is an eigenstate at the starting endpoint. Its charge cannot cross under any duration of that isolated contact. At the calibrated pulse, the two singular values of the material crossing map are one and zero on both the 99-state contact and its positive written repair. Any lossless link acting on both encoded head inputs differs by at least one in operator norm. This excludes a direct one-pulse identification. A route using more states and operations would need its own executable transfer test. [Microscopic contact check](measure.py)

From the repository root, use a fresh output file.

```sh
.venv/bin/python -B -W error model/candidate_laws/interval_echo/openwave/photon_identity/energy_receivers/measure.py --output /tmp/written-field-energy.json
cmp /tmp/written-field-energy.json model/candidate_laws/interval_echo/openwave/photon_identity/energy_receivers/results.json
```
