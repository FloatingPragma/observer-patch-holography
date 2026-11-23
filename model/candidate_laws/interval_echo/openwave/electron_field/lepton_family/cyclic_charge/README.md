# Does the native threefold rotation split a charged level?

The proper threefold rotation gives three **equal-energy** Fourier modes for each intrinsic level of a single mobile charge on the first symmetric cell of the 20-cell federation. For both signs, the balance $Q=(E_0+E_1+E_2)/(\sqrt{E_0}+\sqrt{E_1}+\sqrt{E_2})^2$ is exactly **1/3**. Koide balance would require **2/3**. This geometric attachment supplies no free transition between charged background orientations.

## Physical rotation and conservation

The selected cell has label 0. Its three external contacts are at ports 1, 3 and 10. Its proper rotation cycles those three ports and preserves every physical port degree. Start with the archived sharp charged configuration, then rotate both its charge address and all 30 oriented heads. The three configurations belong to distinct reachable components. During a material hop, one field head and the signed charge change together. A face operation preserves the head count at each port. Thus each component has a fixed background vector $b_p=\text{heads at }p - s\,\mathbf{1}_{p=q}$, where $q$ is the charge address and $s$ its sign. The sum of its entries is $30-s$. An invariant vector would have entries constant on four orbits of length three and a sum divisible by three. Neither charge sign admits such a vector. The computation checks every background and every allowed transition in all three reachable components.

| Charge | Configurations per orientation | Directed material and face moves across all three | Moves between orientations | Lowest character energies, in native units |
| --- | ---: | ---: | ---: | --- |
| Vacancy (−) | 5,790 | 153,000 | 0 | 4.3437596815, 4.3437596815, 4.3437596815 |
| Excess (+) | 1,954 | 44,460 | 0 | 5.4500056053, 5.4500056053, 5.4500056053 |

In the rotated configuration bases the integer contact columns and degree potentials coincide entry by entry. Write $J_r$ for the inclusion of orientation $r$ into the whole space and $H_0$ for one complete orientation block. All three relations $H J_r=J_r H_0$ are exact integer matrix identities. Consequently, for $\omega=\exp(2\pi i/3)$, the isometry $U_k=\sum_r\omega^{-kr}J_r/\sqrt 3$ satisfies $H U_k=U_k H_0$ for every $k=0,1,2$. Every intrinsic eigenvalue, including the positive ground energy, occurs at the same value in all three character spaces. For a matched level the energy ratios are exactly 1, 1, 1. The result even allows an externally prepared coherent superposition of the distinct Gauss references. The face and material contacts conserve each orientation separately for all free evolution times.

## Broken-degree control

Cell 1 of the same physical federation has no nonidentity degree-preserving proper rotation. Keep the same sharp sources and local contact moves and use that cell’s actual degrees. The Gauss sectors remain disjoint, while the integer Hamiltonian fails to commute with the rotation. In the vacancy sector the three orientation ground energies become 5.063591, 3.985858 and 5.494361. In the excess sector they become 4.002839, 2.442386 and 4.373977. A coherent combination can acquire relative phases under this control. Its Fourier character is no longer a conserved quantum number, and each oriented component retains its population. Neither experiment identifies the three charged lepton species. The tests concern one mobile charge with the stated local fields on a finite isolated cell.

## Reproduce

From the repository root with its Python requirements installed, choose an unused output path.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/lepton_family/cyclic_charge/measure.py --output /tmp/echo-cyclic-charge.json
```

The executable verifies the archived charged-sector and head-rotation sources against their SHA-256 manifest. It regenerates every contact directly from the archived patch, checks closure, and compares the complete integer operators. [Design](DESIGN.md) · [Executable](measure.py) · [Receipt](results.json)
