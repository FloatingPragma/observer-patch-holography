# Does an ordinary load move inward under the written energy source?

At the primary time four, **the added ordinary-load response shifts outward by 0.069105 graph hops** when the correlation on the written source is erased. A separate ordinary source placed on the same blank board shifts its added-load response inward by 0.159062 hops at that time. Both launches, one from each side, give these signs. Earlier and later appointments give different signs. [Design and primary sign test](PROTOCOL.md) · [All native port occupations](results.json)

The apparatus is the same closed assembly of two 12-port owners used for the free source experiment. A sharp 01 write on their shared seam (1,14) gives a paired memory and outgoing load. Exchanging memories between independent copies raises its ordinary energy by 17/8. The two ensembles have identical separate memory and outgoing-load distributions. We add an ordinary sharp unit at port 2 or 13, the two opposite-owner ports farthest from the source. Each is three hops away and the extra unit adds 50.5 native energy units to every branch of both ensembles. After preparation all seams run freely. [Original source, energy and one-load control](../README.md) · [Complete generator and readout](measure.py)

Let $n_x(E,2)$ be the simulated mean occupation of physical port $x$ with erased source and added unit, and $n_x(E,1)$ its source-only counterpart. Define $n_x(P,j)$ similarly for the paired source. The local response due to the source change is

$$D_x=[n_x(E,2)-n_x(E,1)]-[n_x(P,2)-n_x(P,1)].$$

The radial shift is $\sum_x r_xD_x$, where $r_x$ is the measured graph distance to the nearer source endpoint. Each one-unit insertion raises total occupancy by exactly one; summing $D_x$ gives zero. The two units are indistinguishable, so this read includes redistribution of the original load and scattering from the source hardware. Putting a distinguishable spectator on an independently evolved copy gives zero connected density. At time four the changes in the four radius shells, from zero to three, are +0.073764, −0.053052, −0.184294 and +0.163581. An increase at the source endpoints alone would conceal the outward first moment. [All shells and both launch directions](results.json)

| Free exposure | Erased minus paired at native contrast | Erased minus paired at zero contrast | Ordinary unit source versus absent, all seams blank |
| ---: | ---: | ---: | ---: |
| 0 | 0 | 0 | 0 |
| 1 | −0.040593 | −0.000463 | +0.011819 |
| 4 | **+0.069105** | −0.047041 | **−0.159062** |
| 16 | −0.229165 | −0.229088 | +0.186209 |
| 64 | −0.343614 | −0.074699 | +0.174132 |

For the ordinary comparison a single ordinary source unit starts with equal probability at either shared edge endpoint. The control subtracts its source-only density and the flight of a far unit on an otherwise empty board. The added-load response therefore measures the native interaction between ordinary matter in this apparatus. The blank binary two-load Hamiltonian includes the hard occupancy constraint and an energy reduction when the two loads occupy neighbouring sites. The erased comparison additionally couples both units to the written memory. Native degree-reference costs enter throughout. The zero-contrast control retains exchange, the same source stocks and the same measurements. [Generator and ordinary control](measure.py)

A retrospective read sampled $2^{32}$ copies for each actual binary port probability. At time four the two native launches gave +0.069138 and +0.069196. The ordinary controls gave −0.159015 and −0.159108. A simultaneous 99% half-width of 0.005734 hops covers the 2,880 underlying mean readings after propagating their individual Hoeffding bounds to the radial contrasts. This uncertainty concerns simulated read sampling. [Read design](SHARP_READ.md) · [Underlying counts and measured shifts](sharp_read.json)

## Onset and time average

The integer Hamiltonians also permit an exact radial expansion. The coefficients through fourth order vanish. For either launch the erased case starts as $-1445t^6/2304+O(t^8)$ and the ordinary case as $-85t^6/576+O(t^8)$. Both give an inward change at sufficiently short positive time. The written coefficient's magnitude is 17/4 times the ordinary coefficient, although its added energy of 17/8 is only 1/24 of the ordinary endpoint source energy of 51. Dividing each response by its added source energy gives a factor of 102 between these specific experiments. The references are physically different assemblies, one with a retained written interaction and the other with blank memory. A common response per unit energy does not describe this instrument. [Exact expansion including higher orders](time_mean.json)

The numerical infinite-time radial shifts on this finite board are −0.135608 for the written source and +0.147007 for ordinary matter. The calculation retains coherences within equal-energy subspaces. The two fixed grouping tolerances reproduce the same numbers; for the native two-load source the smallest separation between groups exceeds $2.00\times10^{-5}$ and the largest spread inside a group is below $4.0\times10^{-13}$. The time-mean and onset checks were chosen after seeing the finite-time results. These are local density responses of directly interacting matter on a closed finite board. A gravitational acceleration, physical distance scale or transfer to rays, clocks and galactic tracers remains open. [Time-mean protocol](TIME_MEAN.md) · [Spectral-projector calculation](time_mean.py)

The written sector has $4\binom{24}{2}=1{,}104$ states, including every placement of the two units and all four magnitude-one memories. The one-unit comparison has 96 states. Native column checks verify that both sectors are invariant. A regression independently reconstructs every two-unit generator column from the full 99-state matrix on each physical seam. Direct spectral evolution and a separate sparse exponential differ by at most $1.07\times10^{-11}$ in state norm over the five clocks. [Physical projection and detector tests](test_probe.py)

Run from the repository root with `PYTHONPATH=code:.` and one BLAS thread.

```sh
prefix=model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.free_probe
python -m "$prefix.measure"
python -m "$prefix.time_mean"
python -m "$prefix.sharp_read"
python -m pytest -q model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary
```
