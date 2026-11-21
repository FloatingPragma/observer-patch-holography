# Global image of the registered bound-pair currents

The twelve Interval Echo currents generate **U(3) × SO(3)** on the two retained triplets of the prepared pair. This fixes the faithful connected image of this one-cell instrument. It also identifies a missing premise for attaching the Standard Model global group. A full rotation of the weak triplet is invisible on the entire pair register. A weak doublet would distinguish the two lifts of that rotation. This is an internal current experiment; physical spacetime spin and propagation of a flux around a defect require separate instruments. [Frozen protocol](PROTOCOL.md) · [Raw coefficients, words and detector outcomes](result.json)

## From twelve addressed services to a kernel

The registered source determines two real three-dimensional band bases. The positive rotation projector in the current formula keeps the first band and annihilates the second. On the second band each Hermitian current therefore contains only the imaginary circulation on the neighbour ring. Its anti-Hermitian generator is real skew-symmetric. On the first band the anti-Hermitian generators lie in u(3). The projected real ranks are 9 and 3, and the joint rank is 12, with smallest singular value 0.36676. The twelve generators therefore span the whole direct sum. Products of their flows give the connected image U(3) × SO(3). On the 60 unused pair modes the registered operations act as the identity. [Native-current construction](../../candidate_laws/interval_echo/CURRENTS.md) · [Rank and block checks](__main__.py)

Calibrate the primitive central generator to act as +1 on the first triplet and 0 on the second. It equals half the sum of the twelve native currents. Let z denote its unit-circle phase. In these coordinates the covering-product action sends (A, B, z) to the block matrix diag(zA, Ad B), with A in SU(3) and B in SU(2). For the whole operator to be the identity, B must be either sign of the identity, A must equal z⁻¹ times the identity, and z³ must equal one. With ω = exp(2πi/3), the kernel is exactly

$$
\{(\omega^a I_3,(-1)^b I_2,\omega^{-a}) : a=0,1,2,\ b=0,1\}.
$$

The weak centre alone is in this kernel. The derived subgroup of the faithful image is SU(3) × SO(3), while the derived subgroup of S(U(3) × U(2)) is SU(3) × SU(2). This difference holds independently of the convention for the abelian phase. A larger group could still act with an undetectable kernel on this pair. The gauge paper constructs a different carrier by a spin lift and trace redistribution, or by a charged triplet and doublet under an additional port-phase construction. Realizing such a response on independently accessible matter and transporting it through the federation are open premises. [Gauge theorem, matrix current and global form](../../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex)

## What a patch can distinguish

Three generators made from the address currents provide a centre turn of SU(3), a full rotation of the weak triplet, and a sixth turn of the primitive centre. All pulses execute in sequence, including reverse address order in each symmetric step. At 512 and 2,048 steps, the largest operator errors over all 36 combinations were 0.00004295 and 0.00000269 respectively. Both runs resolved the same six identities. The weak half-turn sits operator distance two from identity. Completing the turn gives identity in the integrated action. The calculation uses the registered ideal guarded pulse at each address. The original compiler establishes its realization by native gates; these long words use exact pulse arguments and supply the schedule. [Executable](__main__.py) · [Address selection and port coefficients](result.json)

For a local measurement, we used the same antipodal edge-pair source and the same 30-seam-plus-quiet detector as the registered current experiment. Only the intervening sequence of current pulses changed. Ten seam outcomes were predicted to become more likely under the diagonal comparison. All 31 counts per setting were saved. Each arm used 8,192 independent simulated preparations.

| Word (SU(3) cube roots, SU(2) signs, circle sixth turns) | Total variation from sham | Counts on the ten seams |
| --- | ---: | ---: |
| Sham (0,0,0) | 0 | 5 |
| Weak centre alone (0,1,0) | 0.000000086 | 0 |
| Native compensation (1,1,4) | 0.000000058 | 2 |
| Conventional diagonal word in this primitive normalization (1,1,1) | 0.796034 | 6,556 |

In this calibration the last word gives opposite phases to the two triplets. The terminal patch resolves their interference. The tiny calculated differences in the two identity arms come from finite pulse slicing. Calling the abelian charge of this pair a physical hypercharge would require an independent calibration against other matter. These words take place inside one complete cell using the pair’s currents. The source-to-loop equality required for a physical global quotient remains open. [All terminal distributions and counts](result.json)

From the RER repository root:

```sh
python -B -m model.experiments.interval_echo_gauge
python -m pytest -q model/experiments/interval_echo_gauge
```
