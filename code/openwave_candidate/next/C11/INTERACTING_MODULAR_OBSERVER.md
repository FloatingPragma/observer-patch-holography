# Modular observer on the actual interacting source

The input is `canonical_states.npz:selected_protocol_poststate` from the selected C10 construction. Its one 384-state interacting H, actual global N=3 reference, .37 source-hop interval and .29 MIX event are reused. No older one-load modular output, Gibbs preparation or depolarizing floor is transferred. The new driver calls existing `faithful_spectrum`/`modular_flow` and `GroupedEnergyInstrument.continue_branch` APIs.

The declared M2 observer corner consists of two full occupation/record configurations: mask 7 with F=(5,1,5,1,5,1), and mask 11 with F=(1,0,0,0,1,5). This is a supplied coherent global apparatus, not a local-port algebra or partial trace. Its inside probability is 0.0119141715753 and complement probability 0.988085828425. Conditional eigenvalues are 0.041499758393 and 0.958500241607, giving dimensionless modular frequency 3.139682209194. The full source has numerical rank 12 at threshold 1e-10 and is not faithful on M384. The corner is numerically faithful; no floor or pseudoinverse is inserted.

The apparatus was selected after inspecting this known simulated output. A finite lexicographic basis-pair search required both raw diagonals above 1e-8, conditional minimum eigenvalue above 1e-4 and raw offdiagonal magnitude above 1e-4; it selected the disclosed configurations. Execution fixes that pair. The modular frequency is conditional on this supplied apparatus, not a blind physical prediction or a uniquely selected observer. State knowledge, tomography and sample/preparation costs are not derived.

The standard four-dimensional Hilbert-Schmidt/GNS representation uses Delta=(rho^-1)^T tensor rho and Omega=vec(sqrt(rho)). These are algebraic representations, not an additional physical reference preparation. The convention sigma_s(A)=rho^(is)A rho^(-is) has the minus-i KMS shift. Each density eigenvector's largest component is positive real, fixing the declared apparatus eigenframe. This phase convention does not change the state.

The actual observer apparatus first distinguishes inside/outside, then reads A0 followed by As inside. Four ordered sign pairs plus the full outside outcome form one complete instrument. Their full subnormalized poststates and all complete native continuations are archived. No complement is discarded. The sequential real correlation equals the real part of Tr(rho A0 As); the complex algebraic correlation is not a classical joint probability or a measurement of subjective duration.

As is computed from the frozen premeasurement conditional state. A pure observed branch is not faithful on the same M2 algebra; its old modular generator is not silently reused as that branch's state flow. Subsequent native evolution uses the actual full interacting H and MIX, with .37/.29 supplied again, and retains every native record and no-commit complement.

The corner is not invariant under H: coherent source time .37 leaks 0.3622350216 outside it. Even compressed H differs from every scalar multiple of K=-log(rho), up to scalar identity, by residual 0.5004137013. The modular parameter therefore has no established autonomous native-clock, physical-time or subjective/experienced-time conversion. The apparatus and both parameter choices are declared resources. Propagation and faithfulness checks are numerical, without interval enclosures. No physical species, vacuum, mass, gravity or official row award follows.

From the public bundle root:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python next/C11/interacting_modular_observer.py --write
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python next/C11/interacting_modular_observer.py --check
python -m pytest -q next/C11/test_interacting_modular_observer.py
```
