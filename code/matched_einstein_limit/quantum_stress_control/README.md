# One Maxwell quantum and the optical response of its box

The [existing charged Wilson–Maxwell experiment](../../../model/candidate_laws/interval_echo/openwave/photon_identity/quantum_transfer/) supplies a quantized field, its free Hamiltonian and its charged four-circuit preparation. On fresh instances of the same box, a source pointer of zero followed by a number-one filter prepares a single quantum. Dark controls start in the vacuum. Independent measurement ensembles address the same preparations for field stress and optical propagation.

| Box side | One-quantum spectral energy | Measured spatial stress trace at launch | Maximum optical change from adding the photon | Maximum change from the clock control |
|---:|---:|---:|---:|---:|
| 4 | 2.243025 | 2.24295 ± 0.00107 | 1.11 × 10⁻¹⁶ | 0.47737 |
| 6 | 2.245429 | 2.24353 ± 0.00190 | 2.50 × 10⁻¹⁶ | 0.51343 |
| 8 | 2.245684 | 2.24331 ± 0.00286 | 2.78 × 10⁻¹⁶ | 0.44878 |

Energy and all three integrated pressures use the dimensionless box Hamiltonian's units. The uncertainty is one simulated sampling standard error. On side 8 the exact integrated pressures along the three axes are 0.62448, 0.62448 and 0.99673. There are 24 gauge-invariant signed optical channels on each side. A 20% increase in the free clock rate gives the last column, with preparations and detector coefficients held fixed. This is a sensitivity control.

## Stress instrument

At unit edge lengths, the Hamiltonian is half the sum of electric flux squares on edges and magnetic flux squares on faces. A homogeneous logarithmic strain $s_i$ along direction $i$ gives the coefficient $\exp(2s_i-\sum_j s_j)$ for flux along $i$, at fixed canonical flux and face circulation. The integrated Maxwell pressure is the negative metric derivative $S_i=-\partial H/\partial s_i=H-X_i$. Here $X_i$ is the sum of electric squares on parallel edges and magnetic squares on normal faces. It follows that $\sum_i S_i=H$. The strain convention defines a pressure instrument on the supplied uniform box. The free evolution uses the original unit edge lengths throughout.

For the normalized one-quantum packet $\psi$, the excess square of the quadrature $Q(v)=a(v)+a^\dagger(v)$ over vacuum is $2|v^\dagger\psi|^2$. Every local edge or face uses the box's transverse modes and its own normalization. For an electric quadrature an overall factor of $i$ can be dropped from its squared-outcome distribution in these number states. The instrument samples the spectrum of the total-occupation-capped quadrature. At one quantum, a fraction $u=|v^\dagger\psi|^2/\|v\|^2$ occupies the measured mode. The orthogonal component carries a spectator, reducing the available occupation cutoff by one. Vacuum subtraction uses independent draws. A second instrument samples the frequency of the occupied free mode. Both instruments give the same energy within their sampling errors. The spatial stress trace here is the sum of three pressures. The full Lorentzian Maxwell stress tensor has a vanishing four-dimensional trace in the uniform metric convention.

## Optical instrument and its exact limit

An ordered Wilson circuit around the source face applies a small signed kick. The detector reads the quadrature of an independently addressed face at a later time. Its first derivative with respect to kick strength is

$$
\chi_{rs}(t)=i\langle[Q_r(t),Q_s]\rangle
             =-2\,\mathrm{Im}\langle v_r(t),v_s\rangle,
\qquad v_r(t)_k=e^{i\omega_k t}v_{rk}.
$$

The free evolution rotates each mode vector by the same state-independent phase. The canonical commutator is a multiple of the identity. The parent's cap of three preserves this commutator on the vacuum and all states supported at total occupation at most two, since the one-operator intermediate states remain inside the cap. This proves equality for every pair of loop quadratures and every wait on these preparations. The simulations evaluate the two state expectations separately before checking this identity.

Finite signed kicks use the parent's ordered edge phases. On side 6 the exact difference between the two finite-amplitude slopes falls from $1.75\times10^{-5}$ at strength 0.12 to $1.10\times10^{-6}$ at 0.06 and $6.90\times10^{-8}$ at 0.03. The upper occupation boundary modifies finite drives. Extrapolation to zero kick and the exact commutator distinguish this effect from a change in the propagation law.

This control establishes a same-preparation gap. The charged-circuit photon changes local square statistics and carries integrated stress. The box's optical propagation kernel stays fixed by its specified Hamiltonian. An observer using only that kernel obtains the same inferred geometry in both preparations. Other field statistics, including the existing photon absorption instrument, remain sensitive to the state. This calculation treats the Wilson–Maxwell sector as a separate quantized comparator. Demonstrating gravity in the registered material requires its own co-evolving records of field stress and independent geometry. Its cochain geometry is prescribed in the box constructor.

## Reproduction

Run from the repository root with the Python environment used by the parent photon experiment.

```sh
.venv/bin/python code/matched_einstein_limit/quantum_stress_control/measure.py --output /tmp/maxwell-stress-control.json
.venv/bin/python -m pytest -q code/matched_einstein_limit/quantum_stress_control/test_control.py
```

The [protocol](protocol.json) fixes box sizes, preparations, measurement schedules, controls and enrollments. The [packet](results.json) includes original cochain hashes, the parent code hash, spectrometer counts, summed local square and fourth-power records, finite receiver eigenvalues and their counts, and all weak kernels. Counts are aggregate multinomial draws from the finite simulated instruments. Every time and every face is expressed in the supplied box coordinates.
