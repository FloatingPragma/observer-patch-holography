# Internal transport on every written transition

Arbitrary finite-dimensional unitary actions on the existing material hops and field-ring flips cannot reproduce the charged-lepton hierarchy in either fully written cell sector. This includes coherent spin or family rotations, finite gauge representations and noncommuting transport around closed paths. Even allowing a separately dressed neutral ground, any three positive charge gaps in the vacancy sector have Koide statistic below **0.518250**. In the excess sector the bound is **0.466141**. The Koide comparison value is $2/3$. [Scalar certificates](../family_projector/results.json) · [Connection replay](results.json)

## Uniform bounds in every internal dimension

A vertex $i$ specifies the charge port and sharp head word. Every edge joining configurations $i,j$ has positive weight $t_{ij}$, equal to 4 for a permitted material hop or ½ for a ring flip. The native scalar port cost is $V_i=8(6-d_i)$. Take any finite-dimensional internal space, assign any unitary $U_{ij}$ to each edge and its inverse to the reverse edge. The covariant energy form on a normalized state is

$$
\langle\psi,H_U\psi\rangle=
\sum_{\{i,j\}}t_{ij}\|\psi_i-U_{ij}\psi_j\|^2
+\sum_i V_i\|\psi_i\|^2.
$$

Set $f_i=\|\psi_i\|$. The reverse triangle inequality gives $\|\psi_i-U_{ij}\psi_j\|^2\ge(f_i-f_j)^2$ for each edge separately, and $\sum_i f_i^2=1$. The charged ground energy is therefore bounded below by the **scalar** ground energy of the same graph. This bound holds for arbitrary holonomies and entanglement between the internal state and the charge or heads.

In the other direction, $\|\psi_i-U_{ij}\psi_j\|^2\le 2\|\psi_i\|^2+2\|\psi_j\|^2$. Every charged energy is at most $\max_i(V_i+2\sum_jt_{ij})$. The replay reproduces the earlier integer positive-vector lower certificates and the exact scalar row caps, including the witness hashes, denominators and strictly positive margins. The two estimates depend only on the undressed graph and native scalar weights. Restricting the charged trial states to a further internal symmetry subspace preserves both inequalities. [Certificate keys and exact cap calculation](measure.py)

| Charge sector | All absolute charged energies | Maximum/minimum bound | Koide upper bound |
| --- | --- | ---: | ---: |
| Vacancy | 4.343 < E ≤ 55 | 12.66406 | 0.4745850 |
| Excess | 5.450 < E ≤ 53 | 9.72478 | 0.4475344 |

## Releasing the neutral reference

The matched gateway’s neutral sector has only positive ring terms. A single sharp head word with any unit internal state has trial energy equal to ½ times its number of flippable faces. Its value is independent of every neutral unitary link, including assignments unrelated to the charged links. The least-cost audited trial has three flippable faces in the 191-word vacancy background and two in the 188-word excess background. A new constraint that projects out these neutral trials would require a new comparison. With the specified neutral states the dressed ground lies in $[0,1.5]$ or $[0,1]$, respectively.

Subtract the same neutral ground $v$ from all three energies in a sector. Given charged lower and upper bounds $L,U$, the largest possible gap ratio is bounded by $(U-v)/(L-v)$, increasing throughout the allowed interval of $v$. For positive gaps with largest/smallest ratio $R$, the Koide statistic is at most $(R+2)/(\sqrt R+2)^2$. [Proof of the triplet bound](../family_projector/README.md#whole-spectrum-koide-check)

| Charge sector | Matched neutral ground | Largest gap ratio | Koide upper bound on gaps |
| --- | --- | ---: | ---: |
| Vacancy | 0 ≤ v ≤ 1.5 | 18.81815 | 0.518250 |
| Excess | 0 ≤ v ≤ 1 | 11.68540 | 0.466141 |

Even pooling opposite signs with independent matched neutral backgrounds gives a common envelope from strictly above 2.843 through 55 native energy units. Multiplying all gaps by any common positive energy conversion preserves the ratio exclusion. It applies to energy gaps from the specified neutral ground. Splittings between two excited states can be arbitrarily different observables. [Matched trial head words and pooled bounds](results.json)

## Noncommuting SU(2) on real cycles

The first control assigns a local SU(2) frame $G_i$ to each sharp charge and field state and uses $U_{ij}=G_iG_j^*$. Exact integer matrix comparison in real spinor coordinates verifies that this full operator is a unitary conjugate of the scalar Hamiltonian tensored with the spinor identity. Every scalar eigenvalue then occurs with twice its original complex multiplicity.

The curved control chooses a breadth-first tree of the actual configuration graph. Tree links carry the identity. Every remaining material link carries $i\sigma_x$, and every remaining ring link carries $i\sigma_z$, using the inverse for reverse traversal. There are fundamental cycles of both kinds in both charge sectors. Referring them to the same root gives group commutator $-I$. The receipt records complete closed paths together with each charge port and head word.

| Charge | Actual matter and ring edges | Pure-gauge min, max | Curved SU(2) min, max |
| --- | --- | --- | --- |
| Vacancy | 8,266 and 17,234 | 4.343760, 40.929696 | 5.737574, 39.261269 |
| Excess | 2,565 and 4,845 | 5.450006, 40.632994 | 6.814622, 39.534511 |

Equality of the scalar lower bound would require the positive scalar ground-amplitude profile and a parallel internal unit vector on every edge. Such a vector would be fixed by every rooted holonomy. Here the commutator is $-I$, which has no fixed unit vector. The curved ground must therefore rise, in agreement with both computations.

These results cover internal unitary transports on the given configuration graphs with unchanged scalar port costs and positive edge weights. A coupling that changes contact strengths, supplies an internal onsite Hamiltonian, permits new material sectors or changes the physical graph defines a different spectral problem. [Frozen design](DESIGN.md) · [Executable and archived inputs](measure.py)

From this directory, run the repository interpreter with `-B -W error measure.py --output /tmp/covariant-field.json`.
