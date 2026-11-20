# The flat charged band with variable seam couplings

The opposite-charge conversion model has a flat lowest charged band for any static positive rank-one
term per external seam, with the one-cell charged operator fixed. With both endpoint factors nonzero
on every seam, the exact multiplicities are 30, 120, and 480 on the 20-, 80-, and 320-cell carriers.
[Conversion model](README.md)

## Lowest-band dimension

Work in units of the source coefficient $J$. Write $Q_0=D+A$ for the unchanged charged
operator on the 12 ports of one complete cell. Its lowest eigenvalue is $\lambda_0=5-\sqrt{5}$ with
a three-dimensional eigenspace $E_0$. Total occupation is conserved. In the one-charge sector every
on-port conversion term vanishes, whatever its strength. The admissible modified
one-body quadratic form is

$$
\langle z,Qz\rangle = \sum_{c=1}^{C}\langle z_c,Q_0z_c\rangle
 + \sum_{e=(i,j)} |\alpha_e z_i+\beta_e z_j|^2,
\qquad \alpha_e,\beta_e\in\mathbb C.
$$

Positive strengths can be absorbed into the two coefficients. The uniform source uses $\alpha_e=\beta_e=1$.
The difference from $\lambda_0 I$ is a sum of positive semidefinite forms. Its kernel consists of the
vectors in $E_0^{\oplus C}$ obeying the external seam equations $\alpha_e z_i+\beta_e z_j=0$. There
are $3C/2$ external seams, so at least $3C/2$ independent lowest-band directions survive even if some
coefficients vanish. Two adjacent cells provide six amplitudes and meet at most five external seams.
A nonzero lowest-band vector supported on their union therefore exists for every member of this family.
[One-body source](source.py) · [Uniform-source certificate](subspaces.py)

The actual carriers make the count exact when both coefficients on each seam are nonzero. Consider a
linear dependence among the seam equations on $E_0^{\oplus C}$. On a cell with three independent
port evaluation rows, the three incident dependence coefficients vanish. At every other cell the
unique pair of opposite rows forces its two incident dependence coefficients to be proportional by a
nonzero factor, and the third vanishes. Make a graph with the seam equations as vertices and these
proportionalities as edges. The carrier check confirms that each component contains a coefficient
forced to zero. Zero then propagates through the entire component for any choice of the nonzero
endpoint factors, including complex phases. The seam equations are independent and the dimension is
exactly $3C-3C/2$. When strengths change with the endpoint ratio on each seam held fixed, the
subspace itself is unchanged.

## Observer consequence and controls

For every preparation entirely in this lowest band, free evolution is multiplication by
$e^{-iJ\lambda_0t}$. A fixed port or cell detector has the same outcome probabilities at every wait.
The result concerns this band of the static one-charge generator. Excited bands and interactions
that act during the wait are separate dynamical questions.

A first control on the same carrier replaces each external square by the rank-two endpoint cost
$|z_i|^2+|z_j|^2$. At the original energy $J\lambda_0$, every cell's three external
amplitudes must vanish. That leaves zero directions at an independent tripod and one at a dependent
tripod. The exact dimensions at that energy are 0, 24, and 104. A second control turns off one
external square entirely. By independence of the remaining seam equations, the original lowest band
gains exactly one direction, giving 31, 121, and 481.

## Consequence for opposite-charge conversion

There is also a large family of stationary nonconverting pairs. Take both charge signs to have the
same modified one-body operator, as in the original source, and let $E$ be its lowest band of dimension
$m$. The charges are distinguishable. Within their amplitude space, use the subspace
$\operatorname{Sym}^2 E$ of tensors invariant under interchange of their two addresses. Only
coincident amplitudes activate the on-port contact. Their values on any one cell factor
through $\operatorname{Sym}^2 E_0$, of dimension six. Thus the coincidence readout has rank at most
$6C$.

With nonzero factors on each seam, this bound improves to $6C-3C/2$. The one-charge seam equation
makes the two endpoint evaluations proportional on $E$, so their squared evaluations on a symmetric
tensor obey one external equality per seam. These equalities are independent on the direct sum of
the local symmetric tensor spaces. Three independent linear port forms have independent squares, and
at a dependent tripod the opposite rows square to the same form, independent of the third square.
The same component argument forces every dependency coefficient to zero. The coincidence map
therefore has a kernel of dimension at least $m(m+1)/2-9C/2$.

| Cells | Lowest one-charge dimension | Lower bound on address-symmetric nonconverting pair directions |
| ---: | ---: | ---: |
| 20 | 30 | 375 |
| 80 | 120 | 6,900 |
| 320 | 480 | 114,000 |

Every vector in this kernel has charged energy $2J\lambda_0$. Its coincident amplitudes and neutral
detector probabilities stay zero at all free waits. As a preparation control with the converter
enabled, choose a normalized $u\in E$ and prepare $u\otimes u$. Its coincident profile $u_i^2$ is
nonzero. The coefficient of $t^2$ in the probability of a neutral pair is $J^2\sum_i|u_i|^4>0$ at
the original uniform conversion strength. The same source and neutral-pair readout therefore
separate the two kinds of preparation. The dimensions above are lower bounds on an identified subspace of the
finite two-charge model. [Contact and detector definitions](README.md#local-law-and-source-map)

## Reproduction

The verifier checks the spectral polynomial and multiplicity of the one-cell operator, the carrier
structure, the local row relations, and the propagation of zero coefficients. It reduces $\sqrt{5}$
to 45 modulo 101 to certify nonzero minors and checks each opposite-row identity in integer
coefficients. The spectral and detector implications follow from the proofs above.
[Integer verifier](flat_band.py) · [Certificate test](test_flat_band.py)

```sh
python -B -m model.scalar_channel.radiative_repair.flat_band
python -B -m unittest model.scalar_channel.radiative_repair.test_flat_band
```
