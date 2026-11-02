# Actual spatial rotation channels of the finite particle spectrum

This package attaches the **actual binary cubic rotation representations**
to every local charged and neutral spectral line certified in FA85 and FA88.
It uses the unchanged centered open cube, `c=1/8`, `U=1`, and its unique
`N=8` ground. It does not refit energies or identify measured species.
The calculation resolves an essential particle-dictionary ambiguity:
the source's hidden internal SU(2) is not its spatial spin representation.

## The spatial action is read from the source

Write a cube vertex as `x in {0,1}^3`. For each proper cubic rotation `R`,
rotate about the center, `p_R(x)=[R(2x-1)+1]/2`, and choose either SU(2)
lift `U(q)=w I-i(x sigma_x+y sigma_y+z sigma_z)`. The actual one-particle
action in the original Weyl frame is `D(q)=P_R tensor U(q)`.
The link `i c sigma_j/2` rotates covariantly, and the contact term
`U n_up n_down` and centering `-UN/2` are onsite spin scalars. Thus its
full CAR lift commutes with the complete interacting Hamiltonian. At the
central group element the lift is `Gamma(-I)=(-1)^N`; no number sector
is discarded in this assertion.

The 48 unit quaternions are the eight signed coordinate units, the sixteen
`(+-1,+-1,+-1,+-1)/2`, and the twenty-four vectors with two nonzero
coordinates `+-1/sqrt(2)`. The code verifies closure, the spatial product
law, the source onsite frames and signed hopping covariance exactly in
`Q(sqrt(2))`. No rounded matrix decides a representation multiplicity.

In FA85's hidden frame `d_x=R_x c_x`, where
`R_x=(i sigma_x)^x1(i sigma_y)^x2(i sigma_z)^x3`, this becomes

```
R_(p(x)) U(q) R_x^* = s_q(x) W_q,
D_hidden(q) = S_q tensor W_q.
```

The signs `s_q(0)=1` agree with the exact signed automorphisms used in
the spectral certificate. The unique ground is invariant under those
signed automorphisms (FA85's trivial ground block) and is a hidden-spin
singlet (FA88's complete spectrum). Its spatial representation is A1.
The factors `S_q` and `W_q` may individually be projective; their product
is the actual representation. Their adjoint action is unambiguous.

## A complete exact character calculation

Let `a(q)` be the sign of the permutation of the four unoriented body
diagonals, `t=tr R`, and `w=q_0`. The characters, in order, are

```
A1=1, A2=a, E=(t^2+tr R^2)/2-1-a*t, T1=t, T2=a*t,
G1=2w, G2=2aw, H=8w^3-4w.
```

These are characters of concrete representations: scalars, the diagonal
permutation sign, diagonal traceless quadratics, vectors, off-diagonal
quadratics, the spinor, its A2 twist, and the symmetric cube of the spinor.
Exact character inner products equal the identity and the squared
dimensions sum to 48, proving irreducibility and completeness. Restriction
of any SU(2) spin follows from `chi_(n+1)=2w chi_n-chi_(n-1)`.
For example spin 1 restricts to T1, spin 3/2 to H, spin 2 to E+T2,
and spin 7/2 to G1+G2+H. The converse is not unique: a cubic label is
not by itself a continuum spin assignment. This standard use of cubic
representations in lattice spectroscopy is discussed by
[Basak et al., Phys. Rev. D 72, 094506 (2005)](https://doi.org/10.1103/PhysRevD.72.094506);
the finite calculation here is self-contained.

## All 228 charged lines and explicit symmetry filters

Let A be the signed cube adjacency in FA85, with `A^2=3I`. The commutant
of all 48 signed cube automorphisms is exactly `span{I,A}`: the integer
commutator system has rank 62 modulo 1009, while the displayed independent
commuting matrices give the matching upper bound. A modular nonzero minor
is a rational nonzero minor. The same conclusion therefore holds over C.
For a charged line, covariance, hidden-spin invariance and ground
invariance give the full residue `B tensor I_2`. Its spatial rank is four
by FA85. Positivity and the common diagonal weight w then force

```
B = w (I + epsilon A/sqrt(3)), epsilon in {-1,+1}.
```

The certified off-diagonal sign determines epsilon: 114 lines use each
sign. For every proper rotation `tr(A S_q) tr(W_q)=0` exactly. Therefore
either band has half the character of the full sixteen-dimensional
one-particle representation, namely **G1+G2+H**, of dimension eight.
Removal gives the dual representation, with the same real characters;
addition gives the same result. These labels apply to every visible
charged line, not to an unsupported description of the whole eigenspace.

The emitted Gaussian-rational matrices are the explicit central filters
`Q_alpha=(d_alpha/48) sum_q chi_alpha(q) D_hidden(q)`.
The code verifies Hermiticity, idempotence, pairwise orthogonality,
completeness, and commutation with `A tensor I_2` by integer products.
They specify coherent symmetry-filtered CAR probes. Since the residue is
scalar within each occupied band, the fractions of the **complete probe
trace**, summed over all sites and both spin components, are 1/4, 1/4,
and 1/2 for G1, G2, and H. These are not probabilities for an arbitrary
single-site detector. No native execution of those filtered probes is
claimed here.

## All 135 spin and 135 density lines

For neutral bilinears the signs cancel. The site action is ordinary cube
permutation. Its commutant has dimension four (integer commutator rank
60 modulo 1009), with orthogonal Walsh projectors

```
P_l(x,y)=(1/8) sum_(|k|=l) (-1)^(k.x+k.y), l=0,1,2,3.
```

Their ranks are 1,3,3,1 and representations A1,T1,T2,A2. The uniform
spin and density probes annihilate this spin-singlet, fixed-N ground,
so `P_0` has no positive-energy residue. Each certified spatial residue
has rank one or three; positivity in this commutant forces exactly one
remaining projector, `B=(8w/rank(P_l))P_l`. Rank one selects l=3. At rank
three, the certified sign of `B(000,001)` distinguishes l=1 from l=2.
Every retained residue interval is checked against the resulting exact
matrix. Both channels have 57 degree-one, 54 degree-two and 24 degree-three
positive lines. Each retains its one exactly dark candidate separately.

The three hidden-spin components transform by `Ad W_q`, whose exact
character is **A2+E**, not T1. Density is a spin scalar. Consequently:

| Walsh degree | Density | All three hidden-spin components |
| --- | --- | --- |
| 1 | T1 | T1 + 2 T2 |
| 2 | T2 | 2 T1 + T2 |
| 3 | A2 | A1 + E |

The lowest full-Fock excitation, at source-clock energy
`[0.01502125417612,0.01502125417614]`, is thus **A1+E** under actual
spatial rotations. The complete gap certificate shows that these three
visible states exhaust that first eigenspace. It is not a spatial vector
merely because its hidden spin is one. The first density line at
`[0.81180294910720,0.81180294910723]` is T1. The first charged line remains
`[0.36168271180597,0.36168271180599]`, with G1+G2+H visible channels.
The machine-readable output retains every line's original energy enclosure,
weight enclosure, representation and complete-probe trace fractions.

## Scope, reproducibility and physical identification

All energies remain source-clock gaps of one finite open cube. Cubic
rotations supply genuine spatial quantum numbers, but neither boost
symmetry nor relativistic dispersion, continuum spin, mass calibration,
gauge charge, a stable asymptotic particle or an empirical species follows.
The finite result constrains such identifications and provides concrete
symmetry-filtered observables for testing them. It closes a finite
dictionary input, not the physical particle identification requirement.

Run `python -W error -m pytest -q model/cubic_particle_channels`, then
`python -m model.cubic_particle_channels --write` and `--check`.
The cold replay recomputes the group, commutants, central projectors and
classification. FA85/88 supply their already-certified spectral intervals;
their expensive eigensolver witnesses are not silently rerun or replaced.
Direct source hashes and inherited custody are reported separately.
