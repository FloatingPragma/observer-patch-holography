# Conditional particle roles from actual chiral source currents

This construction connects the declared left-Weyl source to its massive
background fields without adding an independent mirror family. The
[sitewise source encoding](../chiral_source/PROOF.md) supplies the actual
native currents on `F(C^2 tensor R15)`. The
[spectral regulator](../chiral_regulator/PROOF.md) supplies `h_L(p)=-sigma.p`
on its complete finite momentum menu, with a stated nonlocal position
kernel. The [Higgs construction](../chiral_higgs/PROOF.md) supplies both a
quantum scalar interaction and the separately specified neutral-background
benchmark. The claims below concern that benchmark. They do not replace
the quantum scalar with its expectation value inside an interacting vacuum.

## A canonical change of field variables

Use the fifteen native subsets and the exact symmetric Yukawa tensors in
`physical_attachment.chiral.representation`. Let `M=v sum_c y_c C_c,1`,
with c equal to up, down or charged-lepton and the neutral Higgs component
indexed by 1. For the exact calculations v and all y are positive real
parameters. The singular-value statement extends to complex coefficients.
The physical Higgs convention here is `H=(0,v)`; comparison with a convention
`H=(0,v_EW/sqrt(2))` requires `v=v_EW/sqrt(2)`.

In spin-major order the pair matrix is `Delta=epsilon tensor M`, with
`epsilon=[[0,1],[-1,0]]`. The Hamiltonian contains
`(1/2)c^T Delta c + h.c.`. For the Nambu vector `(c(p),c^dagger(-p))`,

`B(p) = [[h_L(p), -Delta*], [Delta, -h_L(-p)^T]]`.

A star in this formula means entrywise conjugation. The antisymmetry of
Delta makes B Hermitian. Its sixty rows duplicate a thirty-mode canonical
algebra; they are not sixty independent species.

Each nonzero neutral tensor pairs one weak-doublet index l with a singlet
conjugate index r. Write `C_lr=s_lr`, with `s_lr=+1` or `-1`. Define

`Psi_L,s(p)=c_l,s(p)`,

`Psi_R,s(p)=-s_lr sum_t epsilon_st c_r,t^dagger(-p)`.

The unpaired neutral doublet component is retained. In matrix form
`Psi=A c+B c^dagger`; the executable checks exactly

`AA^dagger+BB^dagger=I`, `AB^T+BA^T=0`,

and the inverse canonical identity. This is a bijective change of canonical
variables, not deletion of states. It gives seven four-component Dirac
fields, counting three colors of each quark type and one charged lepton,
plus one two-component neutral left field: `7*4+2=30`. All finite occupation
states remain. Particle-hole conjugation is also compatible with the
symmetric antiperiodic momentum menu used by the regulator.

With `S=[A B]`, the stronger reducing identity is executed:

`S B(p)=h_physical(p) S`.

Each paired block is exactly

`h_c(p)=[[-sigma.p, v*y_c I2], [v*y_c I2, +sigma.p]]`.

Consequently its two energy signs are `+/-sqrt(|p|^2+|v*y_c|^2)`.
The neutral block is `-sigma.p`. The finite antiperiodic menu excludes zero
momentum; `abs(v*y_c)` is the rest parameter of this background continuum
formula, rather than a detector record at an absent finite momentum.
A separate sixteen-dimensional four-mode
Fock calculation verifies the pair-Hamiltonian identity for both native
Yukawa signs, including occupied inputs and all CAR identities. These
checks detect a wrong pair sign that a list of squared masses would miss.

A massive Dirac field assembled from two different gauge representations
is compatible with chiral weak interactions after a neutral background is
specified. This is distinct from placing the same R15 representation on
both chiralities and adding a gauge-invariant beta coin.

## The physical currents are inherited, not relabelled

For every native Hermitian current T, its spin-lifted one-body matrix is
`G=I2 tensor T`. In the physical variables the number-preserving current is

`G_physical=A G A^dagger-B G^T B^dagger`.

The anomalous pair terms vanish for all twelve subgroup currents. The
normal-ordering constants `Tr(G B^T conjugate(B))` vanish for the complete
family and are explicitly checked. Electric charge `Q=T3+Y` has the same
value on both halves of each massive field: `2/3` for the up-type channels,
`-1/3` for down-type channels and `-1` for the charged-lepton channel.
The neutral field has charge zero. Color acts equally on the physical left
and right quark components. All three weak generators have support only on
the left components. The raising current connects each up/down color pair
and the neutral/charged-lepton pair, with unit matrix elements for both
spin components. Thus these conditional roles match more than a dimension
count: the source supplies their actual chiral current action.

The native subset map is computed directly from the tensors, with stable
`native-chiral-pair-*` IDs in the receipt. Its conjugate-right subsets
represent left-handed charge-conjugate fields. Antiparticle charges reverse
under the physical particle-hole interpretation. Color labels in this
finite field dictionary are representation components; they do not assert
that isolated colored asymptotic particles exist in a confining theory.
The Standard Model's left doublets, singlets and Higgs convention are given
in the [PDG electroweak review, section 10.1](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf).

## What this spectrum determines

| Conditional role | Number of Dirac fields | Charge | Rest energy in source units |
| --- | ---: | ---: | --- |
| Up-type quark | 3 color components | 2/3 | `abs(v*y_u)` |
| Down-type quark | 3 color components | -1/3 | `abs(v*y_d)` |
| Charged lepton | 1 | -1 | `abs(v*y_e)` |
| Neutral left lepton | one Weyl field | 0 | `0` |

This is one family. Nothing selects an electron rather than a muon or tau,
nor an up rather than charm or top generation. A probability of such an
identification needs a likelihood, prior or discriminating prediction;
none is supplied by this algebraic dictionary. The receipt therefore
reports conditional field roles and leaves measured species unset.

Changing `(y_u,y_d,y_e)` from `(1,1,1)` to `(2,3,5)` at v=1 preserves every
current and chirality identity while changing every corresponding mass
ratio. These exact controls exclude mass selection by these constraints
alone. The native dark-module idle generator is zero and cannot secretly
supply the three coefficients. An allowed energy conversion u_E gives
`m_c c^2=u_E abs(v*y_c) GeV`, with both the coefficient ancestry and the
energy standard retained. No experimental mass is consumed.

The one-family massless neutral sector cannot reproduce neutrino mass
splittings; additional flavor content and mass dynamics would be required.
No nonzero interacting self-energy, confining spectrum, physical Higgs
vacuum, gauge-boson pole or scattering prediction is established by the
neutral-background diagonalization. A common interacting state and
controlled gauge/Higgs continuum are separate physical requirements.
