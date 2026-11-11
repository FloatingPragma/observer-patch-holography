# What a centre-class census determines

The [exact calculation](species_identifiability.py) solves an inverse problem
for the component menu in
[the five-level exterior code](../defect_particle_classes/PROOF.md).
Its [receipt](species_identifiability_results.json) contains every matrix
column, a rational nullspace basis, attained nonnegative bounds and
indistinguishable populations. The component dictionary is declared. No map
from native face-flux events to physical particle counts is supplied.

## The component maps

For a five-bit occupation mask, let `n3` and `n2` count its first three and
last two bits. The observed centre class is

```
k = (3 n2 - 2 n3) mod 6.
```

Index the masks by `j=0,...,31` and put `A[k,j]=1` if mask `j` has class `k`,
and zero otherwise. A nonnegative component population `x` produces the
six-vector `c=A x`. Each column has one unit entry, so total count is retained.
The matrix has six disjoint nonempty rows and rank six over the rationals.
For each class choose its first column `p(k)`. The vectors
`e_j-e_p(k)`, for every other column `j` in class `k`, form a kernel basis.
They lie in the kernel, are independent because each has its own nonpivot
coordinate, and span every vector whose sum within each class vanishes.

The full menu includes exterior degrees zero and five, labelled
`anti-nu_R` and `nu_R` in the declared dictionary. The restricted menu omits
exactly these two masks. This is the degree restriction of the charged
exterior grammar in [ExteriorSelection.lean](../../Lean/Screen/ExteriorSelection.lean),
whose premise-register identifier is PR-59. It is a 30-component menu before
the additional chiral/anomaly mask selection. Neither that restriction nor
its selection theorem proves completeness of the 32-component extension.

The boson extension adds two Higgs components and their two antiparticle
components, all of class three, and eight gluon, three weak gauge and one
abelian gauge internal components, all of class zero. These are internal
representation components. Their number is distinct from a physical
helicity count.

| Menu | Class multiplicities, `k=0,...,5` | Matrix size | Rank | Nullity |
| --- | --- | --- | --- | --- |
| Full exterior | `(4,6,6,4,6,6)` | `6 x 32` | 6 | 26 |
| Restricted exterior | `(2,6,6,4,6,6)` | `6 x 30` | 6 | 24 |
| Full exterior plus bosons | `(16,6,6,8,6,6)` | `6 x 48` | 6 | 42 |
| Restricted exterior plus bosons | `(14,6,6,8,6,6)` | `6 x 46` | 6 | 40 |

## Identifiability and sharp bounds

For any real coefficient vector `w`, the aggregate `w.x` is identifiable
from every admissible census precisely when `w` is constant within each
class. This is equivalent to `w` belonging to the row space of `A`. Necessity
follows from the two nonnegative unit populations `e_i,e_j` in a shared
class: they have equal censuses. Sufficiency follows by summing the common
class coefficients against `c`.

For a particular nonnegative census, define

```
L(c,w) = sum_k c[k] min_{j: k(j)=k} w[j],
U(c,w) = sum_k c[k] max_{j: k(j)=k} w[j].
```

Then `L <= w.x <= U` for every `x>=0` with `A x=c`. Both bounds are
attained by assigning each class's entire population to a minimizing or
maximizing component. Their convex combinations attain every intermediate
value. This also proves sharp integer-count bounds when `c` is integral.
At a fixed census, only classes with positive `c[k]` constrain
identifiability. No uniform within-class distribution is assumed.

Write `q=c1+c2+c4+c5`. In both fermion-only menus the aggregate labelled
quark is exactly `q`, and the aggregate labelled lepton is exactly `c0+c3`.
Up versus down is unconstrained between zero and `q`. In the full menu,
neutrino and charged-lepton counts each range from zero to `c0+c3`, and
the right-neutrino count ranges from zero to `c0`. In the restricted menu,
the neutrino interval is `[0,c3]` and the charged-lepton interval is
`[c0,c0+c3]`.

With the boson extension the quark-labelled count is exactly `q`,
while the lepton-labelled count ranges over `[0,c0+c3]`. The fermion count
ranges over `[q,q+c0+c3]`; the boson count over `[0,c0+c3]`. Gauge counts
range over `[0,c0]` and Higgs counts over `[0,c3]`. These bounds refer to
each aggregate separately; arbitrary combinations of their endpoints need
not coexist. For positive total count `N`, divide by `N` to obtain
fractions. If `q>0`, the quark fraction *among fermions* lies in
`[q/N,1]` in the boson extension. If `q=0`, that conditional fraction is
zero when the fermion population is positive and undefined when it vanishes.

Examples retained as full nonnegative vectors exchange `u_R` and `d_R`
inside class four, `e_R` and `nu_R` inside class zero, or `e_L` and `nu_L`
inside class three. In the extended menu, an `e_R` population and an
abelian gauge population have the same class-zero census while their
fermion fractions differ by one. An `e_L` population and a Higgs component
give the corresponding class-three lepton ambiguity.

Equal component weights give `24:8`, or `75%:25%`, in the full fermion
menu and `24:6`, or `80%:20%`, in the restricted menu. These are counting
controls. Different weights between classes change the quark/lepton
fraction over its entire unit interval. Multiplicity supplies neither
these weights nor a kinetic law. Conversely, a supplied code-quantum
census does identify its quark aggregate; within-class ambiguity does not
invalidate that conditional statement.

## The smallest additional label

For a specified collection of linear aggregates, associate to each component
the tuple of its aggregate coefficients. Within each class, two components
can share a deterministic readout label precisely when those tuples agree.
If class `k` contains `m_k` distinct tuples, any label sufficient for all the
aggregates needs at least `m_k` values: two nonnegative unit populations
with distinct tuples otherwise remain indistinguishable. Assigning the
tuple index within each class attains this bound. Because the original
class is retained, the minimum common alphabet has `max_k m_k` values.
The claim concerns deterministic component labels whose complete joint
class/label census is measured, without hidden priors or auxiliary probes.

For quark, lepton, Higgs and gauge totals in either boson extension, one
fermion/boson bit in classes zero and three suffices. Together with the
original class, this gives eight distinct population sums. The quark
classes need no additional distinction. For every aggregate in the
receipt, including up, down, charged lepton, neutrino and right neutrino,
three label values suffice and are necessary. Class three contains the
three distinct signatures represented by `e_L`, `nu_L` and `H`, so a
binary label cannot suffice. A fixed binary encoding therefore uses two
bits. The joint class/label matrices have rank fourteen for the full boson
extension and thirteen for the restricted extension. In either
fermion-only menu the full aggregate collection needs a binary label.
These are exact minimum information requirements for the declared
dictionary; no physical measurement implementing those labels is proved.

## Native events and an absent decoder

The [native population source](PROOF.md) counts two oriented flux faces
per triangular cell. A face of class zero is a zero-flux face. It is not
an observed electron, neutrino or gauge quantum. Face flux, code occupation
and thermal particle occupation are distinct sample spaces. The initial
zero-flux preparation alone cannot choose any of the different species
populations with the same class-zero code census.

Using the inverse problem on run counts requires a source-bound
preparation/readout instrument establishing what constitutes one counted
quantum, its component label and class, its detection efficiency and its
normalization. It must retain all outcomes, including zero-class events,
and preserve the relevant original probabilities. A rule that labels a
face by its class supplies no such instrument. The executable interface
rejects native flux-face records as inputs to species inference.

## A thermal comparison with declared assumptions

The comparison uses `T=1 TeV`, one equilibrium temperature, zero chemical
potentials and an electroweak-symmetric ideal gas. All masses, interaction
corrections and thermal masses are set to zero. It includes three Standard
Model generations with active Weyl neutrinos and their antiparticles,
without populated right-handed neutrinos. This temperature and thermal
preparation are comparison inputs. The native source supplies neither.

The equilibrium thermodynamic convention and the difference between the
ideal energy count and the interacting Standard Model equation of state
are described by [Saikawa and Shirai, JCAP 05 (2018) 035,
sections 3 and 3.4](https://arxiv.org/abs/1803.01038).
Their interacting calculation differs from the ideal `106.75` value even
above the electroweak crossover. The comparison here makes no precision
claim about that plasma and uses no fitted cosmological abundance.

One massless Weyl particle or antiparticle has one physical helicity.
Every one of the 30 restricted components contributes once per generation,
giving 90 fermionic helicity states. There is no extra factor of two from
the code's declared spin doublet. The twelve gauge internal components have
two transverse polarizations, and the complex Higgs doublet contributes
four scalar particle/antiparticle components: 28 bosonic states. Counting
both those scalar components and a further real-component factor would
double count the Higgs.

For zero chemical potential, expand the Bose and Fermi distribution
functions in exponentials. Integrating `x^r exp(-l x)` gives
`r!/l^(r+1)`. The alternating sum has ratio `1-2^(-r)` to the Bose sum,
hence `3/4` for number density (`r=2`) and `7/8` for energy density
(`r=3`). In units `hbar=c=k_B=1`,

```
n   = [g_B + (3/4) g_F] zeta(3) T^3 / pi^2,
rho = [g_B + (7/8) g_F] pi^2 T^4 / 30.
```

Thus the number coefficient is `191/2` and the energy coefficient is
`427/4 = 106.75`. Normalizing to the corresponding total gives:

| Aggregate | Number fraction | Energy fraction |
| --- | --- | --- |
| Quarks and antiquarks | `108/191` | `252/427` |
| Leptons and antileptons | `27/191` | `63/427` |
| Higgs components | `8/191` | `16/427` |
| Gauge polarizations | `48/191` | `96/427` |

The quark fraction among fermions is `4/5` under either weighting. The
normalized class vectors are `(57,27,27,26,27,27)/191` for number and
`(117,63,63,58,63,63)/427` for energy. Given an admitted number-count
decoder, the six class sums and their linear combinations could test the
number benchmark. An energy benchmark additionally needs energy-weighted
records and the same thermal preparation. A class-only census cannot
separately test the leptonic, Higgs and gauge shares. No time evolution,
expansion, freeze-out or simulated cosmological history follows from this
comparison.

## Reproduction

```
python3 -m model.native_abundance_epochs.species_identifiability --write
python3 -m model.native_abundance_epochs.species_identifiability --check
python3 -m pytest -q model/native_abundance_epochs model/defect_particle_classes
```

The independent tests use rational elimination, direct bit enumeration,
attained extremizers and adversarial domain inputs. They preserve the
separate attempted-event population source and its scientific receipt.
