# Topological mass textures and their spectra

This package determines continuum and original-lattice fermion spectra for
explicitly supplied native mass textures. It retains all eight Weyl valleys
and all sixteen folded rows from
[native_mass_channels](../native_mass_channels/README.md). A separate finite
theorem retains the original contact interaction and controls the splitting
of the two localized core orbitals, including their filled negative-energy
sea and all virtual states.

The mass field, its amplitude and its profile are premises. The results
identify modes trapped by that field. They do not derive the field's
formation, its own rest energy or an observed particle assignment.

## Conditional spectra

Here `v>0` is the continuum velocity in the source clock, `m>0` a supplied
asymptotic spin-mass magnitude, `kappa>0` a supplied gradient, and
`r=|m_vector|` for a uniform spin mass. Each signed multiplicity applies
separately to the positive and negative energy shown.

| Supplied quadratic field | Energy spectrum | Multiplicity and domain |
| --- | --- | --- |
| Uniform scalar and spin masses | `E=+/-sqrt(v^2 p^2+(m0+r)^2)` and `+/-sqrt(v^2 p^2+(m0-r)^2)` | Two four-component Dirac copies per gap branch; coincident energies add. |
| Uniform `m0=r>0` | One zero-gap branch and one gap `2r` | Two massless Dirac copies and two massive Dirac copies on the same quadratic source family and clock. |
| Bounded radial hedgehog `m_vector(x)=m x/|x|` | Two zero modes; `E=+/-m sqrt(1-1/n^2)`, `n>=2`; essential spectrum outside `(-m,m)` | Every discrete subgap level is classified; each signed excited level has multiplicity `2n^2`. The core is singular at one measure-zero point. |
| Smooth radial core `m x/sqrt(|x|^2+ell^2)`, `ell>0` | Exactly two zero modes; essential thresholds `+/-m` | The excited Coulomb levels above are not transferred to this profile. |
| Original lattice bounded hedgehog | Each fixed continuum Coulomb gap level is approached as `a->0` by at least its continuum multiplicity of discrete lattice levels | All eight valleys are sampled into the original two-spin lattice. Extra levels and exact bounded-lattice zero modes are not classified. |
| Linear continuum texture `kappa x` | Two zero modes; `E=+/-sqrt(2v kappa q)`, `q>=1` | Each signed level has multiplicity `4q^2+2`; this texture is unbounded and confining. |
| Original lattice linear texture | Two exact Bessel zero modes; `E=+/-sqrt(b_i^2+b_j^2+b_k^2)` | The `b_i` are positive simple eigenvalues of the exact one-dimensional Jacobi operator, with `b_0=0`. Every coordinate, sign and coincident sum is retained. |

The uniform tuning `m0=r` closes the bulk gap. The hedgehog index theorem
uses `m0=0`, `m>0` and an open essential gap, so these results do not
construct one regime with both massless bulk fields and the stated
protected localized zero modes.

The bounded degree `+1` hedgehog has chiral index `+2`; reversing its mass
direction gives index `-2`. The grading and its sign are calculated from
the original folded matrices. The protected quantity is this chiral index
under continuous chiral Fredholm deformations. Preservation by arbitrary
native repair branches is not assumed.

For the bounded constant-magnitude radial profile, the first positive
energies in units of `m` are `sqrt(3)/2`, `2sqrt(2)/3`, `sqrt(15)/4`,
`2sqrt(6)/5` and `sqrt(35)/6`, with multiplicities 8, 18, 32, 50 and 72.
The profile fixes these energies; topology fixes neither their values nor
a physical mass unit.

The full continuum and lattice proofs are in [PROOF.md](PROOF.md).
For the bounded hedgehog, all-valley sampling has second-order kinetic
consistency on smooth compact fields away from the core, and its finite
Gram matrices converge to their continuum values. The lattice's exact
essential gap and a graph-core approximation then give the spectral
inclusion in the table. This result supplies no global spectral convergence
rate or exclusion of additional lattice levels.

The six native lattice levels at `v=kappa=a=1` have exact rational
enclosures:

| Positive one-dimensional level | Certified interval |
| --- | --- |
| `b_1` | `[1.283497330929, 1.283497330930]` |
| `b_2` | `[2.130865638776, 2.130865638777]` |
| `b_3` | `[3.085057131288, 3.085057131289]` |
| `b_4` | `[4.063229009991, 4.063229009992]` |
| `b_5` | `[5.050373826697, 5.050373826698]` |
| `b_6` | `[6.041883198534, 6.041883198535]` |

These are enclosures for the infinite Jacobi operator. Exact residuals
retain both omitted boundary rows; a tail-resolvent bound and exact Sturm
counts certify that each interval contains the stated level. The stored
vectors are untrusted proposals. A finite eigensolver value is not used
as a certificate. No finite-spacing equality with the continuum oscillator
levels or uniform continuum spectral convergence rate is asserted.

## A controlled interacting finite core

The 27-site cube `{-1,0,1}^3`, with `v=kappa=a=1`, has two zero orbitals,
26 negative levels and 26 positive levels. Filling the negative sea gives
a four-state low subspace inside the complete `2^54` Fock carrier. For
the unchanged centered contact with coupling `U>0`, its exact compression
is

```
C+(U/8)[n1*n2-(n1+n2)/2].
```

Let `w=27U/4` and `epsilon=w^2/(6/5-2w)`. When `4w<6/5` and
`U/16>2epsilon`, the full interacting ground is an exact Kramers doublet.
The two particle-hole-related charged levels have excitation gap

```
U/16-2epsilon <= gap <= U/16+2epsilon.
```

At `U=1/65536`, a conservative decimal enclosure is
`[9.3599070381e-7, 9.7135792901e-7]` in source energy units. The
[finite core proof](CORE_PROOF.md) bounds the entire Fock complement;
it does not execute a `2^54` matrix. Independent numerical controls retain
all 54 one-body cube rows and all 64 Fock states of a three-site chain.
The finite cube and continuum radial hedgehog have different profiles and
domains, so their excitation energies are not identified.

## Replay and scope

```sh
python -W error -m pytest -q model/topological_defect_spectrum
python -W error -m model.topological_defect_spectrum --write
python -W error -m model.topological_defect_spectrum --check
```

The receipt uses a fresh Python process and binds both proofs, every
owning Python file, the rational spectral witness and the complete
inherited native mass source manifest. It verifies every inherited digest
without rerunning the parent's expensive centered-ground proof. Tests
compare independent Jordan-Wigner and differential realizations, Coulomb
radial equations, original lattice links, full finite Fock controls and
negative spectral-certificate cases. Sampling tests also compare literal
Weyl hopping and staggered masses against all sixteen continuum rows at
three spacings, retaining the alternating cross-valley Gram terms.

The unresolved physical premises are source selection and preparation
of the texture, calibrated parameters, autonomous defect formation and
rest energy, repair protection, physical charges/species, and a controlled
interacting relativistic continuum with a common Einstein stress tensor.
The receipt accepts no completion gate and adds no physical comparison.
