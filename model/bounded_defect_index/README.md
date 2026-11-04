# Exact localized zeros of the bounded lattice hedgehog

The supplied bounded radial mass from
[topological_defect_spectrum](../topological_defect_spectrum/README.md)
has at least two exact normalizable zero modes on the original infinite
lattice, at every positive spacing. Their localization has an explicit
bound that remains finite in physical distance as the spacing decreases.
The proof uses the same source hopping and mass profile; it adds no
texture profile or envelope ansatz.

The operator on the complete two-spin lattice carrier is

```
(h0 psi)(n)=iv/(2a) sum_j sigma_j[psi(n+e_j)-psi(n-e_j)],
B_n=eta*m sum_j(-1)^n_j(n_j/|n|)sigma_j  for n!=0,
B_0=0,                 H=h0+B,
a,v,m>0,               eta=+1 or -1.
```

All eight folded valleys are retained. The site-centered mass field,
its amplitude and the source clock remain specified inputs. No claim
is made that the interacting source selects this texture.

## The index and the admitted perturbations

Spatial inversion `R_inv psi(n)=psi(-n)` anticommutes with H. Its
even-to-odd block has complex Fredholm index two. The proof continuously
reduces the hopping to zero while the essential gap stays open. Every
nonzero inversion pair then has an invertible mass block, while the
origin supplies two even kernel rows. The index therefore forces at
least two exact inversion-even zero modes of the original operator.

This grading is nonlocal spatial inversion. Its index and the parent's
continuum internal Clifford index have different orientation behavior:

| Texture orientation | Continuum internal index | Original-lattice inversion index |
| --- | --- | --- |
| Hedgehog | `+2` | `+2` |
| Antihedgehog | `-2` | `+2` |

The exact kernel may contain further modes; this package proves a lower
bound, not its exhaustion. The same lower bound survives any compact
bounded self-adjoint V satisfying `R_inv V R_inv=-V`, and also any
bounded self-adjoint V with that anticommutation and `||V||<m`.
These are precise Hamiltonian perturbation classes. Preservation by
the native repair instrument is not inferred.

A compact scalar potential at the fixed site removes both atomic
zero levels without closing the essential gap. A uniform scalar term
also moves every original zero eigenvector to its added energy. These
countercontrols show why the anticommutation hypothesis cannot be
replaced by an essential-gap assumption alone.

## Explicit localization

For the unperturbed nearest-neighbor source define

```
R=max(1,ceil(12v/(am))),      q=1+am/(12v).
```

Every normalized exact zero mode satisfies

```
sum_{|n|_infinity>=2R+L} |psi(n)|^2 <= q^(-2L),
integer L>=0.
```

For `am<=12v`, the right side is at most
`exp[-m(aL)/(12v)]`, and `2aR<=24v/m+2a`. Thus the stated physical
reference radius and probability decay length are bounded as `a->0`.
At `a=v=m=1`, the explicit values are `R=12` and `q=13/12`.

The same estimate applies to a self-adjoint perturbation whose two
matrix indices lie in a finite core cube, after enlarging R to contain
that cube. The anticommutation condition is needed to guarantee zero
modes under that perturbation. Localization is not extended to arbitrary
nonlocal compact perturbations.

[PROOF.md](PROOF.md) derives the essential spectrum, index, exterior
compressed-square gap, weighted resolvent and tail bound. The exterior
calculation squares the compressed operator and retains its boundary
terms. Finite tests independently compare the source matrix with the
parent's literal sampling action, retain the full inversion spaces,
vary the hopping, apply strong compact odd perturbations, check scalar
countercontrols, and evaluate actual zero-vector tails. Those controls
are not an infinite-volume extrapolation.

## Replay and physical scope

```sh
python -W error -m pytest -q model/bounded_defect_index
python -W error -m model.bounded_defect_index --write
python -W error -m model.bounded_defect_index --check
```

The cold receipt binds the proof, all owning Python files and the
complete parent source manifest. It verifies every inherited digest
without rerunning the parent's spectrum proof. Its exact checks concern
the original Clifford and edge identities, the atomic index and rational
localization bounds; independent finite numerical controls remain tests.

This is a statement about a supplied quadratic source operator. It does
not establish autonomous defect formation, a many-body interacting
particle spectrum, native repair protection, physical species or charge
assignments, calibrated masses, or Einstein dynamics. No completion gate
or physical comparison is accepted.
