# Charged composite states invisible to single-fermion probes

A local three-fermion probe reaches **138 certified charged spectral lines**
outside the centered cube's entire linear-fermion response. Every weight is
strictly positive. The smallest gap is enclosed by
`[0.40246851577010,0.40246851577013]` in the source's energy unit. Direct
source calculations independently check the total response weight and its
energy moment. These finite excitations have hidden spin three halves;
physical particle identity requires a moving, localized state and a common
physical clock.

The complete [line table](lines.csv) and [receipt](results.json) retain every
candidate in the specified reducing block. No weight cutoff or empirical
mass input selects a line.

## The source and the composite probe

The source is the open eight-site Weyl/contact cube with its exact hidden-spin
frame and centered generator

```
K = dGamma(A tensor I_2)/16 + sum_x n_xup n_xdown - N/2,
A=A^T, A^2=3I,             c=1/8, U=1.
```

The reference is the unique `N=8`, hidden-spin-zero, eta-zero ground. The
[centered-spectrum construction](../centered_particle_spectrum/README.md)
supplies its integer ground witness, the exact signed spatial symmetries,
and the ground-selection theorem. Its imported Lieb uniqueness theorem
and full-number-sector bound are explicit mathematical inputs. This
packet checks the source hypotheses and sector-bound arithmetic and
reverifies the integer ground witness. It reports inherited receipt
custody separately; it does not rerun the earlier all-number calculation.

At vertex zero the native neighbors are vertices `1,2,4`. Define

```
C_0 = (S_1^+ + S_2^+ + S_4^+) c_0down,
S_y^+ = c_yup^* c_ydown,                ||C_0|| <= 3.
```

This is a cubic, number-changing observable of the same CAR algebra.
On the singlet reference, its output has `Nup=5,Ndown=2`, spin weight
`3/2` and eta weight `-1/2`. The exact total `S+` and `eta-` annihilate
the output. Therefore it has total hidden spin `3/2` and eta `1/2`.
The code verifies both annihilation identities on the entire exact
ground reducing space and verifies that `K` intertwines those ladders.
The local linear CAR probes of the same singlet have total hidden spin
`1/2`. Their complete dynamical cyclic spaces are orthogonal to this
response because the generator commutes with hidden spin. Degeneracy of
energies would not invalidate that orthogonality. Hidden spin is an exact
internal symmetry here; it is not identified with Lorentz spin.

The origin's six spatial permutations preserve the sum over neighbors and
the CAR operator at zero. Their signed Fock lifts give a complete
262-dimensional invariant orbit space in the 1568-dimensional `(5,2)`
occupation sector. In that orbit space the exact kernel projectors of
`S-S+` and `eta+eta-` use the respective nonzero roots `{5,12}` and
`{3,8}`. Their product is a weighted self-adjoint projector of rank
**138**. Polynomial idempotence, exact trace, source commutation and a
rank witness modulo 1009 prove complete coverage of this reducing space.
No approximate eigensolver selects the subspace. The actual probe image
of all 32 ground-block columns lies in it, so its omitted spectral weight
is exactly zero.

## Exact finite spectrum and retained weights

The proposal generator computes a full generalized eigensystem and rounds
its coefficients and energies to integers with denominators `2^60` and
`2^56`. The cold verifier uses only those integers. If `Y` denotes its
full orbit-coordinate columns, `D` the diagonal energy proposal and
`R=KY-YD`, exact integer products and outward rational square roots give

```
delta >= ||Y*Y-I|| < 1,
epsilon >= ||R||/sqrt(1-delta).
```

The basis is complete and invertible. Its similarity transform is
`Y^-1 K Y = D+Y^-1 R`, whose perturbation norm is at most epsilon.
Resolvent bounds and continuity along `D+tY^-1R` place exactly one
eigenvalue in each disjoint radius-epsilon disc. All 138 discs are
disjoint. For a normalized column with residual bound `r_i` and
distance `gamma_i` from the other discs, a phase of its exact eigenvector
differs by at most `2r_i/gamma_i`. This proof uses neither an ordered
Hermitian perturbation inequality for a nonunitary similarity nor a
floating residual as a bound.

The transition overlaps are integer products of the retained ground and
charged columns with the actual cubic CAR map. Their normalization gets
outward rational enclosures. If the ground and charged vector errors
are `eta_g` and `eta_i`, the amplitude error is at most
`3(eta_g+eta_i)`, using the independent operator bound on `C_0`.
Squaring the resulting nonnegative intervals gives every residue.
Every lower bound is positive, including the weakest line, whose weight
is approximately `3.84e-13`. All 138 states are visible to this probe.

| Quantity | Certified interval in source energy units |
| --- | --- |
| Smallest excitation | `[0.40246851577010,0.40246851577013]` |
| Largest excitation | `[2.97475650816793,2.97475650816795]` |
| Lowest-line weight | `[0.18203143278621,0.18203143278821]` |
| Total weight | `[0.74354609240642,0.74354609241091]` |
| First energy moment | `[0.39064379045340,0.39064379049794]` |
| Mean energy conditional on this response | `[0.52537938729092,0.52537938735399]` |

These are removal-sector excitations relative to the selected ground.
The local highest-weight probe does not count full spatial/spin
multiplet degeneracies. The calculation does not exhaust every charged
state invisible to linear probes, and its 138 poles are not 138 species.

## Independent source energy readout

The two spectral sums obey

```
sum_i w_i = <G,C_0^* C_0 G>,
sum_i w_i Delta_i = <G,C_0^* [K,C_0] G>.
```

An independent calculation acts with the original sparse integer source
and cubic CAR matrices on all 4900 ground-sector and 1568 charged-sector
configuration rows. It uses no composite eigenvalue or eigenvector.
The exact rational expectation on the normalized ground proposal is
enclosed using its certified vector error. The two source norm bounds,
from exact absolute row sums, are `11/2` and `71/16`; the observable
bounds are `9` for `C_0^*C_0` and `9(11/2+71/16)` for
`C_0^*[K,C_0]`. For normalized vector error eta, each expectation changes
by at most twice its operator norm times eta. The spectral and direct
intervals overlap and are intersected in the displayed values. Division
by the separately certified positive total weight gives the conditional
mean; no spectral line is dropped or renormalized out of the response.

The contraction `C_0/3` has the complete mathematical instrument with
complement `(I-C_0^*C_0/9)^(1/2)`. Its success probability is enclosed by
`[0.08261623248960,0.08261623249011]`; the complementary outcome is
retained. This establishes a finite observable and its source spectral
law. It is not a compiled or executed native preparation and detector.
The clock is the supplied generator `K` with `hbar_model=1`, with no
conversion to seconds, electronvolts or a relativistic rest mass.

The finite source, repulsion, centering and probe are specified without
fitting observed masses. Autonomous defect protection, motion, physical
charge and spin, continuum pole limits and the shared gravity conditions
are separate requirements on a common realization.

## Reproduce

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/charged_composite_search
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.charged_composite_search --check
```

`--propose` creates a candidate integer witness and admits it only after
the exact certificate passes. `--write` produces the receipt and complete
line table in a fresh process. `--check` uses the retained integers and
replays the source, coverage, residues and independent energy identities.
