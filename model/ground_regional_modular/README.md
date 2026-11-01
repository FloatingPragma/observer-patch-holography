# Regional entropy and modular flow of the actual particle ground

This constructs the reduced density and modular generator of the invariant
ground `P_ground/5` already used in the
[full ground-family spectrum](../spin_resolved_spectrum/README.md). The
reference is the original `2 x 2 x 2` Weyl/contact cube at `c=1/8,U=1`,
not an independently chosen Gibbs state. Exact signed CAR traces retain
all regional occupation states, including the zero eigenspace. The result
supplies a finite same-reference entropy first law and identifies its
central correction. It does not identify a geometric cap or discharge
the common Einstein premises.

## A formula for every spatial region

In the source's hidden spin frame, the highest ground vector is the
spin-up Slater determinant filling the four negative spatial orbitals.
Write its spatial projector as `P_-=(I-A/sqrt(3))/2`. The entire ground
space is its spin-two multiplet. Averaging the highest vector over common
spin rotations therefore gives exactly `P_ground/5`.

This average can be established without assuming a spin integration
formula. A common rotated spinor has up/down squared amplitudes `t,1-t`
and a relative phase. Uniform phase integration deletes unequal down-spin
counts. In the symmetric k-spin subspace, the remaining diagonal weight
at down count l is

```
binomial(k,l) integral_0^1 t^(k-l) (1-t)^l dt = 1/(k+1).
```

The beta integral follows by repeated integration by parts. The result
is the identity on the symmetric spin-k/2 space divided by `k+1`.

Let `C_R` be the compression of `P_-` to any set R of whole spatial
sites, with natural-orbital eigenvalues `lambda_j` in `[0,1]`. A Slater
Schmidt decomposition across R and its complement writes the polarized
regional density as independent occupied/empty natural orbitals. For an
occupied subset I it has weight

```
w_I = product_(j in I) lambda_j product_(j not in I) (1-lambda_j).
```

One way to derive this decomposition is to singular-value decompose the
restriction of the occupied one-particle orbitals to R. The resulting
orthogonal occupied orbitals have the form
`sqrt(lambda_j) a_j^* + sqrt(1-lambda_j) b_j^*`; expansion of their Slater
product and the exterior trace give the displayed weights. Eigenvalues
zero or one fix the corresponding occupancy. Empty regional natural
orbitals are added when R is larger than the occupied one-particle space.

Partial trace commutes with the common spin rotation. Consequently, in
each I block the exact mixed-ground density is `w_I P_sym,I/(|I|+1)`.
Its nonzero eigenvalue is `w_I/(|I|+1)` with multiplicity `|I|+1`.
Every other spin/occupation state has eigenvalue zero. No probability is
discarded or conditionally renormalized. In particular

```
S(rho_R) = sum_j h_2(lambda_j) + sum_I w_I log(|I|+1).
```

All logarithms are natural. When all natural occupations are strictly
between zero and one, the modular generator on the positive support is

```
K_R = -sum_j log(1-lambda_j) I
      + sum_j log((1-lambda_j)/lambda_j) N_j + log(N_R+1).
```

The occupation operators here act on the support with at most one fermion
per natural orbital and symmetric occupied spins. This is not an identity
on the entire regional CAR Fock space. Occupations fixed by endpoint
eigenvalues are omitted from the logarithmic expression and handled by
their fixed blocks. The spectral formula itself covers those endpoints.
For the entire cube it gives the five eigenvalues `1/5`, as required.

## Exact one-site and two-site results

Every site has density `diag(1/2,1/4,1/4,0)` in vacuum/up/down/pair order.
Every pair of distinct sites has rank eight in its full sixteen-state
Fock space. The positive support consists of the vacuum, all four
one-particle states, and the three two-particle spin-triplet states.

For adjacent sites, `C_R` has eigenvalues `(1+1/sqrt(3))/2` and
`(1-1/sqrt(3))/2`. The density spectrum is

| Eigenvalue | Multiplicity |
| --- | ---: |
| `1/6` | 1 |
| `(2-sqrt(3))/12` | 2 |
| `(2+sqrt(3))/12` | 2 |
| `1/18` | 3 |
| `0` | 8 |

Let `H_R` be the actual internal regional Weyl/contact Hamiltonian. The
contact term annihilates this support. Its hopping has supported
one-particle eigenvalues `+/-1/16` and vanishes on the supported vacuum
and triplet. Hence the exact modular identity is

```
K_R = log(6) I + beta_R H_R + log(N_R+1),
beta_R = 16 log(2+sqrt(3)).
```

The code checks every one of the cube's twelve signed bonds directly
against the five exact source ground vectors. This coefficient is a
regional modular/energy conversion for this reference, not a physical
temperature or an independently calibrated time unit. For nonadjacent
sites the spectrum is `1/4` once, `1/8` four times, `1/12` three times,
and zero eight times; `K_R=log(4)I+log(N_R+1)` on the support.
All sixteen nonadjacent pairs are also checked, including signed
nonconsecutive CAR restrictions.

The hidden frame is related to the original local spin frame by the
source's sitewise spin rotations. Conjugating the regional density,
Hamiltonian and support by that same local unitary preserves these
identities and spectra. It does not identify hidden spin with Lorentz spin.

## The same-ground first law, its center and its domain

For every Hermitian trace-zero variation delta supported in the positive
regional support, functional calculus gives

```
DS_rho_R(delta) = Tr(delta K_R)
               = beta_R Tr(delta H_R) + Tr(delta log(N_R+1))
```

on an adjacent pair. Commutation of delta with rho is not required.
On the supported number-preserving observable algebra
`B(H_0) direct_sum B(H_1) direct_sum B(H_triplet)`, `log(N_R+1)` is central.
It is not central in the full supported matrix algebra or the full CAR
algebra. This identifies an actual central finite correction rather than
assigning the whole logarithm to bare regional energy.

More explicitly, write a number-block state as `rho=direct_sum p_n rho_n`.
Its entropy is `H(p)+sum_n p_n S(rho_n)`. At the reference, the triplet
block is the maximally mixed state, and the spin multiplicity contribution
is exactly the `log(n+1)` term already displayed. A formal definition of
a bulk pairing by dividing `beta_R Tr(delta H_R)` by `2*pi` still does not
prove that this pairing is a geometric diamond charge. The existing
Einstein theorem stack additionally needs its actual geometric kernel,
area/edge interpretation, common stress, scaling and normalization.

The smallest positive eigenvalue is `(2-sqrt(3))/12 > 1/48`; the strict
bound follows from `sqrt(3)<7/4`. If a density sigma has the same support,
delta=sigma-rho and `0<=epsilon<1`, the positive path remains at least
`(1-epsilon)/48` on that support. The logarithm divided-difference bound
then gives the explicit relative-entropy remainder

```
0 <= D(rho+epsilon delta || rho)
   <= 24 epsilon^2 ||delta||_HS^2/(1-epsilon).
```

The logarithm defines supported modular automorphisms
`sigma_t(X)=rho_R^(it) X rho_R^(-it)=exp(-it K_R) X exp(it K_R)`.
Tensoring the supported operator identity with a reference identity
preserves it; no joint entropy first law for an arbitrary
correlated extension is inferred. The construction does not define a
faithful modular group on all sixteen Fock states.

The kernel cannot be ignored for arbitrary detector disturbances. If
sigma is supported orthogonally to rho, then

```
S((1-epsilon)rho+epsilon sigma)
 = (1-epsilon)S(rho)+epsilon S(sigma)+h_2(epsilon).
```

The entropy derivative at zero diverges as `-log(epsilon)`. Thus the
supported first law cannot be extended to such a variation by deleting
the zero eigenvalues. This explicit countercontrol is why the result is a
finite same-ground theorem input, not a completed G1-G6 attachment.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.ground_regional_modular --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/ground_regional_modular
```

The exact matrices use rational coefficients of `1,sqrt(3)`. Separate
floating full-Fock restrictions, logarithms and entropy remainders test
the implementation; the displayed derivations and exact source traces
carry the claims. No new Lean proof, physical particle assignment,
Lorentz metric, Newton constant or Einstein-register inhabitant is asserted.
