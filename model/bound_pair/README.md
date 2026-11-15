# One bound excitation in the Weyl-contact action

CP-03/04/05/07, issues [#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139),
[#1140](https://github.com/FloatingPragma/observer-patch-holography/issues/1140),
[#1141](https://github.com/FloatingPragma/observer-patch-holography/issues/1141) and
[#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143).

This construction studies the **same supplied local action** as
[`interacting_field`](../interacting_field/README.md),
`H = dGamma(h) + U sum_x n_x,0 n_x,1`, with
`h(k) = -c sum_j sigma_j sin(k_j)`, `c=1/8`, `U=1`.
The periodic boxes used for momentum resolution change the boundary condition
of that package's open boxes; no equivalence of finite boundary observations
is assumed. Spatial lattice spacing, CAR dictionary and action are supplied.
The infinite-lattice integrals hold spacing fixed. They are not a spatial
continuum limit or a physical vacuum construction.

Repulsively bound lattice pairs are established physics; see
[Winkler et al., Nature 441, 853–856 (2006)](https://www.nature.com/articles/nature04918).
The contribution here is the explicit attachment, with error and scope checks,
to this repository's spin-dependent Weyl-contact action. No discovery of a
new elementary particle is claimed.

## Bound band and its reference

Work on the complete antisymmetric two-fermion Hilbert space. Put `T=dGamma(h)`
and let `P` project onto the on-site spin singlets. The Fourier symbol gives
`||T|| <= b = 2 sqrt(3) c`, and a single hop leaves the on-site-pair space,
so `PTP=0`. Min–max applied to that space and its orthogonal complement gives
exactly one upper-band eigenvalue per site, each at least `U`, while all
remaining eigenvalues are at most `b`. Thus `U>b` gives a uniform separating
gap at least `U-b`. For the declared parameters this exceeds `0.56698` in
model energy units. This is spectral protection under the closed action,
not autonomous error correction against arbitrary source operations.

At total momentum `K`, the free fiber is
`B_K(k) = h(k) tensor I + I tensor h(K-k)` and the contact vector is the
normalized constant spin singlet `d_K`. The interaction is `U|d_K><d_K|`.
For `E>b`, the scalar resolvent equation

```text
g(E,K) = average_k <s|(E-B_K(k))^-1|s>,     U g(E,K) = 1
Z(K) = 1 / [U^2 average_k <s|(E-B_K(k))^-2|s>]
psi_K = U sqrt(Z(K)) (E-T_K)^-1 d_K
```

has the unique isolated upper pole. `Z` is the contact spectral weight;
all other spectral weight remains a distinct output, rather than being
renormalized away. `fiber.py` evaluates the four spin-projector weights
and independently checks the resulting eigenstate against the full matrices.

At `K_eta=(pi,pi,pi)`, `h(K_eta-k)=h(k)`. The total spin kills a singlet,
so the pair has exactly `E=U`, `Z=1`. On every even torus of side at least
four, `r(k)=sqrt(2) sin(k_x) d(k)` is a normalized, orthogonal, physical
antisymmetric two-fermion state with `Tr=0` and `Pr=0`. Thus it has energy
zero. The two columns `(r,d)` span an exactly invariant `M_2` algebra with
`H=diag(0,U)`. Both states have particle number two and the same total
momentum. Their relative phase measures the energy difference `U`; no
vacuum–pair coherence or absolute vacuum energy is assumed.

[`particle_clock`](../particle_clock/README.md) attaches this finite algebra
to an explicitly heralded native source corner. The isometry defining that
attachment and the selected thermal reference are supplied. It is not an
implemented encoding of the entire spatial field into the 66-state source.

## Curvature and what its mass means

At `K_eta`, the band has isotropic curvature `kappa = d²E/dK_j²`.
Second-order perturbation of the full resolvent, with the singlet and triplet
blocks retained, gives

```text
kappa = (2 c²/U) average_k cos²(k_j)
        [1 + 4 c² (sum_l sin²(k_l) - sin²(k_j))
                 / (U² - 4 c² sum_l sin²(k_l))].
```

This formula is for the infinite Brillouin-zone integral. A quadrature is
reported only as a diagnostic. `mass.py` encloses it using exact rational
moments `m_n = binomial(2n,n)/4^n`. With `a=c²/U²`, expand the positive
denominator as a geometric series. The exact lower sum through order `N` is

```text
R_N = 1 + 16 a sum_(n=0)^N (4a)^n
             sum_(i+j+k=n) multinomial(n;i,j,k)
                 (m_i-m_(i+1)) m_(j+1) m_k,
R = U kappa / c²,
0 <= R-R_N <= 4a (12a)^(N+1)/(1-12a).
```

The tail bound uses `sum sin² <= 3` and
`average cos²(k_x)(sin²(k_y)+sin²(k_z))=1/2`.
It is valid for `12a<1`, which is the same gap regime. The zeroth
term and elementary denominator bound alone give `17/16 <= R <= 14/13` at the
declared parameters. The retained order-16 rational interval narrows this
to approximately `R=1.06904777858`, `kappa=0.0167038715404` and kinetic
coordinate mass `1/kappa=59.8663608`.

The value that the same `c` and energy difference `U` would assign through
a relativistic relation is `U/c²=64`. Their disagreement is a proved
property of this candidate, not a numerical fit. Under a common time
conversion `t_phys=tau t` and length conversion `x_phys=ell x`,
`E_phys=E/tau`, `kappa_phys=ell² kappa/tau` and
`c_phys=ell c/tau`; the ratio `E kappa/c²` is unchanged. Unit calibration
cannot remove this mismatch. The value is consequently neither an SI mass
prediction nor evidence for a relativistic particle of an identified species.

The interaction and kinetic energies also refer to this same normalized
state: `E_int=U Z(K)` and `E_kin=<psi_K|T_K|psi_K>`, with their sum equal
to the pole energy. Hellmann--Feynman gives `dE/dU=Z(K)` at fixed `c,K`.
At eta the complete energy difference is the contact cost `U`; away from
eta both terms contribute. Independent finite differences check this
response against the eigenstate and its contact residue. This derives the
response to the declared action coupling; it does not derive `U` from
native repair maintenance costs or count coherent apparatus work as energy.

Tests compare complete small-fiber diagonalization, normalized physical
exchange, the rational interval, independent quadrature and a finite
difference of the secular pole. The finite difference uses continuous `K`
as a quadrature parameter; it is not mislabeled as an untwisted finite-torus
momentum. Removing `U` destroys the isolated-pair premise.
