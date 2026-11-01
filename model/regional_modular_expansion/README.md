# A controlled interacting regional modular generator

This CP-09 / [#1145](https://github.com/FloatingPragma/observer-patch-holography/issues/1145)
component constructs the regional modular Hamiltonian of the **same**
complete-Fock Weyl/contact action used by the particle calculations, with
an exact polynomial and an outward rational error. It turns the existing
numerical regional mismatch into a controlled positive approximation and
an analytic exclusion of the bare local Gibbs identification.

The retained field is the open `2 x 2 x 1` plaquette, `c=1/8`, `U=1`:
eight CAR modes and all 256 occupation states. The region contains the
first two sites, or modes `(0,1,2,3)`, with all sixteen regional states.
No number sector, boundary bond, thermal eigenvalue or detector outcome
is removed. The action, spatial lattice and inverse temperature remain
declared simulator parameters. The existing source implements its field
action, but this package does **not** execute preparation of its full
plaquette Gibbs reference. That would exceed the smaller, already executed
[two-site Gibbs preparation](../field_thermal/README.md).

## Exact construction and uniform remainder

Write `E(X)=Tr_B(X)/dim(B)` for the complete signed CAR partial trace.
It is a unital completely positive map. A scalar shift of H changes only
the scalar part of its logarithm, so use `Hc=H-I`. The exact maximum
absolute row sum is `B=3`, proving `||Hc|| <= 3`. No chemical potential
is inserted into the field generator.

Define `X(beta)=E(exp(-beta Hc))`. Its positive normalized state is
`rho_A=X/Tr_A(X)`. If `M_n=E(Hc^n)`, the noncommuting power-series
logarithm has coefficients

```
-log X = beta C1 + beta^2 C2 + beta^3 C3 + ...,
C1 = M1,
C2 = (M1^2-M2)/2,
C3 = M3/6-(M1 M2+M2 M1)/4+M1^3/3.
```

The implementation preserves every multiplication order and computes
`C1,...,C8` with exact Gaussian-integer sparse arithmetic followed by
rational matrix arithmetic. Its signed partial trace is independently
checked against the existing full-Fock CAR restriction, including an odd
operator on a nonconsecutive region. The actual 256-dimensional field
Hamiltonian is also checked independently.

For complex `z` with `x=|z| B < log(2)`,

```
||E(exp(-z Hc))-I|| <= exp(x)-1 < 1.
```

The norm-convergent logarithm series therefore applies. In every
coefficient, replacing each matrix moment by its norm bound gives the
positive scalar majorant

```
f(x)=-log(2-exp(x))=sum_(n>=1) a_n x^n.
```

It follows directly by the triangle inequality, with no commuting
assumption, that

```
||-log X(beta)-P_p(beta)|| <= R_p(|beta| B),
P_p=sum_(n=1)^p beta^n C_n,
R_p(x)=f(x)-sum_(n=1)^p a_n x^n.
```

This bound is independent of the field and regional dimensions once B
is fixed. B itself can grow with volume; this is **not** a uniform
thermodynamic-limit or fixed-temperature locality theorem. A different
region and Hamiltonian require their own moments and B.

The rational checker uses the smaller domain `0 <= beta B <= 1/2`.
It bounds the exponential by 32 Taylor terms plus a geometric tail and
the logarithm by 128 terms plus a geometric tail. All arithmetic and
rounding are upward exact rational operations. It does not use floating
eigenvalues, an unbounded numerical series, or an empirical tolerance to
prove the remainder.

At the retained `beta=1/32`, the eighth-order operator error is at most

```
9582815241821 / 5000000000000000000000 < 1.917e-9.
```

The exact regional modular Hamiltonian is `K_A=-log X+log(Tr_A X) I`.
The normalization scalar must be retained in statements about K itself.
For **every** Hermitian trace-zero variation delta, including variations
that do not commute with the reference, it cancels and gives

```
|DS_rho_A(delta)-Tr(delta P8)| <= 1.917e-9 ||delta||_1.
```

Tensoring the operator estimate with an untouched reference identity
does not enlarge its norm. This is a finite same-action modular and
entropy-first-law attachment, rather than an identification made by
renaming a local energy operator.

## The first boundary correction and a certified exclusion

`C1` is the bare regional Hamiltonian minus its scalar trace mean, and
`C2=-(99/512) I`. The first nonscalar boundary correction occurs at cubic
order. In the regional occupation-qubit order, highest bit first,

```
C3 - Tr(C3) I/16 =
   (IIIZ+IIZI+IZII+ZIII)/4096
 - (IIZZ+ZZII)/3072
 + (IXXI+IYYI-XZZX-YZZY)/98304.
```

It does not commute with `C1`. Deleting U removes the Z and ZZ terms but
retains the hopping correction. Deleting both crossing bonds makes
**every** nonscalar correction through order eight vanish; exactly, the
disconnected Gibbs product has only the bare regional logarithm plus a
scalar at every order. These controls locate the effects of interaction
and boundary tracing, rather than attributing every mismatch to U.

The regional full-minus-empty diagonal difference of `P8-beta C1`,
together with twice the rational remainder, is strictly negative at
`beta=1/32`. Consequently, for **every** scalar a,

```
||K_A-beta C1-a I|| > 2.7e-8.
```

This is a rational separation from the bare local Gibbs generator on
this declared field, region and temperature. It excludes neither a
corrected modular generator (which is constructed here) nor other
temperatures, algebras, preparations or continuum limits. Independent
full-spectrum logarithms check the bound at four temperatures; those
floating replays are diagnostics, separate from the exact certificate.

## Gravity premise boundary

The construction closes an explicit finite corrected-generator and
first-law-error map for the interacting regional algebra. It does not
identify the polynomial with a geometric cap or null charge, derive a
Lorentz stress tensor or its null translations, select a Newton constant,
or prove fixed-proper-volume entropy/area stationarity. There is no
shrinking-ball comparison, `EinsteinRegisterData` inhabitant or accepted
Einstein/particle completion gate. The high-temperature domain is stated
rather than extrapolated to the physical vacuum.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.regional_modular_expansion --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/regional_modular_expansion
```

Intentional receipt regeneration uses `--write`. Source hashes bind the
complete field action, exact arithmetic, proof implementation and tests.
