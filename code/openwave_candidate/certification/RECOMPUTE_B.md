# Independent definition-only recomputation B

This is a local mathematical recomputation of the supplied scalar symbol. It
is neither external certification nor an empirical prediction.

## Input isolation and reproduction

The only external task-content file read was `DEFINITION_ONLY.md`. No candidate
dispersion implementation, other recomputation, OPH source, target value, or web
source was read. The implementation reads its own file solely to record its
SHA-256 hash. After the run, its own `recompute_b.json` was inspected to prepare
this note. Expected coefficients were not supplied or used to choose weights
or normalization (`targetvalues_not_provided: true`).

Reproduce from the repository root:

```sh
python3 openwave-tuning/candidate/audits/dispersion_oracle.py --output openwave-tuning/certification/recompute_b.json
```

The relative paths in this command assume the repository root; both arguments
may instead be absolute. The script locates its definition packet relative to
its own location. It requires Python, SymPy, and mpmath.

## Derivation

Write `r2=x^2+y^2+z^2`, `L2=1+phi^2`, and use the coordinate-plane pairs
`(x,y), (y,z), (z,x)`. Eliminating the independent signs gives the exact moment

```text
M_(2m) = 4/L2^m * sum_pairs sum_(j=0)^m binomial(2m,2j)
                                  * first^(2m-2j) * phi^(2j) * second^(2j).
```

The implementation calculates this expression for degrees two, four, and six,
then applies the Taylor multiplier `(-1)^(m+1)/(2*(2m)!)`. An independent exact
sum over the twelve explicit unit vectors checks each generated moment.

Rotationally invariant coefficients are calculated by exact spherical
integration of Cartesian monomials:
`mean(x^(2i)y^(2j)z^(2l))=(2i-1)!!(2j-1)!!(2l-1)!!/(2i+2j+2l+1)!!`.
The homogenized sum of the packet's degree-six Legendre polynomial is
normalized at a supplied vertex. The normalized polynomial has zero spherical
mean, vanishing Cartesian Laplacian, and value one at all twelve vertices.
These are checked by exact algebra.

The homogeneous terms, with `p=k*n=(x,y,z)`, are

```text
omega_squared = H2 + a^2 H4 + a^4 H6 + O(a^6 |p|^8)
H2 = r2
H4 = -r2^2/20
H6 = (x^6+y^6+z^6)/900
   + (5-sqrt(5))*(x^4*y^2+y^4*z^2+z^4*x^2)/1200
   + (5+sqrt(5))*(x^2*y^4+y^2*z^4+z^2*x^4)/1200.
```

The Legendre-sum value at a vertex is `132/25`. With exactly the requested
normalization,

```text
omega_squared/k^2 = 1 - (a*k)^2/20
                    + (1/840 + (2/7875)*I6(n))*(a*k)^4
                    + O((a*k)^6).
```

The first nonconstant angular part is degree six. The dimensionless ratios
`c6_I6/c6_isotropic`, `c2*c6_isotropic/c4^2`, and `c2*c6_I6/c4^2` are respectively
`16/75`, `10/21`, and `32/315`. At fixed `k,n`, replacing `a` by `a/2` multiplies
the quartic term by `1/4`, both sextic terms by `1/16`, and the eighth-order
remainder bound by `1/64`; these ratios are unchanged.

## Independent numerical and negative checks

The direct numerical path independently constructs the original unit vectors
and evaluates the original cosine symbol at 110 decimal digits. It also checks
the equal `sin^2(t/2)` representation. Five directions, five decreasing momenta,
and two lengths produce 50 comparisons. The direct Legendre sum is separately
compared with the exact normalized angular polynomial.

For every row, the absolute degree-six truncation error is bounded by
`sum_u |a*k*(u.n)|^8/(2*a^2*8!)`, itself at most `a^6*k^8/6720`. This follows
from Taylor's theorem at degree seven and the bound on the eighth derivative
of cosine, so verification does not rely solely on subtracting nearby values.
All comparisons pass; errors decrease with momentum and show the expected
eighth-order behavior.

The negative control removes one original port. The traceless quadratic
tensor then has squared Frobenius norm `1/24`, whereas the original tensor's
value is zero. Exact isotropy testing rejects the altered operator, including
after scalar rescaling to restore the original total weight. All recorded
checks pass. The full exact polynomials, hashes, numerical rows, and Boolean
checks are saved in `recompute_b.json`.
