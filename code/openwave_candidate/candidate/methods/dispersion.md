# Definition-only local mathematical recompute

The sole scientific input was `certification/DEFINITION_ONLY.md`. No OPH
source, other result files, target values, experimental data, web results,
or other agents' output was read. The final file-read audit consists of that
packet and this recomputation's own module, tests, note, and JSON artifact.
Installed Python, SymPy, and mpmath libraries are ordinary implementation
dependencies, not additional scientific inputs. Target values were not provided.

Write `r2=x²+y²+z²`, and let `Md=Σu(u·(x,y,z))^d`. The module constructs
the twelve ports and expands these sums over the exact field Q(√5).
The normalized spherical monomial integral used to derive invariant parts is
`<x^(2p)y^(2q)z^(2r)>=(2p−1)!!(2q−1)!!(2r−1)!!/(2p+2q+2r+1)!!`;
an odd coordinate power has zero integral.

The resulting moments are `M2=4 r2`, `M4=(12/5) r2²`, and

```
M6 = (8/5)(x^6+y^6+z^6)
   + (6−6√5/5)(x^4 y^2+y^4 z^2+z^4 x^2)
   + (6+6√5/5)(x^2 y^4+y^2 z^4+z^2 x^4).
```

The homogeneous contributions to the symbol are
`Q2=M2/4`, `Q4=−a² M4/48`, and `Q6=a⁴ M6/1440`.
Their invariant parts are respectively `r2`, `−a² r2²/20`, and
`a⁴ r2³/840`. In particular the first angular contribution is degree six.
Exact differentiation shows `H6=M6−(12/7)r2³` is harmonic.

Generating `P6(t)=(231t^6−315t^4+105t²−5)/16` and homogenizing its port
sum gives `ΣP6=(231/16)H6`. Its value at every normalized port is
`132/25`, derived by exact substitution. Thus `I6=(175/64)H6` and
`M6=(12/7)r2³+(64/175)I6`, where the homogeneous I6 polynomial is
evaluated on a unit vector when used as an angular function. Its full exact
monomial coefficients are included in the JSON. It has zero spherical mean
and equals one at all twelve ports.

Consequently, for unit n,

```
omega_squared = k² − (a²/20)k⁴
              + a⁴[1/840 + (2/7875)I6(n)]k⁶ + O(a⁶k⁸).
```

With the explicit factors of a stripped off, the angular/invariant sextic
ratio is `16/75`; the invariant and angular sextic coefficients divided by
the square of the quartic coefficient are respectively `10/21` and `32/315`.
At fixed k and n, replacing a by a/2 multiplies the quadratic, quartic,
sextic, and eighth-order bound by `1`, `1/4`, `1/16`, and `1/64`.
These ratios and the angular normalization are unchanged.

The numerical check reconstructs the ports independently in mpmath at
80 decimal digits. It evaluates the original cosine sum for four directions,
two lengths, and five decreasing momenta. It also evaluates the equivalent
squared-sine expression to check loss of accuracy from subtracting cosines.
For every real argument Taylor's theorem gives the analytic bound

```
|omega_squared − Q2(k n) − Q4(k n) − Q6(k n)|
    <= B8 = a^6 k^8 Σu(u·n)^8 / (2·8!).
```

All selected arguments have magnitude at most 0.4. Alternation further gives
`−B8 <= error <= −B8+B10`, where
`B10=a^8 k^10 Σu(u·n)^10/(2·10!)`. The result records each direct value,
truncation, error, both analytical bounds, and the error/bound ratio. The
comparisons allow only `1e-70` absolute numerical tolerance. This combines
analytical remainder control with numerical validation, rather than relying
only on a small difference of two floating point values. The residual falls
by more than 200 for each momentum halving, as consistent with order eight.

Two deliberately incorrect operators are constructed: one omits one port,
and the other doubles one port's weight. The same exact quadratic isotropy
test rejects both. The main calculation uses exactly unit weights throughout;
no tuning against a target is performed.

Run `python candidate/oph_openwave/dispersion.py` to emit JSON, and
`python -m unittest discover -s candidate/tests -p test_dispersion.py` to
run the verification suite. Acceptance here is internal mathematical
consistency only. This is not external certification or an empirical
prediction. No physical identification of a field or length is supplied or
inferred, and `full_criterion_met` remains false.
