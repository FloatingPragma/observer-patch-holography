# A Gaussian history action from native fluctuations

This packet derives a Gaussian history law from explicitly independent copies
of the existing finite native process. Its temporal covariance and quadratic
action follow from native transition counts. Gaussian noise and oscillator
coefficients are not assigned to a single carrier. Independence, stationary
preparation, the binary total-two sector and the uniform thirty-seam scheduler
are premises of this product construction.

It also predicts an accessible native measurement: signed currents at one port
have successive-attempt correlation **−1/20**, including zero currents. The shared primitive runtime
executes three such records, retaining all 27 outcomes and their complete
terminal quantum states. No measurement of the initial hidden load is needed.

## Native covariance and the limiting action

Let `x_n` be the twelve binary loads, with total two, after native attempt `n`.
The stationary law is uniform on its 66 configurations. Write
`z_n=x_n-(1/6,...,1/6)` and let `L` be the source carrier Laplacian. Exact
native enumeration, checked on every configuration, gives

```
E[z_(n+1) | z_n] = A z_n,       A = I-L/60,
C = E[z_n z_n^T] = (5/33)(I-J/12),
E[z_0 z_n^T] = C A^n.
```

The last identity follows by conditional expectation and is checked against
direct propagation on the finite configuration kernel for lags zero through
four. `J` denotes the all-ones matrix. The conserved total-load mode is removed;
the remaining sector has eleven dimensions. This number is not a claim about
physical spatial dimensions.

For `M` independent copies of the **whole stationary history**, set
`F_M(n)=M^(-1/2) sum_r z_n^(r)`. Every fixed finite set of these variables has a
Gaussian limit, with the same covariance `C A^|n-m|`. The characteristic bound
below proves this directly for arbitrary fixed linear combinations. Independence
at each single time alone would not suffice; the joint histories are independent.

The Gaussian limit is the stationary vector autoregression

```
F_(n+1) = A F_n + eta_n,        Cov(eta_n)=D=C-A C A^T.
```

Here independent Gaussian innovations are a representation of the derived
limit. A single native carrier has state-dependent, discrete innovations. On
the mean-free sector the source graph has `0<L<=10 I`, so `0<A<I` and `D>0`.
The finite-history density, with respect to Lebesgue measure on that sector,
is proportional to `exp(-S)` with

```
S(f_0,...,f_N) = (1/2) f_0^T C^+ f_0
              + (1/2) sum_(n=0..N-1)
                  (f_(n+1)-A f_n)^T D^+ (f_(n+1)-A f_n).
```

The superscript `+` denotes inversion on the mean-free sector, zero on the
conserved direction. The initial term is necessary for the stationary law;
omitting it describes a conditional or improper initial-state problem.
The implementation builds this precision and compares its inverse with the
native history covariance, with residual below `4e-14` in the retained run.

The action determines relative coefficients in declared native load units and
unit native-attempt spacing. It is an action for the Gaussian fluctuation
history, not an action selected for the laboratory oscillator apparatus.
`D^+` can couple distant ports. Neither spatial locality nor a relativistic
continuum follows from this fixed twelve-port construction.

## A finite bound instead of assuming Gaussianity

For any centered scalar linear functional `Z` of a finite one-carrier history,
let `v=E[Z^2]` and `b=E[|Z|^3]`. Finite carrier states imply `b<infinity`.
For `M` independent histories and Fourier coordinate `theta`, Taylor's remainder
and telescoping products give

```
|E exp(i theta sum_r Z_r/sqrt(M)) - exp(-theta^2 v/2)|
    <= |theta|^3 b/(6 sqrt(M)) + theta^4 v^2/(8 M).
```

To see the bound, put `a=theta^2 v/2`. The single-carrier characteristic
function at scale `1/sqrt(M)` differs from `1-a/M` by at most
`|theta|^3 b/(6 M^(3/2))`; `exp(-a/M)` differs from `1-a/M` by at most
`a^2/(2 M^2)`. Both characteristic factors and exponential factors have
modulus at most one, so the error of their `M`th powers is at most `M` times
their single-factor difference. Applying this to every linear functional
gives the finite-dimensional Gaussian limit. No convergence in total variation
to a continuous Gaussian is claimed for the discrete finite-copy law.

This is an analytic proof with exact finite moments, not a new Lean CLT.
The source response and conservation identities reuse
[native_fields.py](../common_source/native_fields.py). The positive mode
multipliers meet the premises of the existing
[bond-reflection Gram identity](../../Lean/SourceLaw/ReflectionPositivity.lean),
with mode covariance `c0=5/33` and `r` an eigenvalue of `A` on the mean-free
sector. The supplied time-register action in that file is not silently
identified with this derived precision.

## What the observer executes and predicts

For a selected seam `(a,b)` and committed native current `j`, continuity gives

```
delta x_port = -(1_port=a - 1_port=b) j.
```

Thus the signed change is computed entirely from the admitted seam/current
record. The program coarsens all actual native C Kraus branches to outcomes
`-1,0,+1` at port zero, starting from `I66/66`. It executes three successive
attempts. Hidden seam and coin alternatives are probability-summed; outcomes
of zero probability remain in the complete dictionary. Every terminal density
is compared with exact native path counts. A checkpoint after the first read
is serialized, independently replayed and continued.

The exact one-carrier predictions are

| Quantity | Native value |
| --- | --- |
| Variance of the signed current per native attempt | `5/198` |
| Covariance of signed currents on successive native attempts | `-1/792` |
| Their correlation | `-1/20` |
| Covariance of signed currents two native attempts apart | `-161/142560` |

For the fixed three-record contrast `Z=delta x_0-2 delta x_1+delta x_2`,
the full finite law has scores `-4,...,4`. Its variance is `11359/71280`,
its absolute third moment is `105643/356400`, and its fourth cumulant is
`918641807/1693612800`, which is nonzero. The contrast is therefore not
Gaussian at one carrier.

The product prediction includes every outcome by exact factorization, avoiding
a dense `66^M` state simulation. At 64 independent copies and `theta=1`, its
characteristic function is approximately **0.92373845**, compared with
**0.92341296** in the Gaussian limit. The analytic error bound is **0.00622496**;
the evaluated difference is about **0.00032549**. The fourth cumulant is exactly
the single-copy value divided by `M`. Perfectly correlated copies instead give
variance `M v` under the same normalization: replication without independence
does not establish this limit.

These are predictions of the specified native construction. There is no new
laboratory dataset in this packet. The executed native channel semantics do
not emit source pulse words; float matrix residuals are diagnostics rather
than outward numerical or physical service certificates. The Gaussian limit
does not select C against E, derive coherent preparation, supply an observer's
physical clock or establish an interacting quantum field theory. Its new
content is a controlled history-action attachment and complete accessible
current predictions from the same native transition law.

From the repository root:

```sh
python -m model.native_fluctuations --check model/native_fluctuations/results.json
python -m pytest -q -W error model/native_fluctuations
```
