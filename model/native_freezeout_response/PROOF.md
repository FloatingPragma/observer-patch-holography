# Complete native absorption and source-count freezeout

The source is the unchanged 1,296-state native population kernel. Every
attempt keeps its seam, coin, applied or rejected twist and no-operation
record. The complete 7,776-record law is bound through
`model/native_thermal_admission/results.json`. No thermal sampler, cooling
schedule or physical particle decoder is introduced here.

## All terminal classes and lifetimes

Let `C_0,...,C_45` be the 46 closed classes classified by the exact native
positive-edge graph. They contain 66 states in total. Write `tau` for the
first entry into their union, counting native attempts from zero. For
each starting state define

`h_c(s)=Pr_s(X_tau in C_c)`, `t(s)=E_s tau`,
`g(s)=E_s tau²`, `u_c(s)=E_s[tau 1_(X_tau in C_c)]`.

Their boundary values are `h_c=1_(C_c)` and `t=g=u_c=0` on the closed set.
On every transient state the exact backward equations are

`h_c=K h_c`, `t=1+K t`, `g=2t-1+K g`, `u_c=h_c+K u_c`.

The source orders states by load square, register cost, total load and
flux. Every transition to another block lowers that order. A transient
permutation block has at most six states; rational elimination solves
all 46 harmonic functions and the lifetime equations on those literal
blocks. The source absorption bound
`Pr_s(tau>n)<=25^(-floor(n/32))` makes every moment here finite and proves
uniqueness of these bounded harmonic and Poisson solutions. The receipt
checks the equations on every original state, nonnegativity, `sum_c h_c=1`,
`sum_c u_c=t` and `g>=t²`, using exact rational arithmetic.

For the specified independent uniform loads and zero initial flux,

`E tau = 5557709/1125000`,
`E tau² = 177070889/2812500`,
`Var(tau) = 48793770721319/1265625000000`.

The maximum mean over all 1,296 input states is exactly 15 attempts.
The complete terminal-class probabilities and conditional mean lifetimes
are recorded, including zero-probability classes. These are lifetime
statistics of the native relaxation protocol, not particle decay widths.

## Full class distributions and stationary limit

The complete finite population is `p_n=p_0 K^n`. The receipt retains its
entire 1,296-state vector at every stated cut, rather than a renormalized
occupied subset. The limiting state has class weight `p_0 h_c`, uniform
within each closed class. Singleton classes are fixed; the three-state
classes have kernel `(1/2)I+(1/6)11^T` and mix to uniform.

For the native uniform zero-flux preparation, the limiting oriented-face
class distribution in order `0,...,5` is exactly

`(12817/13824, 499/27648, 25/1728, 1/128, 25/1728, 499/27648)`.

It reproduces the persistent occupied fraction `1007/13824` and supplies
its six-class resolution. Opposite faces carry `k` and `-k`; averaging
both faces makes the class distribution even for every individual source
state. The oriented face can carry an odd contrast under asymmetric
input, so this distinction is retained explicitly. Neither observable
is identified with baryon number or a Standard Model particle family.

Once the source reaches a closed class its face class is fixed. Coupling
the finite trajectory to its eventual class therefore gives

`TV(class(p_n),class(p_infinity)) <= Pr(tau>n)`.

For the full state, split at `m=floor(n/2)`. The histories absorbed by
that time have at least `n-m` remaining mixing attempts, with exact
three-state worst-case distance `(2/3)2^(-(n-m))`. The other histories
have probability at most the native survival bound. Thus

`TV(p_n,p_infinity) <= 25^(-floor(m/32)) + (2/3)2^(-(n-m))`.

An operational source-count freezeout criterion is the first declared
attempt budget whose survival upper bound is at most a tolerance epsilon.
For example 96 attempts guarantee face-class error at most `1/15625`,
below `10^-3`, for every initial native law. The receipt supplies exact
budgets for tolerances through `10^-12`. This is cessation of native
class changes, without expansion, temperature or cosmological time.

## Input and boundary response

Absorption is linear in preparation: for
`p_epsilon=(1-epsilon)p+epsilon q`, the terminal class law is
`(1-epsilon)(p h)+epsilon(q h)` exactly. Stochasticity gives
`TV(p h,q h)<=TV(p,q)`. The universal coefficient one is sharp, attained
by two different absorbing states. At fixed total load not divisible by
three, there is a single recurrent class, so its limiting class weight
has zero preparation sensitivity. Loads divisible by three retain six
independent frozen flux possibilities.

The literal asymmetric input `(0,1,5,1)` has terminal oriented contrast
`Pr(k=1)-Pr(k=5)=-91/384`. The receipt records its exact derivative
relative to the native preparation. This response is driven by the
specified input; no spontaneous asymmetry is inferred. For any supplied
terminal payoff vector `b`, the full boundary readout is `h_b=H b`,
with `||h_b-h_d||_infinity<=||b-d||_infinity` because every row of `H` is a
probability vector. All boundary outcomes, including zero classes, are
retained. This is the finite-record push/pull convention of
`Lean/Dynamics/ProtectedCharge.lean`; the conjugation-even special case
uses the preserved symmetry in `Lean/Dynamics/ChargeConjugation.lean`.

A separate executed boundary control starts at `(0,1,5,0)`. After eight
native attempts it resets to `(0,1,5,1)` with probability `1/10`, retaining
the public reset/no-reset bit. Both states have total load six. This
supplied preparation operation is not part of the native repair kernel;
it can inject cost and changes no claim about native dissipation. The
complete post-reset population uses the same kernel. Its eventual law
and all finite cuts are calculated without dropping the reset failure.
Each post-reset cut stores both unnormalized record-conditioned state
vectors. They are propagated independently with the native kernel, have
the original masses `9/10` and `1/10`, and sum exactly to the separately
propagated marginal state. Thus the retained-bit joint law is executed,
including its correlation with every terminal state.
Any relaxation budget for this control starts at the final reset, not
at the earlier preparation or a physical cosmic epoch.
