# Complete coupled-loop spectra and observer records

Use the exact source and complete Gauss basis in
[`openwave_coupled_haar_source`](../openwave_coupled_haar_source/PROOF.md).
The squared couplings are declared regulator choices, with spacing one.
No measured mass or coupling is fitted.

## Full-Haar eigenvalue enclosure

Write `A=P_N H P_N`, `D=Q_N H Q_N` and `B=P_N H Q_N`.
The native electric and magnetic identities give `D>=d_N` and
`||B||²<=4 beta²`. For a test energy `x<d_N`, the Schur complement obeys

`A-x-4 beta²/(d_N-x) I <= A-x-B(D-x)^-1 B* <= A-x`.

The tail has no negative direction. The negative index of the full
operator therefore equals that of the Schur complement. Suppose the
first four eigenvalues of `A` have certified intervals `[lk,uk]`, and put
`U=u3<d_N` and `delta=4 beta²/(d_N-U)`. The variational upper bounds and
the displayed form inequality give, for these four levels,

`max(0,lk-delta) <= Ek(H) <= uk`.

The exact complete-frame residual proof is the consumer in
[`fixed_witness.py`](../openwave_colour_string/fixed_witness.py).
Every frame column and eigenvalue is a declared dyadic proof coordinate.
Exact integer Gram and residual products are recomputed against the actual
algebraic source matrix, including its rounding bound. No stored error
estimate is accepted. The finite frame is complete, while the separate
Schur inequality controls all omitted representations.

There is a sharper low-spectrum certificate when the first levels are
distinct. Three independent fixed dyadic trial vectors have their residuals
evaluated in the full source, including the native boundary contribution.
The spectral theorem places a spectral point within each residual-radius
interval about its trial energy. If these three intervals are disjoint and
lie below the independently certified fourth-level floor, each contains
exactly the corresponding first, second or third eigenvalue. No additional
low state can fit below that floor. The program checks every separation
before using these tighter intervals. The degenerate independent-loop
control retains the coarse enclosure instead.

## Ground preparation and tighter energy

The normalized trial has dyadic coefficients in the exact ordered
spin-network basis. Each algebraic operator entry has a rational square.
Integer square-root bounds enclose it before rounding to a dyadic matrix;
the maximum row error bounds the operator error. Exact integer products
give the finite trial residual and Rayleigh interval. Only states with
`2j1=N` or `2j2=N` can leave the projection in one action, so the full
residual is at most the finite residual plus
`2 beta ||P_boundary trial||`. This includes the infinite source tail.

Let `mu` be its declared approximate Rayleigh value, `r` the full residual
bound and `L1` the independently proved first-excitation floor. If
`mu<L1`, the ground-projector distance is at most `r/(L1-mu)`. The code
requires this bound to be below one hundredth. For an exact trial Rayleigh
value `m` in `[mlo,mhi]`, its variance is no greater than the residual
about `mu`. Temple's inequality therefore gives

`max(0,mlo-r²/(L1-mhi)) <= E0 <= mhi`,

with the positive denominator checked explicitly. The ground is unique by
positivity improvement of the full elliptic source. The reported gap
uses this ground enclosure and the full-source first-excitation interval.

## Connected response and complete outcomes

Both Wilson multipliers have norm at most one. A projector distance
`epsilon` gives expectation error at most `2 epsilon`. Operator-rounding
errors are bounded separately by exact row sums. `W1` raises only the
first loop spin and `W2` only the second. Their omitted one-step vectors
have disjoint labels. Thus the mixed moment of a supported trial includes
its complete source value when computed from the two finite products.
The bounds for the individual second moments also include the nonzero
outside component. Interval subtraction gives the connected covariance
and a strictly positive variance.

Set the declared readout response to `eta=1/2`. The four joint effects are

`E_(s,t)=(I+s eta W1)(I+t eta W2)/4`, for `s,t=+-1`.

They commute, have floor `1/16`, and sum exactly to the identity on the
whole physical Haar space. Positive square roots define a complete
instrument with an arbitrary external quantum reference retained. Every
outcome probability is enclosed; no selected successful branch is used.
The recorded sign covariance is `eta² Cov(W1,W2)`.

This measurement has apparatus work. For the ground state the increase of
source energy is nonnegative. Writing each Kraus factor as
`f_s(W)=sqrt((1+s eta W)/2)` gives
`sum_s (f_s')²=eta²/[4(1-eta² W²)]` and `sum_s f_s f_s'=0`.
The cross gradients cancel after summing all outcomes. On one fundamental
SU(2) Wilson loop, the sum of its four link-gradient norms is
`1-W²<=1` in the stated Casimir normalization. The complete nonselective
source-energy gain is consequently at most
`alpha eta²/[2(1-eta²)]=alpha/6`. This finite positive budget is an
apparatus input, not a claim of an energy-free readout.

Removing the shared-link electric coupling gives
`4 alpha(C1²+C2²)+beta(2-W1-W2)`. Its unique ground factors into the two
one-loop grounds, so its connected Wilson covariance vanishes exactly.
The separately executed control must enclose zero. The shared source
must have a strictly positive covariance after every error is included.

For `v=(W1-<W1>)Omega`, the ground-state double-commutator identity gives
`<v,(H-E0)v> <= alpha`. Its variance lower bound therefore supplies a
finite energy upper bound for the normalized nonzero physical vector.

All energy units and correlations refer to the specified finite graph
with its full Haar Hilbert space. These inequalities do not identify a
physical glueball, fix a GeV conversion, establish a uniform-volume gap,
or construct a four-dimensional Yang--Mills continuum.
