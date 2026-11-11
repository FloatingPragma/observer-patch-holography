# Native abundance and persistent defects on a closed cell

## Construction and complete records

The declared cell has three ports, oriented seams `(0,1)`, `(1,2)`, `(2,0)`, and two triangular faces with opposite
orientations. Their fluxes are `k` and `-k mod 6`. The cell is a closed finite complex; it is not an induced cell of
the icosahedral federation. The local repair and register operations are imported directly from
`model/defect_integer_lift/model.py`.

A state is `(x0,x1,x2,k)`, with each load in `{0,...,5}` and `k in Z6`. On the selected seam `(a,b)`, the native
repair balances `xa+xb` according to a coin. Write `j=xa-new_xa` and `d=xa-xb`. The register applies `k -> k+j`
when its pair-cost change is at most `5 j(d-j)`. The pair cost is `C(k)+C(-k)`, where
`C=(0,5,8,10,8,5)`. Rejected twists retain the load repair. Zero currents retain the attempted event.

The scheduler gives each of the six seam/coin records probability `1/6` per attempt. The exact kernel merges
coincident target states by adding these record probabilities. It has 1,296 states and 7,776 complete attempt
records. The preparation has independent uniform loads in `{0,...,5}` and zero initial flux. No parameter is
selected using a particle abundance or cosmological observation.

## Conservation and dissipation

Every transition preserves the integer load sum and the total face class modulo six. Let `Q=sum_i xi^2` and
`R=C(k)+C(-k)`. The balancing identity is

```
Q_after - Q_before = -2 j(d-j),       j(d-j) >= 0.
```

For an applied twist, `R_after-R_before <= 5 j(d-j)` by the native acceptance rule. For a rejected twist this
difference is zero. Thus `E=5Q+2R` is nonincreasing on every complete record. This is a native finite Lyapunov
quantity. It is not a calibrated physical energy or a thermodynamic temperature.

For fixed total load `M=3q+r`, `0<=r<3`, the least possible load square is
`Qmin=(3-r)q^2+r(q+1)^2`. On this domain `Q-Qmin<=16`. Until all load contrasts are at most one, at least one of
three seams has contrast at least two. Either coin then lowers `Q` by at least two. At most eight such decreases
are needed. The exact population therefore satisfies

```
Pr(loads unbalanced after n attempts) <= sum_{a=0}^7 binom(n,a) (1/3)^a (2/3)^(n-a),
```

with the sum truncated at `a=n`. This bound follows by domination by Bernoulli trials of probability `1/3`; it
does not assume independent successive states. Every stored epoch is checked against it in exact arithmetic.

## Abundance law and its limit

The complete population law is `p_(n+1)(s') = sum_s p_n(s) K(s,s')`. The face share of class `a` is
`sum_s p_n(s) ([k=a]+[-k=a])/2`. The six shares sum to one, including zero flux. A class-only transition matrix
cannot reproduce this law for arbitrary preparations: states with the same zero flux but different load
contrasts have different birth probabilities. The joint load/register state is retained throughout the
calculation.

Balanced equal loads permit no further current, so any class on such a state is persistent. Balanced unequal
loads permit only class-cost-decreasing twists; a suitable seam and coin reduce the class distance from zero.
Zero flux on balanced loads is closed, although load permutations can continue. Define the terminal set as
balanced states with either zero flux or all three loads equal. Every state has a positive-probability route to
this set: at most eight strict load descents followed by at most three class descents. Since there are finitely
many states and each selected record has probability `1/6`, absorption occurs with probability one. A uniform,
loose bound on survival after `11b` attempts is `(1-6^(-11))^b`.

Exact backward propagation of the terminal indicator strengthens this result on every one of the 1,296 states.
The minimum probability of absorption within three attempts is `1/108`. Within 32 attempts the minimum is
`39899059579819656509585/41451359947637504606208`, exceeding `24/25`. The Markov property therefore gives

```
Pr(outside the terminal set after n attempts) <= (1/25)^floor(n/32)
```

for every initial distribution, including distributions correlated across loads and flux. This is a uniform
operational relaxation bound. It uses the complete attempted-event kernel and no fitted rate.

The eventual occupied-face probability `h(s)` is the unique bounded solution of `h=Kh`, with boundary value one
on nonzero frozen classes and zero on the zero-flux terminal states. The implementation solves this system in
exact rational arithmetic. It orders blocks by load square and then register cost; equal-block transitions are
solved by elimination. Every state is checked against the complete backward equation. Averaging over the
declared preparation gives exactly

```
lim_(n->infinity) expected occupied face share = 1007/13824 = 0.0728443287037037...
```

At any finite epoch, the probability of reaching a nonzero frozen state is a lower bound on this limit. Adding
the probability outside the terminal set gives an upper bound. At 128 attempts the interval width is
`5.852613853690619e-10`. These are exact rational probability bounds, not sampling errors.

## Physical scope

An epoch counts attempted native events, including rejected twists and no-operation events. The source has no
expansion law, physical seconds, thermal preparation, temperature readout or particle dictionary. Its class
populations are finite simulation-generated records. They supply a quantitative population and persistence
law for this declared source. They do not establish a thermal relic abundance, nucleosynthesis or a cosmological
freeze-out temperature, and no transfer to the full federation is claimed.
