# Exact record-instrument discrimination

## Carrier, repair and recorded channels

Let X be the 78 nonnegative integer configurations on twelve ports with
total load two. The independently specified compatibility condition on seam
`e=(0,1)` is `|x0-x1|<=1`. For a fair coin `c` the native update is

```text
d = x0-x1
j = [d-(2c-1)(|d| mod 2)]/2
F_c(x) = x+j(unit_1-unit_0).
```

This balances the endpoints by floor and ceiling and conserves total load.
The two candidates differ in their quantum records:

```text
C: L_(c,j) = (1/sqrt(2)) sum_{x: j(x,c)=j} |F_c(x)><x|
E: M_(c,d) = (1/sqrt(2)) sum_{x: x0-x1=d} |F_c(x)><x|.
```

For fixed current, the endpoint transformation is an injective translation.
For fixed contrast, the unchanged spectators and endpoint sum reconstruct
the original endpoints uniquely. Each input basis column has two entries of
amplitude `1/sqrt(2)` in the joint record/source carrier. Distinct input
columns have disjoint joint supports. Consequently

```text
sum_r K_r^* K_r = I
```

for each candidate on all 78 input states. The Kraus construction proves
complete positivity; the displayed identity proves trace preservation,
including for sources entangled with a reference system. The implementation
checks every column exactly. It also retains the subnormalized density of
each record branch.

These are the C/E channels specified in the
[record-law construction](../record_laws/README.md).
The calculation here adds a fixed discrimination experiment and its complete
count decision rule. It does not introduce a distinct general channel theorem.

## Interference witness with fixed record meaning

Write `x={2,3}`, `z={1,2}` and `w={0,2}` for the corresponding binary
occupation states. In each notation the two listed ports carry one load each.
For the pure preparations set

```text
|phi_theta> = (|x>+exp(i theta)|z>)/sqrt(2),
P_plus = |phi_0><phi_0|,   P_minus = I-P_plus.
```

The fixed phases are zero, pi and pi/2. A fourth preparation is the
incoherent density `(|x><x|+|z><z|)/2`. All have the same initial classical
load probabilities. The preparations and projectors are specified
independently of the candidate being tested and of any comparison data.

For C the nonzero recorded source densities after one repair are

```text
(c,j)=(0,0):   (1/2)|phi_theta><phi_theta|
(c,j)=(1,0):   (1/4)|x><x|
(c,j)=(1,-1):  (1/4)|w><w|.
```

For E, the contrast label distinguishes x from z. Its nonzero branches are

```text
(c,d)=(0,0):   (1/4)|x><x|
(c,d)=(0,-1):  (1/4)|z><z|
(c,d)=(1,0):   (1/4)|x><x|
(c,d)=(1,-1):  (1/4)|w><w|.
```

The common readout probability is the sum of
`Tr(P_plus rho_r P_plus)` over every recorded branch. Therefore

```text
p_C(theta) = 3/8+(1/4)cos(theta),
p_E(theta) = 3/8.
```

The incoherent preparation gives `3/8` under both channels. For the pure
settings, the complete common-interface distribution is:

| Coin, current | C: plus | C: minus | E: plus | E: minus |
| --- | --- | --- | --- | --- |
| `0,0` | `(1+cos(theta))/4` | `(1-cos(theta))/4` | `1/4` | `1/4` |
| `1,0` | `1/8` | `1/8` | `1/8` | `1/8` |
| `1,-1` | `0` | `1/4` | `0` | `1/4` |
| Every other coin/current pair | `0` | `0` | `0` | `0` |

The diagonal control has the same distribution as the table with
`cos(theta)=0`. In particular, current `+1` has zero probability for every
setting. The comparison retains all twelve common bins, including these
zeros. E's additional contrast records are summed as probabilities for the
common interface; classical grouping cannot restore their lost coherence.

The final readout applies `P_y rho_r P_y` within each retained record block.
Since `P_plus+P_minus=I`, summing its two probabilities preserves the
probability of each old record exactly. The source measurement never couples
different classical record values. Its interference concerns x and z
within the same C branch. Neither erasure of a committed record nor
postselection is required.

## Fixed statistical contract

There are four settings and 8,192 fresh independent attempted trials at each
setting. No optional stopping is permitted. Each trial contributes to one
of the twelve common record/outcome bins or to abort. The plus indicator
equals one for any completed plus outcome and zero for completed minus
outcomes and aborts. Every empirical frequency uses all 8,192 attempts.

At each setting score twelve joint indicators, the plus indicator and the
abort indicator. This gives `m=4*(12+2)=56` Bernoulli statistics. Independence
between these statistics is unnecessary. For each statistic, independence
of its fresh trials gives the two-sided Hoeffding bound

```text
Pr(|observed_frequency - actual_probability| > b) <= 2 exp(-2 n b^2).
```

Use `n=8192` and `b=1/32`, so `2 n b^2=16`. A union bound yields

```text
Pr(any of the 56 intervals misses its actual probability) <= 112 exp(-16).
```

The software avoids a rounded exponential in the certificate. With
`S_32=sum_{k=0}^{32} 16^k/k!`, positivity of the exponential series gives
`S_32<exp(16)` and hence

```text
112 exp(-16) < 112/S_32 < 0.000012606 < delta=0.001.
```

All quantities in this certificate use rational arithmetic. The reported
observation interval is `[frequency-b,frequency+b]` clipped to `[0,1]`.

The apparatus contract separately requires each scored actual probability
to differ from the ideal candidate probability by at most `eta=1/64`.
Preparation, repair, storage and readout contribute to this allowance. The
budget applies to each aggregate probability as well as each joint bin;
it is not inferred by multiplying or summing unverified component estimates.
The ideal abort probability is zero, so the same contract bounds actual
abort probability by eta.

The candidate's allowed interval is `[p-eta,p+eta]` clipped to `[0,1]`.
A disjoint observation interval rejects that candidate under the contract.
If the contract holds for a candidate, rejection can occur only when at
least one simultaneous observation interval fails. This proves the stated
familywise bound. Both candidates use the same observation intervals, so
testing both requires no additional union bound over candidate names.

The budget must be established independently of the comparison counts.
If the calibration certificate itself fails with probability at most beta,
the combined error bound is at most `112/S_32+beta`. The fixed sampling
certificate alone supplies no calibration guarantee. Likewise, acquiring
the declared number of counts does not establish independence or verify
the physical preparation and readout identification.

## Missing data, failed services and scope

An abort retains an attempted trial for which no complete valid record and
readout are available. It is never dropped or redistributed among completed
outcomes. If the lower observation bound for its probability exceeds eta,
the result is `apparatus_budget_exceeded`. An ideal-rule attribution then
lacks its service premise.

The decision rule requires complete counts for all settings. An unfilled
template is invalid data. With no supplied dataset it emits `no_data` and
no candidate verdict. Synthetic draws exercise both discrimination
directions but are labeled software controls. External files remain
unverified; the program cannot authenticate acquisition, calibration or
provenance from a JSON label. Every result sets empirical confirmation to
false.

The experiment orders preparation, one repair and one measurement. It makes
no duration or rate prediction and requires no universal clock or conversion
from host runtime to physical time. Coherent preparation, record storage and
readout are supplied apparatus services that an independent realization
must implement. Both candidate instruments obey ordinary quantum mechanics.
A successful conditional discrimination identifies their different recorded
channels, without establishing continuum physics or a complete OPH model.

## Separate electromagnetic comparison

The [endpoint audit](endpoint_audit.py) imports the pinned interval from the
[P-closure implementation](../../code/P_derivation/README.md). Its incomplete
`thomson_structured_running` output has inverse alpha near `136.994835177413`.
The [NIST CODATA 2022 table](https://physics.nist.gov/cuu/pdf/wall_2022.pdf)
reports `137.035999177(21)`. The exact endpoint difference corresponds to
about `300.388218` parts per million of the published central value.

The audit uses a comparison band of five reported standard uncertainties.
It treats that band as a declared allowance, not a rigorous confidence
interval. The numerical source enclosure is disjoint from it. Missing
physical transport has no supplied remainder bound, so no theory sigma
significance follows. The inconsistent hypothesis is the extra assumption
that the incomplete map supplies the entire physical endpoint with zero
missing correction. The retrospective audit does not refit P, establish a
new prediction, or decide OPH as a whole.
