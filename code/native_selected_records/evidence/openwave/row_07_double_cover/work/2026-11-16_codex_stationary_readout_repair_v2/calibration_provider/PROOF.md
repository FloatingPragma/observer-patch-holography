# Continuing-state calibration theorem and cost

Assume an admitted signed unit sigma=±1 occupies one public target vertex and
the target's remaining two entries are zero, with an arbitrary fixed spectator
background elsewhere. Assume the complete target triangle is addressable,
source coins are IID fair, and the native deterministic endpoint rule is
stationary. Admission's larger unknown-input proof is inherited explicitly;
the calibrator itself is not a proof of its starting premise.

For an addressed oriented edge(a,b), the observed current changes the tracked
unit vector by(-j,+j). Thus its next position is decoded solely from its known
previous position, sign and observed current. Invalid or splitting records
are rejected rather than projected back into the unit domain. Native repair
of a signed unit against zero either holds or moves it, each with probability
1/2. This remains true for a negative unit, with the raw coin roles exchanged.
The compiler retains those raw coin/current records.

Routing between any two target positions needs one addressed triangle edge.
Repeated native attempts reach the requested destination with conditional
failure probability2^-C after cap C. The position is reset only by actual
current-confirmed moves, never by state assignment. Each probe is made after
such a route has actually succeeded, so its labelled input position is an
observed state of the same continuing particle. The source is not reinitialized
between columns.

There are nine input-position/edge groups. A group observes both coin values
unless its first C probe coins are all equal: probability2^(1-C). Intervening
routing is adapted to prior records; the future IID probe coins remain fair
at these stopping times. There are at most9C probe events and9C+1 routing
phases including final return. A union bound therefore gives

```
Pr(core calibration timeout) <= (9C+1+18)/2^C = (9C+19)/2^C.
core endpoint calls <= C(9C+1)+9C = 9C²+10C.
```

The optional blank coupon check adds failure at most2/2^C and at most C
endpoint calls. At C=64 the combined failure bound is597/2^64 and the maximum
number of endpoint opportunities is37568. These counts include every failed
route attempt and every duplicate probe coin. Timeouts are not falsifications
of the asymptotic statement and are not successful calibrations.

An uncapped two-coin coupon group takes3 expected probes. Hence expected probe
count is at most27; every necessary route takes2 expected endpoint events.
There are at most one route per probe plus final return, giving expected
core cost at most27+2(27+1)=83. Blank adds3, giving86. Early caps and termination
cannot increase this bound. Conditional expected cost given success is at
most86/(1-597/2^64). This is a theorem under IID fair coins, not a statistical
claim that the finitely seeded cross-input replay trials are independent.

Only target edges and the blank edge execute after admission. On target
events, subtracting the tracked sigma unit from the full state leaves the
spectator background unchanged. The blank edge is a declared vacuum
reference, so its events leave that background unchanged as well. Consequently
the final return restores the entire post-admission vector. Visiting the
shared core/target vertex0 changes its signal occupancy temporarily, never
the spectator component. All other reservoir/marker/copy coordinates remain
unchanged event by event.

For a stationary deterministic source conditioned on input position, addressed
edge and observed coin, the two observed coin columns determine that finite
endpoint table. Repeated observations must agree; inconsistency is retained
as a failed calibration. This is not a universal identification theorem for
arbitrary hidden-memory or noisy laws. Fair coin weights are a declared source
premise, not inferred merely by having seen both coin values.

Let t index declared edge tickets, with probabilities w_t, and let T_tc be
the observed one-hot transition map for coin c. The measured mean is

```
P = sum_t,c (w_t/2) T_tc + w_blank I.
```

Current-resolved maps keep the observed current labels. Their weighted Gram
sum is checked exactly from the table; no completeness assertion is inserted
by setting a result to I. The proposed reversal maps(p,t,c) to
(observed_output,t,c xor[current!=0]); permutation and involution are also
checked rather than assumed. The inactive fair-coin reference is identity;
folding to its declared coherent coin vector is an explicit memory convention.
These identities determine a conditional coherent implementation given the
stated register services. Classical measurements alone do not establish its
coherence or relative column phases.

All source constants derived here are operational calibration inputs to the
later consumer. The derived mean cannot count as an independent mass/phase
prediction of itself. Native versus always-swap probe control can alter P
with the same triangle and reader. Conservation, integrality and elementary
normalization identities receive no emergence credit. Further physical use,
complete quantum instruments, clocks and implementation error belong to the
separately registered common-source protocol.
