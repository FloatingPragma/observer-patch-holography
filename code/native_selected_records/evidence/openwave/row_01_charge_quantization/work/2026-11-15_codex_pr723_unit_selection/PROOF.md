# Dynamical selection and record-only observation

This proof builds on E-33's finite balancing path. The new conclusion combines it with a continuing fair scheduler, quantitative convergence, and the executed observer. Integer occupancy is a supplied alphabet. An electric unit and a topological configuration class are not derived here.

## 1. Local repair and finite state domain

For seam(a,b), d=x_a−x_b and coin c, literal repair reports j=(d−(2c−1)(|d|mod2))/2 and x'=x−j e_a+j e_b. Equivalently, the endpoint pair becomes floor((x_a+x_b)/2),ceil((x_a+x_b)/2), with their order selected by c. This holds for negative sums too. TotalQ is conserved and neither endpoint leaves the previous closed interval[min(x_a,x_b),max(x_a,x_b)]. Thus each finite initial integer vector stays in a finite box. No finite global load cutoff is required.

Let V(x)=sum x_i². A repair decreases V by (d²−r²)/2, where r=|d|mod2. This is a nonnegative even integer; it is at least2 when |d|≥2. For |d|=1 the two coin outcomes are hold and exchange.

## 2. All integer inputs reach adjacent levels

On a connected graph of n vertices, assume max x−min x≥2. Take a shortest path from a maximum M to a vertex of value≤M−2, stopping at the first such vertex. All internal values are M or M−1. Move an M along this path by permitted exchanges with M−1; an encountered M already supplies the moving maximum without an event. The last edge has contrast≥2 and a repair strictly lowers V. This uses at most n−1 events. This is E-33's reachability argument; it handles local gradients that have no immediately large edge contrast.

At fixed Q=n b+k with0≤k<n, the minimum Vmin=(n−k)b²+k(b+1)² is attained exactly at vectors with entries b,b+1. It is a lower bound because a pair of values differing by≥2 can be moved one unit toward each other and reduce sum of squares; iterating reaches this balanced form. The balanced class is closed under every native event. Integer V decreases by at least2 along each constructed path, so finitely many such paths reach it. Within the balanced class adjacent swaps move the k upper-level entries along any connected graph, giving one communicating class.

For the uniform IID scheduler on m seams, each specified edge/coin occurs with probability1/(2m). Let L=n−1, p=(2m)^(-L), K=(V(initial)−Vmin)/2. From every unbalanced state there is probability at leastp of decreasing V within L events; other events cannot increase V. The conditional waiting time for each decrease is dominated in L-event blocks by a geometric(p) law. There are at mostK decreases before balance. Therefore entry is almost sure and E[tau]≤LK/p. This finite bound is deliberately loose; it is a proof of termination, not the measured mean. K is a nonnegative integer because sum x_i² and Q have the same parity. ForK=0, tau=0.

For the actual12-port graph, m=30 and every vertex has degree5. AtQ=2 a double2e_v waits until an incident edge is addressed; both coins split it into two ones. Its exact entry probability per event is5/30=1/6, survival(5/6)^t and mean6 events. Every binary total-two configuration stays binary. Swap/hold only permutes the original vector and hence retains its double forever.

For0≤Q≤n, balanced entries are0/1 (Q=n gives all ones). ForQ outside this range they are b,b+1. In particular all-two is fixed. It would be false to call the latter a derivation of absolute0/1 occupancy for arbitrary density. WhenQ=0, the signed source converges to zero. This is not stable pair creation.

## 3. Operational reading without hidden state

On one probed seam the first current is j0, with coin c0 and initial difference d. Repair leaves difference s r, s=2c0−1. Repeating c0 gives zero current. At the first opposite coin the current is s r. Hence d=2j0+jopp. The reader uses independent fair coins, does not command their values, and retains all repeats. With per-edge capC, failure probability is2^(-(C−1)). Eleven sequential tree probes have completion probability(1−2^(-(C−1)))^11. Every failed probe produces an explicit inconclusive result.

Before probing a tree edge(a,b), let Delta be the accumulated incidence update from earlier probe currents. The initial pre-probe difference is 2j0+jopp−Delta_a+Delta_b. Integrate these eleven differences with r_0=0. Connectivity gives exactly the pre-probe state modulo a uniform shift. The reader takes only the probe tape and cap; neither Q nor the preparation is an input. A uniform shift changes no seam contrast or current law, so no absolute-background measurement is claimed.

The body tape gives each state change. Transport the terminal relative state backwards by those recorded changes. This recovers every body state up to the same constant offset, so max−min, first balanced entry, and postentry escape are directly determined by observer records. Initial states are retained separately only for audit. The same incidence identity applies to the control body. Its terminal probes still run the actual native law, so the same reader remains valid. A native probe can split a control double during measurement; reconstructed body states all refer to times before that intervention.

## 4. Exact evidence and limits

The finite artifact enumerates all nonnegative statesQ0..4 and every native edge/coin event, checks closed binary support and directed reachability, then executes fixed-seed histories from double, quadruple, signed and shifted inputs with complete operational readings. The independent checker uses floor/ceiling allocation rather than the producer current formula. Finite sampled histories supplement the proof; their horizon does not prove all-input convergence.

No part of this theorem turns supplied integer loads into an emergent real-valued charge quantization. It is a positive answer to the narrower dynamical unit-selection route offered in PR723, conditional on acceptance of the new source and observable dictionary. Fractional-state selection and a non-bookkeeping topological charge remain open alternatives.
