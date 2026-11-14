# Complete native attempt-cap law

This addendum qualifies both field-feedback executors. The position
semigroups in `FIELD_FEEDBACK.json` and `RECIPROCAL_FEEDBACK.json` describe
uncapped dynamics. Their 96 and 9,216 position outcomes are complete for
those uncapped semigroups. An executor capped at 4,096 scheduled native
events instead has a nonzero additional `cap_exhausted` output. The earlier
position-law calculations must not be identified with its complete
instrument without that outcome and its probability.

## Exact finite instrument

Let Q be the position generator, Lambda(x) the total rate of all scheduled
seams, and A=Q+diag(Lambda). Then A_yx is the rate of native events taking
x to y. Its diagonal includes every native zero-current hold. Consequently
sum_y A_yx=Lambda(x).

For cap C, use count levels0 through C. At levels k<C the diagonal generator
block is -diag(Lambda), and the block from k to k+1 is A. Level C is
absorbing, retaining its terminal position and the cap flag. The initial
state is the configured source at level0. Every column of this augmented
generator sums to zero. Its exponential is therefore a complete positive,
normalized probability law. Conclusive position probabilities sum over
levels0 through C-1; the remaining mass is the cap output.

The fine instrument additionally retains each seam/coin mark and ordered
waiting time. Its event density is exp(-Lambda dt) r_e/2, with a final
no-event survival factor on conclusive outputs. It stops immediately after
event C on cap outputs, exactly as the executors do. A holding event
increments the count just as a moving event does. No normalization on
conclusive paths is used.

`capped_law.py` executes the augmented generators for caps1 and2, both
carriers, and kappa0,1,2. It compares each full count/position distribution
with independent analytic integrals over zero, one and two ordered waiting
times. At kappa0 the surviving count probability also equals the exact
Poisson sum with mean M*T or2M*T. Separate actual cap1/2 source runs are
replayed and retain the failure flag. These small development caps expose
the failure mass instead of relying on a numerically invisible tail.

## A rigorous nonzero bound at cap4096

For a connected unweighted graph, let v=G(delta_x-delta_y) for adjacent x,y.
The discrete maximum principle puts every v(z) between v(x) and v(y),
whose difference is the effective resistance R_xy. The single-edge unit
flow gives R_xy<=1. Centering this range yields

    |U(y,z)-U(x,z)| <= kappa R_xy/2 <= kappa/2.

The finite local-field approximation adds at most2*kappa*epsilon_K to
one-hop energy differences. With beta=4, kappa<=2 and K=512, every incident
seam attempt rate is therefore less than

    exp(2+8 epsilon_K) < 8.

The last inequality is certified by rational bounds e<11/4 and
exp(d)<=1/(1-d), with d=8 epsilon_K<1. Inactive seam rates are exactly one.
Thus Lambda<=288-6+6*8=330 for one copy and Lambda<=660 for two copies.
These bounds concern the exact-arithmetic declared scheduler; actual
floating rates are separately checked with a wide numerical margin.

The scheduled event count up to time T is dominated by a Poisson variable
with mean Lambda_max*T, even though the state-dependent clock includes
holds. At T=1/2 its mean is at most165 or330. For C=4096, the Chernoff bound
is the explicit positive number

    P(cap_exhausted) <= exp(C-mu) (mu/C)^C.

`CAPPED_LAW.json` retains this exact expression and numerical logarithms.
It is extremely small, never zero. Mapping all cap paths to their explicit
flag gives a total-variation difference from the corresponding uncapped
terminal-position law of at most this bound. It also bounds any expectation
difference for a detector with range[0,1]. General floating-library error
remains outside the exact-arithmetic certificate.
