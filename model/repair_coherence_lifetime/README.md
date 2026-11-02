# Complete-record coherence lifetimes and faithful finite references

This package supplies quantitative rates for the bare binary endpoint repair
instrument classified in `native_fixed_algebra/OCCUPATION.md`. It retains both
C/E laws, every endpoint coin, raw record and CAR sign. On the actual uniform
twelve-port selector, fixed-number coherence contracts by at least `5/6` per
repair tick; the bound outside the extremal cat block is `11/12`. These are
repair rates, with no identification as Hamiltonian energies or particle masses.

The theorem also controls classical record-adaptive selection with a uniform
cut bound. It does not exclude quantum memories, coherent recovery, a changed
repair law, or every OPH refinement. It supplies finite inputs for issues
#1138 and #1139 and a specific refinement control relevant to #1142.

## Exact endpoint identity

Let `G=(V,E)` be a finite connected graph, `q_e>0`, `sum q_e=1`, with fair
endpoint coins. For two binary configurations `S,T`, write
`D_x=1_S(x)-1_T(x)`. A Kraus branch is a signed partial injection. For a
matrix unit `|S><T|`, it is nonzero exactly when its two columns give the same
coin and raw record. Its absolute coefficient is then `q_e/2`, independently
of the CAR sign. The four local configurations, ordered `00,10,01,11`, have
the following common-record probabilities under C:

```
  1   1/2 1/2 1
  1/2 1   0   1/2
  1/2 0   1   1/2
  1   1/2 1/2 1
```

Consequently the complete-record coefficient mass is exactly

```
r_C(S,T) = 1 - (1/2) sum_{ab} q_ab |D_a-D_b|.                 (1)
```

E additionally measures local occupancy parity, so its mass is
`r_E=1-sum_{ab}q_ab 1[D_a != D_b] <= r_C`. Reserved abort slots have zero
support. The code checks the source's full column inventory and the local
record-overlap table; no numerical rank or selected trajectory enters (1).

Let `lambda(q)` be the minimum weight of a nonempty proper graph cut. With
`A={D=1}` and `B={D=-1}`, the exact discrete coarea identity is

```
sum q_ab |D_a-D_b| = q(boundary A) + q(boundary B).            (2)
```

Except for identical configurations and the full/empty pair, at least one of
these is a nonempty proper cut. Thus `r<=1-lambda/2`. If `|S|=|T|` and `S!=T`,
both sets are nonempty and proper, giving `r<=1-lambda`.

Each branch preserves particle number and is injective on its support. It
therefore sends an off-diagonal matrix unit to another off-diagonal unit and
cannot move a nonextremal pair into the full/empty block. Define the coherence
norm as the sum of the absolute off-diagonal entries, omitting that block.
Triangle inequalities and (1) give contraction of this norm. An untouched
reference is allowed by replacing each absolute entry with its trace norm.
Distinct classical records are separate blocks and their norms are summed.
Thus retaining the complete classical tape cannot evade this bound.
The blocks are unnormalized, so this is a probability-weighted bound.
Conditioning on a rare record can preserve normalized coherence; it cannot
remove the record's probability cost or justify a per-trajectory decay claim.

For record-adaptive `q_{e|h}` and repair probabilities `p_h`, the same argument
applies separately to every history block. If
`p_h lambda(q_h)>=ell` at every step, contraction is at most `1-ell` within a
fixed-number sector and `1-ell/2` outside the extremal block. The bound permits
record-controlled phase changes and occupation permutations that preserve
the stated sectors. It does not cover coherent operations that create fresh
source coherence. Neither a fixed joint stationary reference nor autonomous
clock selection follows from this adaptive statement.

## Actual finite graph and source checks

All 4094 proper cuts of the thirty-edge native graph are enumerated with exact
rational weights. The minimum is `5/30=1/6`. A second enumeration of all
`3^12=531441` ternary difference fields independently attains both coarea
floors. Both signs and both laws are checked on all `4096` input addresses,
thirty seams and two coins: `983040` native columns.

A two-configuration equal superposition has initial coherence norm one.
After 26 full repair ticks its retained fixed-number coherence is bounded by
`(5/6)^26 < 0.01`. With repair probability `p`, use `(1-p/6)^t`, not the full
repair rate. The extremal full/empty cat remains an explicit zero-decay
countercontrol; it is outside the fixed-total-number observer class and has
the one-port fragility proved in the occupation package.

For distinct one-particle sites `x,y`, every common-record surviving branch
holds both configurations. Hence under uniform C repair

```
Phi_*(|x><y|) = (5/6) |x><y|.                               (3)
```

This is an attained decay rate, not a bound inferred from a small eigenvalue.

## Faithful finite-sector convergence

For a fixed selector, diagonal entries follow symmetric exclusion: an
occupied-empty edge swaps with probability `q_e/2`. Each edge heat bath is a
self-adjoint positive projection on diagonal functions; their mixture has
spectrum in `[0,1]`. The uniform law in each fixed-number sector is invariant.

We explicitly import the proved weighted-graph spectral-gap theorem of
Caputo, Liggett and Richthammer, *Proof of Aldous' spectral gap conjecture*,
Theorem 1.1 and its exclusion corollary (4.1), printed page 19:
<https://arxiv.org/abs/0906.1238>. Its premises here are a finite connected
undirected graph, symmetric positive edge rates `q_e/2`, and fixed particle
number `1,...,|V|-1`. The exclusion process is a quotient of interchange; the
one-particle eigenfunctions lift by summing over occupied sites. These give
both spectral inequalities and therefore equality of the gaps. The imported
theorem is not claimed as a new proof in this repository.

Exact integer arithmetic gives the native adjacency characteristic polynomial

```
(x-5)(x+1)^5(x^2-5)^3.
```

The exclusion gap in every nontrivial number sector is consequently
`g=(5-sqrt(5))/60 > 11/240`. Here `Phi_*` is the nonselective source channel:
the classical tape is traced out for the stationary-state assertion. The
joint source and tape need not forget initial populations. With
`d=binomial(12,N)`, the quantum source state obeys

```
||Phi_*^t(rho) - I_N/d||_1
  <= sqrt(d-1)(1-g)^t + (d-1)(5/6)^t.                       (4)
```

The diagonal term follows from the reversible `L2` bound and Cauchy--Schwarz;
the off-diagonal term uses positivity,
`sum_{i!=j}|rho_ij| <= (sum_i sqrt(rho_ii))^2-1 <= d-1`.
For an untouched reference, apply the classical transition-column bound to
the positive diagonal reference blocks and the block trace-norm version of
the coherence inequality. The limit is `I_N/d tensor rho_reference`, with
the same dimension-only bound. The calculator conservatively replaces
`sqrt(d-1)` with `d`, uses `11/240`, and caps the result at the trivial trace
distance two. The state is faithful on its declared sector. It does not
select a distribution across number sectors or the centered Hubbard ground.

## A specified diffusive refinement cannot retain one-particle coherence

Consider a periodic cubic nearest-neighbor graph with `V=L^3`, `L>=3`, lattice
spacing `a`, degree six and `3V` edges. Under uniform C repair a one-particle
off-diagonal unit has exact multiplier `1-2/V`, by the same holding argument
as (3). The diagonal generator with a Poisson event clock of rate `nu` is
`-nu L_graph/(6V)`. If its physical diffusion coefficient is held fixed at
`D=nu a^2/(6V)>0`, then

```
Phi_t(|x><y|) = exp(-2 nu t/V)|x><y|
             = exp(-12 D t/a^2)|x><y|,   x != y.             (5)
```

Here `Phi_t` traces out the classical tape. Let `Pi` pinch onto the position
diagonal, acting trivially on an untouched reference. Equation (5) gives the uniform bound
`||(Phi_t tensor id)(rho)-(Phi_t tensor id)(Pi rho)||_1
<=2 exp(-12Dt/a^2)` for every source/reference density. Thus the conclusion
does not rely on following one vanishing matrix entry as the volume grows:
the reduced one-particle output loses its initial position coherence.

The lattice diagonal diffusion symbol approaches `-D|p|^2` because
`2 sum_j(1-cos(a p_j))/a^2 -> |p|^2`. In contrast, (5) vanishes at every
positive macroscopic time as `a->0`. This specific bare, classically recorded
repair scaling cannot supply reversible one-particle propagation on the same
coherent matrix units. This assertion does not address a coherently retained
environment, a recovery code, another clock scaling, or the separately
constructed Hamiltonian Weyl fields. It is not a whole-program obstruction.

## Reproduction

```
python -m pytest -q -W error model/repair_coherence_lifetime
python -m model.repair_coherence_lifetime --write
python -m model.repair_coherence_lifetime --check
```

The new receipt binds actual source bytes through the existing many-body
source closure. It does not assert that every older release receipt has been
regenerated. The endpoint table, cut witnesses and rate calculator are
recomputed independently of those older numerical receipts.
