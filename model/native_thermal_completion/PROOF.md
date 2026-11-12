# Minimal reversible support and a complete cooling experiment

The state space and integer cost are those of the exact native closed
three-port/two-face cell. A state is (x0,x1,x2,k), with loads in 0,...,5
and k modulo six. There are 1,296 states. The unchanged native kernel K
and cost E=5 sum(x_i²)+2[C(k)+C(-k)] are imported literally. The construction
below is a separately declared bath law; it does not describe native
repair as a thermal process.

## Reverse support and detailed balance

For distinct x,y, set s_xy=max(K_xy,K_yx). This is the smallest symmetric
support containing every native nonidentity transition. Positive full-support
detailed balance requires a reverse edge for each forward edge, so retaining
the native support while demanding reversibility forces these additions.
It does not uniquely fix their rates. The chosen symmetric proposal and
acceptance rule are explicit design choices.

There are 4,104 native directed nonidentity edges and 3,924 added reverse
edges. The resulting undirected support has 26 connected components:
one complete conserved-load sector for each load sum 1,...,14, and twelve
singleton flux states at total loads zero and fifteen. No transition changes
the total load. Cost differences along any edge are even integers.

For 0<q<=1 put beta=-log(q)/2 and Lambda=2 max_x sum_(y!=x)s_xy=12.
At x, propose y with probability s_xy/Lambda and accept with probability
a_xy=min(1,q^((E_y-E_x)/2)). Retain the accepted record, the rejected
record returning to x, and the idle record of probability
1-sum_y s_xy/Lambda. These define 17,352 original records per attempt,
including zero-weight rejection records at q=1. Every row sums to one
exactly; its total self-loop probability is at least one half.

On every edge,

```
P_xy/P_yx = q^((E_y-E_x)/2) = exp[-beta(E_y-E_x)].
```

Consequently the Gibbs law proportional to exp(-beta E) on each connected
component satisfies detailed balance and stationarity. Each nonsingleton
component is irreducible and aperiodic, hence has this unique conditional
equilibrium. Component weights are conserved and free. They are not selected
by the Gibbs formula. The executed experiment fixes them to the weights of
the native independent-uniform-load, zero-flux preparation. Relative energy
floors are removed inside each component before taking exact rational powers
of q; this does not change its normalized law.

The producer checks all edges and all stationarity equations in rational
arithmetic at q=1,3/4,1/2,1/4. Its receipt contains every proposal edge,
cost difference and idle weight, sufficient to reconstruct all accepted and
rejected probabilities exactly. No large lattice or cosmic clock is used.

## Conjugation and thermodynamic records

Load complementation x_i->5-x_i and k->-k sends E to
E+375-50M, where M is the conserved load sum. Thus edge cost differences,
symmetric proposals and acceptance probabilities are equivariant under the
same native conjugation. The existing finite Markov conjugation theorem
applies to this explicitly verified kernel. A symmetric preparation has
zero mean conjugation-odd class imbalance at every attempt. These classes
are not baryon or lepton numbers.

For a fixed beta during one attempt, let p'=pP and Q=sum_xy p_x P_xy(E_y-E_x)
be heat into the system. With natural-log entropy S, detailed balance gives

```
S(p')-S(p)-beta Q
 = sum_xy p_x P_xy log[p_x P_xy/(p'_y P_yx)] >= 0.
```

The right side is relative entropy between two normalized joint laws.
The identity includes rejection and idle events. It is valid on the active
components even when the other component weights vanish. For a sequence of
temperature settings the identity is applied at each actual bath setting;
the system cost is fixed, so no change of its Hamiltonian is hidden as heat.
External bath-control resources are not assigned a physical work value.

For the fully recorded law, reverse an accepted x-to-y event to the accepted
y-to-x event. A rejected proposal retains its proposal label at the same
state, and an idle mark remains idle. Conditional on a self-loop, the
relative rejection/idle probabilities agree in the two joint laws. Thus
retaining these marks gives the same relative-entropy identity as summing
them into the diagonal transition probability.

## Executed cooling and error budget

The initial law is the declared q=1 componentwise equilibrium with the
specified component weights. The protocol uses q=3/4,1/2,1/4,1/2,3/4,1,
with either 32 or 256 attempts at each setting. These two fixed schedules
retain complete state laws, class shares, mean cost, heat, entropy change
and relative-entropy production. Finite-time populations are compared with
the exact component Gibbs law; they are not replaced by equilibrium values.

Exact transition probabilities are converted once to binary64. Let eta be
the maximum row sum of the exact coefficient-rounding errors and d the
maximum number of incoming terms. Under ordinary binary64 products and
sums without overflow or underflow, gamma=2du/(1-2du), u=2^-53, bounds the
positive sparse multiplication error. Put delta=eta+gamma(1+eta). If the
input has L1 error e, one step has error at most (1+e)(1+delta)-1.
The code computes eta and the initial error as exact fractions and propagates
this bound through every stage, rounding the bound upward to a multiple
of 2^-80 after each stage. It verifies the non-underflow domain and
never renormalizes the calculated probabilities. Entropy and heat displays
are numerical evaluations; their sign law is the exact identity above.

The native attempt count, declared cost and bath beta have no physical
second, GeV or temperature calibration. This construction supplies a
thermal completion and explicit finite cooling records. It does not prove
that the native repair law, a Standard Model gas or the early Universe
instantiates this added bath.
