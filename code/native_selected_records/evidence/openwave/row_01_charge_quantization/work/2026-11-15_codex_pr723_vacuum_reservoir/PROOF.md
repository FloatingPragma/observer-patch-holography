# Selected units and transport on a native vacuum-reservoir family

## 1. All admitted loads settle to absolute single units

Choose public B≥2 before preparation. Put arbitrary x∈[-B,B]^12 in the first
twelve-port copy and zero on B−1 added copies. Connect corresponding ports in
successive copies. The graph is finite and connected, with N=12B and m=42B−12.
The conserved total obeys |Q|≤12B=N. This bound follows from the declared input
domain; observers do not need to know Q to operate or read the experiment.

Native repair replaces an addressed pair by floor and ceiling of its mean.
It preserves the total and initial min/max box. Its decrease of V=Σx_v² is
(d²−r²)/2, r=|d| mod2. This is nonnegative and at least2 for |d|≥2. As in E-33,
a shortest path from a maximum to a value at least two lower allows exchanges
of adjacent levels, followed by a strict V decrease, in at most N−1 events.
Consequently every state reaches the unique closed balanced class, with entries
floor(Q/N),ceil(Q/N). Under IID uniform seams and fair coins the probability of
such a decreasing word in N−1 attempts is at least (2m)^(-(N−1)). Finiteness of
the initial box and descent give almost-sure finite entry. Explicitly, letting
Q=Nb+k, Vmin=(N−k)b²+k(b+1)² and K=(V(initial)−Vmin)/2,
E[tau]≤(N−1)K(2m)^(N−1). This intentionally loose bound is not a fitted runtime.

Since -N≤Q≤N, the selected class consists of {0,1} when Q≥0 and {-1,0} when
Q≤0. Q=0 selects vacuum, Q=±N selects the uniform ±1 state. Thus this family
handles every admitted initial first-cell vector, including arbitrary signed
loads and all-double input. Increasing B extends the input domain and native
vacuum capacity together. The claim does not cover arbitrary preloaded added
ports. An all-two entire graph has Q=2N and remains fixed, as it should.

Swap/hold only permutes the initial vector. Every main preparation contains a
magnitude B≥2; that magnitude persists indefinitely. The current observer can
distinguish this from selected absolute units on the same graph and readout.
Integrality and conservation alone therefore do not imply the target property.

## 2. Absolute operational readout with a declared vacuum datum

On a probed edge the initial contrast is d=2j_first+j_first_opposite. Repeated
identical coins after the first repair contribute zero current until that first
opposite coin. Correct the measured contrast by the accumulated incidence of
earlier probe records and integrate an actual graph spanning tree. This gives
the terminal pre-probe vector r, equal to the true vector plus one unknown
constant. Reverse every recorded body current to obtain r_initial. The known
INITIAL reference x_12=0 fixes the constant: x_initial=r_initial−r_initial[12].
Propagating the same body currents forward gives every absolute body state.
No assumption that the reference remains empty is used. It may end occupied.

The probability one capped edge probe fails is 2^(-(64−1)); a union bound gives
P(any failed edge)≤(N−1)2^-63. Such failure is explicitly inconclusive. The
reader cannot authenticate a secretly shifted vacuum preparation: that is a
supplied preparation premise, not a dynamical theorem or an inferred charge.
The zero reference fixes the origin only. The disappearance of magnitudes ≥2
is independently produced by native repair and defeated by the control.

## 3. Native transport after selection

In either selected sign sector, a same-level pair produces j=0. A mixed pair
holds or exchanges with probability1/2 each, and an exchange transfers exactly
one unit. Thus the selected process is symmetric exclusion with one edge
attempt per native event. If D(x) graph edges join unequal levels, then
P(j≠0 | x)=D(x)/(2m), and j² equals the nonzero-event indicator in source units.
This statement is about a source-derived state class, not all integer states:
actual larger currents are possible before selection and under swap/hold.

On any connected graph the class with k upper-level ports is irreducible by
adjacent exchanges. Its transition matrix is symmetric; the uniform law on
k-subsets is stationary. For an addressed edge, exactly k(N−k)/(N(N−1)) of
uniform configurations have either specified mixed orientation. A fair native
coin moves in half those cases. Therefore

    P(j=+1)=P(j=-1)=k(N-k)/(2N(N-1)),
    P(j!=0)=E[j^2]=k(N-k)/(N(N-1)).

This is an exact stationary prediction, not an assumption about finite burn-in.
The executable exact census checks it in every available total0,1,2 state of
both graphs, including column sums, from literal native events.

For the held-out trajectories, X_t=1(j_t≠0) has conditional mean D(x_t)/(2m)
after entry. The centered variables have conditional range length one. The
martingale Hoeffding bound is P(|Σ(X_t−p_t)|≥a)≤2exp(-2a²/T), for fixed T.
Preparation and reconstruction do not alter this statement: the controller is
IID and blind to loads, and the terminal reader reconstructs the earlier x_t
exactly from its records. The main interpretation requires the independently
checked post-burn unit-state gate. If any path is not settled, that gate fails
and the packet cannot use this bound as if its hypothesis held.

One source-record magnitude u is calibrated from a fixed earlier window; no
external electric-charge value enters. Later amplitudes and occurrence counts
are distinct tests. A change in u cannot repair a failure of the event-count
prediction. This construction supplies an operational unit-selection component,
not a physical SI charge magnitude or a new topological invariant.
