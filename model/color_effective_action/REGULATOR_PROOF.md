# A uniform regulator limit for the implemented quantum rotor

The [quantum-flux source](../interacting_extension/chiral_lattice/rotor_response.py)
implements an integer-flux cutoff K on a closed L-edge graph, for L>=3.
Its U(1) coordinate completion has its own Haar rotor and supplied
electric, hopping and loop coefficients. It is a different configuration
from the native sixty-state raising-link Higgs source. The exact relation
to that native source is the conditional coordinate attachment in
[the quantum-link construction](../particle_gauge_band/QUANTUM_LINK_PROOF.md).
No reducing-subspace identification with the full nonabelian rotor is used.

This construction removes K at the specified coefficients and clock,
retaining the complete local resolvent and its effective action. The error
bound is independent of the largest electric eigenvalue and of L. This
avoids a global Hamiltonian-norm bound that diverges with the regulator.

## The source family and its limiting operator

For each source charge q in {-4,-3,1,2,6}, the exact Gauss basis is
|x,ell>, with 0<=x<L and ell an integer. The endpoint electric fields are

    E_e(x,ell)=ell+q 1_(e<x).

They obey every Gauss constraint including the static charge at x=0.
At the source defaults the operator is

    H=D+V,
    D(x,ell)=[x(ell+q)²+(L-x)ell²]/2,
    V=-(T+T*)-(W+W*)/2,                                    (1)

where T sends (x,ell) to ((x+1) mod L,ell+q 1_(x=L-1)) and W increases
ell by one. On the full countable space both T and W are unitary.
Thus D>=0, ||V||<=3 and H is selfadjoint on Dom D with H>=-3.
The finite source is exactly P_K H P_K on |ell|<=K, including its flux
boundary. `source_replay` compares every finite matrix column, initial
state, electric label and Gauss row with the literal source implementation.

For mu>3, R_0=(mu+D)^(-1) has norm at most 1/mu. The convergent series

    (mu+H)^(-1)=sum_(n=0)^infinity (-R_0 V)^n R_0             (2)

defines the actual source resolvent. A truncation through n=N has operator
error at most

    r_N=(3/mu)^(N+1)/(mu-3).                                (3)

Every finite term is evaluated exactly on finitely supported source inputs.
The diagonal factors preserve support. Each V insertion changes loop flux
by at most s=max(1,|q|).

## Uniform comparison with the finite regulator

Let P_in select all positions with |ell|<=K0<=K and put
r=ceil((K-K0+1)/s). A term of degree less than r cannot leave the
retained cutoff. A term beginning and ending in P_in must make at least
2r insertions to leave it and return. The corresponding terms of the full
and finite series are identical. Bounding both remaining tails gives

    ||[(mu+H)^(-1)-(mu+H_K)^(-1)P_K] P_in||
      <= 2(3/mu)^r/(mu-3),

    ||P_in[(mu+H)^(-1)-(mu+H_K)^(-1)]P_in||
      <= 2(3/mu)^(2r)/(mu-3).                               (4)

The finite inverse in the second line is understood on the retained
subspace. These bounds hold with arbitrary untouched references and are
uniform in every larger cutoff. They also hold at every L>=3, since the
bounded part and its flux step have the same norms. They imply strong
resolvent convergence on the full countable Hilbert space: finite-support
vectors form a dense set, and all resolvents have the common bound
1/(mu-3). This fixes the source limit through the actual operator (1).

There is no assertion that its dimension, state or coupling is a physical
three-dimensional color field. K removes a gauge-coordinate regulator on
the specified finite spatial graph. A limit in spatial lattice spacing
has separate hypotheses.

## A convergent complete coarse action

For the four-edge source let P select all four positions at loop flux zero.
The complete local resolvent G=P(mu+H)^(-1)P is positive. The same Schur
minimization as in the [finite action proof](PROOF.md) gives

    K_eff(mu)=G^(-1),
    G >= g_* I,       g_*=1/(mu+3q²/2+3).                   (5)

Indeed PHP has diagonal at most 3q²/2 and bounded interaction norm at most
three. The Schur action is no larger than mu+PHP; inversion proves the
floor. This uses all discarded flux sectors through their resolvent.

The exact rational series (2), compressed to these four inputs, gives
G_N with ||G_N-G||<=r_N. When r_N<g_*,

    ||G_N^(-1)-K_eff(mu)||
      <= r_N/[g_*(g_*-r_N)].                                (6)

The receipt gives all sixteen action entries and resolvent entries for
each native charge at mu=8 and N=24, rounded outward after exact rational
arithmetic. Matrix inversion is only four-dimensional and exact. No
finite-cutoff eigenvalue fit, suppression of boundary outcomes or guessed
running parameter enters this regulator limit.

Complete energy measurements on coherent superpositions of these four
coarse basis vectors give the same resolvent by the score construction in
the finite action proof. Because H>=-3, (mu-3)/(mu+E) is a bounded score
in [0,1], including every energy outcome. The exact source matrices also
permit direct specification readout. Preparation, energy measurement and
coherent-control implementation remain the declared experiment interfaces.

The conclusion is a controlled action correspondence for the implemented
U(1) rotor completion at its fixed graph and clock. Its additional flux
operators are retained in the full four-by-four action, so no closure onto
one electric/magnetic ratio is presumed. This proves neither equality to
the original raising-link Higgs dynamics nor a nonabelian QCD matching
constant. A spatially refining nonabelian source, its state and the
laboratory scheme dictionary are not supplied by this rotor limit.
