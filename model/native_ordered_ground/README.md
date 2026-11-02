# Native preparation of the selected ordered spin ground

Every finite bipartite Heisenberg sample with a positive rational staggered
selection field has a unique ground vector that can be computed and prepared
to any requested precision using the existing native port interface. This
removes a supplied-amplitude assumption from the finite inputs to FA92 and
FA93. An eight-site example is actually executed, including its spin pointer,
on the full native success and failure instrument. The thermodynamic statement
uses ordered iterated limits; no effective convergence rate to a specified
infinite-volume marginal is claimed.

## Exact finite Hamiltonian and a terminating ground algorithm

On a finite bipartite graph with N vertices and E edges (multiplicities counted),
let J,B be positive rationals and eps_x be its bipartite signs. The specified
selection Hamiltonian is

```
H_B = J sum_edges (S_x.S_y-1/4) - B sum_x eps_x S_x^1.
```

Conjugation by the product of Z on the odd sublattice gives a real matrix H'
with diagonal -J/2 per unlike-spin edge, off-diagonal -J/2 for its swap, and
off-diagonal -B/2 for each individual bit flip. The common logical spin frame
of FA93 has (Xbar,Ybar,Zbar)=(2Sx,-2Sy,-2Sz); it preserves the exchange and
the X selection field. The Marshall transformation is explicitly undone in
the preparation. Arbitrary nonbipartite graphs are rejected.

Put s=JE+BN, delta=B/(2s), alpha=delta^N, D=2^N, and M=I-H'/s.
Every diagonal of M is at least one, every single-flip entry is delta, and
every row sum is at most two. Also -s<=H'<=BN/2, so M is positive definite.
Hypercube paths padded with diagonal steps prove (M^N)_ij>=alpha for every
i,j. Perron--Frobenius gives a unique strictly positive normalized vector v,
with eigenvalue 1<=lambda<=2. Its eigen-equation gives
v_min/v_max>=alpha/2^N. The stochastic h-transform of M^N consequently has
every entry at least alpha^2/4^N. Its Doeblin constant is

```
kappa = alpha^2/2^N,
gap(H_B) >= s*kappa/N,
|| normalized(M^(mN) 1)-v || <= 2 sqrt(D) (1-kappa)^m.       (1)
```

For the last inequality the normalized all-ones vector overlaps v by at
least 1/sqrt(D); normalizing the power iterate costs at most twice the
orthogonal/parallel ratio. All other eigenvalues of M are positive, so the
Doeblin contraction bounds their ratios directly. Concavity of the Nth root
gives the displayed Hamiltonian gap. The conservative rational choice
m>=ceil(2D/(kappa*epsilon)) guarantees vector error<=epsilon, since
(1-kappa)^m<=1/(1+m*kappa). Integer matrix multiplication and gcd reduction
therefore terminate without a ground-state oracle or an unproved gap.
This is a computability statement, not an efficient algorithm.

## Rational angle construction and the native interface

For positive integer power-iterate amplitudes, each binary tree node needs
theta=-atan2(sqrt(b),sqrt(a)), where a,b are exact subtree squared norms.
After normalizing a+b=1, write

```
x = sqrt(b)/(1+sqrt(a)),
t = x/(1+sqrt(1+x*x)),
theta = -4 atan(t),             0<=t<=sqrt(2)-1<1/2.
```

`angles.py` encloses each square root by integer square roots and atan by its
alternating rational series, and refines until a dyadic angle has error at
most 2^-b. The exact rational Walsh transform supplies the Gray-code
uniformly controlled rotations. A level has disjoint conditional blocks,
so its error is the maximum node error; N levels cost at most N*2^-b.
The tree uses 4(2^N-1)-3 pulses. Marshall phases use two per odd site and
one-spin encoding uses five per site. The optional X pointer uses 25.

FA93's explicitly specified dyadic integrated-port controller implements
each ideal pulse to isometry error<=29*2^-p. It derives its commands from
the exact twelve six-level source generators, including the pi coefficients
in the CNOT and Marshall words. Thus for P pulses the entire prepared pure
state has vector/isometry error

```
epsilon_ground + N*2^-b + 29P*2^-p.                           (2)
```

The unhalved trace error is at most twice (2), including complement exports.
Blank registers, load/export service, the port precision family and the
specified schedule remain explicit apparatus inputs. No added physical
noise/memory law, unlimited binary64 precision, or autonomous stationary
observer is inferred. Repairs occur at the completed logical operation
boundaries of FA93; an unprotected partial pulse is not asserted invariant.

## Executed eight-site certificate and energy readout

The retained example is the open 2x2x2 cube (12 edges), J=1, B=1/4.
Its untrusted eigensolver proposal is retained as `witness.npz`; replay
does not regenerate it. Every entry is a literal integer times 2^-48.
`certificate.py` checks all 256 eigenvectors by exact rational arithmetic.
If V*V differs from I by norm<=d<1 and HV-VE has Frobenius norm<=r, the
polar unitary Q satisfies ||Q-V||<=d and ||HQ-QE||<=r+2sd. Weyl's theorem
therefore encloses the entire spectrum. The isolated first eigenvalue and
the first column's exact residual then enclose the true ground vector.
Strictly positive amplitudes fix its sign using Perron--Frobenius.

The resulting enclosures in source units are

```
ground energy       [-8.2399123998819,-8.2399123998817]
first excited       [-7.147627722206,-7.147627722184]
finite-field gap    [ 1.092284677676, 1.092284677698]
staggered spin/site [ 0.335475070930, 0.335475070931].
```

These are H_B selection energies, not relativistic particle rest masses.
The 1,090-pulse preparation and X pointer run on all 131,072 logical rows of
the 17-bit work register through the actual six-level generators. Every
first-failure weight is retained. The target pointer-zero probability is
0.8354750709305182; its success contribution is 0.8354750709305093, with total
native first-failure weight about 1.73e-30. The certified complete isometry
error is less than 5e-13. No outcome is selected or success-renormalized.

The finite example uses FA89's checked binary64 angle proposals for speed;
its exact node/rounding bounds, literal port controls and complete ground
certificate prove its error independently. The rational angle algorithm
above supplies the separate unbounded-precision existence construction.

## Joining the ordered state and source dynamics

Take even three-dimensional periodic boxes of side at least 4, so the graph
is bipartite. First take L to infinity at fixed B>0, then B down to zero,
with preparation errors tending to zero. Compactness gives the same
field-selected limit states used in FA92. FA92's DLS/KLS order bound and
Koma's maximal-staggered-magnetization selection apply to this family; the
preparation changes local expectations by at most twice its vector error.
For a local contraction A supported on s_A sites the finite B ground obeys

```
<A*[H_0,A]> >= -B s_A.
```

Indeed H_B is its ground Hamiltonian and the selection-field commutator
has norm<=B s_A. Each site touches at most six bonds, and each exchange
term has norm<=J. For a prepared vector within epsilon of that ground,
the lower bound becomes -B s_A-24J s_A epsilon. FA90's protected source
Hamiltonian differs from H_0 by norm<=3J/E on its selected block, adding
at most 6J/E. All three defects vanish in the stated limits. Finite B
states need not be stationary for the zero-field source Hamiltonian.

FA93 can run finite local histories from each of these prepared states;
its comparison bounds add the input trace error from (2). The result
constructs the finite states used to select the ordered ground. It does
not provide a computable box-size bound for approximating a particular
infinite marginal, nor identify every ground subsequence with the selected
ordered state. FA97's normalized packet limits must retain this distinction.
No massive species, Lorentz dispersion, interacting spatial continuum or
geometric area/stationarity premise follows solely from this preparation.

## Verification and custody

Run `python -W error -m pytest model/native_ordered_ground -q`, then
`python -m model.native_ordered_ground --write` and `--check`. Replay is in
a clean subprocess, verifies the retained proposal, and executes all pulses.
The receipt binds the current executing sources, exact generator/control
data and proof. Prior broad receipts are historical evidence, separately
identified; they are not accepted merely because a new source hash matches.
