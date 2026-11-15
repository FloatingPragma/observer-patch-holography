# Quantum mixing of the native C/E repairs

For independent uniform seam attempts on the existing twelve-port carrier,
both native repairs have quantum spectral gap `(5-sqrt(5))/60`. Their slowest
binary coherences decay by `5/6` per event for C and `7/10` for E. The binary
66-state channel is entanglement breaking after 134 events by an explicit
sufficient certificate. The full 78-state channel converges to the same
binary stationary state but is **not entanglement breaking at any finite
event count**.

These are statements about the source channel after discarding every seam,
coin and repair record. The [recorded source implementation](source_discrimination/DERIVATION.md)
fixes the underlying C/E instruments, including arbitrary untouched
references. Discarding records is a further channel operation, not an
identification with the observer's growing classical tape. The independent
seam schedule here is also distinct from the two-event E/F feedback
controller. Event counts do not supply a physical clock, and added N2/N4
protected quantum registers are outside this source algebra.

## Exact sectors on any weighted connected graph

Let `G` be a finite connected undirected graph with `n >= 2` vertices,
positive edge weights `w_ab`, and `W = sum_edges w_ab`. Each event chooses
edge `ab` with probability `w_ab/W` and a fresh fair placement coin `c`.
Fix an orientation only to label the current. On the binary fixed-load
space

```text
H_k = span{|x>: x in {0,1}^n, sum x = k},
1 <= k <= n-1,                 D = binomial(n,k),
d_x = x_a-x_b,                 j_x(c) = [d_x-(2c-1)|d_x|]/2,
F_c(x) = x-j_x(c)e_a+j_x(c)e_b.
```

C records `(ab,c,j)`; E records `(ab,c,d)`. Write `Phi_C` and `Phi_E` for
their record-discarded channels. A matrix unit `|x><y|` contributes in a
branch exactly when its two records agree. For C this means `j_x=j_y`;
for E it means `d_x=d_y`, which also implies equal currents. In either
case `F_c(x)-F_c(y)=x-y`. Thus signed difference

```text
z = x-y in {-1,0,1}^n
```

is invariant. If `z` has `r` positive and `r` negative entries, the common
occupied set `S` lies in `Z={v:z_v=0}` and has size `k-r`. This labels every
matrix unit in the sector once. Define the positive exclusion Laplacian
on these shared particles by

```text
(L_excl f)(S) = sum_{ab in E, a,b in Z} w_ab [f(S)-f(S^{ab})],

K_C(z) = sum_edges w_ab |z_a-z_b|,
K_E(z) = 2 sum_edges w_ab 1[z_a != z_b].
```

Here `S^{ab}` exchanges occupations of `a,b`; exchanges of equal
occupations contribute zero. Directly comparing the two coin branches gives

```text
Phi_L restricted to z = I - [L_excl(G[Z],k-r)+K_L(z) I]/(2W),   L=C,E.  (1)
```

Indeed, an edge inside `Z` swaps a shared particle with probability one
half. On all other edges the matrix unit either stays put or its records
disagree: the C loss fraction is `|z_a-z_b|/2`, and the E loss fraction is
`1[z_a!=z_b]`. No surviving branch moves a nonzero difference site.
This proves (1), including its normalization, without diagonalizing a
large superoperator.

The induced graph need not be connected. Its constant vector is still a
zero mode, so the largest sector eigenvalue is exactly `1-K_L(z)/(2W)`.
All binary eigenvalues lie in `[0,1]`: a single-edge binary C channel has
Hilbert--Schmidt eigenvalues among `0,1/2,1`, and E among `0,1`. This follows
by separating unchanged matrix units, killed matrix units, and the two
swapped shared-particle units. Their weighted average is therefore a
self-adjoint positive semidefinite contraction. In particular the sector
maximum is also its spectral radius.

## The full binary quantum gap equals the classical gap

The discrete-channel gap here means one minus the largest modulus of a
nonstationary eigenvalue.
Let `lambda_1(G)>0` be the first positive eigenvalue of the weighted graph
Laplacian, with quadratic form `sum_edges w_ab(f_a-f_b)^2`. The diagonal
sector `z=0` in (1) is the symmetric `k`-particle exclusion chain. Its
Laplacian gap equals `lambda_1(G)` by the weighted exclusion consequence
of the Aldous gap theorem: Caputo, Liggett and Richthammer,
[Section 4.1.1, equation (4.1)](https://www.stat.berkeley.edu/~aldous/260-FMIE/Papers/caputo.pdf).

For every nonzero sector, `sum z=0`, `sum z_v^2=2r`, and

```text
K_E(z) >= K_C(z)
       >= (1/2) sum_edges w_ab(z_a-z_b)^2
       >= (lambda_1(G)/2) sum_v z_v^2
        = r lambda_1(G) >= lambda_1(G).                         (2)
```

The first inequality uses `|z_a-z_b|<=2`; the second uses
`|h|>=h^2/2` for `h in {0,+/-1,+/-2}`. Thus no off-diagonal sector has a
slower mode than the diagonal exclusion gap. Equality is already attained
in the diagonal sector. Consequently

```text
gamma_C = gamma_E = lambda_1(G)/(2W),
tau = I_D/D,                  Fix(Phi_L) = span{I_D},
||Phi_L^t - Dep_D||_{2->2} = (1-gamma)^t,
Dep_D(X) = Tr(X) I_D/D.                                          (3)
```

Both channels are unital on this binary space. Their fixed observable
algebra is just the scalars, and they have no other eigenvalues on the
unit circle. This rules out an unchanging nontrivial quantum subsystem of
this forgotten-source channel.

Equation (3) gives a bound uniform over every source input and every
untouched reference `R`. With trace distance `dist(rho,sigma)=||rho-sigma||_1/2`,

```text
dist((Phi_L^t tensor id_R)(rho), tau tensor rho_R)
    <= min{1, (D/2)(1-gamma)^t}.                                (4)
```

For a pure input, its Schmidt rank is at most `D`. The output difference
therefore acts on a space of dimension at most `D^2`, where
`||Y||_1 <= D||Y||_2`. Tensoring the Hilbert--Schmidt operator in (3) with
the identity preserves its `2->2` norm, and a pure density matrix has
Hilbert--Schmidt norm one. This proves (4) for pure inputs; convexity proves
it for mixtures, with their respective reference marginals. No bound on
the reference's original dimension is required.

## The retained twelve-port carrier

For its thirty unweighted seams, `n=12`, `W=30`, all vertex degrees are
five, and the [existing carrier spectrum](../common_source/native_fields.py)
has `lambda_1=5-sqrt(5)`. At total load two the binary dimension is `D=66`.
Write

```text
gamma = (5-sqrt(5))/60,          r_mix = 1-gamma.
```

The sharp off-diagonal rates follow from (1). If `r=1`, with opposite
difference sites `u,v`, then `K_C=10` regardless of adjacency, while
`K_E=2(10-1[uv is an edge])`. For `r=2`, write `e_pp`, `e_mm` and `e_pm`
for edges within the two positive sites, within the two negative sites,
and between these pairs. Then

```text
K_C = 20-2(e_pp+e_mm) >= 16,
K_E = 2[20-2(e_pp+e_mm)-e_pm] >= 24.
```

The bounds use `e_pp<=1`, `e_mm<=1`, `e_pm<=4`. Therefore the off-diagonal
spectral radii are exactly

```text
r_C = 5/6,                     r_E = 7/10.                       (5)
```

These rates have an operational witness. Choose adjacent `u,v`, write
`|uw>` for the binary configuration occupied at `u,w`, and set

```text
A_uv = sum_{w not in {u,v}} |uw><vw|,
rho_+/- = (1/20) sum_{w not in {u,v}}
             (|uw> +/- |vw>)(<uw| +/- <vw|).
```

The constant sector vector `A_uv` has eigenvalue `r_L`. The two density
matrices have orthogonal supports and identical load populations, so

```text
dist(Phi_L^t(rho_+), Phi_L^t(rho_-)) = r_L^t,
dist(Phi_C^t(rho_+), Phi_E^t(rho_+)) = (r_C^t-r_E^t)/2.            (6)
```

The first integer event count reducing the first distinguishability to
at most one half is four for C and two for E. The second distance is
`1/15` at one event and peaks at four events, at `12259/101250`.
For the latter, the increment from `t` to `t+1` changes sign when
`(25/21)^t` crosses `9/5`, and `(25/21)^3 < 9/5 < (25/21)^4`.
The binary effect `(I+A_uv+A_uv*)/2` attains this C/E probability
difference: `A_uv+A_uv*` has eigenvalues `-1,0,1`, so the effect is valid.
This is an analytic apparatus-level preparation and measurement, without
a new executed source pulse word. These numbers concern the uniform seam
schedule, not the existing discriminator's fixed-seam experiment.

All load populations, classical load histories and diagonal geometry
observables obey the same exclusion process under both candidates when
only these native events and load observations are used. Off-diagonal
sectors never feed populations, even for coherent initial inputs. Such
classical geometry observations cannot choose between C and E; (6) uses
coherent preparation and readout. Their common overall gap therefore
does not imply operational equality of their quantum channels.

## A finite binary entanglement-breaking certificate

Let `J_t=(id tensor Phi_L^t)(|Omega_D><Omega_D|)` be the normalized Choi
state, with `|Omega_D>=D^(-1/2) sum_i |i,i>`. For an orthonormal
Hilbert--Schmidt eigenbasis and all nonconstant eigenvalues `lambda_i`,

```text
||J_t-I_(D^2)/D^2||_2^2 = D^(-2) sum_{i>0} lambda_i^(2t).         (7)
```

This follows from the invariance of Hilbert--Schmidt norm under Choi
reshuffling; the identity eigenoperator gives the maximally mixed term.
[Gurvits and Barnum, Corollary 3](https://arxiv.org/pdf/quant-ph/0204159)
give separability whenever a normalized bipartite state of total dimension
`d` is within squared Hilbert--Schmidt distance `1/[d(d-1)]` of `I_d/d`.
Taking `d=D^2`, (7) therefore proves entanglement breaking if

```text
sum_{i>0} lambda_i^(2t) <= 1/(D^2-1).                            (8)
```

Here we use the equivalence between separable Choi state and
entanglement-breaking channel, established in
[Horodecki, Shor and Ruskai, Theorem 4](https://arxiv.org/abs/quant-ph/0302031).
Entanglement breaking means separable output with every untouched
reference, including inputs not used in any finite test.

There are `D-1=65` nonconstant diagonal modes and `D(D-1)=4290`
off-diagonal modes. Equations (3) and (5) make the following sufficient:

```text
65 r_mix^(2t) + 4290 r_L^(2t) <= 1/4355.                         (9)
```

For both candidates (9) first holds at `t=134`. The exact rational
enclosure `559/250 < sqrt(5) < 2237/1000` bounds
`14309/15000 < r_mix < 19079/20000`: the upper bound certifies (9) at
134, and the lower bound fails it at 133. Monotonicity settles all other
integers. Treating all 4355 modes with the slower `r_mix` would instead
certify 178 events. These are sufficient-certificate times, not a claim
that 134 is the channels' earliest entanglement-breaking time. Subsequent
channel composition preserves the entanglement-breaking property.

## All 78 states: transient collisions and the same gap

For general `n>=3` at total load two, let `P` project onto the binary
subspace and `Q` onto the collision states `|2u>`. Binary states remain
binary. Touching the occupied site of a collision turns it into a binary
pair, which never becomes a collision again.

Moreover, a surviving matrix unit involving a collision either stays
unchanged or enters the binary operator space. If a collision survives,
the chosen edge did not touch it, so its current is zero. Record equality
then forces the other side to stay unchanged too. Thus, ordering binary
matrix units first, the entire channel has the form

```text
Phi_full = [ Phi_binary   B ],
           [     0       T ]
```

where `T` is diagonal with entries in `[0,mu]`, and

```text
mu = 1-d_min/W,             d_min = min_u sum_v w_uv.
```

The upper endpoint is attained by a collision population at a vertex of
minimum degree. The spectrum is the union of the binary spectrum and
these transient diagonal entries. The Rayleigh vector
`e_u-(1/n)1` gives `lambda_1 <= n d_min/(n-1) <= 2d_min`, so
`mu <= 1-gamma`. The full total-two gap is therefore still exactly
`gamma=lambda_1/(2W)`. Its unique stationary density matrix is
`tau_P=P/binomial(n,2)`, and its Heisenberg fixed observables are scalars.
For the twelve-port carrier, `mu=5/6`. The full channel is not unital or
Hilbert--Schmidt self-adjoint; the binary norm argument must not simply
be applied to dimension 78.

There is nevertheless a reference-uniform finite forgetting bound.
After `s` events the collision probability satisfies `q_s<=mu^s`.
For the resulting joint state `sigma`, replace its collision part by

```text
R(sigma) = P sigma P + tau_P tensor Tr_source(Q sigma Q).
```

This preserves the reference marginal. Positivity of the two-by-two
`P,Q` operator block gives `||P sigma Q||_1<=sqrt(q_s(1-q_s))`.
Thus `dist(sigma,R(sigma))<=sqrt(q_s(1-q_s))+q_s`.
Contractivity for the next `t-s` events followed by (4) proves, for any
integer `0<=s<=t`,

```text
dist((Phi_full^t tensor id)(rho), tau_P tensor rho_R)
 <= min{1, mu^(s/2)+mu^s+(D/2)(1-gamma)^(t-s)},
D = binomial(n,2).                                               (10)
```

For example `s=floor(t/2)` proves uniform convergence. This simple bound
is not asserted to be sharp.

## The full carrier is never entanglement breaking in finite time

On the twelve-port carrier choose adjacent `u,v`, set
`alpha=|2u>`, `beta=|uv>`, and entangle with a reference qubit:

```text
|Psi> = (|alpha>|0> + |beta>|1>)/sqrt(2).
```

The transient diagonal coefficient of `|alpha><beta|` is

```text
a_C = [25-4/2]/30 = 23/30,       a_E = [25-4]/30 = 7/10.          (11)
```

There are 25 edges avoiding `u`; of these, four touch `v`. C retains one
coin on each of those four, while E retains neither. All other edges
avoiding `u` retain both coins. An edge touching `u` cannot leave this
matrix unit unchanged. By the triangular structure, after `t` events its
coefficient is exactly `a_L^t`, with no return from the binary space.

The output on the `|1>` reference branch comes from a binary input and
has exactly zero collision population. Partial transpose on the reference
therefore has the principal block on `|alpha,1>,|beta,0>`

```text
[     0       a_L^t/2 ],
[ a_L^t/2       b_t   ],          b_t >= 0.
```

Its determinant is `-a_L^(2t)/4<0` for every finite `t`. Separable states
have positive partial transpose, as proved by
[Peres](https://arxiv.org/abs/quant-ph/9604005), so this is an entangled
output at every finite time. The full 78-state powers are consequently
never entanglement breaking. This is consistent with (10): the surviving
entanglement becomes arbitrarily small, and the limiting replacement
state has support only on the binary subspace. It does not provide a
stable quantum memory or a protected subsystem.

## Executable correspondence and scope

The existing source-discrimination executor uses the endpoint seam `(0,1)`.
The native channel on another seam is its exact port relabelling; this
calculation does not claim thirty separately executed source pulse tapes.

[quantum_mixing.py](quantum_mixing.py) assembles every signed-difference
block directly from the native primitive and actual record equality. It
separately constructs (1), checks exact equality on all 4356 binary matrix
units, and checks the diagonal transient structure on all 6084 full
matrix units. [quantum_memory.py](quantum_memory.py) checks the rational
entanglement-breaking bound and the finite-time witness coefficients.
[Small weighted-graph controls](test_quantum_mixing_graphs.py) independently
compare floor/ceiling Kraus outputs with (1) for every nontrivial load on
four graphs of two through five vertices. A disconnected two-edge control
has an exact fixed coherence, exhibiting why connectedness is required.
The [retained receipt](results/quantum_mixing.json) binds these calculations
to the existing thirty-seam carrier. Exhaustive checks establish that
the executable channel is the one in the proof; they do not replace the
weighted exclusion theorem or the separable-ball theorem.

The result settles the quantum gap, eventual entanglement breaking on the
binary shell and its failure on the full carrier for these two specified
forgotten-source channels; the exact shell entanglement-breaking index is
not determined. It does not discard records from the actual observer,
derive coherent code services, calibrate the event schedule in physical
time, or decide the full OPH programme's realizability.
