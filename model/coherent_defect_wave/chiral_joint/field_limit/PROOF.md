# Controlled low-energy interaction on the quantum source

The source operator is the complete one-link, two-cell quantum lepton/Higgs
Hamiltonian in [the actual operator package](../two_cell/hamiltonian.py).
Both local Yukawa coefficients and the scalar-number coefficient are one.
The link-number coefficient kappa is a declared parameter. Every one of the
614400 active basis states is retained. The unit-parameter ground energies
are not transferred to another kappa. The periodic path observations below
use the separately specified [full periodic family](../periodic_family/PROOF.md).
No scalar or link operator is replaced by its expectation.

The construction gives a finite interacting effective Hamiltonian, a dressed
physical state map and complete-record error bounds. It does not identify
that Hamiltonian with a relativistic field or a measured particle. The
[existing spin reduction](../../../bipartite_spin_limit/README.md) supplies
a reusable spectral-graph method on a different source. The derivation here
uses the actual quantum link and the full local scalar interaction.

## 1. An exact dressed local parity

For each cell x define

`C_x=(-1)^N_matter,x product_(incident edges e) (-1)^N_link,e`.

Every actual link raising matrix anticommutes with its excitation parity.
A path hop changes matter parity at its two endpoints and link parity on
every traversed edge. An interior vertex sees two changed incident edges;
each endpoint sees one changed edge and one changed matter occupation.
Thus every C_x commutes with every complete path hop and its adjoint.
Local Yukawa terms change matter number by two, and scalar/link energies
are diagonal, so `[H,C_x]=0` on the full source. All these operators commute
with the actual gauge currents. The argument covers the three-axis periodic
path family, including paths across a seam and nonzero background matter.

`star_parity_checks()` checks all 24576 actual open hopping columns and
all 320 complete local matter/scalar columns, with independent path-incidence
controls. This is an operator identity rather than a vacuum-cyclic test.

Let P project onto every link blank. On this subspace C_x is ordinary local
matter parity. Consequently the exact compressed propagator `P exp(-itH) P`,
the exact Feshbach operator where defined, and a symmetry-compatible dressed
low-energy Hamiltonian preserve each site's matter parity. A matrix element
between two blank-link states differing by movement of a single fermion is
zero at every time. A single-particle hopping Hamiltonian cannot describe
that particular compressed domain.

This does not forbid charged motion in the complete quantum source. Its
link registers can retain the transported charge, as in the certified
[propagation experiment](../two_cell/PROPAGATION_PROOF.md). A charged field
must retain a justified dressing and reference instead of discarding that
record. Gauge-neutral odd composites also require examination of this extra
local parity; a global CAR relation alone does not settle their mobility.

## 2. A nontrivial virtual interaction from the actual Hamiltonian

For the one-link source put Q=I-P. Because a primitive raising link moves
only between blank and excited, the exact block form is

`H = [[A,B^dagger],[B,kappa I+D]]`.

A and D contain both full local quantum Yukawa interactions and the scalar
energies. Their norms are at most 18: each local Yukawa has eight unit-norm
terms and each local scalar-number projector has norm one. On the active
five link channels,

`B=sum_(a,b) |a,b><blank| tensor A_ab`,

`A_ab=(i/2) sum_s (sigma_3)_ss c_left,s,a^dagger c_right,s,b`.

Each `||A_ab||<=1`, so `B^dagger B=sum A_ab^dagger A_ab<=5 I` and
`||B||<=sqrt(5)<3`. These bounds hold on every occupation and scalar state.
The link projection and all intermediate states are retained in the bounds.

For real z with `kappa>18+|z|`, exact elimination gives

`F(z)=A-z-B^dagger(kappa I+D-z)^(-1)B`.

Its leading operator is `A-z-B^dagger B/kappa`, with error at most

`5*(18+|z|)/(kappa*(kappa-18-|z|))`.

This is the resolvent identity with the actual full off-diagonal operator.
It does not replace a full residual by a compressed zero.

The virtual term is genuinely interacting. On the four occupation inputs
with the native spin-up e^c mode empty or full at each endpoint and every
other matter mode empty, the diagonal values of B^dagger B are

| Left occupation | Right occupation | Value |
| --- | --- | --- |
| 0 | 0 | 0 |
| 1 | 0 | 0 |
| 0 | 1 | 1/4 |
| 1 | 1 | 0 |

The connected difference for `-B^dagger B/kappa` is `1/(4 kappa)`. It would
vanish for an additive one-body diagonal energy. The Pauli-blocked hopping
term supplies this quartic response. Yukawa terms have not been removed;
they belong to A and D and to the controlled remainder.

The charged test inputs require their retained static compensators to be
physical. Supply at each cell a compensator charge opposite to its e^c
occupation, retain the four settings, and declare zero compensator energy.
The references carry opposite lepton number and spin projection and matching
fermion parity, so the total preparation can preserve those charges.
All source-plus-apparatus Gauss constraints hold. The construction uses one
common source Hamiltonian and its direct-sum reference domain. These are
energy expectation controls, not an identification of charged asymptotic
particles.

## 3. Constructive exact low-energy graph

Use the conservative bounds a=18 and b=3. Let

`g=kappa-2a`, `g>=4b`, `r=2b/g`, `lambda=4b^2/g^2`.

The Sylvester operator `S(X)=(kappa I+D)X-XA` has inverse norm at most 1/g.
For example, the integral of `exp[-t(kappa I+D)] Y exp(tA)` gives that bound.
The Riccati equation for an invariant graph is

`S(X)=-B+X B^dagger X`.

Starting with X_0=0, define

`X_(n+1)=S^(-1)(-B+X_n B^dagger X_n)`.

On the operator ball `||X||<=r<=1/2`, the map sends the ball to itself and
has Lipschitz constant lambda<=1/4. Thus it converges to a unique X there,
with `||X-X_n||<=r lambda^n`. The graph is invariant under the complete H.
Varying B continuously from zero preserves the spectral separation: the
lower band lies below a+b and the upper band above kappa-a-b. The graph
therefore spans the complete lower band, of dimension rank(P), rather than
selected eigenvectors. This argument also holds within any full conserved
physical block. Its source symmetries are inherited by uniqueness.

Let `W_X=(I;X)` and `J_X=W_X(I+X^dagger X)^(-1/2)`. Then J_X is an isometry,
and `K=J_X^dagger H J_X` is a self-adjoint effective Hamiltonian satisfying
`H J_X=J_X K`. No frozen-link approximation appears in this exact definition.
The same construction from X_n gives J_n and K_n. The inverse-square-root
integral representation implies

`||J_n-J_X|| <= (1+2r) r lambda^n`.

With `h=kappa+a+b>=||H||`, a conservative full residual bound is

`||H J_n-J_n K_n|| <=4h(1+2r) r lambda^n = epsilon_n`.

Duhamel then bounds the state error by `|t| epsilon_n`; every complete effect
and the total variation of any common finite outcome distribution have the
same bound, capped at one. For the latter, take the supremum over unions of
retained effects. The initial state is J_n rho J_n^dagger and the reference
comparison is the same rho evolved with K_n. Bare P inputs require their
additional preparation/dressing error. The exact finite isometry admits a
Gauss-preserving unitary extension on its source and retained controller;
an efficient local gate compiler is not asserted.

This is a source-specific realization of
[EffectiveQuantumComparison.lean](../../../../Lean/QFT/EffectiveQuantumComparison.lean).
The full generator residual and the common reference are supplied, rather
than assumed small. The general spectral method is consistent with the
[Schrieffer-Wolff construction of Bravyi, DiVincenzo and Loss](https://arxiv.org/abs/1105.0675);
its source constants and convergence estimates are derived here.

## 4. Exact finite physical execution and a positive interaction bound

At kappa=128 the code uses the complete neutral Gauss block `(L,2Sz)=(0,2)`.
It has fifteen states, of which six have a blank link and nine an excited
link. The physical embedding and its positive rational Gram matrix come
from the [full Gauss and spectrum source](../two_cell/SPECTRUM_PROOF.md).
The actual phase D of that embedding is retained. No physical state in this
block is dropped from H or from the residual.

Three exact rational Sylvester/Riccati iterations produce X_3. To avoid
assuming Euclidean coordinates for a nonorthogonal physical basis, the code
exports the reduced metric and Hermitian form

`G_3=W_3^dagger G W_3`, `Kform_3=W_3^dagger G H W_3`.

The coordinate generator is `G_3^(-1) Kform_3`; it is self-adjoint for G_3.
Its canonical operator is obtained by the positive metric square root. For
`R_3=H W_3-W_3 G_3^(-1) Kform_3`, the exact rational number

`Tr(G_3^(-1) R_3^dagger G R_3)`

bounds the squared canonical operator residual. Its outward square-root
bound is below `10^-12`. Thus all dressed inputs in this six-dimensional
reference have complete-record dynamical error below `10^-12 |t|`, while
all fifteen original physical states participate in the residual. This is
an interacting quantum source with controlled finite approximation. The
block is not the global ground sector, and its level differences are not
relativistic masses.

There is also a parameter-domain result on the full source. Write
`y_bound=(2ar+b r^2)/kappa`. The exact graph equation gives
`||X+B/kappa||<=y_bound`. The canonical generator satisfies

`||K-(A-B^dagger B/kappa)|| <= b*y_bound+r^2*(a+b*r)=delta_K`.

Consequently the four compensated dressed inputs from section 2 have
connected energy at least `1/(4 kappa)-4 delta_K`. The executable checks this
lower bound is strictly positive at kappa=65536. This is a rigorous
nonadditive finite energy response, with both local Yukawa operators kept.
It neither uses measured masses nor proves that the interaction survives a
physical continuum limit. Kappa=128, 65536 and the unit-parameter spectrum
are distinct explicitly specified members of one operator family.

## 5. Raw periodic paths and a bounded alternative

In a native irreducible block of dimension d, a length-r path has

`W_ab |blank,...,blank> = sum_(r-1 internal indices) |link word>`.

All `d^(r-1)` words are distinct and orthogonal. Thus the fixed-endpoint
squared norm is `d^(r-1)`, and the full transporter acting on one internal
endpoint vector has squared norm `d^r`. The code verifies actual native
labels for d=3,6,1,2,3 and path lengths one through four. This is a property
of the declared raising matrices, not a heuristic estimate.

In the raw periodic source, even neighboring cells are connected by a
forward path of length N-1 around the torus. Its coefficient magnitude is
`pi/(2Na sin(pi/N))`, approaching `1/(2a)` at fixed spacing. On blank links
the corresponding squared response is that coefficient squared times
`d^(N-1)`. A source fermion can be paired with a retained conjugate internal
reference to make the input gauge neutral; orthogonality of those reference
states preserves the same norm. The d>1 raw path therefore supplies no
volume-uniform blank-reference local response bound. Orthogonality of its
excited path records prevents the cancellations of an identity-background
Fourier symbol. This rejects that specific free-reference transfer; it does
not exclude a different correlated interacting reference.

A permitted architecture alternative is explicit: replace U_ab by
`U_ab/sqrt(d)` in each irreducible block. On endpoint-internal space times
its blank/excited link, define

`V=sum_ab |a><b| tensor U_ab/sqrt(d)`.

Then `V^dagger V=I_d tensor P_blank`; V is a partial isometry of norm one.
Every path product is contractive and has norm one on a normalized internal
state with all links blank. Native current Ward identities survive the
block-scalar normalization. The actual raw family is not modified here.
If identity is substituted into the original unnormalized U, this change
also changes the free comparison coefficients; that source identification
cannot be kept silently.

For the normalized alternative, orthogonality of distinct directed paths
and `sum_(r=1)^(N-1) csc^2(pi r/N)=(N^2-1)/3` give the exact three-axis
one-fermion hopping response

`||H_hop psi_blank||^2 = pi^2*(N^2-1)/(4 N^2 a^2) <=pi^2/(4a^2)`.

The code checks the normalized matrix identity in every native block size
and the momentum sums for N=2,4,6,8. The normalization removes the stated
exponential path amplification. It does not turn a nilpotent quantum link
into a classical unitary, remove star parity, provide a magnetic plaquette
response, or prove a local interacting relativistic limit. Direct quantum
spectroscopy or a controlled dressed comparison is required on the chosen
normalized construction as well.

## 6. Physical interpretation and the next quantitative consumer

The established full-source graph construction supplies a controlled finite
interaction and state map. The periodic family supplies three compatible
translations, with its own quantum reference. Combining them requires a
common family and error estimates; their states and clocks are not joined
by a matching coefficient name. Growing-volume local estimates must control
actual path norms and retained quantum references. The local spin-reduction
theorem cannot be imported across these long-range path operators without
its own locality hypotheses.

The source has no explicit magnetic plaquette term. Its absence alone does
not prove that no magnetic response can be induced by matter. An induced
term needs its actual gauge-invariant correlation or effective-operator
calculation, with the eliminated sector retained in the error. A chosen
magnetic term is also an admissible architecture change, with new source
ancestry. The general existence of finite quantum-link constructions in
[Chandrasekharan and Wiese](https://arxiv.org/abs/hep-lat/9609042) supplies no
continuum theorem for these particular link matrices and reference.

Spin classification requires the physical spatial/little-group action on
the same excitation channel. Internal sigma matrices, conserved spin
projection and Kramers doublets have their separate finite meanings. Mass
identification accepts the relevant observable: a pole for an appropriate
stable gapped sector, an infrared-dressed charged threshold with inclusive
records, a resonance for unstable particles, or a specified renormalized
confined-quark parameter. A sharp isolated pole is not required of every
observed charged particle. None of these identifications follows from the
finite virtual-interaction coefficient alone.
