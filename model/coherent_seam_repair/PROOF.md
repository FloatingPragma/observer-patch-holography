# Finite coherent repair and observer readout

## Source and observation contract

`coherent-seam-repair-v1` is a declared finite program. Its carrier is the pinned
public twelve-port, thirty-seam, twenty-face icosahedron obtained through
`model/common_source/slow_geometry.py`. Its complex state and instrument
framework, coherent record, address reflection, preparation filter and optional
phase probes are architecture choices. No measured particle mass or fitted
coupling enters this package.

The scalar one-load state space is `C^12 tensor C^2 tensor C^60`: port, parity
memory and oriented seam address. One clock step is the unitary W defined below.
Preparation produces a reference on this walk-register algebra. Optional probes
are separately typed interventions; they are absent from idle evolution. The
five-level version adds an encoded fibre and a declared static link background.
Its stationary reference is recalculated for its own walk.

The committed observer's coin/current records and its two-slot tape define a
different program. Matching a first-event population kernel does not identify
its later histories, state or clock with those of the coherent program. Neither
construction's finite results are automatically premises of the separately
supplied Wilson gauge–matter field family.

## Local reversible repair

On a seam define `|+>=(|10>+|01>)/sqrt(2)` and
`|->=(|10>-|01>)/sqrt(2)`. The pure agreement pump has the complete local Kraus
operators, in source order `00,10,01,11,20,02`,

```
K0=|00><00|+|11><11|+|+><+|,
K1=|+><-|, K2=|11><20|, K3=|11><02|.
```

Their squared effects sum to I. On a blank record the reversible implementation
maps `|->|0>` to `|+>|1>`, `|20>|0>` to `|11>|2>` and `|02>|0>` to `|11>|3>`.
It exchanges each displayed input/output pair and fixes the orthogonal
complement. The three differences `v_j=input_j-output_j` are mutually
orthogonal with squared norm two, so `U=I-sum_j v_j v_j*` obeys
`U=U*=U^(-1)`. Its dimension is 24 and it conserves total load. Its four
blank-record blocks are exactly the Kraus operators above.

Committing or discarding each record and supplying a fresh blank gives the
irreversible pump. Reusing the coherent record gives a different evolution:
one repeated local operation is an echo because U squared is I. A classical
archive alone does not restore the coherence removed by a measurement.

`common.py` exposes both local operators. `program.py` executes the coherent
one-load rule; the double-load result is the exact local completion, without a
claim of a complete many-load coherent federation.

## Walk and reference

For every oriented seam `k=(a,b)`, U_k exchanges its antisymmetric source/parity-0
subspace with its symmetric source/parity-1 subspace. Let A append the uniform
coherent superposition of all sixty addresses, and set

```
Q=sum_k U_k tensor |k><k|, R=2AA*-I, W=QR.
```

Q and R are Hermitian involutions. W is unitary on dimension 1440. Averaging
the two orientations cancels cross-parity terms, giving the exact discriminant

```
A*QA = diag(P,P_sign),
P=I-L/60, P_sign=I-(Degree+Adjacency)/60.
```

Here P is the actual native one-load population kernel. From `WA=QA` and
`WQA=2QA(A*QA)-A`, induction gives `A*W^t A=T_t(A*QA)`. A discriminant
eigenvalue lambda has phases satisfying `cos(omega)=lambda`. This is the
reflection-walk spectral mapping used by the repository's
[coherent continuation](../observer/COHERENT_CONTINUATION.md) and the
[abstract quantum-walk theorem](https://arxiv.org/abs/1506.06457).

On the blank-parity block the carrier eigenvalues are
`1`, `11/12+sqrt(5)/60`, `11/12-sqrt(5)/60`, `9/10`, with multiplicities
`1,3,3,5`. The full walk has +1 multiplicity 38, -1 multiplicity 1356 and 46
oscillatory dimensions. The A-cyclic subspace has dimension 47; the blank-parity
source-cyclic subspace has dimension 23. Auxiliary record modes are included
in these counts and have no particle identification.

Starting at one source port, the first blank-record population update agrees
with the native rule. The second return probabilities differ:
the coherent observer gives `643/900`, while the committed Markov law gives
`101/120`. The callable implementation reproduces this difference exactly.

For a unitary W and `0<alpha<1`, prepare a maximally mixed walk register and
an ancilla `( |0>+sqrt(alpha)|1> )/sqrt(1+alpha)`. Controlled W followed by a
complete ancilla readout gives

```
K_success=(I+alpha W)/(1+alpha),
K_failure=sqrt(alpha)(I-W)/(1+alpha).
```

Their squared effects sum to I. Conditioned on the recorded success, the
reference is

`rho_W=[(1+alpha²)I+alpha(W+W*)]/Z`,
`Z=(1+alpha²)d+2alpha Re tr(W)`.

It commutes with W, and its eigenvalues are bounded below by
`(1-alpha)²/Z`. Each address reflector has trace 22, hence
`tr W=(2-60)22=-1276`. At the declared preparation choice `alpha=1/2`,

`rho_W=[(5/4)I+(W+W*)/2]/524`,

with eigenvalue lower bound `1/2096`, success probability `131/810`, and
failure probability `(679/810)^k` after k independent fresh restarts. A bounded
retry protocol reports its final failure flag. Preparation costs include
controlled W, source reinitialization and the retained failure records.

The reference is nontracial. At address `k=(a,b)`,
`rho_(a,0,k),(b,0,k)=-29/31440`, so an ordinary register projector does not
commute with it. The modular action is nonconstant on a readable effect.
Its generator is a function of `W+W*`; opposite walk phases receive the same
weight. This does not identify modular time with the walk clock or prove
geometric cap normalization.

The reference is explicitly reduced to the walk-register algebra. If a basis
label used to initialize the maximally mixed ensemble is retained and read,
conditioning on it changes the state. Stationarity of the joint walk-plus-
preparation archive is not asserted. Experimental pure preparations such as
`A|port,internal>` are separate, named interventions within this same program.

## Optional apparatus

The apparatus consists of the complete detectors `{P_i,I-P_i}` for the owned
register basis and the optional quarter-turns `V_i=I+(i-1)P_i`, together with
finite W steps. `program.py` returns both unnormalized detector branches and
their labels. Detector operations can disturb the source; they are excluded
from the idle lifetime calculation.

Every W diagonal is nonzero, with magnitude at least `29/60`. Every nonzero
entry has magnitude at least `1/60`. `apparatus.py` checks all 2,073,600 ordered
register pairs and constructs a directed support path of length at most four.

For the complexification of the real Hermitian history-effect span, let
`Ad_V(E)=V*EV`. Applying `(I-Ad_(V_i²))/2` and
`(I-Ad_(V_j²))/2` isolates the i,j pair. Applying
`(I+i Ad_(V_i))/2` then isolates the ordered matrix unit coefficient E_ij.
This is a linear reconstruction identity, not a physical outcome with a
complex or negative probability. Each term uses at most five quarter-turns.

Starting with a basis projector, a W pullback produces nonzero coefficients
along its support. If E_ij is available in the complexified span, pulling it
back by W and isolating a pair propagates j along the support because W_ii is
nonzero. The checked paths therefore give every matrix unit and the full
2,073,600-dimensional real Hermitian span, at depth at most 24 unitary settings
before the final detector. The coefficient amplification bound is recorded;
the result does not promise efficient finite-sample tomography.

This supplies complex tomography of the 1440-dimensional observer with the
declared optional phases. It does not implement every positive-operator-valued
measurement or every modular unitary. The gauge-coupled support splits into
charge blocks, so this dimension count does not transfer to the 7200-dimensional
version. Complete infinitesimal gauge response is also distinct from the
larger apparatus control algebra; its operational quotient requires proof.

## Encoded holonomy

`holonomy.py` reconstructs the proper rotations from the pinned carrier's exact
frame, lifts each to two unit quaternions, and compiles words from a pair of
source-ordered generators and their inverses. The declared code is a real
triplet plus the Galois-conjugate spinor doublet. Its central lift acts as
`Z=diag(I3,-I2)`, which is not scalar on the five-level code.

The source executes typed arrows `(source chart, gate, destination chart)` on
their five-level fibres. Each arrow has its recorded inverse. The code is
transformed and re-presented; no unknown quantum state is copied, and no
superposed chart register is unitarily reset. Exact checks cover 14,400 group
products, 480 generator representation products, 120 unitary matrices and
2,820 basis-state transitions. The longest compiled word has nine gates. A
fixed 1128-pulse schedule visits all lifts and returns its code and chart.

For the declared geometric port response `X_p=(P_p/2) direct_sum (-Q_p/2)`,
the exact Gram annihilator has positive roots `5/2,1/2,1/10`, hence rank twelve.
The image is the traceless block algebra `s(u(3)+u(2))`. Each implemented lift
belongs to its connected `SU(3) x SU(2)` subgroup and implements the same
recorded carrier action. The [A2 holonomy bridge](../../Lean/Screen/A2HolonomyBridge.lean)
permits a cover; it does not require an `A5 -> PU(5)` section. The existing
[finite port cover](../../Lean/Screen/SL2F5PortCover.lean) is an abstract finite
theorem; identification of its implementation with the concrete producer is
a separate formalization claim.

These are local same-history witnesses for the declared two-chart controller.
They do not establish a growing federation's accepted overlaps, refinement,
full scheduler covariance, or complete response for all apparatus operations.
The continuous twelve-direction response interface is declared separately
from the finite transport gates.

## Gauge-coupled walk

Attach a five-level fibre to the same moving source. A unitary link G from a
to b defines isometries

`m(h)=(|a>h-|b>Gh)/sqrt(2)`,
`p(h)=(|a>h+|b>Gh)/sqrt(2)`.

The local reversible operation exchanges `m(h)|0>` and `p(h)|1>` and fixes
the complement. Reversing the orientation uses G*, reverses m and preserves p
up to the fibre change. Averaging orientations gives

`diag(I-L_G/60, I-(Degree+A_G)/60)`.

The same Q/R construction gives a unitary W_G on dimension 7200. For block
unitary links, `Q5=diag(2,2,2,-3,-3)` is an exactly conserved internal charge.
Under a passive fibre recharting `V_a`, links transform as `G_ba -> V_b G_ba V_a*`
and the entire walk transforms covariantly. If the link alphabet is restricted
to the 120 cover values, local frames from that cover preserve it. A common
global continuous block center also preserves every link. Arbitrary independent
continuous local frames require admitting the corresponding enlarged link orbit.

Each controlled reflector has negative rank five and trace 110. Thus
`tr W_G=-6380`. The complete filter applied to this W_G prepares

`rho_G=[(5/4)I+(W_G+W_G*)/2]/2620`,

with eigenvalue lower bound `1/10480` and success probability `131/810`.
The scalar walk's reference is not substituted for this coupled reference.

For a concrete source signal, apply the actual compiled central word `[0,0]`
to select Z on one carrier edge, with identity on the other edges. The two
adjacent triangular Wilson loops have trace one, compared with trace five
for flat links. The exact adjacency trace changes from 600 to 552, giving
`delta tr(P_G³)=-1/4500`. The embedded three-step return amplitude of the
prepared doublet endpoint changes by `-1/6750`.

The callable observer executes these three steps and then returns both branches
of the detector onto the prepared vector `A|endpoint,internal3>`. Its return
probability changes by `-9353/91125000`; the complement probability changes by
the opposite amount. This is a complete simulation-generated signal from
static holonomy in one finite source and walk clock. The pure experiment
preparation is explicitly distinct from the stationary ensemble preparation.

## Protection and physical boundaries

The coupled walk conserves the static block charge Q5 under idle evolution.
This is an additional classical label, not a protected arbitrary quantum code,
physical spin, chirality or particle identity. The central flux cannot be
removed by a vertex recharting, but one allowed link replacement erases it.
No dynamical or topological lifetime is inferred from gauge invariance alone.

Each 12-dimensional component of the frustrated doublet connection Laplacian
is positive on the finite carrier and has determinant 13,132,800. The full
24-dimensional doublet determinant is its square. For any N-vertex graph with one negative edge,
the constant trial vector has Rayleigh quotient `4/N`. Thus
`lambda_min(L_signed)<=4/N`, and the corresponding complete-attempt walk phase
is at most `arccos(1-2/(M N))`. A finite carrier phase gap is not a volume-uniform
physical mass gap. A spatial scaling and clock conversion require their own
controlled family.

The committed two-setting phase pump is a separate comparison: it preserves
native populations for every density input, has a faithful stationary density
and a 144-dimensional history-effect span. Its coherence norm is bounded by
`3/4+sqrt(2)/12` per attempt, approximately 0.014235 per thirty-attempt sweep.
The graph-family bound excludes the desired undamped infrared wave for that
unchanged committed rule. The reversible construction changes its memory and
continuation law rather than inheriting its claimed propagation.

The package supplies no dynamical gauge field, backreaction, many-particle
scattering, spin-statistics theorem, controlled relativistic continuum,
Standard Model identification or GeV mass prediction. Its finite reference,
apparatus and gauge-signal statements have the domains specified above.

## Finite resources and replay

The literal scalar W entries have denominators dividing 60. Its reference
entries have common denominator 62,880. After t steps, an unnormalized density
admits common denominator `62,880*60^(2t)`; register-projector outcomes add no
denominator before conditional normalization. Exact rational storage grows
with the finite horizon. Optional quarter-turns require both real and imaginary
numerators. A classical emulator's memory is larger than `log(1440)`.

`python -m model.coherent_seam_repair --check` recomputes the local identities,
carrier spectra, controller words, apparatus support and actual observer
signals. The receipt hashes the public implementation and pinned source files.
No private workspace file, fitted mass table or research-note narrative is an
input. Passing this finite replay does not certify the stated physical
extensions.
