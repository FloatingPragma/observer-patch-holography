# Coherent source codes, entropy selection and charged clocks

## Question and answer

The explicit massive walk in `code/m1_operational_clocks` has a readable
inertial clock, but supplying its matrices is not a source realization.
The relevant source already exists in `code/source_selection_model`: its
accessible processor is M6, its record centre has twelve primitive ports,
and its complete port response is u(3) plus so(3). The present construction
uses that response. It does not append an arbitrary quantum gate set.

There are three distinct conclusions to establish. Classical record
agreement permits complete loss of noncentral phase, and the declared
maximum-randomness transition problem selects that loss. Agreement for two
complementary noncentral read experiments forces a coherent code channel.
Finally, the existing response supplies quantum computation on such codes,
including the clock circuit, with its durations and errors retained.

The extra physical requirement is stated below as coherent code gluing.
It strengthens the proposed classical RG law; it is not installed as a
fourth core axiom. The construction concerns an explicit source family,
not selection of a unique representation, particle mass or physical energy.

## 1. The response is an active processor

Write the twelve existing generators as

```
D_p = diag(hat(v_p/2)+i Q_p, hat(sigma(v_p)/2)),
Q_p = v_p v_p^T/2 - r^2 I/6 + I/12.
```

The v_p are the exact icosahedral calculation frame; sigma conjugates
sqrt(5). They are reconstructed from the retained ordered source response,
not chosen from a desired gate matrix. The parent proves that their real
span is u(3) plus so(3), is commutator closed and has observable rank twelve
on M6. In particular the relative phase between the two triples is visible.
Restricting the accessible algebra to M3 plus M3 would erase that direction
and is not the source used here.

The frame identities sum v_p=0 and sum v_p v_p^T=4 r^2 I give

```
J := sum_p D_p = diag(i I3, 0).
```

Thus a uniform drive of the twelve ports supplies the relative phase J.
The coefficient vector is exactly twelve ones. For any desired skew
generator H in the first u(3) block, its source controls are uniquely
reconstructed by

```
G_pq = -Tr(D_p D_q),
b_p = -Tr(D_p H),          h = G^(-1) b.
```

G is positive definite. This is an exact control reconstruction in
Q(sqrt(5)), not a fitted Hamiltonian. It supplies both real and imaginary
two-level rotations and all three diagonal phase controls.

Use the first two processor coordinates as a qubit code. Its SU(2)
rotations are available in u(3). Use coordinates 0,1,2,3 as the pair code
00,01,10,11. The same source pulse has restriction

```
exp(-theta J)|pair = diag(e^(-i theta), e^(-i theta), e^(-i theta), 1)
                  ~ diag(1,1,1,e^(i theta)),
```

where ~ removes one global phase on the code, not a relative phase.
At theta=pi this is controlled-Z. Acting on |+>|+> gives concurrence one;
at general theta the concurrence is |sin(theta/2)|. Single-qubit rotations
and this entangler supply ordinary finite quantum circuits once codes can
be grouped, moved and separated coherently. No coincidence between a coin
index and one of the twelve central ports is asserted.

More specifically, a two-mode SU(2) mixing acts on the occupation code
00,01,10,11 as diag(1,V,1), V in SU(2). This is already a native u(3)
operation on coordinates 0,1,2. Individual occupation phases are also
native, up to a common phase: diag(1,e^(i alpha),1,e^(i alpha)) is the
restriction of e^(i alpha) diag(e^(-i alpha),1,e^(-i alpha),1,1,1).
Adjacent two-mode rotations and all final diagonal phases therefore
implement the entire tetrahedral coin, with no dropped matrix entries.

The logical pair is a code used by a composite experiment. It is not a
replacement twelve-port carrier. Unrestricted two-qubit control has Lie
algebra su(4), of dimension fifteen; it cannot be identified with a
twelve-dimensional primitive response. The source distinction below is
essential to this construction.

## 2. Proper codes and complete physical instruments

A code request is a fixed isometry V:C^d -> C^6, with d=2 or d=4.
Its ordered computational basis is part of the calibrated interface.
The local calibration/preparation interfaces are explicit interventions,
not ontic entropy optimizers. In particular the reset/seed maps on the
processor have Kraus operators |b><j|, j=0,...,5, with b=0 or 1.
Their adjoint products sum to I6 and their output is |b><b| for every input.
They are rank-one, irreversible operations. The sender reset in a code move
is the b=0 map. Local destructive binary readback uses effects E and I-E
and a fixed reset output; its retained outcome is never postselected.
These finite instruments extend the source's local preparation/readback
interface. A pure blank is not claimed to follow from the unconstrained
tracial A3 minimum.

Packing |a>|b> into |2a+b> and unpacking are fixed coordinate identifications.
They do not permit an arbitrary input-dependent change of quantum basis.
Changing a source chart transports V, the state and every effect together;
it does not execute a free active gate. All nontrivial rotations in the
clock circuit use the response pulses of section 1.

Let P=VV*. The full source instrument first distinguishes P from 1-P.
On the valid code it applies the fixed code-service channel T and writes its
output into a receiver whose previous contents are discarded. The full
CPTP family below is the transition-selection domain; T is its selected
service, not an arbitrary per-instruction gate argument. In the exact
coherent completion T is identity, and in the classical completion it is
dephasing. The program cannot choose an extra unitary through this slot.
The sender is reset. On the complement it emits a retained failure flag and a fixed
output. For the identity code channel an explicit decoder into a code plus
failure space has Kraus operators

```
K_ok = sum_(j<d) |j>_out <j|_processor,
K_j  = |fail>_out <j|_processor,       d <= j < 6.
```

Their adjoint products sum to I6. Packing two qubits uses the analogous
tensor-code projector and fixed isometry into the four-dimensional code.
The previous receiver state is traced out in both cases. Source and
receiver central records retain the public success/failure outcome and
the immutable writer version. No successful branch is selected after
execution; the failure outcome belongs to the channel and the evidence.

These channels preserve arbitrary reference entanglement on their stated
code when T is the identity. They are not injective on the entire physical
algebra: |j><k| with j in P and k in 1-P is in their kernel. Their receiver
overwrite supplies another explicit kernel. An inverse logical move into
a new blank receiver recovers the code, not the overwritten physical state.
No unknown quantum state is copied into a retained classical record.

### Response completeness, including controllers

The primitive reversible port responses act uniformly on the twelve
central sectors and are exactly the source D responses and their proper
chart transports. Record-dependent raw blockwise unitaries are not added.
Classically controlled quantum instructions enter through the proper-code
instrument, even when their requested angle is zero. Internal destructive
readback and repair retain the parent source's nontrivial kernels.

Consider a putative additional reversible response of a complete carrier,
with every other input, cache and controller held fixed. Until the first
non-response operation touches its quantum state, its processor has only
undergone source response unitaries and its buffers are unchanged. If the
first touch accesses the processor code, the proper-code operation kills
a conjugate nonzero code/complement matrix unit. If it accesses a private
buffer, the mandatory load overwrites the unknown processor and kills a
traceless processor perturbation tensored with a fixed buffer state.
Destructive readback or reset likewise loses a nonzero perturbation. No
later linear instrument can recover it, even from the entire collection
of outputs and flags. If neither occurs, the connected reversible response
is a word in the original closed Lie algebra. Classical sector permutations commute with the
uniform private response, and the admitted proper rechartings normalize it.

For randomized read requests a common linear kernel need not survive the
classical mixture. The stronger channel argument is therefore needed. Work inside a fixed
primitive central block, where the quantum algebra is a full matrix algebra. Each
Kraus operator of a proper-code export, processor-overwriting buffer load,
or destructive read has rank strictly below the full carrier input
dimension (tensor unchanged buffers and append the retained flag). If a
later CPTP recovery made the total response a unitary U on that full input,
every composite Kraus operator L_j K_i would have to be proportional to U:
the unitary channel has rank-one Choi state. But rank(L_j K_i) cannot exceed
rank(K_i), so all these proportionality constants would be zero, contradicting
trace preservation. More generally one nonzero rank-deficient Kraus
operator already suffices: a trace-preserving recovery cannot annihilate
all of its subsequent branches. A nonzero irreversible branch among otherwise reversible
branches gives the same contradiction. Fixed auxiliary inputs can be
purified before this argument; appending a flag or keeping the entire
transcript cannot repair the lost rank. Thus randomization, feedback and
external controllers cannot turn the proper-code service into an omitted
full-carrier reversible direction. The claim is about CPTP reversibility,
not mere injectivity of a state-tomography map.

At a fixed finite cutoff, the same distinction survives channel limits.
In quantum dimension D, a Kraus family of rank at most D-1 has fidelity
with any full-space unitary Choi state at most (D-1)/D, by
|Tr(U* K)|^2 <= rank(K) Tr(K* K). A finite admitted program is a sum of
untouched branches proportional to original response unitaries and
rank-deficient touched branches. Convergence to a unitary forces the latter
weight to zero. The original U(3) times SO(3) response group is compact, and
a unitary is an extreme channel, so a unitary limit of the remaining
mixtures still belongs to that group. Weakening an irreversible operation
or taking a limit cannot silently supply a new full-carrier response.
This argument is at fixed cutoff; it does not identify a growing encoded
logical algebra with one primitive response tangent.

The processor is active, rather than an unused guard. It loads every private
buffer, supplies every native gate and executes readback/repair. Nevertheless
the exact load/overwrite rule is part of this *named* source grammar. A
hardware design allowing nondestructive access to all buffers must repeat
the response-completeness analysis; it cannot cite this proof unchanged.

The restriction is substantive and testable. Admitting a direct adjustable
buffer phase with no processor overwrite adds a public reversible
direction. Admitting arbitrary central-record-controlled raw responses
likewise adds blockwise directions. Those are different source grammars.
The construction does not assert response completeness for them.

## 3. A finite net and its all-level source completion

Here are the model-membership obligations, rather than just a circuit
placed beside a spherical picture. Let S_n be the parent's normalized
midpoint refinement of the icosahedral sphere. Replace each support vertex
v by a nonempty finite fibre F_v of carrier charts. A collection is a nerve
simplex exactly when its distinct support labels form a simplex of S_n.
The resulting nerve N has a simplicial projection b:N -> S_n. Selecting
one carrier in each fibre gives a section s. Hence

```
z=s_*[S_n],                 b_* z=[S_n].
```

The bridge has degree one. Fibres may be enlarged to accommodate any fixed
finite program without changing the carrier boundary packet. Programs may
be hosted at any chosen support vertex; that choice is transported with
the experiment under presentation change. Each primitive carrier still
has its own twelve-port boundary, independently of the federation nerve.

For each carrier i use a private centre C^12_i, active processor Q_i=M6
and finitely many **privately owned** matrix buffers B_(i,k). Define

```
A_i      = C^12_i tensor M6_i tensor (tensor_k B_(i,k)),
A_global = tensor_i A_i.
```

The centre of A_i is exactly C^12_i. Every static pair and higher seam is
the parent's scalar algebra; restrictions are normalized traces and their
cocycles commute. Include the finite unions and intersections of the
primitive patches needed by the experiment, with the corresponding tensor
subalgebras. This is an isotone net; disjoint factor regions commute. The
union observer is not declared a primitive twelve-port carrier.

A *process link* is distinct from a static overlap restriction. It is a
registered joint CP update between named private input/output interfaces,
with a departure, flight and arrival in the additional RG grading. It
moves a code into a different receiver factor and overwrites that factor;
it never identifies the two tensor factors as one shared matrix algebra.
All such links used by one computation are hosted in one support fibre.
Each fibre has the same cofinal capability, at arbitrary operational bulk
positions. Support labels are not positions in the additional affine E.
This is permitted by the axiom's typed update and interface data; it is a
new declared transport grammar, not an update inherited from the parent's
noncommunicating scalar-seam grammar. Scalar *static* seams alone do not
forbid newly specified causal joint updates. The tensor-factor net and
actual update channels together specify the model.

Every quantum buffer access, including a sender departure and receiver
arrival, uses its owner's active processor as the interface engine. Loading
a buffer overwrites that processor's old state. Access to a code already in
the processor uses the proper-code projection of section 2. Arbitrary full
M6 readback is destructive, as in the parent. Direct buffer phases and
remote access through another carrier's processor are not operations of
this grammar. Thus the first non-response quantum operation on a carrier
kills either a full processor direction or a code/complement coherence.
There is no dormant spectator added as an artificial guard: this same
processor is the engine for every exported code and every native gate.

The ownership rule is necessary for the completeness proof. If a matrix
buffer belonged to two carriers and its other owner could rotate it while
the first owner's entire processor was unchanged, the first carrier would
acquire a new full-algebra reversible direction. Merely calling that a
composite experiment would not remove the direction. That tempting shared-
buffer completion is excluded and its rank-thirteen control is tested.

The response Hilbert space remains the six-dimensional processor space;
the represented action on A_i is D on Q_i and identity on the buffers.
There is no unobserved response spectator: the same Q_i performs the
entangling phase, rotations, readback and repair. Its full M6 effects
retain the parent's rank-twelve process tomography, including the
cross-triple probes. The source's actual Cayley chart words implement the
proper carrier automorphisms on this same processor. Code isometries and
probe requests transform with those chart words. The buffer's abstract
logical labels are transported with its interface, not independently used
to supply the missing endogenous response.

Give every physical register, buffer access and event an address in the
affine response space E and a nonnegative duration. These addresses and
durations implement RG's additional physical attachment; A1's spherical
support is not being reidentified with bulk position. A buffer has one
physical location at each access event. Algebra membership does not grant
an instantaneous remote read: the typed read, write, controller and report
interfaces include the corresponding flights. A buffer belongs to one carrier. A flight transfers its code to a different
private factor, resetting the sender and overwriting the receiver; it does
not grant both carriers persistent access to the same quantum factor.
Local operations have zero displacement; a flight of displacement v costs
||v||/c; waits preserve the stated record. Sequential time adds, parallel
time takes a maximum. Every payload and metadata dependency is in this
same ledger. In particular, a nonlocal effect on the union observer is not
declared a zero-time local operation at one endpoint.

The classical service copies a central symbol into a buffer in the stated
basis, overwriting its previous state, and copies the decoded symbol into
the receiver's centre. The standard Kraus operators
e_p tensor |p><b| sum in adjoint product to the identity. Finite classical
permutations, preparation, retained versions and discard give the parent
NOT/Toffoli compiler. Unknown noncentral inputs are not cloned. All inputs,
intermediate interventions, failures and receiver overwrites are retained.
Since fibres and buffer lists are cofinal, every finite classical program
has a finite realization. Thus the classical source has the complete RG
service, rather than merely some successful M1 routing examples.

Refinement uses the parent's normalized midpoint subdivision, declared
simplicial coarsening f:S_(n+1) -> S_n and controlled mesh estimate. Each
old support vertex persists. Retain every old carrier (v,j) and its private
processor and buffer factors once. New carriers at old vertices enlarge their fibre;
new vertices start with a nonempty fibre. Under nerve coarsening a retained
carrier maps to itself, a new slot over an old vertex maps to its designated
slot zero, and a carrier over a new vertex v maps to slot zero over f(v).
This is simplicial: the distinct labels of a simplex map under f to a
coarse simplex. It satisfies b_n n = f b_(n+1) on every vertex, hence on the
whole nerve. Choose slot zero as the section at every stage. The section
cycles then push forward exactly by the parent's degree-one chain identity.
The choices are transported under presentation equivalence, not asserted
canonical in an unlabeled carrier.

The private ownership rule makes the algebra refinement literal. Retained
carrier factors, including each buffer, embed once by tensoring with
identities on new independent factors; states restrict by partial trace.
All static seams remain scalar. Restriction diamonds therefore commute;
there is no inclusion of one old matrix into two disjoint fine factors.
A coarse support observer is represented by the union of its connected
coarsening fibre; the parent's tagged paths supply its support routing.
The response is unchanged on retained processors and copied as a *law* to
new processors, never as an unknown quantum state.

Refining a flight adds co-located or collinear intermediate receiver slots
within the same support fibre, with the same endpoints and total graded
time. Identity transfers compose to identity; classical dephasing composes
to the same dephasing. Equation (5) supplies the noisy case. Intermediate
interventions act on the retained physical register at its actual location,
and all added processing is charged. For any finite program and positive
overhead tolerance choose the local service times small enough that the
sum over the *whole program* is below that tolerance. This is a capability
of the declared source family, not a bound on a fixed finite-strength device.

Regulators comprise a support level, finite per-fibre inventories, a finite
program and positive service-time bounds. Unions of inventories, common
support subdivision and smaller bounds give a directed common refinement.
Old data are retained once. The support cover, cross-triple response probes,
classical NOT/Toffoli instructions, destructive local tomography, repair
reset, mismatch comparison, retained writer versions and checkpoints are
all available. Deleting a response direction removes its independent
ordered probe; deleting a code buffer removes its chosen transport; deleting
its version merges two distinguishable intervention histories; deleting a
support face breaks the displayed oriented cycle. None of the source laws
mentions a mass, coin, empirical value or desired clock signal.

### Actual A3 objects and complete constraints

For the unconstrained ontic object take a global density on A_global and
all its restrictions. This declares a joint extension in this model; it
does not infer one from generic overlap consistency. The finite constraint
grammar is the complete Hermitian matrix-coordinate grammar of this
algebra, its partial traces, the displayed primitive CP instruments and
their compositions. Product local tomography spans the joint linear
observables; it does not supply an instantaneous joint measurement.

Choose the normalized full matrix/counting trace as the faithful reference,
and an injective observer cover containing the joint observer, with
strictly positive exact weights. All scored references are restrictions of
that same reference. Relative entropy is nonnegative and the joint term
vanishes only at the reference, so the unique unconstrained selected state
is the product tracial state. Its central port weights are 1/12. Adding
independent factors and taking the stated partial traces preserves this
selected family exactly. No general commutation of optimization and
refinement is asserted.

Prepared experimental inputs are declared local interventions and their
actual conditioning data. Their subsequent statistics follow the actual
CP instruments. They are not claimed to be the unconstrained ontic
minimizer. The finite transition object is a second, separately typed A3
problem: normalized Choi states of the code instrument, the input marginal
constraint, the source's accepted read-agreement constraints and the
faithful reference I/d^2. The next section solves that problem. No statement
about its optimizer is transferred from the stationary ontic state.

## 4. What agreement and maximum randomness select

Let d>=2, Z={|j>} be a computational basis of C^d and X its discrete
Fourier basis. For a trace-preserving code channel T put

```
J_T=(id tensor T)(|Omega><Omega|),   |Omega>=sum_j |j,j>/sqrt(d),
F_Z=(1/d) sum_j <j|T(|j><j|)|j>,
F_X=(1/d) sum_j <x_j|T(|x_j><x_j|)|x_j>.
```

These are two declared prepare/transport/read experiments, not a test that
a later clock happened to work. On the Choi space their effects are

```
P_Z=sum_j |j,j><j,j|,
P_X=sum_j |bar(x_j),x_j><bar(x_j),x_j|.
```

Direct Fourier summation gives P_Z P_X=P_X P_Z=|Omega><Omega|.
Consequently

```
F_e := <Omega|J_T|Omega> >= F_Z+F_X-1.                 (1)
```

If both read tables are exact, positivity puts J_T in the intersection of
the two projector ranges, which is one-dimensional. Normalization gives
J_T=|Omega><Omega|: the channel is identity, including on inputs entangled
with any reference. A2 meaning-preservation for this pair of declared
experiments therefore fixes the coherent transfer; no matrix entry of a
desired walk was inserted in the agreement constraints.

One basis is insufficient. Exact Z agreement makes every Kraus operator
diagonal, since each basis input has a pure unchanged output. Thus

```
T(|j><k|)=G_jk |j><k|,     G>=0, G_jj=1.             (2)
```

Conversely every such correlation matrix defines a CPTP channel. This is
the entire classical-agreement family, not just two selected examples.
Its Choi state is supported on span{|j,j>}. Relative entropy to I/d^2 is
minimized uniquely by the uniform state on that d-dimensional support:
G=I, the completely dephasing channel. All classical basis records survive
and all off-diagonal phase information is lost. With both bases the
feasible face is the singleton identity channel, so the same A3 rule
selects coherent transfer. The reference, trace convention and objective
have not changed.

### Exact noisy optimizer

Impose F_Z>=1-e_Z and F_X>=1-e_X, where
0<=e_Z,e_X<=1-1/d. Twirl a candidate Choi state by the commuting Bell
stabilizers. This preserves its input marginal and both constraints and
cannot decrease entropy. The Bell basis is
|Omega_ab>=(I tensor X^a Z^b)|Omega>. In that basis P_Z selects a=0
and P_X selects b=0. Subadditivity of classical entropy gives the unique
maximizer

```
p_ab = p_a q_b,
p_0=1-e_Z,  p_(a!=0)=e_Z/(d-1),
q_0=1-e_X,  q_(b!=0)=e_X/(d-1).                      (3)
```

Each marginal maximizes entropy with its zero outcome probability at least
the stated bound; below the uniform threshold the bound is saturated.
Strict concavity gives uniqueness, including equality cases by restriction
to the feasible face. For looser bounds replace the corresponding error
by 1-1/d. The result is a Weyl channel, not a channel chosen for clock
success. Orthogonality of the Bell states and a triangle bound give exactly

```
(1/2)||T-id||_diamond = 1-p_00 = e_Z+e_X-e_Z e_X.     (4)
```

The lower bound is the maximally entangled test; the upper bound writes
T=p_00 id+(1-p_00)N. For an arbitrary feasible channel rather than the A3
optimizer, (1) only gives the conservative bound
(1/2)||T-id||_diamond <= min(1,d sqrt(e_Z+e_X)).
Using the sharper formula (4) without the optimization argument would be
incorrect.

For exact Z agreement put
e_X(t)=(d-1)(1-exp(-lambda t))/d. The selected channel is

```
T_t = exp(-lambda t) id + (1-exp(-lambda t)) Delta_Z,
T_t T_s = T_(t+s).                                  (5)
```

This proves compatibility with subdividing a timed flight. A fixed error
per substep would instead destroy coherence under unlimited subdivision.
The two-error family likewise composes through its two commuting Weyl
dephasing semigroups. The rates are supplied agreement data, not values
derived from an ontic entropy minimum.

### A complete classical countermodel

Use the source of sections 2--3 with only the classical read table in its
accepted code-transport data. Its A3 transition optimizer is the dephasing
channel just derived. It retains the same complete response, endogenous
charts, full support/refinement family, classical RG compiler, ontic
reference and state-selection problem. Local quantum probes and response
pulses remain available; the model is not secretly a commutative algebra.

Every carrier-to-carrier quantum path in this grammar crosses a dephasing
service. With independent input and receiver preparations, any protocol
using these services, classical records, local operations, feedback and
stopping is a measure-and-prepare channel from the unknown input to the
receiver: condition on the complete classical transcript and compose the
local branches. The same argument applies to limits of finite stopping
protocols in fixed finite input/output dimension, since the separable
Choi set is closed. Pre-shared entanglement is not silently supplied by
the tracial product preparation or by the classical service.

For any such channel the Choi state is separable and F_e<=1/d. Indeed
|<Omega|u tensor v>|^2<=1/d for unit vectors, and mixtures preserve the
bound. Equation (1) implies the operational diagnostic

```
(F_Z+F_X)/2 <= (1+1/d)/2.                            (6)
```

Thus the four qubit prepare/read tests have average fidelity at most 3/4
in this entire transport class (the general qutrit bound is 2/3). They cannot realize coherent transfer, whose
average is one. A receiver cannot evade the bound by retaining classical
flags, waiting, making more copies of those flags or using a different
decoder. Abstention is an output failure, not removed from the denominator.
The obstruction is to remote unknown-state transport, not to all local
quantum interference inside an individual processor.

## 5. The additional principle

**Coherent code gluing (CCG).** In the already graded RG service, finite
tensor products of the specified qubit codes admit moving, waiting,
co-located packing into the specified pair code, and unpacking. Receivers,
storage, overwritten contents and failure outcomes are counted. For each
such code interface, both its fixed computational and Fourier read tables
belong to the accepted shared data and are preserved for every declared
input and context. The same complete clock grading, intermediate-writer
semantics and whole-program overhead requirement apply as for RG.

CCG concerns the transmission of two incompatible descriptions of the
same quantum information. It does not name a coin, mass, target state,
preferred radius, population or successful clock trajectory. It does not
permit copying an unknown code or changing its meaning by selecting a
different output basis for free. On the exact finite channel domain,
section 4 reduces its coherence clause to two finite read families and
proves its equivalence to identity code transport. The source construction
shows consistency with the stated primitive response and complete algebra
interfaces. The classical source above shows that the clause is not a
consequence of A1--A3 plus classical RG.

This is a precise additional physical requirement, not a claim that the
two-basis fidelity theorem itself is new. Its necessity and sufficiency
here concern the specified code service. They do not select the unique
microscopic OPH source or identify Floquet accounting with physical energy.

## 6. An actual source compiler for the clock

Use one ordinary occupation qubit per walk mode, including vacuum. For two
masses there are sixteen modes at every spatial cell. A pair of modes is
packed into the fixed four-state occupation code of section 1, operated on
by a native six-level response pulse, and unpacked into fresh receiver
slots. Both pack and unpack include their failure outcome and overwrite.
On the code the composition is the required channel even when the modes
are entangled with all remaining registers. Outside it the same physical
instrument is completely positive and irreversible, as already specified.

The QR elimination of the entire four-channel coin has six adjacent SU(2)
rotations and four final phases. The diagonal is retained. The source uses
these ten instructions for each coin. A macrostep for each mass consists
of the first-block coin, eight parallel register flights of displacement
plus/minus a(1,1,1), a(1,-1,-1), a(-1,1,-1), a(-1,-1,1), the second-block
coin, and four two-mode mass rotations. Their one-particle matrix is exactly

```
exp(-i m tau beta) diag(F C, C F*),       tau=sqrt(3)a/c.
```

Every mode flies; an empty input is a quantum register, not a discarded
branch. For two masses the ledger has 48 pulses, 96 pack/unpack events,
32 departure/arrival events and 16 register flights per cell per macrostep.
Independent full-basis periodic executions check the complete matrix, not
only the low band or a successful prepared input. The detector uses local
mass-beam rotations, occupation reads and the specified classical acceptance
weight. Its outcome is a bounded effect on the vacuum/one-particle sector.
Weighted acceptance can equivalently be implemented as a local destructive
two-outcome CP instrument; no nonlocal instantaneous read is introduced.

**Precise sector.** Tensor occupation-register permutation agrees with the
fermionic walk on vacuum and the one-particle sector, where this clock
lives. It is not fermionic second quantization on all occupation sectors:
ordinary SWAP leaves |11> unchanged, whereas fermionic SWAP gives -|11>.
This explicit countercontrol prevents transfer of the parent's all-Fock CAR
read theorem. Its one-particle evolution, packet, fringe, positive read
and accounting estimates do transfer. The classical RG causal/count cone
also persists. No all-Fock or interacting-matter source realization is
claimed by this compiler.

### Control time and finite hardware

A native two-mode pulse or phase can be executed with integrated operator
norm at most pi. Let the processor drive norm be bounded by Omega and let
each code event take eta>0. Choose positive h and

```
Omega = 96 pi/h,          eta = h/256.
```

Forty-eight pulses then cost at most h/2, and all 128 code events cost h/2.
Pad the parallel cells to equal duration. The actual macrostep is

```
tau_wall = tau+h = (1+kappa) tau,       h=kappa tau.
```

This is an explicit schedule, not an uncharged retiming convention. At its
actual timestamps, all parent comparison bounds apply with native tick n,
and the reference parameters become

```
v_wall = (c/3)/(1+kappa),  m_j,wall = m_j/(1+kappa),
u_wall = u/(1+kappa),     gamma = (1-u_wall^2/v_wall^2)^(-1/2).
```

Consequently the moving clock phase is Delta m_wall times t_wall/gamma.
The maximum record-flight speed is still c. At fixed spacing the required
Omega is finite. For fixed positive kappa it grows as 1/a; at bounded Omega
this particular sequential pulse compiler has nonvanishing macrostep
control time and its propagation speed tends to zero. A3 does not supply
free unlimited-strength control. CCG's cofinal-overhead clause states the
needed capability explicitly. No lower bound on every conceivable quantum
hardware implementation is inferred from this compiler's pulse count.

### Finite preparation and report

There is a finite source program for the actual truncated sampled packet.
Start with blank mode registers and one excitation in a root register.
Enclose its compact preparation support in a complete dyadic cube of N^3
leaf cells, side length Na. Each octree node uses seven native two-mode
rotations to split its excitation amplitude among eight children according
to their target squared-norm sums. Zero-weight children remain allocated
vacuum registers. Fly the eight children to the child cube centres in
parallel and repeat. At leaves use fifteen rotations and sixteen phases to
prepare the sixteen internal mass/spin amplitudes. This constructs arbitrary
specified one-particle amplitudes, so in particular the parent's normalized
smoothly truncated sampled packet. It does not assume a successful clock
state exists without a preparation program.

There are N^3-1 split pulses, 31 N^3 leaf pulses, and
8(N^3-1)/7 flights. The flight depth, rather than their total number, bounds
elapsed preparation time:

```
T_prep <= sqrt(3) a (N-1)/(2c)
          + (7 log2 N+31)(pi/Omega+2 eta) + (9 log2 N+66) eta.
```

There is one active processor per cell or tree node. At each tree level
it must perform eight departures **serially**; the eight child processors
then receive concurrently, costing one arrival time on each path. Initial
blank preparation serially loads and resets at most 32 buffers per
processor (64 events), and preparing/storing the one root excitation costs
two more events. Those are the 9 log2 N+66 events in the last term.
The gate pack/unpack events are already in the middle term. The resulting
path has 7 log2 N+31 pulses and 23 log2 N+128 non-pulse events. The emitted
integer schedule is checked separately from the much larger flight time,
so an omitted tiny local operation cannot hide within a floating tolerance.
For support inside radius P+sigma, choose the least dyadic Na >= 2(P+sigma+a). The flight
term is less than 2 sqrt(3)(P+sigma+a)/c, independent of leaf count. All
branches and their native events count as resources. The small executable
octrees include entire zero subtrees; they are checks of this construction,
not executions of the enormous analytic witness.

For n macrosteps allocate a finite workspace containing the entire native
speed-c forward cone of the compact preparation and the detector region.
Two banks of sixteen mode buffers per cell suffice for overwrite-safe
moves; each cell has an active processor. The tree uses at most sixteen
additional buffers and one processor per internal node. Immutable classical
success/failure records require at most one twelve-symbol record slot per
native event; allocate these separately as central slots of auxiliary
carriers. A conservative finite bound is

```
B = ceil(log2(2+N_buffers+N_prep_pulses+N_prep_flights+n N_workspace_cells)),
N_records <= (1+8B) [4 N_buffers + 4 N_prep_pulses + 4 N_prep_flights
                      + 256 n N_workspace_cells + 100 N_workspace_cells].
```

Eight B-bit identifiers/counters per event suffice for its finite indexed
program, endpoints and retained version. Analogue pulse settings are the
fixed instruction table, not arbitrarily precise unknown records. The bound
includes blank/seed preparation, all pulses, code events, flights, local
readout and record versions; it is not a claim of optimal storage. Program
parameters and local counter schedules are distributed before the blank
preparation; their classical compilation and distribution take finite
charged time. No quantum coherence is present during that initialization.
An adaptive replacement must charge its new dependencies by the same rule.
The counter schedule is an external control of the experiment, not an
independently derived physical time standard.

### Serial physical readout, including no-click outputs

After the native clock duration, freeze the walk and apply eight phase
pulses followed by eight native two-mode rotations between corresponding
mass modes. The first output amplitude in spin mode j is
(psi_(1,j)+exp(-i theta) psi_(2,j))/sqrt(2). There is still one processor at
each cell: load each of the eight first-output buffers separately, then
perform its destructive weighted occupation read. These are sixteen buffer
events in addition to the 32 pack/unpack events for the sixteen pulses.

On vacuum plus the one-particle sector, let |0> denote vacuum and |j> the
currently read mode after this basis change. With detector weight 0<=w<=1,
the actual destructive instrument is

```
C_j = sqrt(w) |0><j|,                 click,
M_j = sqrt(1-w) |0><j|,               no click, consumed mode,
O_j = I-|j><j|,                      no click, excitation elsewhere.
```

C_j* C_j+M_j* M_j+O_j* O_j=I. Keep both no-click branches; their probability
is not discarded or used to renormalize the next read. This is the reduced
payload map of a proper load followed by a local destructive binary effect
read on M6, with the processor reset and the classical outcome retained.
Distinct measured modes are orthogonal. Propagating the read effects
backwards through all the no-click instruments therefore gives exactly

```
E_clock = (w/2) sum_j (|1,j>+exp(i theta)|2,j>)
                           (<1,j|+exp(-i theta)<2,j|).
```

Its eigenvalues lie in [0,1] on this sector and its vacuum block is zero.
Thus the actual serial source read implements the parent detector, including
spectator entanglement, rather than merely reading a fitted phase. Across
cells the sum is still an effect on the global at-most-one-particle sector;
at most one click occurs. The executable controls reconstruct the entire
17-dimensional effect and every first-click probability independently.

The positive measurement duration is consequently

```
T_read = 16 pi/Omega + 48 eta.
```

Label the sample with the last macrostep endpoint and retain this measurement
latency separately. Outcomes become available after their read events.
Reports can then travel to a collector; for the witness below, collection
from the radius-600 detector costs at most sqrt(3)600/c using its containing
cube. This delays availability of the report without changing the recorded
local outcome. No instantaneous distributed aggregate is asserted.

## 7. Noise that preserves the signal can still fail the energy test

The A3-selected exact-Z phase channel has the same rate form on qubit and
pair codes. For the robust extension, permit its fixed-basis dissipative
semigroup on stored private registers as well as transported codes. This
adds no reversible response: at nonzero rate its Kraus family contains
nonzero proper basis projectors, so the preceding rank argument applies.
It does not permit an adjustable reversible buffer phase. The memory
semigroup acts continuously during storage/pulses; it is not implemented
by infinitely many positive-duration load events. Take disjoint local groups of modes and restrict to exactly one
excitation across the whole apparatus. For exposure t let s=exp(-lambda t).
If the excitation labels lie in the same group, their coherence is multiplied
by s; in different groups it is multiplied by s^2, since both groups can
record which alternative occurred. Thus the *entire* channel on this sector is

```
T = s^2 id + (s-s^2) Delta_group + (1-s) Delta_mode.
```

The coefficients are nonnegative and sum to one. This proves

```
(1/2)||T-id||_diamond <= 1-exp(-2 lambda t) <= 2 lambda t.          (7)
```

There is no population multiplier. Vacuum registers are included and have
identical states on alternatives in which they remain empty. Changing the
grouping, adding vacuum registers and interleaving number-preserving source
gates preserves the bound by telescoping and contractivity. The Trotter
limit covers phase noise during bounded pulses. Preparation, storage,
flights and readout therefore add at most 2 lambda T_exposure to any outcome
probability, with no conditioning on an error-free trajectory.

Destructive readout can remove the particle, so its intermediate payload
space also contains vacuum. The linear bound in (7) extends to this whole
space, including vacuum/particle coherence and reference entanglement.
Write P_j=|j><j|, P_g=sum_(j in g) P_j and P=sum_j P_j. The generator of
independent fixed-basis group dephasing, with rate lambda, is

```
L(rho) = lambda [sum_j P_j rho P_j + sum_g P_g rho P_g - P rho - rho P].
```

Each of the first two completely positive maps has diamond norm at most
one, because its adjoint sends I to P<=I. Left and right multiplication by
P each have norm at most one. Hence ||L||_diamond<=4 lambda; integrating
the CPTP semigroup gives (1/2)||exp(tL)-id||_diamond<=2 lambda t, independently
of the number of empty registers. Interleaving the complete read instruments
and retaining their outcome flags preserves this bound by contractivity.
This supplies the readout step without incorrectly applying the exact-N=1
convex-mixture identity to a larger Hilbert space.

During initial blanking the payload need not yet lie in the at-most-one-
particle space. That does not invoke the sector bound: each local reset
erases its incoming state, and subsequent phase noise fixes the prepared
computational vacuum. After all blanks and the basis-state seed are ready,
the entire payload is exactly in the declared sector. Their positive time
is retained in the exposure upper bound as a conservative excess; no
population-independent estimate is assumed for arbitrary uninitialized
many-particle states.

The qualification exact-Z is necessary. With independent occupation flips
of probability p on M blank registers, the probability of leaving vacuum is
1-(1-p)^M. Its dependence on population is real. Applying (7) to the general
two-error optimizer with nonzero bit-flip rate would be false. The general
A3 optimizer is still solved by (3); it simply does not have this stronger
one-particle stability property.

Contrast alone is also insufficient for energy. During a wait with
individual mode dephasing, the sector channel is
s^2 rho+(1-s^2) Delta_mode(rho). The parent's *same fixed* positive observable
E_a=|i log U|/tau obeys

```
0 <= E_a <= pi/tau,
Tr(E_a Delta_mode(rho)) >= asin(1/sqrt(6))/(2 tau).
```

Every diagonal component is a single-site particle, so the latter follows
from the parent's all-internal-state lower bound. For any fixed positive
wait time, uniformly bounded accounting energy requires lambda(a)=O(a).
A constant nonzero rate has divergent cost even if the measured fringe is
still large. Conversely (7) and the spectral range give

```
|Tr E_a(rho_noisy-rho_ideal)| <= (pi/tau) 2 lambda T_exposure.       (8)
```

For bounded total exposure, lambda=o(a) makes this extra cost vanish. For a finite workspace, use the compression of the parent's fixed E_a to
that workspace's one-particle subspace. Before readout the native forward
cone lies inside it, so the embedded state and its accounting expectation
are exactly those of the infinite-lattice experiment. Compression preserves
the upper spectral bound and all the stated single-site expectations; a
periodic logarithm is not silently substituted for the fixed observable.
These energy statements apply before the destructive final measurement. They do
not confuse mode energy with the control apparatus's energy. Positive
timing overhead multiplies the corresponding wall-time accounting by
1/(1+kappa) if that convention is used; no laboratory identification follows.

### Informative finite source-clock witness

Use the parent's fully bounded analytic witness:
c=3, a=10^-9, masses 100 and 100.1, reference velocity .6, sigma=10,
preparation cutoff P=80, detector scale 100 and extent 600. Set kappa=.01,
Omega=96 pi/(kappa tau), eta=kappa tau/256. The finite tree and workspace
sizes, all register counts and the drive bound are computed in the receipt.
The source construction is exact under CCG, lambda=0. As a separate robust
extension, choose the *declared* phase rate 10^-14 in these units; it is not
inferred from A3 or the axioms. It can be embedded in the family
lambda(a)=10^4 a^2 after fixing the displayed length/time units.

The actual preparation takes less than 105 time units, the full moving
cycle less than 80, and the charged read duration is included. Their noise
error is below 4*10^-12. The parent's outward interval calculation bounds
all other probability errors by .035153 and gives an actual noise-free
swing greater than .6753. The complete source-clock swing is consequently
**greater than .6752**, with extra accounting-energy cost below .022. Its
beat frequency is .08/1.01 in wall time. The output is informative after
preparation, flight, native control and readout costs have all been charged.
The enormous finite witness is an analytic existence construction, not a
claim that its lattice or hardware has been simulated or is practical.

For the alternative classical source, run the same detector instructions,
replacing every code service by its A3-selected dephasing channel. The
phase-stage code accesses erase the inter-mass matrix entries before the
mass mixers. Each mixer then sends either input occupation to its measured
output with probability 1/2. Propagating the actual serial read instruments
backwards gives E_classical=(w/2) P_one, independent of theta. The executable
comparison checks this full effect, not only its value on a chosen state.
For the equal two-mass superposition the coherent detector probabilities
at theta=0, pi/2, pi are 1, 1/2, 0; the classical completion gives 1/2 at
all three settings. Thus this source implementation loses the same clock
diagnostic when its code services are dephased.

The all-program statement is separate: averaging the four complementary
unknown-input tests cannot exceed 3/4. Native local interference is still
possible in that model, so neither conclusion says that every clock of
every classical-gluing source is impossible.

## 8. What is closed and what is not promoted

The missing gate supply has been replaced by an exact source-response
construction. The missing noncentral transport condition has a finite
operational criterion, a unique solved A3 transition problem, a complete
classical countermodel and an explicit sufficient source completion.
The clock has an executed compiler, an analytic preparation program, a
positive time budget, and a population-independent phase-noise estimate
with a separate ultraviolet-energy test. None of these depend on fitting
a desired continuum curve or retaining only successful histories.

CCG remains an additional proposed information-transport principle. A1--A3
and classical RG do not force it; both completions above obey the stated
source grammar. This construction establishes consistency and sufficiency
in the named source family, not unique microscopic selection by the core
axioms. Particle masses, empirical energy/time calibration and arbitrary
many-particle fermion dynamics are not promoted by a one-particle clock.
The finite M1 population and radius remain regulator choices under the
already established causal/count replacement results.

## Attribution

The complementary-fidelity inequality is due to Hofmann,
[Complementary classical fidelities](https://arxiv.org/abs/quant-ph/0409083).
The measure-and-prepare characterization is standard; see Horodecki, Shor
and Ruskai, [General entanglement breaking channels](https://arxiv.org/abs/quant-ph/0302031).
The source-specific claims require the response reconstruction, complete
source grammar, A3 optimization, clock compilation and resource analysis;
the quantum-information identities alone do not establish them.
