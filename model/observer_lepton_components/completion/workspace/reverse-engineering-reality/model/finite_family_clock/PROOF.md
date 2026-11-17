# Family records and finite Hamiltonian spectra

The executable calculations distinguish a family measurement effect, a
unitary tick, a neutral Higgs background and a quantum scalar memory. Each
has its own operator, preparation and clock. Their finite identities do
not supply a transfer to the primary Wilson source.

## One-tick family effect

The neutral family register is the regular representation of the cyclic
group of order three. Write `R^3=I`, with eigenvalues
`1, exp(2*pi*i/3), exp(4*pi*i/3)`. Its three characters are mathematical
sectors of the declared register. Calling them observed generations
requires the particle and preparation dictionary.

A common neutral register and three independently occupiable families are
different carriers. For a one-particle space K of dimension d,
`Fock(K) tensor C^3` has dimension `3*2^d`, while
`Fock(K tensor C^3)` has dimension `2^(3d)`. For six lepton modes these
are 192 and 262144. The three-family memory declares the latter CAR
construction; its multiplicity does not follow from the three characters.

Each seam object has twelve states, indexed by load (tail, head, idle),
record (zero or one), and address orientation. The exact single-object
matrix is `B=V Q R2`, where `Q` is the native seam repair, `R2` exchanges
the two addresses and `V` changes the sign of record one at the reversed
address. `native_event_check()` compares all six local columns at every
one of the sixty native oriented addresses, including an off-seam port,
against the actual native `Observer.repair` implementation. All 360
columns agree exactly. The two-address `R2` is a declared replacement for
the native sixty-address reflection.

Let `p` be the product preparation with each load amplitude
`sqrt(1/3)|tail>+sqrt(2/3)|idle>`, blank records and equal address
amplitudes. Let `Phi` be the normalized sum of the three orthogonal flags
`C idle x E agree`, `C repair+ x E idle`, and `C repair- x E idle`.
For the declared controlled family-shift tick `U(theta)`, exact compression
gives

```
M=<Phi|U|p>=(sqrt(2)/3) C(theta),
C(theta)=a I+b(exp(i theta)R+exp(-i theta)R^-1),
a=1/sqrt(6), b=1/sqrt(12).
```

The three eigenvalues of `C` are
`lambda_k=a+2b cos(theta+2*pi*k/3)`. On a chamber where every eigenvalue
is nonnegative, `sum lambda_k=3a` and `sum lambda_k^2=3a^2+6b^2=1`.
For the declared squared-eigenvalue mass reading,

```
Q = sum(lambda_k^2)/(sum(lambda_k))^2 = 2/3.
```

The positivity hypothesis is essential. It holds at `theta=2/9`, and
through `|theta|<=pi/12` modulo `2*pi/3`. Taking absolute square roots
outside this chamber changes the denominator.

The analyzer's complete binary effect pair is `A=M^dag M` and `I-A`.
At `theta=2/9`, the success probabilities in character order are about
`0.209693556`, `0.000060300559`, `0.012468366`. The maximally mixed family
success probability is exactly `2/27`. Keeping the selected amplitude
while omitting the complement is a conditioned experiment. The actual
complement can have several Kraus operators; its effect alone does not
specify the full conditional state.

The independent flags have probabilities `1/12` and `1/6`; their exclusive
or has probability `2/9`. This probability does not fix the coefficient
of a Hamiltonian coupling to the flag. The family phase is an architecture
input unless that coupling and its time standard are specified.

## The actual tick and its clock

In family character `k`, put `zeta=exp(i theta) exp(2*pi*i*k/3)`. The full
tick reduces exactly to `(D(zeta)B) tensor B`, with `D` multiplying the
two oriented repair records by `zeta` and `zeta^-1`. Direct symbolic
calculation gives

```
det(tI-D(zeta)B)
  = (t-1)^2(t+1)^2(t^2+1)^2(t^2+zeta)(zeta*t^2-1)/zeta,
det(tI-B)=(t^4-1)^3.
```

The actual Floquet spectrum has multiple branches in every character
sector. The `arccos(1-10^-6 lambda_k^2)` spectrum in `fr4_ratios.py` is a
separately specified reading; it is not the spectrum of this matrix.

For every nonzero `zeta`, the literal two-tick compression is

```
<Phi|U(theta)^2|p> = sqrt(3)/9 I3.
```

`family.py` proves this with twelve-state symbolic matrices. The
independent 432-state numerical reduction reproduces the one-, two-,
three- and four-tick amplitudes. In particular, the second tick cannot
equal `M^2`, and its three squared amplitudes are identical. A one-tick
effect does not establish a persistent Hamiltonian spectrum.

## An unconditional squared-effect clock

A permitted declared extension does realize `A=M^dag M` as an effective
Hamiltonian. This theorem applies to any finite unitary `U`, normalized
ancilla preparation `p` and normalized analyzer `Phi`; it uses no target
masses.

On ancilla and family let

```
P = U^dag (|Phi><Phi| tensor I) U,
A = (<p| tensor I) P (|p> tensor I).
```

`P` is a projector and `0<=A<=I`. Couple a fresh ancilla for duration
`delta` under the Hermitian generator `P`, and retain its output. The
family marginal channel is

```
E_delta(rho) = Tr_anc exp(-i delta P)
                  (|p><p| tensor rho) exp(i delta P).
```

Every ancilla outcome is retained in the complete record. There is no
success conditioning. The first derivative is exactly
`E'_0(rho)=-i[A,rho]`. For arbitrary reference extensions, the commutator
superoperator with a norm-one generator has diamond norm at most two.
The integral Taylor remainder of each unitary conjugation therefore has
diamond norm at most `2 delta^2`. Preparation and partial trace are
channels and cannot increase this norm. Comparing with
`V_delta(rho)=exp(-i delta A)rho exp(i delta A)` gives

```
||E_delta-V_delta||_diamond <= 4 delta^2.
```

Telescoping `N` channel products, using their diamond norm one, proves

```
||E_(t/N)^N - V_t||_diamond <= min(2,4 t^2/N).
```

The bound includes a retained quantum reference on the family input. It
concerns the family/reference marginal, rather than equality of all
outgoing ancilla states with an idle environment. For this register,
`A=(2/9)C(theta)^2`. The executable computes the exact finite collision
channel as a Schur multiplier: in the family eigenbasis,

```
u_k=U_k^dag Phi,
p_k(delta)=p+(exp(-i delta)-1)u_k<u_k,p>,
rho_jk -> <p_k(delta),p_j(delta)>^N rho_jk.
```

At `t=1,N=10000` the largest coherence-multiplier error is approximately
`9.16e-6`, within the uniform bound `0.0004`. The generator, inverse tick,
projector phase gate, fresh ancillas and schedule are explicit architecture
choices. This theorem supplies an executable clock extension. It does not
identify `A` with the original FR tick's logarithm or the primary Wilson
Hamiltonian.

## Three-family quantum scalar memory

`memory.py` uses the actual native lepton Yukawa tensors, CAR
sequence compiler and five-state quantum scalar matrices. It adds three
family copies with real positive coefficients `y_f`. There is one shared
scalar memory, scalar energy `omega>0`, and no spatial hopping. The full
local carrier has `5*2^18` states. All 320 columns with only family zero
active agree exactly with the source's complete one-family operator.

Graph closure of the vacuum preparation yields a thirteen-state invariant
star: a center with diagonal zero and twelve leaves with diagonal omega.
Its squared coupling norm is `S=4 sum_f y_f^2`. The charged creator
`T_f=c^dag_(f,spin0,e^c)` maps this entire component into an eleven-state
invariant star with squared norm `S_f=S-2y_f^2`. The executable generates
the graph and proves both invariance and this rectangular map exactly.

Let

```
D0=sqrt(omega^2+4S), Df=sqrt(omega^2+4S-8y_f^2),
E0=(omega-D0)/2, Ef=(omega-Df)/2.
```

The reference is the normalized `E0` eigenvector in the vacuum-generated
component. Its density is stationary under the same full local `H`.
The exact normalized ground vector is the star center plus its bright
leaf vector divided by `E0-omega`, normalized by its computed norm. A
Householder on the difference between the blank and this vector prepares
it deterministically and acts as identity on the orthogonal complement.
The executable verifies this algebraic preparation on the exact rational
coefficient control, retains a one-dimensional idle reference, and
asserts no global ground-state property. A charged
preparation carries the explicitly opposite-charge apparatus reservoir
with idle Hamiltonian zero. Omitting that reservoir would not produce a
gauge-invariant charged experiment.

The native twelve-current tensor identities hold independently in each
family copy, hence for their sum with arbitrary diagonal coefficients.
The empty fermion/scalar state is invariant. The selected reference lies
in the span of that state and its image under H, so it is invariant too.
The Householder between two invariant vectors commutes with the currents.
The executable checks that the actual native creator column is a weak
and color singlet with hypercharge one. A two-state apparatus with charge
`diag(0,-1)` and zero Hamiltonian has a transition of charge minus one;
the compensated creator therefore obeys all twelve Ward identities.

The complete target spectrum consists of `Ef`, `omega-Ef` and omega with
multiplicity nine. The elementary transfer has zero overlap with every
dark omega eigenstate. It has exactly two energy lines,

```
m_f = (D0-Df)/2 = 4y_f^2/(D0+Df),
M_f = (D0+Df)/2.
```

There are no omitted nonzero dark residues. If
`N(E,S)^2=1/(1+S/(E-omega)^2)`, their elementary residues are

```
R_-=N(E0,S)^2 N(Ef,Sf)^2
        [1+Sf/((E0-omega)(Ef-omega))]^2,
R_+=N(E0,S)^2 N(omega-Ef,Sf)^2
        [1+Sf/((E0-omega)(-Ef))]^2.
```

Both are positive when `y_f>0`. Their sum equals
`N(E0,S)^2[1+Sf/(E0-omega)^2]=<T_f^dag T_f>`, retaining the complete
elementary spectral weight. The elementary operation is a contraction.
Its complementary no-operation effect is exactly `I-T_f^dag T_f`, a
projector. The receipt retains its probability, and the two energy
weights plus that probability sum to one. Composite residues are squared
operator matrix elements; they are not declared probabilities for an
unscaled composite preparation. The composite is independently compiled
from the native scalar and CAR formula
`J_f=-y_f sum_(k,j) C^k_(e^c,j) phi_k^dag c_(f,spin1,j)`.
Scalar number commutes with the creator, so
the actual composite `J_f=[H_Y,T_f]` obeys on every retained input column

```
H_target T_f - T_f H_source = J_f.
```

Consequently its two residues are exactly `m_f^2 R_-` and `M_f^2 R_+`.
The source, reference and Hamiltonian clock `exp(-itH)` are the same for
the elementary and composite spectra. A finite stationary transfer is
proved; a relativistic dispersion, physical rest mass and GeV conversion
are outside this local theorem.

## Exact exclusion of one proposed dressing

Set `y_f=y lambda_f(2/9)` and keep any positive `y/omega`. For a heavier
and lighter family,

```
m_h/m_l = (y_h/y_l)^2 (D0+Dl)/(D0+Dh) > (y_h/y_l)^2,
```

because `Dh<Dl`. The inequality is exact and holds at every coupling,
including beyond perturbation theory. The weak-coupling limit gives the
squared-circulant ratios `206.770315972729...` and
`3477.47283710460...`. Every positive-coupling ratio exceeds these limits. Thus this local
memory self-energy cannot produce a negative correction to the
squared-circulant ratios. At `y=omega` the ratios are
about `208.1468625` and `3984.108802`. No phase or coupling is fitted by
these calculations.

The neutral-background source has a different exact relation,
`m=|vY|`, with `H=(0,v)` and hence `v=v_EW/sqrt(2)`. Choosing `Y=y C`
there gives ratios `14.3795103` and `58.9701012`; choosing `Y=y C^2`
gives the squared-circulant ratios. These are different declared Yukawa
operators. A finite quantum memory, a specified scalar background and
the collision clock cannot exchange spectra by sharing the symbol `C`.

## Source custody and replay

The receipt validates the immutable coherent-seam manifest and the native
quantum-Higgs receipt, including every inherited input pin. It binds this
package's Python files, tests and prose, the native representation data
and the transitive project Python import closure. All runtime paths are
repository relative. The minimal seam-register primitives are vendored
with their attribution in `PROVENANCE.md`; their correctness is checked
against the actual native event, rather than an external worker path.

`python -m model.finite_family_clock --write` regenerates the receipt.
`--check` validates pins and independently recomputes the complete result.
Exact strings, symbolic identities, scope flags and custody hashes compare
exactly. Only binary64 diagnostics have absolute tolerance `1e-10`.
Mutation tests reject altered exact claims and changed diagnostic spectra.
Neither the receipt nor its hashes transfer the source domains.
