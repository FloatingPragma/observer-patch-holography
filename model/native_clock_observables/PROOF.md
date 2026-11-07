# Native clock observables

## Source and allowed protocol

The source is the fixed weighted two-occupation carrier H66 in
`local_repair_family.modular_clock.parameters`: twelve ports, weights
(25,25,8,...,8), thirty seams and the complete tagged C/E instruments.
Its faithful stationary density is

    rho_* = diag(theta_i theta_j)/7505,  i<j.

The probe is the existing clock seam (1,7). The observer protocol uses the
existing nine-qubit proper-code cap rotation and quadrature words in
`native_modular_reference.cap_local_clock` and `phase_readout`, followed by
one native E repair on this seam. These coherent services and their code
load/export resources are declared inputs of the named coherent clock member.
They are not consequences of passive C/E repair. The composition is an explicit
observer protocol; the repair-only `native_program` grammar does not itself
encode the coherent-word steps. Its `source_word` admission checks every word;
the owning full-512 compiler executes them because the generic small-word
interpreter has a smaller execution limit.

The result concerns this fixed carrier. The independent product capacity
embedding for repair/hold commands supplies no encoding or resource theorem
for these coherent services. Intermediate compiler pulses are not proved
spatially local. Locality below means membership of the completed two-port
occupation algebra.

## Three accessible local observables

For every spectator k distinct from 1 and 7, write |+k> for occupation {1,k}
and |-k> for {7,k}. There are ten orthogonal pairs. On their direct sum Q define

    X = sum_k (|+k><-k| + |-k><+k|),
    Y = sum_k (-i|+k><-k| + i|-k><+k|),
    Z = sum_k (|+k><+k| - |-k><-k|).

Each vanishes on P0=I-Q. The ranks of Q and P0 are 20 and 46, respectively;
P0 includes both empty and doubly occupied endpoints with every spectator
configuration. Exact full-carrier identities give X²=Y²=Z²=Q and [X,Y]=2iZ.
Let J=sum_k(|+k><-k|-|-k><+k|). The existing ideal cap and quadrature words are

    U_X=P0+(Q+J)/sqrt(2),   U_Y=P0+(Q-iX)/sqrt(2),   U_Z=I.

They satisfy U_A† Z U_A=A. Hence the ideal projectors are
P_A,±=(Q±A)/2 and P_A,0=P0, without conditioning on Q.

For a fixed native E command the complete interpreter alphabet has 25 slots:
one global hold, eight C slots and sixteen E slots. Nineteen have zero maps.
The six live branches are indexed by the retained coin, raw flag, incoming
contrast delta and current j:

| Coin | Raw | delta | j | Probability |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 8/33 |
| 0 | 4 | -1 | 0 | 8/33 |
| 0 | 5 | +1 | +1 | 8/33 |
| 1 | 0 | 0 | 0 | 25/33 |
| 1 | 4 | +1 | 0 | 25/33 |
| 1 | 6 | -1 | -1 | 25/33 |

Probability here is the coin factor in the branch map, not its unconditional
probability. Write K_h=sqrt(p_coin)V_h for the native partial injection. The
raw analyzed instrument is K_h U_A rho U_A† K_h†, with all flags retained.
Summing its effects at a given contrast gives exactly P_A,delta. Each live
zero-contrast V_h equals P0 on all 46 dimensions: the flags reveal neither
empty versus double occupancy nor any spectator. All 66×66 operator identities
are checked, so the result holds with every finite untouched reference.

The native update is destructive: on an odd pair, forgetting coin and contrast
replaces its occupation by the biased placement (25/33,8/33). It is not a
Lüders update for X or Y merely because its effects are their projectors.

## Record-controlled Lüders instrument

Let S=U_X²=P0+J. For each raw record apply S if j is nonzero, then U_A†.
This uses two existing positive cap rotations for S and the literal inverse
word for U_A†. There is no added primitive or Hamiltonian. Directly on the
whole carrier,

    U_A† S^(1[j!=0]) V_h U_A = epsilon_h P_A,delta,

where epsilon_h=-1 only for coin 1, raw 6; all other live branches have sign +1.
The sign is constant across all spectators. Thus every fine corrected map is

    rho -> p_coin P_A,delta rho P_A,delta.

All zero addresses remain zero, and both zero-contrast branches preserve
arbitrary coherence throughout P0. Coarse-graining only the coin gives the
three-outcome Lüders instrument. Retaining the coin gives the same instrument
with its independent coin refinement. The full matrix identity implies its
tensor-reference extension and arbitrary subsequent protocol composition.
No outcome, padding address or failure is normalized away.

## Literal service execution and error

The four words U_X, U_X inverse, U_Y and U_Y inverse are executed as their
unmodified full-512 binary64 arrays. Their compiled source gates are existing
six-level pulses. Exact rational certificates bound the returned-array operator error on
all 512 addresses, including the 446 addresses outside H66. Each host checks
its recomputed residual against the fixed bound 10^-8. The receipt serializes
that uniform bound, rather than platform-dependent residual digits. A separate
multiplication bound compares the returned array with the exact product of
the stored binary64 primitive matrices. This uses IEEE binary64 round-to-nearest
arithmetic. Put u=2^-53 and gamma_n=nu/(1-nu). Each local complex dot has at
most four entries, so gamma_32 bounds its real arithmetic; each of the ten
512-dimensional block products uses gamma_4096. The checked primitive norm
is at most 1+g, g=2/10^11. For N<=12000, the inequality

    N[g+64 gamma_32(1+g)] + 10*1024 gamma_4096 + 2^-900 < 1/4

bounds all intermediate Frobenius norms by 32 and every subsequent propagation
gain by 2. Local multiplication error is at most
64 gamma_32(1+g), while a block product contributes at most
1024 gamma_4096. Summing and propagating gives

    e_arithmetic <= 128 N gamma_32(1+g) + 20480 gamma_4096 + 2^-900.

The absolute term covers gradual-underflow or flushed-subnormal arithmetic;
all other entries stay far below overflow by the displayed norm bound. The
reported e is the returned-array error plus this product error. This binds
the exact product of the stored primitive matrices. It does not certify each
matrix exponential against an unrounded analytic exponential or model physical
hardware errors. The finite implementation is the specified numerical source.
The diagnostic density evaluator also uses floating arithmetic, including
coin-weight multiplication; its rounding slack is reported separately and is
not included in the ideal rational Z-instrument's zero error. The complete
address export is applied as a completely positive trace-preserving map; its
abort output retains the untouched reference after tracing the escaped address.
`compiled_output` materializes every native slot and every address-export
abort on each branch.

For each word, the inherited proper-code theorem bounds the total mass f of
all first primitive failures. The bound includes their recorded pulse ordinal
and complement outcome, although these extremely long lists are represented
by the bound rather than individually allocated matrices. The code checks
all distinct executed primitive pulses against the theorem's local leakage
and norm gates. For an arbitrary source/reference input, the complete
half-diamond error is at most

    delta = [e(2+e)+f]/2,       gain <= (1+e)^2+f.

The address-export contribution is included in e(2+e) by contractivity. The
failure contribution f/2 is separate. No polar normalization or survivor
conditioning is applied to the executed arrays. Consecutive stages obey
D_next <= gain_next D_previous + delta_next. Classically controlled correction
uses the maximum stage bound over the retained raw labels. The worst branch
has one analyzer, two cap rotations and one inverse analyzer. The receipt
counts every native pulse and every proper-code load/export call. These are
mathematical implementation bounds, with no physical hardware error model.

The positive effects of the literal successful analyzer have perturbations
at most a=e_X(2+e_X) and b=e_Y(2+e_Y). Consequently their commutator norm is at
least 1/2-2(a+b+ab)>0. This is an executed noncommuting finite menu, even with
the stated compiler errors; it is not an assertion of all complex effects.

## Reference, uncertainty and clock

On the full normalized rho_*, the spectator weights sum to 97, giving

    <X>=<Y>=0,  <Z>=1649/7505,  <Q>=3201/7505,
    Var(X)=Var(Y)=3201/7505,
    |<[X,Y]>|²/4 = (1649/7505)².

Thus the source supplies a density matrix, actual Hermitian observables and
operational measurement instruments for the finite theorem
`EventAlgebra.Robertson.finite_state_robertson_commutator`. The strict margin
is 301088/2253001. The bound is evaluated without postselecting the odd sector.
The calculation does not identify position and momentum, a canonical
commutation relation, or a relativistic field.

The same existing source clock has K=-log(rho_*) and U_s=rho_*^(-is)=exp(isK).
Every probe pair has population ratio 25/8. Therefore, with Delta=log(25/8),

    U_s† X U_s = cos(s Delta) X - sin(s Delta) Y,
    U_s† Y U_s = sin(s Delta) X + cos(s Delta) Y.

These are exact ideal identities for real s; the existing compiled clock
certificates cover the declared parameters s=1/2 and s=1. They are cited as
such, without claiming an arbitrary-parameter compiler theorem. This K on
the working H66 carrier differs from the positive record-transfer generator
H=-log(S) reconstructed from stationary scalar histories. Their source is
shared, but their states, Hilbert interpretations and generators are not
identified. No physical time calibration is derived here.

The native mixing theorem prepares the reduced source/reference within half
trace distance 542159/10^12 after 1600 fixed-law unread attempts. The complete
preparation archive is traced in this replacement statement. A subsequent
fixed protocol has error at most its gain times this preparation bound plus
its complete compiler error. Adaptive conditioning on that archive is outside
this statement; archive correlations are not erased by rhetoric about an
unread marginal.

The finite measurement probabilities and state updates follow from the
specified quantum source model. The universal operational Born bridge also
requires implementation-independent probabilities for equal effects and
joint realization of every coexistent complex effect pair. This finite menu
does not produce those global premises. It supplies actual finite instruments
and uncertainty inputs on the declared coherent observer, without deriving
quantum primitives from a solely classical repair law.
