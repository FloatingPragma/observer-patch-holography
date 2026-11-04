# Complete radial readout with Hamiltonian drift

The source is the full finite-CAR radial Hamiltonian of
[radial_interacting_core](../radial_interacting_core/PROOF.md), on an open
odd cube containing the origin and its six neighbors. No number sector is
discarded. The radial mass, centered contact coupling U, spacing a and
hopping velocity v are supplied. The measurement is the same complete
origin-occupation instrument and recorded integer repair of
[radial_geometry_join](../radial_geometry_join/PROOF.md). The exact repeated
reference instrument is [radial_reused_readout](../radial_reused_readout/PROOF.md).
The reference comparison in this proof does not turn off terms in the
actual source.

## A local full-Fock norm

At the origin let

```
D0=(n_up+n_down-1)^2,        P0=2D0-I,
W0=U(n_up-1/2)(n_down-1/2)=(U/4)P0.
```

D0 is a projector, whereas P0 is the sign-parity observable. Denote the
origin-incident hopping by Hstar. Its one-body off-diagonal block V has
two columns and

```
V*V=b^2 I2,                 b^2=3v^2/(2a^2).
```

This follows by summing the six actual `±iv sigma_j/(2a)` edge blocks.
The two normalized columns define orthogonal remote CAR modes f1,f2.
Without deleting any other mode, the full hopping operator is unitarily
equivalent to

```
Hstar=b sum_(j=1,2)(c_j* f_j+f_j* c_j) tensor I_spectators.
```

Each even two-mode hopping has spectrum `{-b,0,0,b}`; the two terms
commute. Thus `||Hstar||=2b=sqrt6 v/a` on the complete Fock space. This
is a norm reduction by a CAR basis change, not a four-mode truncation of H.
Every incident hop flips P0, while every other source term commutes with
D0. Therefore

```
{Hstar,P0}=0,
||[H,D0]||=||Hstar P0||=sqrt6 v/a,
||[H,P0]||=2sqrt6 v/a.                                  (1)
```

The origin radial mass is zero. Put

```
Hfar=H-Hstar-W0,    B=Hstar+W0,
beta^2=6v^2/a^2+U^2/16.
```

Order the origin CAR modes before the remote modes. Every term of Hfar
is even and supported away from the origin, so it is the identity on the
origin factor. Anticommutation in (1) gives

```
B^2=Hstar^2+(U^2/16)I,       ||B||=beta.                  (2)
```

Neither the radial mass nor the number of sites occurs in (1) or (2).
The one-dimensional chain controls have `6` replaced by `2`. The exact
16-row CAR normal form and the actual 54-row source block are checked
independently in `algebra.py`. This volume independence concerns readout
error; the parent's interacting spectral-state preparation bound has its
own finite-volume restrictions.

## Interleaved complete native events

First define an event process without a continuous-pulse assumption. Between
complete native readout events, evolve by the actual `exp(-i H dt)` for
arbitrary nonnegative intervals with total T. Each event is the parent's
complete instrument, including its native failures, fresh probe registers,+nine coins and exact classical integer-repair map. All outcomes, prefixes
and halted suffixes are retained. Every native readout event has complete
half-diamond error below delta=10^-9 relative to its ideal instrument.
An untouched reference is permitted throughout.

Replace each complete native event by its ideal instrument. Contractivity
and telescoping cost at most N delta in half-diamond norm. In the ideal
process, replace H by `Hrest=H-Hstar`. Duhamel on each free interval gives
operator error at most `dt ||Hstar||`; unitary-channel half-diamond error
is no larger than that operator error. Complete instruments are contractions,
so summing the actual elapsed intervals costs at most `T sqrt6 v/a`.
The ideal Kraus operators

```
K0=(7I-6D0)/sqrt50,       K1=i(I+6D0)/sqrt50
```

commute with Hrest. Its common final unitary can therefore be moved past
the entire ideal frozen-reference instrument. The exact event-process bound is

```
half_diamond(actual,
  Ad(exp(-i Hrest T)) after ideal_reused_instrument)
 <= min(1,N delta+T sqrt6 v/a).                           (3)
```

This compares the complete field, reference and archive, rather than only
one parity probability. The final source state is not compared to an
unevolved field. The event model admits actual drift between every event;
it does not assign an internal duration to an abstract native event.

## Drift during timed local controls

For a continuous-pulse implementation, additionally specify a common source
time and a local apparatus realization of the native control word. On a
joint source/apparatus Hilbert space its Hamiltonian is

```
H tensor I_apparatus+C(t),
```

where C(t) acts only on the origin and apparatus. Its zero-drift complete
instrument must be the certified native word, including import/export,
all complement outcomes, fresh coins and retained records. Every elapsed
control, transfer and controller interval is included in T. Any physical
port encoding of the source factor must preserve the remote tensor factor
and extend B with norm at most beta. An arbitrary unbounded drift on native
complement states does not meet this interface. The source's angle-only
API supplies no physical pulse duration or bandwidth by itself.

Compare this process to the same local controls with Hfar in place of H.
Equation (2) and time-dependent Duhamel give an operator error at most
`beta dt` on every continuous interval, independent of control strength.
All intermediate complete measurements, native integer repairs and failure
continuations are the same contractions on both sides. Their records are
retained. A common fixed terminal time T handles declared halting: pad the
halted branch with field evolution and no further readout until T, on both
sides. The same local failure continuation and source clock are used. The accumulated complete half-diamond error is at most beta T.

Since Hfar commutes with every allowed control, it factors from the
zero-drift apparatus process. Comparing that process to the ideal repeated
instrument costs at most N delta, using the parent's complete control
bound. Thus

```
half_diamond(actual_timed,
  Ad(exp(-i Hfar T)) after ideal_reused_instrument)
 <= min(1,beta T+N delta).                               (4)
```

The diamond norm without the factor one half is at most twice (4).
In particular the half trace-norm bound holds for every source input and
arbitrary untouched reference. Hstar and W0 are omitted only from the
comparison generator; the actual timed process includes the full H.
No Hamiltonian compensation is used. Equation (4) is conditional on the
stated timed local realization, rather than a derivation of hardware speed
or an unspecified extension of the source to erased port states.

## A positive-time one-preparation classifier

Use the actual source ground mixture or either unique charged singlet at
`a=v=1, m=64, U=1/65536`, on the radius-one cube. Their parent parity bounds
give wrong-parity probability at most

```
e=188263256928/480690783803197.
```

Take N=32, score every failure and unused halted slot zero, and classify as
charged at pointer count at least 16. This is the parent's total scoring
rule; no survival branch is conditioned upon. The frozen ideal classifier
has error less than `e+10^-6`, by the retained exact binomial and rational
Hoeffding certificate. Its records are unaffected by the common final
Hfar unitary. Equation (4) therefore gives

```
error < e+10^-6+32*10^-9+beta T.
```

For the explicit positive interval `T=1/10^6`, (2) gives `beta<5/2` and

```
error < e+10^-6+32*10^-9+1/400000 < 1/2500.              (5)
```

All inequalities are rationally certified in `bounds.py`; floating matrix
calculations are independent controls. The classification uses one supplied
field preparation. It does not resolve the sign of charge, preserve the
initial spectral sector, or infer an ensemble parity expectation from one
field. Under actual drift, the exact record law need not be the frozen
latent-parity mixture. Equation (4) controls its difference. Each native
load coordinate lies in `[0,2]`, so its mean differs from
`[(1-q)M0^N+qM1^N]x0` by at most twice the right side of (4).

## Commands, records and energy

`commands.py` retains the same 40 parity pulses and nine four-pulse coin
words per round. It assigns each of the 2,432 pulses duration
`T/(2*2432)` and each round a nonpulse budget `T/(2*32)`. The latter must
cover transfer, measurement, integer-repair and controller delays. Each
actual binary64 integrated native port coefficient is represented by its
exact rational value; dividing by the positive pulse duration gives an
explicit finite rate. Integrating that rate recovers the unchanged source
coefficient exactly. The certificate records every one-round command and
the reindexing rule for 32 fresh apparatus blocks. It retains 736 fresh
logical bits, 32 integer repair attempts and the whole finite archive.
Meeting the displayed rate and latency budget is an apparatus input.

For the ideal event model, free H evolution conserves field energy, and
each readout changes it by `-(18/25)<Hstar>` in its actual premeasurement
state. Thus the total measurement work obeys
`|W|<=N(18/25)sqrt6 v/a`; replacing those premeasurement states by the
initial state is unjustified. A continuously driven implementation also
exchanges work with its controls. The frozen formula with `(7/25)^N` is
not asserted for that process, and no control energy is labeled heat.

The proof does not include uncontrolled field-repair operations that move
origin parity. The recorded integer work-register repairs are the admitted
parent law and are all retained. Source preparation, texture selection,
physical clock calibration, long-time protection, continuum field theory
and an Einstein limit are not inferred from the finite-time readout bound.
