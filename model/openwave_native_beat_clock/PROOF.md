# Free native beats with retained charge references

## Actual source and transition algebra

The source is the same complete eight-cell PW119 Hamiltonian
`H=H0+lambda T+4 N_ref`, with `lambda=2^-34` and
`||lambda T||<=v=360*2^-34`. It retains the full thirty-mode CAR
reference at cell zero. The preparation begins with the existing physical
ground preparation followed by the complete `{P0,I-P0}` measurement.
P0 requires the true local ground doublet at every cell and blank links.
Here P also includes the reference-vacuum projector. Every preparation
failure and untouched environment is retained.

For r=e^c or L, let B_r be the actual native source-to-reference removal
pump and S_r its selected local atomic charged-energy projector. The
existing exact local resolvents supply a strictly positive weight w_r.
The pump intertwines all twelve joint source/reference currents and all
three spin generators. S_r commutes with them. Schur's lemma on the
actual two-dimensional local ground spin therefore gives

```
P B_r* S_r B_r P=w_r P.                              (1)
```

Its scalar is the complete native pole weight, not a fitted coupling.
The e^c spin-resolved residue identity is executed in
[`openwave_energy_resolved_transport`](../openwave_energy_resolved_transport/PROOF.md).
The L weight comes from its complete native resolvent and the same
source covariance. Equation (1) holds on arbitrary states and
entanglements of the eight ground doublets.

Define the finite source operator

```
A_r=S_r B_r P/sqrt(w_r).
```

It is a partial isometry: `A_r* A_r=P`, `A_r A_r*=Q_r`, and `A_r^2=0`.
The latter follows because its input reference is vacuum and its output
reference has one particle. Its coefficients, source domains and target
subspaces are specified by the actual native pump and local energy
projector; this is not a new freely chosen two-level Hamiltonian.

Write A=A_ec and Q=Q_ec. The operators

```
X=A+A*,       Y=i(A-A*),       Z=P-Q
```

obey `X^2=Y^2=Z^2=P+Q` and `[X,Y]=2iZ`, with zero action on the retained
orthogonal complement. A single initial source unitary

```
U_prep=I-P-Q+(P+Q+A-A*)/sqrt(2)                      (2)
```

maps every input psi in P to `(psi+A psi)/sqrt(2)`. It preserves joint
Gauss constraints and total parity. It is an explicitly declared finite
coherent service, computed from the native operators; an efficient
primitive-pulse compilation or physical implementation error is not
claimed. It is applied once. All subsequent evolution uses H without
phase drives, repeated swaps or reference rotations.

If E0 is the neutral local ground energy and e_r the selected charged
local energy, set `Delta_r=e_r-E0`. On P the atomic energy is `8E0`;
on Q it is `8E0+4+Delta_ec`. Therefore, with `Omega=4+Delta_ec`,

```
[H_atomic,A]=Omega A,
i[H_atomic,X]=Omega Y,     i[H_atomic,Y]=-Omega X.    (3)
```

The initial expectation of X is one and of Y zero, independently of
the input spin state or external reference. Atomic free evolution gives
`<X(t)>=cos(Omega t)` and `<Y(t)>=-sin(Omega t)`.
Since `P T P=Q T Q=0` and T preserves reference number, the mean work
of (2), relative to the localized input, is exactly `Omega/2` even for
the full H. The initial coherent preparation is not an autonomous
energy source. A stationary input without such a preparation cannot
self-start: `[rho,H]=0` makes every fixed-observable expectation constant.

## The same full-H bins and a common herald

Let E_n be the actual full-H bin of radius 1/1024 about the neutral
atomic energy, restricted to reference vacuum, and let E_c be the
existing charged bin, restricted to the one-e^c reference family. The
sum `F=E_n+E_c` has one common success flag. It does not reveal which
arm is occupied, so it retains their coherence. Measuring the two arms
as separate classical outcomes would be a different experiment.

On either atomic subspace the full-H residual has norm at most v. If d
is the smaller distance of its atomic energy to the appropriate bin
edge, the spectral theorem gives

```
||(I-F)(P+Q)||<=eta=v/d,
Pr(F)>=1-eta^2.                                    (4)
```

The two reference sectors are orthogonal and invariant, so the maximum
of the two errors suffices. The normalized successful preparation is
within unhalved trace distance `2eta` of the input to F, including its
purifying environment. The nonselective measurement `{F,I-F}` preserves
mean H energy because both effects commute with H. Its failure output
is retained.

After this optional but executed common herald, every frequency in the
finite spectral expansion of the cross-reference observable lies in

```
I_c-I_n=[lower(I_c)-upper(I_n), upper(I_c)-lower(I_n)]. (5)
```

This is an exact frequency-support statement on the actual full H,
not an isolated charged-pole assertion. The support window is independent
of initial P0 spin populations. The prepared coherence has a positive
contrast bound; the theorem does not claim that every arbitrary state
in those bins has nonzero visibility or a unique sinusoidal frequency.

## Autonomous waveform and disturbance bounds

For any bounded observable O, unitary Duhamel comparison gives

```
||exp(itH) O exp(-itH)-exp(itH_atomic) O exp(-itH_atomic)||
 <=2v |t| ||O||.
```

Apply it to the fixed X and Y, each of norm one, and include (4):

```
|<X(t)>-cos(Omega t)| <= min(2,2eta+2v|t|),
|<Y(t)>+sin(Omega t)| <= min(2,2eta+2v|t|).          (6)
```

An additional error epsilon in unhalved trace norm of the normalized
accepted clock state, after the common herald, and a Hermitian generator
disturbance of norm k add `epsilon+2k|t|`. Errors specified before a
conditioning step require their own conditioning bound. These
are stated disturbed constructions; their Hamiltonians need not retain
the original exact frequency window (5).

The four fixed effects `(I+X)/4,(I-X)/4,(I+Y)/4,(I-Y)/4` are positive
and sum to identity on the entire carrier. They assign uniform outcomes
on the orthogonal complement rather than discarding it. The code bounds
their probabilities using rational root intervals, a Machin-series pi
enclosure, exact phase reduction and an outward sine/cosine Taylor kernel.
Large source times therefore do not rely on an uncertified floating
trigonometric evaluation.

The source values give `Omega approximately 4.691656680938395`; the
reference-subtracted gap is approximately `0.691656680938395`.
At one million source-time units the waveform bound still controls over
700000 atomic total-beat cycles. Reference cost four is part of the
measured beat; subtracting it is a declared model energy comparison,
not a continuing compensating phase drive. No laboratory second,
physical Compton frequency or independent absolute rest mass is inferred.
If the initial coherent preparation is omitted and the reference stays
in vacuum, both X and Y are exactly zero at every time: H preserves
reference occupation. This is an unchanged-source control.

## An interchannel beat cancels the common reference cost

The two actual reference families are orthogonal, so
`A_ec* A_L=0`. The native operator

```
D=A_ec A_L*,       D*D=Q_L,       DD*=Q_ec,       D^2=0
```

is consequently a partial isometry between two joint physical source
subspaces. Both have one reference particle and even total parity.
Both obey every joint Gauss constraint; no unphysical superposition of
uncompensated local charges is introduced. The Hamiltonian identity is

```
[H_atomic,D]=(Delta_ec-Delta_L)D.                   (7)
```

The reference cost cancels algebraically because `[N_ref,D]=0`. There
is no changed reference Hamiltonian. A finite initial swap from P into
Q_L, followed by a quarter-turn mix of Q_L and Q_ec, prepares
`(A_L psi+A_ec psi)/sqrt(2)`. Their mean initial energy increase relative
to P is `4+(Delta_ec+Delta_L)/2`; cancellation in the frequency does
not remove that preparation work.

Replace A by D and the two bins by their actual L and e^c bins in
(4)--(6). The resulting free frequency is approximately
`0.594768183046114` in the same source units. The same seven complete
records apply, with one common herald preserving reference-family
coherence. Dephasing or measuring that family removes both beat
quadratures at every subsequent time, because H preserves reference
family and D is off-diagonal between the two families. This is an exact
control. The L pole is a native doublet transition, not an assigned
neutrino or muon mass.

## Complete records and interpretation

Let q be physical ground-preparation success, a the subsequent P0 success,
and h the common-bin herald success. If p_i are the four unnormalized
terminal weights, `sum_i p_i=h`. The seven original-input records are

```
1-q; q(1-a); qa(1-h); qa p_0; qa p_1; qa p_2; qa p_3.
```

They sum exactly to one, with the source copies, quantum reference,
environments and failures retained. No operation depends on elapsed time
after preparation; the readout effects are fixed. This is a native
prepared free beat with state-independent atomic frequency and explicit
full-H errors. It supplies neither passive self-starting, frequency
locking under arbitrary source-parameter changes, an asymptotically
free charged particle, nor a physical time/energy dictionary.
