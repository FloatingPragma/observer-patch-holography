# Complete records of a finite quantum-Higgs transition

The declared finite quantum-Higgs construction has an exactly measurable
neutral transition on its complete 307-state `N_F-2N_H=0` sector. The same
source operator supplies its ground preparation pulse, evolution, occupation
pointer, spectral gap and residue. All 305 dark states are retained. The two
complete occupation-read protocols at times `t` and `2t` recover the gap from
their unconditional joint counts, on an explicit unaliased time window.

These are observer outputs of a chosen architecture. The family coefficient,
number cost, coherent preparation service and program clock are declared
inputs. Their calculation requires no proof that the architecture is unique.
No measured mass or calibrated energy conversion enters the result. The
transition is a finite gauge-invariant scalar/fermion conversion line; no
relativistic particle identity or GeV mass is assigned to it.

The owning executable is [`higgs_records.py`](higgs_records.py). The
[readout contract](READOUT_CONTRACT.md) attaches the complete logical
instruments to declared apparatus or the existing finite source compiler,
including references, outside inputs and implementation errors.

## One actual source carrier and Hamiltonian

The source is the literal finite quantum scalar and full CAR vertex of
[`INTERACTING_CHIRAL_JOIN.md`](../particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md),
with the sector identified in
[`PROOF_PHASE_ATTENUATION.md`](../particle_calibration_selection/PROOF_PHASE_ATTENUATION.md).
The scalar states are its vacuum and conjugate doublet. Eighteen active CAR
modes carry two spins of the three internal charged-Higgs components for
each of three declared family coordinates. For an exact nonzero matrix `Y`,

\[
 W(Y)=\sum_{h,f,g}Y_{fg}B_h^\dagger
 (c_{L_h,f,\uparrow}^\dagger c_{e,g,\downarrow}^\dagger
 -c_{L_h,f,\downarrow}^\dagger c_{e,g,\uparrow}^\dagger),
 \qquad H=W+W^\dagger+\epsilon_FN_F+\epsilon_HN_H.       \tag{1}
\]

The actual source Ward identity gives twelve conserved gauge generators.
The number `Q*=N_F-2N_H` is conserved and is distinct from electric charge.
Spatial hopping, scalar transport and quantum-link interaction terms are
excluded from this declared local Hamiltonian. They require their own
output calculation when included in an architecture.

The complete `Q*=0` sector contains the vacuum `Omega` and every occupied
scalar with two fermions, of dimension `1+2 binomial(18,2)=307`. Define

\[
 v=W\Omega,\quad g^2=\langle v,v\rangle
 =4\operatorname{Tr}(Y^\dagger Y)>0,\quad b=v/g,
 \quad P_0=|\Omega\rangle\langle\Omega|,
 \quad B=|b\rangle\langle b|,
 \quad D=I-P_0-B,
 \quad d=2\epsilon_F+\epsilon_H.                          \tag{2}
\]

The scalar cutoff forbids another creation step after `W Omega`, while
`W†` sends every two-fermion/occupied-scalar input either to the vacuum or
to zero. Consequently the actual full sector operator is

\[
 H=|v\rangle\langle\Omega|+|\Omega\rangle\langle v|
       +d(B+D).
                                                               \tag{3}
\]

Set `Delta=sqrt(d²+4g²)`. Its spectral projectors and energies are

\[
 \Pi_-={1+d/\Delta\over2}P_0+{1-d/\Delta\over2}B
       -{ |v\rangle\langle\Omega|+|\Omega\rangle\langle v|\over\Delta},
 \qquad \Pi_+=P_0+B-\Pi_-,
\]
\[
 H=E_-\Pi_-+E_+\Pi_++dD,\qquad E_\pm=(d\pm\Delta)/2.
                                                               \tag{4}
\]

Their ranks are `1,1,305`. `E_-<min(0,d)`, so `Pi_-` is the ground of
the whole 307-state sector. Equations (3)-(4) are checked against every entry
of the actual occupation-basis source block, rather than introduced as a
separate two-level Hamiltonian. The other sectors of the full
`3*2^18=786432`-state scalar/CAR carrier are not assigned spectra by this
calculation. Their complement must be retained by any larger preparation
or measurement contract; the declared readout contract does this explicitly.

## Complete preparation by the same phase-shifted vertex

The source vertex at coefficient `iY` is
`J=V(iY)=i(W-W†)`. It is Hermitian, gauge invariant and preserves `Q*`.
It also satisfies `J=exp(i pi N_H/2) V(Y) exp(-i pi N_H/2)`, so the
same family coefficient and a scalar-number phase gate supply the pulse.
Choose the declared preparation angle `eta` by

\[
 \alpha={1+d/\Delta\over2}=\cos^2\eta,\qquad
 \beta={1-d/\Delta\over2}=\sin^2\eta,\qquad 0<\eta<\pi/2.
\]

The preparation pulse has exact complete-sector action

\[
 U_{\rm prep}=e^{+i\eta J/g}
 =D+\sqrt\alpha(P_0+B)
    +{\sqrt\beta\over g}
       (|\Omega\rangle\langle v|-|v\rangle\langle\Omega|).
                                                               \tag{5}
\]

It sends the empty vacuum to `sqrt(alpha) Omega-sqrt(beta)b`, whose
density is `Pi_-`, and fixes every dark input. The code verifies unitarity
on all 307 states and the direct equality with the actual phase-shifted
source vertex. No input or environmental record is reset or discarded by
this pulse. The initial empty vacuum is a specified gauge-invariant
preparation.

The number term is paused or refocused during (5). Preparation is a driven
apparatus operation; free evolution under (1) does not prepare its own
ground from a different-energy state. The coherent-service implementation
and its resource/error accounting belong to the readout contract. An
autonomous ground-relaxation mechanism is unnecessary for this declared
experiment and is not asserted.

## The same clock and a complete gauge-invariant pointer

Evolution uses the entire sector propagator

\[
 U(t)=e^{-iE_-t}\Pi_-+e^{-iE_+t}\Pi_++e^{-idt}D.         \tag{6}
\]

The occupation observable is `N_H=I-P0`. It and all three projectors
`P0,B,D` are gauge invariant. The binary measurement has projectors
`P0,N_H`; the refined measurement has `P0,B,D`. Their sums are the full
sector identity. The measurement is complete as an instrument; it is not
claimed to be informationally complete tomography of an arbitrary
307-dimensional density matrix.

For two binary reads separated by (6), the quantum map attached to each
ordered record `(a,b)` is

\[
 \rho\longmapsto K_{ab}\rho K_{ab}^\dagger,
 \qquad K_{ab}=P_bU(t)P_a,
 \qquad\sum_{a,b}K_{ab}^\dagger K_{ab}=I.                \tag{7}
\]

All four records are kept with their original weights. All nine records of
the refined instrument are kept, including its four zero dark/bright
transition maps. The equality is on the full matrix algebra, so (7) also
holds with arbitrary untouched quantum references.

For any admitted initial density, write
`alpha0=Tr(P0 rho)`, `beta0=Tr(B rho)`, `delta0=Tr(D rho)`. Define

\[
 T(t)={4g^2\over\Delta^2}\sin^2(\Delta t/2).
\]

Its exact complete record law is

| First, second record | Unconditional probability |
|---|---|
| empty, empty | `alpha0 (1-T)` |
| empty, occupied | `alpha0 T` |
| occupied, empty | `beta0 T` |
| occupied, occupied | `beta0 (1-T)+delta0` |

Coherences are retained by the quantum output maps; their contributions to
these four pointer probabilities cancel. In the prepared ground,
`alpha0=alpha`, `beta0=beta`, `delta0=0`. The first projective read generally
disturbs that ground, which is included in the law rather than suppressed.

The program clock is the same `U(t)` in both protocols. Repeating one chosen
evolution macro twice realizes `U(2t)`. Preparation and detector services
have their declared evolution paused or explicitly accounted; their compiler
overhead is not silently added to the parameter `t`. Seconds and GeV require
a separate physical unit dictionary.

## A finite interacting spectral observable and its count inverse

In the actual ground, the occupation mean and ordered connected correlator
are

\[
 \langle N_H\rangle=\beta,\qquad
 \langle N_H(t)N_H(0)\rangle-\beta^2
 ={g^2\over\Delta^2}e^{-i\Delta t}.                     \tag{8}
\]

Thus the neutral transition has gap `Delta`, positive residue
`g²/Delta²` and zero dark residue in this particular ground readout.
The zero dark residue is a matrix-element statement; it does not remove
dark states from the source. Equation (8) is a quantum ordered correlator,
which can be complex. It is not the classical probability of two occupied
projective records. Equation (7) gives that probability with its disturbance.

The gap is also recovered directly from complete counts. Let `J(t)` be the
unconditional probability of the record `(empty,occupied)` in (7), with
identical preparation and pointers at `t` and `2t`. If `J(t)>0`,

\[
 R={J(2t)\over4J(t)}=\cos^2(\Delta t/2),\qquad
 \Delta={2\over t}\arccos\sqrt R,
 \quad 0<\Delta t<\pi.                                  \tag{9}
\]

The alias window is part of the protocol. A source upper bound on `Delta`
can set an admissible tick without particle-data input. Equation (9) uses
frequencies of one event among all trials; it does not condition the data
on the first empty result or discard any of the other histories.

For complete probability intervals `J(t) in [p_-,p_+]` and
`J(2t) in [q_-,q_+]`, with `p_->0`, put

\[
 R_-={q_-\over4p_+},\quad R_+={q_+\over4p_-}.
\]

If `0<R_-<=R_+<1`, the exact monotone inverse gives

\[
 {2\over t}\arccos\sqrt{R_+}\le\Delta
       \le {2\over t}\arccos\sqrt{R_-}.                 \tag{10}
\]

Sampling and complete apparatus errors enter the probability intervals
before this inversion. Unresolved denominators or alias bounds are rejected
by the executable. An uncertain time interval requires its own monotone
propagation, rather than treating the chosen program time as a measured
laboratory duration.

## Exact examples and controls against dark-state loss

The declared rational fixture `Y=diag(1,2,3), d=1` uses no empirical inputs.
It gives `g²=56`, energies `(-7,8,1)`, gap15, mean occupation `7/15` and
connected residue `56/225`. At `t=pi/45`, the complete ground record law is

```text
(empty,empty)       1352/3375
(empty,occupied)     448/3375
(occupied,empty)     392/3375
(occupied,occupied) 1183/3375.
```

At twice that time, `J(2t)=1344/3375`; therefore
`J(2t)/(4J(t))=3/4`. Equation (9) returns15 on the certified source window
`Delta*t=pi/3`. This is a finite observable gap in the chosen common clock,
not a fitted particle mass.

The same executable also evaluates the actual coherent analyzer coefficient
at phase1, pass amplitude1 and `d=1`. Its unchanged native branch weights
give

\[
 g^2={425+100\sqrt2\over1944},\qquad
 \Delta^2={911+100\sqrt2\over486}.                       \tag{11}
\]

The ground, residue and full four-record probabilities at the fixed tick
`t=1/8` follow from the same maps. In particular `Delta<2`, so this tick
satisfies the required unaliased window. No rational fixture coefficient
is substituted for the source coefficient in this control.

A physical neutral dark input exists in the same sector. Choose a nonzero
family matrix `X` with `Tr(Y†X)=0`. The vector `W(X)Omega`, after
normalization, is gauge invariant by the same Ward identity and lies in
`D`. Its complete binary record is `(occupied,occupied)` with probability
one. Projecting the carrier onto `P0+B` instead would give output trace
zero. For an equal mixture of that state and the empty vacuum, silently
conditioning on the bright subspace changes the complete law by total
variation `1/2`. The tests verify these controls and retain an entangled
untouched reference through every refined measurement branch.

No conclusion about an interacting relativistic limit follows from a finite
neutral transition. Such an identification must use the same observable in
a spatial theory with its propagation, scattering and physical clock. The
finite result determines the actual outputs of the specified
preparation and pointer experiment, including their complete normalization
and controlled inverse readout.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -W error -m model.particle_observer_outputs.higgs_records
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -W error -m pytest -q model/particle_observer_outputs/test_higgs_records.py
```
