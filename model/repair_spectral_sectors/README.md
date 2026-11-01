# Repaired spectral sectors and exact energy transfer

The six-site field's protected spin is only part of its particle response.
The complete `N=4, S=1` sector splits into **90 repair-active and 15
repair-transparent multiplicities**. The former has no nonconstant
multiplicity observable conserved by both the joined Hamiltonian and every admitted
repair. The latter carries an explicit, nontrivial finite spectrum that
survives repair exactly, through its outside outcome. This positive sector
is inaccessible to simple pair-removal probes of an eta-singlet ground.

These results refine [common spin protection](../common_spin_protection/README.md)
on its unchanged six-site action and complete repair law. They advance
CP-03/#1139 and CP-07/#1143, with an energy-accounting input for CP-04.
No coupling, reference column, outside record or primitive is replaced.
No physical species, absolute mass or localized protected particle is
identified.

## Exact active/dark classification on one declared sector

Use the same centered generator `K=H-N/2`, with `c=1/8`, `U=1`, and the
same hidden spin and eta algebras. On six sites, `N=4` means `eta_z=-1`.
The exact representation count gives

```
N=4, S=1:  (multiplicity90, eta=1) + (multiplicity15, eta=2),
total dimension = (90+15)*3 = 315.
```

In the highest-spin weight, the 120-dimensional occupation space is cut
by `S_+ v=0` and `eta_- v=0` to the active 90-dimensional space. The three
old localized pair columns meet both conditions exactly. In the Hubbard
frame the unnormalized hole is `|up,0>-|0,up>` and the ground is
`S+(4-sqrt(17))D`, with
`S=|up,down>-|down,up>` in the fixed CAR address convention and `D` the
sum of the two doublons. The implementation binds these signs to the
existing analytic columns. The numerical source dictionary uses normalized
columns separately. Its exact ground ratio is
`b/a=-4t/(sqrt(U²+16t²)+U)=4-sqrt(17)` for `t=1/16,U=1`.
The full global frame adds pair-column signs `(-,+,-)`; the receipt checks
all three normalized source columns against these exact seed formulas.

Every complete repair support lies in eta one; its movement operators
preserve that subspace. Its outside operator is `Q=I-P`. Consequently the
instrument preserves `eta²` and `eta_z` on this declared sector. It does
**not** claim full eta raising/lowering covariance of the number-selective
repair. The Hamiltonian has the stronger existing eta covariance.

The native triangle `(0,1,5)` gives the actual local block seams. All 72
C/E raw slots, including zeros, have an exact three-state partial-injection
dictionary, checked against the full-Fock operators. Their effects separate
the three occupied pairs; connected movements generate every corner
matrix unit `|v_i><v_j|`. No projection or zero record is omitted.

The joined Hamiltonian makes this corner cyclic on all 90 active
multiplicities. The exact certificate constructs `16K`, with integer
entries, and the unnormalized seed columns in `Z[sqrt(17)]`. It reduces
`sqrt(17)` to 44 modulo the prime 101; `44²=17 mod101` and the
discriminant 68 is nonzero. There are no rational denominators in these
rank inputs, and the seed norm is nonzero modulo 101. Powers zero through
42 give 129 Krylov columns; a retained **90-by-90 nonzero minor** proves
rank at least 90 in characteristic zero. The exact spin/eta count gives
rank at most 90. A failed pivot or small numerical singular value supplies
neither bound. Independent row elimination verifies the minor; removing
the inter-dimer hopping leaves rank three.

Products `K^a |v_i><v_j| K^b` then span all active multiplicity matrix
units. Thus the common generated algebra on this sector is `M90 tensor
I3`, with commutant `I90 tensor M3`. No nonconstant multiplicity energy
observable is conserved by both laws. This statement does not assert that
every single repair changes every state, or classify the entire Fock
space. In particular, it must not be extended to the dark complement.

## A positive repair-transparent spectral sector

For each all-up two-particle occupation column `|x up,y up>`, set

```
D|x,y> = eta_+ |x up,y up> / 2,  0<=x<y<6.
```

`D*D=I15` exactly. These are the eta-two, spin-one highest-weight columns
at `N=4`. Their spin-lowered partners give all 45 states. The eta symmetry
and the spin-polarized no-doublon reference give the exact intertwiner

```
K D = D (A_pair/16 - I),
```

where `A_pair` is the integer two-fermion nearest-neighbor adjacency of
the open six-site chain. Their overlap with each of the three repair
support columns vanishes as an identity in `Q(sqrt(17))`. All native
outcomes therefore annihilate these states, while the complete outside
outcome is identity. Spin covariance extends this statement to arbitrary
spin/reference states. This is a real spectral sector, although its
protection here is repair transparency rather than a native defect readout.

The exact energies are

```
E_pq = [cos(p*pi/7)+cos(q*pi/7)]/8 - 1,  1<=p<q<=6.
```

For `y=16(K+1)`, the integer characteristic polynomial is

```
y^3 (y^3-7y-7)(y^3-7y+7)
    (y^3-2y^2-y+1)(y^3+2y^2-y-1).
```

The four cubics are squarefree, mutually coprime and nonzero at zero.
There are 13 distinct energies: `K=-1` has multiplicity three; the other
twelve have multiplicity one. Including spin, the dark common commutant
is `M9 + 12 M3`, of complex dimension 189. Active/dark cross terms vanish
in the common commutant because the repair corner is cyclic on the active
space. The complete common commutant on `N=4,S=1` therefore has dimension
**198**, including the active `M3`. This is not an all-Fock maximality claim.

A coherent superposition of mode pairs `(1,2)` and `(1,3)` has nonzero
frequency about `0.05012110849`. The complete matrix experiment evolves
it and an untouched reference for model duration `1.25`, then attempts
E, C and E repair on all three local seams. The complete stopping instrument
terminates on the first native record and otherwise continues its
unnormalized outside branch: all 40 possible first native records and the
final outside record are retained, without postselection. Its phase and full
state/reference correlations survive; every nonzero outcome is outside.
This evaluates the declared ideal matrices, not a new native-pulse
preparation or detector implementation.

There is also an exact access constraint. For any eta-singlet reference,
`eta_-|g>=0`, and every two-annihilator probe commutes with `eta_-`.
Its `N=4` response is therefore lowest weight at `eta_z=-1`, hence eta
one, with zero eta-two component. A complete numerical check retains all
66 two-annihilator probes of the previous supplied ground; both its
eta-singlet residual and dark overlaps obey the stated admission. The
dark pair cannot be renamed the previously visible ground pair pole.
Its explicit eta-raised preparation is a positive mathematical access
control, not an implemented source service or a physical particle claim.

## Exact repair energy transfer and complete response law

For any complete Kraus instrument and a fixed reference generator,

```
W = sum_r L_r* K L_r - K
  = (1/2) sum_r (L_r* [K,L_r] + [L_r*,K] L_r).
```

For the actual unital repair,
`Tr(K W) = -(1/2) sum_r ||[K,L_r]||_HS²`. The latter equality needs
unitality; a complete nonunital amplitude-damping control fails it. For
any normalized input, including correlations with an untouched reference,
the field energy change is exactly `Tr(rho W)`. Orthogonal raw/outside
records remain in the sum. This is energy transferred into the modeled
field; naming it heat or controller work needs a controller Hamiltonian,
which is not supplied here.

On the original two-dimer charged-pair carrier, every native branch holds.
The **complete full-Fock** instrument is exactly the pinching `{P,Q}`, so
`W=-[P,[P,K]]`. The previous leakage identity gives

```
||W||_HS² = 3/128,
spectrum(W) = ±sqrt(3/512) (once each),
              ±sqrt(1/512) (three times each), 0 (248 times),
max_rho |Tr(rho W)| = sqrt(6)/32.
```

The integer Gaussian-dyadic certificate checks these constants against
the actual joined action and full repair. Input states can gain or lose
energy; the transfer is not a positive rest-mass functional.

On the six-site active sector, the receipt retains the complete 90-by-90
matrix `T_nm=sum_r |<n|L_r|m>|²`, including its outside contribution.
Every row and column sums to one. One complete E repair raises the lowest
active energy by approximately **0.06624866367**. All 90 retained
Hamiltonian eigenvectors have a resolved second energy-change moment in
this numerical calculation, the smallest about `6.2010e-7`. These
floating spectral diagnostics have admission bounds; the exact algebra
result above is the independent, basis-free common-conservation statement.
No source energy measurement, stable repaired mass pole, physical mass
calibration or three-dimensional localized defect follows.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.repair_spectral_sectors --check model/repair_spectral_sectors/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/repair_spectral_sectors
```
