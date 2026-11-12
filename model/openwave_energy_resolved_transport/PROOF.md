# Energy-filtered charged arrivals on the same source

## Preparation, source and measured effects

Use the unchanged N=2 Peter–Weyl119 source
`H=H0+lambda T+4N_ref`, with `lambda=2^-34` and `||T||<=tau=360`.
The physical ground preparation, retained environments, static thirty-mode
CAR reference at cell zero, complete localization measurement
`{P0,I-P0}` and elementary removal pump `B=B_ec` are exactly those of
[`interacting_charged_transport`](../coherent_defect_wave/chiral_joint/interacting_charged_transport/PROOF.md).
Here P0 requires every cell in its true local ground doublet and all links
blank. On this source `P0 B* B P0=P0/2`. No interacting ground is silently
replaced by an atomic input.

After pump success, measure the existing full-H energy-bin projector E,
centered on the lowest visible atomic e^c pole and of radius 1/1024.
Its complement is retained. Apply the same clock `U(t)=exp(-itH)`.
The six existing arrival projectors P_e require one specified incident
link nonblank, every other link blank, and the true local ground at
every nonarrival cell. They retain the entire arrival-cell CAR/scalar
space and the reference. Their full complement is retained.

Resolve each P_e further by the actual local atomic charged-energy
projector S_a at the arrival cell. Write `F_e=P_e S_a,y`. The detector is
the complete family

```
{F_e (six routes), P_e-F_e (six routes), I-sum_e P_e}.
```

All these are gauge-invariant effects on the inherited physical sector.
Let a denote the selected atomic total energy, including the seven other
local ground energies and the static reference cost four. The tagged
arrival has exactly one excited link of energy one, hence

```
H0_total F_e=(a+1)F_e,       F_e T F_e=0.             (1)
```

The second identity follows from the matter projectors: every hop has
two distinct matter endpoints, at least one of which is a nonarrival
cell compressed to its true gauge-singlet ground. A single nontrivial
fermion operator has zero singlet compression there. No absence of
nonblank-to-nonblank Peter--Weyl link transitions is assumed.
In particular `F_e H F_e=(a+1)F_e`. The tag is an active
measurement: its conditional output has this energy expectation and
need not remain in the original full-H energy bin. No energy-free
position measurement or automatic work-free transport is assumed.

## The actual pole residue fixes the first-hop coefficient

The native 17-row and 7-row elementary-removal blocks give the lowest
visible pole as the smallest root of

```
x^4-2x^3-212x^2+213x+11114.
```

At that pole, the exact native down-removal residue is twice the
up-removal residue. The executor establishes this by divisibility of
the difference of their rational resolvent numerators by the displayed
polynomial, over the actual neutral-ground algebraic field. The pole
occurs in the raw-lepton-number-two spin-one block and not in its
spin-zero block. This is checked on the actual characteristic polynomials.

The existing local SU(2) Ward identity fixes the remaining tensor
components. If w is the complete pump probability at this pole, then

```
P_g c_s* S_a c_t P_g
 = w delta_st I -(2w/3) sum_i sigma_i,ts S_i.          (2)
```

Indeed the up and down diagonal weights on the spin-up ground are
`2w/3` and `4w/3`. Covariance determines the off-diagonal entries.
Summing the retained reference spin gives `P_g B* S_a B P_g=w I`;
each other atomic pole likewise has its actual scalar weight w_j on
the full ground doublet. These identities also hold on arbitrary
entanglements of the eight ground doublets.

Insert (2) at both ends of the actual hopping term, retaining its
Pauli matrix, coefficient magnitude 1/4, the pump factor 1/sqrt(2),
and the actual e^c Peter–Weyl blank-link Gram. The three exact four-row
matrix contractions yield

```
(F_e T S_a B)* (F_e T S_a B)
 = w^3 [I/12+(2/27) S_0 dot R_axis S_y]              (3)
```

on the two relevant ground spins. The bond-axis rotation and the actual
native cube trial give, with `x_cube=4 epsilon_cube`,

```
m=||F_e T S_a B u||^2=w^3(1/12+x_cube/648)>0.         (4)
```

The source pole weight and cube energy are enclosed by rational root
intervals. No measured mass, fitted residue or separate transport
Hamiltonian enters this coefficient.

## A sharp projector tail without a global charged-band hypothesis

Here a global isolation theorem for the complete charged energy band
is unnecessary. Set `v=lambda tau` and let `g>0` be the distance from
`a+1` to the upper edge of the actual bin. Equation (1) and the
spectral theorem give

```
||F_e E|| <= beta=v/g.                               (5)
```

For any normalized input z, write `y=B P0 z=sum_j x_j`, where x_j is
its exact atomic component of energy a_j. Its norm is at most
`sqrt(w_j)`; index zero is the selected pole a. Let d_j be the distance
from a_j to the appropriate bin edge. The inherited spectral-theorem
bound gives

```
||(I-E)x_0|| <= (v/d_0)||x_0||,
||E x_j|| <= (v/d_j)||x_j||, j!=0.                  (6)
```

Since `[E,H]=0`, applying F_e to the commutator with H and using
`b-a=(a+1)-a=1` gives the exact identity

```
F_e E x_0
 = -lambda F_e T x_0
   +lambda F_e E T x_0 +lambda F_e T(I-E)x_0.        (7)
```

For the other atoms, apply (5) after (6), gaining two factors of v.
Combining every atom therefore bounds

```
||F_e E y +lambda F_e T x_0|| <= R,
R=lambda tau (beta+v/d_0) sqrt(w_0,+)
  +beta v sum_(j!=0) sqrt(w_j,+)/d_j.                (8)
```

Every full-H state outside the retained bin stays in its complementary
outcome. In particular, this proof assumes neither a global gap around
the atomic pole nor the absence of unseen charged states inside E.

The prepared purification is within `eta=eta_g+epsilon_p` of the native
cube trial and an idle environment after the common contraction P0.
Only the selected native edge contributes to `F_e T S_a B P0`, whose
norm is at most 15. Its preparation error is consequently at most
`15 lambda eta`, rather than a global probability error. Equations
(4) and (8) enclose the actual initial tagged amplitude by

```
lambda sqrt(m_-) -R-15lambda eta
 <= ||F_e E B P0 Psi||
 <= lambda sqrt(m_+) +R+15lambda eta.                (9)
```

The lower endpoint is strictly positive for the frozen native source.
Squaring and multiplying by the retained preparation-success lower
bound gives an original-input probability. No intermediate successful
branch is silently normalized.

For each of the six routes, the original-input time-zero probability
lies between `1.04092e-23` and `1.68686e-23`. The exact rational
endpoints and their complete six-route union are retained in the receipt.

## Time-zero dressing and complete time comparisons

The nonzero arrival in (9) is present at t=0, immediately after energy
filtering. Without the filter, `F_e B P0=0` because all links are blank.
With zero hopping, the atomic energy filter leaves them blank as well.
The positive result is thus a virtual tail resolved by the measurements;
it is not a proof of propagation.

A time comparison is bounded independently. Define

```
C=v[1/sqrt(2)+sum_(j!=0) |a_j-a| sqrt(w_j,+)/d_j].
```

For every P0 input, `||(H-a) E B P0 z||<=C`. As E commutes with U,

```
||F_e [U(t)-exp(-ita)] E B P0 z||
 <= min(beta |t| C, sqrt(2) beta).                  (10)
```

This retains the reference phase four and the same absolute atomic
energy a. Combining (9) and (10) gives the timeline intervals. The
probability change from t=0 is bounded by `d(2A_+ +d)`, where d is
the right side of (10) and A_+ is the initial upper amplitude. A positive
time change is not inferred when these intervals include zero change.
At the existing transport comparison time `t=2^-24`, the per-route
probability change is at most `2.24579e-34`, retaining the much larger
time-zero filtered tail. Positive joint arrival therefore cannot be used
by itself as a time-of-flight result.

Let q be preparation success, a0 the P0 success conditional on q, b the
joint pump/bin weight conditional on q, and r_e,s_e the unnormalized
tagged and untagged route weights after the clock. The seventeen records
have probabilities

```
1-q; q(1-a0); q a0/2; q(a0/2-b);
q r_e (six); q s_e (six); q(b-sum_e r_e-sum_e s_e).
```

They sum exactly to one. The complete reference, second source,
environments, energy complement, screened states and other link patterns
are retained. This establishes a joint energy/arrival experiment on the
same Hamiltonian. It does not establish an autonomous free charge,
relativistic mass shell, continuum limit or physical electron dictionary.
