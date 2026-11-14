# Prepared quantum field states and a controlled cycle limit

The source is the pure-gauge sector of the declared fermion/plaquette
process. Every fermion mode is empty. The actual field register has one
coordinate on each of288 carrier edges. Its electric factor is
`exp(-i alpha sum E_e²/2)` and its magnetic factor on face f is
`exp(i beta_f (U_f+U_f*)/2)`. All other face coefficients are zero in the
executed cycle member. This is a specified member of the same field family,
with an explicit preparation, source word, clock and complete instrument.

## Exact field sector and preparation

The first actual face is the oriented triangle with edge coefficients
`A0[0]=1, A0[1]=-1, A0[6]=1`. It satisfies `B A0=0` and `A0* A0=3`.
The map `|z> -> |F=z A0 mod m>` is an isometry for every admitted modulus.
It intertwines the face shift with `z -> z+1` and the electric energy with
`3 center(z)²/2`. The kinetic and magnetic factors therefore preserve this
entire cycle code. No projection after evolution is used.

`native_cycle.py` prepares the pivot coordinate with actual native scalar
words, then copies its value with signed modular additions into the two
other supported edges. The decoder reverses those additions. The complete
physical program acts on the existing field allocation and a fresh clock
wire3203. It includes all288 electric words per substep. Empty matter is a
reducing vacuum of the number-conserving matter factors, so omitting their
identity action in this pure-gauge specialization is exact. This empty
state is not identified with a Dirac sea or an interacting matter vacuum.

The finite initial amplitudes are normalized integers, together with a
literal binary64 phase array. Those integers are fixed by the preparation
API. Every Gray rotation is checked against exact sums of their squared
weights using outward integer trigonometric intervals. The literal phase
word has its own operator bound. `gaussian_bounds.py` separately compares
this exact target with an infinite Gaussian using outward exponential,
algebraic-parameter and Gaussian-tail bounds. Neither the host exponential
nor a numerical residual is assumed to be an exact preparation certificate.

The emitted X and Y programs contain296440 and296441 native pulses. The
full circuit proof is factorized through its exact16-dimensional cycle
code; a dense matrix on3204bits is not claimed. Native execution of the
compressed field/clock word has vector residual below9.35e-14. The full
physical source/preparation/readout bound is at most1.614e-5 in total
variation, including3204blank and5pointer-ready resets. Both C and E native
readout laws retain all32 angle-bin/clock outcomes. Their complete fine
effects and tensor-product probability rule are supplied. Every first
reset failure and native export complement is an output; no division by
success is performed. The rank-one conclusive reset Kraus structure gives
the conditional union bound for coherent, correlated input carriers.

## One clock and coupling family

Write `rho=2/3`, the independently derived coarse field normalization.
For N field substeps in one matter macro, let

`delta=1/N`, `alpha=delta epsilon²/rho`, `beta=delta/epsilon²`.

Then the homogeneous leading field speed per macro is
`N sqrt(rho alpha beta)=1`, and its conditional charge coupling squared is
`epsilon²`. At `N=8, epsilon²=1/3`, this is exactly the proposed Dirac field
setting `alpha=1/16, beta=3/8`. The attached native tick has angle
`-delta/16` per field substep, and relative clock phase `1/8` per macro.
Physical time is `tau0 h` per macro. Processor service duration is separate.
The emitted member fixes h=1. A different spatial scale requires the
corresponding native tick and resource certificate; no spatial refinement
is executed by changing only the register modulus in this package.

On the selected loop, put `p=epsilon z` and `q=theta/epsilon`. Its quadratic
Hamiltonian in cell-rescaled time is

`H=(a p²+q²)/2`, with `a=3/rho=9/2`, `Omega=sqrt(a)=3/sqrt2`.

This is an optical loop frequency `Omega/(tau0 h)`. It is not a
low-momentum photon branch. The weak-angle limit changes the common
coupling: `epsilon -> 0` gives `g_charge² -> 0`. It supplies no
fixed-coupling interacting charged-field limit.

## Quantum Gaussian and the exact Floquet comparison

In canonical coordinates `Q=q/a^(1/4), P=a^(1/4)p`, the quadratic split
step has classical matrix

`M_s=[[1-s²,s],[-s,1]]`, `s=delta Omega`,

and stable angle `omega=2 asin(s/2)` for `s<2`. Its pure invariant
covariance is

`Sigma=(1/(2 kappa)) [[1,s/2],[s/2,1]]`,
`kappa=sqrt(1-s²/4)`.

The positive symplectic square root `S` satisfies `S S*=2 Sigma`.
The preparation is the corresponding squeezed vacuum displaced by
`S(1,0)`. It has a genuine quantum covariance with determinant1/4.
Under the quadratic source, its covariance is fixed and its center evolves
by the actual split matrix. Its one-quantum weight in the squeezed
oscillator basis is `(1/2) exp(-1/2)`. That is a reference oscillator
weight, not a measured photon species assignment.

The global phase of the source is retained. The exact cosine factor is
`exp(i beta) exp(i beta(cos(epsilon q)-1))`; the quadratic factor is
`exp(i beta) exp(-i delta q²/2)`. The phase-sensitive numerical comparison
propagates the complex Gaussian normalization through every square root
and every `exp(i beta)` and performs no output phase alignment.
`reference_phase.py` supplies an optional actual native I-pulse removing
the known total cosine constant. A continuous-Hamiltonian comparison with
that reference includes this additional operation and its error. The
original controlled source phase is not silently discarded.

The finite split step is distinct from continuous evolution under H.
Let `r=(1/2) atanh(s/2)`. The squeeze generator K on a coherent state with
mean occupation `nbar=1/2` obeys
`||K psi|| <= sqrt(nbar²+4 nbar+2)` by the two creation/annihilation terms.
Thus `||(S_hat-I)psi||` is at most r times that constant. Spectral calculus
bounds a phase change by
`|n omega-n s| sqrt(nbar²+2 nbar+1/4)`.
The inequalities

`r <= s/[4(1-s²/4)]`,
`|omega-s| <= s³/[24(1-s²/4)^(3/2)]`

give the explicit separate substep/squeeze bound. Initial Gaussian phase
conventions are fixed before evolution and used in both references.

## Controlled cosine-to-Gaussian theorem

The integer rotor is unitarily represented on a circle. For a normalized
line Gaussian phi, define its periodization on a cell of length
`L=2pi/epsilon` by `P phi(q)=sum_j phi(q+jL)`. Its Fourier coefficients are
`sqrt(epsilon) phi_hat(epsilon z)`. The free electric factor commutes
exactly with periodization. The circle cosine multiplication is periodic.

Periodization is not assumed isometric. If the momentum variance is vp,
its squared norm differs from one by at most

`tau=2 exp(-c)/(1-exp(-3c))`, `c=vp L²/2`.

This follows by summing the Gaussian translation overlaps. The lower
normalization is `sqrt(1-tau)` whenever tau<1.

For any eta>0, Cauchy-Schwarz on the translated cells gives

`||P f||² <= C_eta integral exp(eta q²)|f(q)|² dq`,

where `C_eta <= 1+2 exp(-eta L²/4)/(1-exp(-2 eta L²))`.
The global cosine remainder satisfies

`|exp(i beta(cos(epsilon q)-1))-exp(-i delta q²/2)|
 <= delta epsilon² |q|^4/24`.

Telescoping the unitary source difference along the *quadratic* Gaussian
states therefore bounds the unnormalized output error by

`delta epsilon² sqrt(C_eta)/24
 * sum_j sqrt(E_j[exp(eta q²) q^8])`.

This uses known reference states rather than assuming that the nonlinear
states stay Gaussian. For a Gaussian with mean mu and variance v, put
`D=1-2 eta v>0`. Its weighted eighth moment is

`exp(eta mu²/D)/sqrt(D) * M8(mu/D,v/D)`,

with
`M8(u,w)=u^8+28u^6w+210u^4w²+420u²w³+105w^4`.

Dividing by the initial normalization and adding the explicitly bounded
normalization drift gives the normalized rotor comparison.
The implementation uses eta=1/16 and integer-outward enclosures throughout.

The finite-register comparison uses the separate proved Dyson tail:
if the initial support is `|z|<=K`, then
`2 sum_(j>=m/2-K) (n |beta|)^j/j!` bounds the finite/integer difference.
The exact native initial target differs from normalized infinite Gaussian
samples by its certified amplitude, phase and tail error. These three
terms compose by the triangle inequality. The compilation and register
cutoff estimates retain arbitrary references. The nonlinear estimate is
state-specific: a reference-entangled preparation uses the displayed
weighted moments of its field marginal. The executed pure Gaussian has a
product reference. Any common complete terminal instrument contracts the
resulting trace-distance bound.

At `m=2048, epsilon²=1/192, N=64`, the finite state differs from the
normalized periodized Floquet Gaussian by at most0.00926640. The separate
continuous-time/squeeze term is at most0.01720654. The direct numerical
source comparison gives0.00392769 without phase alignment. These large
register executions are ideal source calculations; the native q16
compiler/readout precision is not transferred to them.

For continuous angle bins, a further quadrature comparison is required.
For a normalized finite Fourier polynomial f, the q-point angle law and
equal-width continuous bins differ in total variation by at most
`(2pi/m) ||f'||_2`, from `||( |f|²)'||_1 <=2||f'||_2` and the cellwise
fundamental theorem of calculus. Projection tails of the reference
Gaussian and periodization tails are added before applying this bound.
No discrete atomic law is equated in total variation with an unbinned
continuous density.

## Preparation, ultraviolet modes and the successor contract

Definite electric F0 has a uniform angle law, at total variation0.59184
from the prepared Gaussian angle law. It does not satisfy the same
weak-angle preparation premise. The Gaussian cycle's exact translation
overlap gives an elementary dressed-transport defect
`||Tg-g||=sqrt(2-2 Re <g,Tg>)`. The executable reports the defect and its
electric energy; a full many-cycle Gauss dressing requires its own map.

`uv_state.py` retains all190 positive coexact modes of the actual L2 face
curvature, rather than removing optical or ultraviolet modes. Its
unwrapped quadratic Floquet Gaussian has a verified pure covariance and
Gauss support. At the proposed fixed setting, all208 plaquette variances
lie between0.305 and0.634. The one-step global Taylor norm bound is5.095
and therefore supplies no useful precision. This is a limitation of that
unrenormalized state bound; it is not a proof excluding a fixed-coupling
Maxwell limit. The95 gauge and3 harmonic coordinates are treated
separately, with no normalized uniform vector on an infinite coordinate.

A quadratic/noncompact successor must declare its real coexact canonical
state and harmonic/charge sectors, a local quadratic function of each
actual face angle, the same routed electric current and Gauss law, native
register preparation/measurement and one clock. The homogeneous spatial
and quantum-state limits must retain nonzero transverse spectral weight.
Replacing the cosine generator is an architecture amendment. None of the
cycle, weak-coupling or fixed-register receipts assigns an official row,
an interacting photon vacuum, a calibrated force, or a matter species.
