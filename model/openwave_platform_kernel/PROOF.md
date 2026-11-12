# Complete native evolution in a headless OpenWave kernel

The declared source is the periodic kinetic-plus-analyzer-pairing member of
model/openwave_kinetic_phase_cycle. The actual native charged Higgs tensor
and analyzer intervals generate its complete 180-row Nambu symbol. Its
ninety physical CAR modes per site, conjugate block, three charged families
and 39 massless internal coefficients are retained. The export expresses the
family coordinates in the native analyzer character frame. This is a
unitary change of coordinates, not a selection of three Dirac matrices.

## Native action and portable stencil

Write the source symbol as H(p)=M+sum_j Gamma_j sin(p_j).
The native construction gives Hermitian M and Gamma, the Clifford identities
{Gamma_i,Gamma_j}=2delta_ij I, and {Gamma_j,M}=0. The exact analyzer
intervals give ||M||<1. Therefore ||H||<2 on every periodic finite cube
at spacing one. The site action uses the central difference
Gamma_j(psi_(x+e_j)-psi_(x-e_j))/(2i), including the opposite-momentum
hole block. The exported sparse matrices contain every nonzero entry.

The Taichi CPU kernel applies this action to all 180 coordinates at every
site. An independent dense NumPy stencil checks its derivative on a fixed
full-carrier complex probe. Reversing one spatial kinetic matrix is a
retained mutation control. This detects an action/sign mismatch rather than
only testing a selected charged branch.

## Integrator bound and numerical records

For Hermitian H let R4(-i dt H) be the fourth-order Taylor polynomial.
For real x,

    |R4(-ix)|^2=1-x^6/72+x^8/576.

Thus the polynomial is a contraction when |x|<=sqrt(8). The integral
Taylor remainder for a unitary exponential gives one-step error at most
dt^5 ||H||^5/120. Telescoping n contractive steps gives

    ||R4(-itH/n)^n-exp(-itH)|| <= t^5 2^5/(120n^4)

in the tested 2t/n<=2 domain. This is an exact-arithmetic integration
bound. The floating residuals and norm drifts are separately reported;
the proof does not absorb hardware roundoff into the analytic enclosure.

The source-defined four-column charged intertwiners are used only to
prepare a localized positive-band addition and independently evaluate its
Fourier evolution. The kernel evolves the complete source. The addition
has original success probability 1/2; its full preparation failure is
retained. The detector contains every site and every native Nambu
coordinate on this one-extra-particle preparation. Its complete effects
sum to identity on the prepared positive band, with the full source
complement assigned its own zero-probability outcome. No physical local
electron detector is inferred from this filtered CAR instrument.

The Fourier reference supplies the actual complete source probabilities.
The numerical solution is normalized back to its success norm only for
comparison. Its unnormalized norm drift is retained. If the operator error
is epsilon, normalization changes the conditional vector by at most
2epsilon. Complete-record total variation is consequently at most
2epsilon, including the unchanged failure branch. Numerical normalization
is not an additional physical channel or a discarded outcome.

## Platform execution and physical domain

The portable engine imports no RER module. It runs under
openwave.xperiments.oph_observer.engine in the pinned OpenWave checkout,
using its installed Taichi dependency. Its mathematical matrix packet is
regenerated from the public RER source and compared by hash. The overlay
does not modify any upstream tracked file or award a model column.

Sides four and six retain11520 and 38880 Nambu coordinates. All three
charged preparations are executed with ten and twenty steps to event
parameter one. Exact source and approximate records, failures, numerical
drifts, full-carrier derivative checks and the sign mutation are retained.
The controlled step-size change tests fourth-order convergence.

The eight sine valleys, massless spectator coefficients, retrospective
phase 2/9 choice and empirical electron scale retain their original domains.
This port supplies no removal of doubling, interacting mass correction,
physical elapsed time, continuum theory or upstream scientific award.
