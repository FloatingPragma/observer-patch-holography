# Complete local energy records on the centered CAR source

The full centered Weyl/contact source has a complete local approximation
to its charged energy readout. Its normalization retains addition, removal
and every clock record. The approximation error depends on observation
time and distance to the cut boundary, rather than the outer source volume.
A specified local ground marginal is sufficient; the cut Hamiltonian's
own ground is not substituted for that marginal.

## Source, ground and finite clock

Use the full-Fock centered source K and its stationary ground state omega
from `LOCAL_MATTER_READOUT.md`. The signed Hubbard chart transports the
field probes and action together. Fix a finite-support unit one-particle
vector f, put A_+=a*(f), A_-=a(f), and use

    tau_t(A)=exp(itK) A exp(-itK).

Both the finite full ground projector and the positive-energy ground GNS
representation are admitted. In the latter, tau_t(A)Omega=exp(itL)A Omega,
L>=0. Define the original positive measures

    mu_q(I)=||1_I(L) A_q Omega||^2,   q in {+,-}.

Their total masses add to one by CAR. A general spatial smearing need not
have separate mass 1/2 in each charge direction. Let n>=2 and Delta>0,
t_l=l Delta for l=0,...,n-1, T=(n-1)Delta. The following finite instrument
uses the same source clock and requires no numerical ground energy:

    K_(q,j)=(1/n) sum_l exp(2 pi i j l/n) tau_(-t_l)(A_q),
       q in {+,-}, j=0,...,n-1.                           (1)

Discrete Fourier orthogonality gives the operator identity

    sum_(q,j) K_(q,j)* K_(q,j)
       =(1/n) sum_l tau_(-t_l)(A_+*A_+ + A_-*A_-)=1.     (2)

It holds on every source input and after adjoining an arbitrary reference.
No success branch is normalized separately. On the ground, the record
probability is

    p_(q,j)= integral F_n(Delta E-2pi j/n) dmu_q(E),
    F_n(theta)=|sum_(l=0)^(n-1) exp(-il theta)/n|^2.       (3)

The energy phase arises from the ground identity in the same dynamics.
The finite phase labels are periodic energy responses, not exact energy
projectors or spatial momentum labels. On an excited input, (1) measures
energy differences through its operator history, rather than the energy
of the post-addition state relative to a chosen ground.

One dilation prepares a uniform time register, applies controlled inverse
source evolution, the complete two-charge transfer isometry, controlled
forward evolution, and the Fourier transform of the time register. The
finite cut described below makes every operation local to one finite CAR
algebra. This uses the declared coherent control and fermionic transfer
services. Their full failure records and approximation errors are retained.
The proof supplies no native execution of a large phase tape and no
preparation of an unspecified infinite ground.

## A complete local approximation with no phase-count penalty

Suppose f is supported on s sites in the radius-r cube. Cut only the source
bonds leaving a radius-R cube, R>r, and replace tau by tau^R in (1).
Every cut clock is an exact automorphism, so the cut family K^R also obeys
(2) exactly. The source locality theorem gives, for either q and all |t|<=T,

    ||tau_t(A_q)-tau_t^R(A_q)|| <= B_R(T),
    B_R(T)=min(2,6cs(2R+1)^2 T 2^(-(R-r))
                                      (exp(48cT)-1)).    (4)

The onsite U drops out in the interaction-picture proof. This estimate
applies to odd CAR fields because the omitted source bonds are even.
Let W,W_R be the Stinespring isometries with all 2n labels. Fourier
orthogonality before taking the norm gives, for every vector psi,

    ||(W-W_R)psi||^2
      =(1/n) sum_(q,l) ||(tau_(-t_l)(A_q)
                         -tau^R_(-t_l)(A_q))psi||^2
      <=2 B_R(T)^2 ||psi||^2.                            (5)

Consequently the complete output half-trace distance is at most
min(1,sqrt(2) B_R(T)), uniformly over arbitrary untouched references.
There is no factor n from summing individual record estimates. The same
bound applies after any classical grouping of the records.

If two unit smearings f,g are compared, the exact CAR norm gives an
additional isometry error sqrt(2)||f-g||. For any actual symmetry alpha_x
commuting with the source time evolution, equation (1) also gives
K_(q,j)(alpha_x f)=alpha_x(K_(q,j)(f)). This is covariance of the operational
packet family. It does not attach open-cube Walsh coefficients to a
translation generator or assume a relativistic dispersion relation.

A cut instrument acts only on the radius-R algebra. Replacing the true
input marginal on that algebra and its accessible reference by a supplied
marginal at half-trace distance delta_R costs at most delta_R. With a
complete implementation error epsilon_service, the total bound is

    eta <= min(1,sqrt(2)B_R(T)+delta_R+epsilon_service).   (6)

A finite source whose ground is used to approximate the infinite state
must provide delta_R. Spatial locality alone does not bound that ground
marginal error. The selected centered ground subsequence gives convergence
for fixed R, without an effective subsequence modulus. For the global
uniform-attempt repair protocol, `LOCAL_MATTER_READOUT.md` applies when the
finite cut readout is an accessible local intervention; its specified
local-record error is added to (6). Distant raw repair values remain in
the source archive and are not compared by this local theorem.

## From measured records to true energy windows

The phase period is P=2pi/Delta. At wrapped energy distance at least d
from a source energy E, the sum of the corresponding phase-record kernel
obeys

    tail(n,Delta,d) <= min(1,1/[2(a-1)]),
    a=n Delta d/(2pi),                                  (7)

when a>1; otherwise use one. The grid-tail proof is the one in
`CUBE_SPECTRUM.md`: the squared Dirichlet kernel is bounded by
1/(4n^2 u^2), and both grid tails sum to at most 1/[2(a-1)].

Let J be a selected group of phase centres and W a real energy set.
Suppose every E in [0,Ecut] outside W has wrapped distance at least d
from every centre in J. Write p_q(J) for the original grouped probability
and M_q for mu_q([0,infinity)). Any proven original ultraviolet tail
bound v_q(Ecut)>=mu_q((Ecut,infinity)) then gives

    mu_q(W) >= max(0,p_q(J)-tail*M_q-v_q(Ecut)).           (8)

A local measured probability interval [p_low,p_high] at complete-output
error eta can replace p_q(J) by p_low-eta. The original second moment
already gives v_q(Ecut)<=M_2,q/Ecut^2; higher source moment bounds can
improve it. The tail is necessary because a local ground response has no
volume-independent finite spectral support. Aliased high energies cannot
be discarded from the observed record by declaration.

Conversely, if every E in a known interval I has wrapped distance at
least d from every excluded centre, then

    p_q(J) >= (1-tail) mu_q(I).                           (9)

This direction needs no ultraviolet truncation. Both inequalities keep
original probabilities, including the opposite charge and all unselected
phase records.

## Quantitative recovery of the original energy measure

A bounded real Lipschitz test g can be read by assigning scores to the
complete phase records. Choose signed phase centres in [-P/2,P/2], then
replace negative centres by zero. Suppose Ecut+d<=P/2. For E in [0,Ecut],
any wrapped phase error smaller than d is then an ordinary signed error
smaller than d; clamping to zero cannot increase its distance from E.
Write G=||g||_infinity and ell for its Lipschitz constant. Equations (3)
and (7) give the original-measure comparison

    |sum_j g(max(0,e_j_signed)) p_(q,j) - integral g dmu_q|
       <= ell*d*M_q + 2G*(tail*M_q+v_q(Ecut)).            (13)

A complete output error eta adds 2G eta. No distribution is conditioned
on the chosen charge or energy record. The zero score on the opposite
charge retains the same probability convention.

For fixed source parameters and packet, choose Ecut tending to infinity,
d tending to zero, Delta small enough for the half-period condition, and
n Delta d tending to infinity. A volume-independent moment bound makes
the ultraviolet term vanish. Equations (4)--(6) then prescribe a finite
cut for each observation. An actual ground-state approximation needs its
local marginal error at that cut. Along the selected ground subsequence
those errors vanish in the corresponding ordered choice. Thus bounded
Lipschitz energy tests are recovered from complete local records without
an isolated spectral pole or a global energy projector. This construction
has no effective selected-ground modulus for a simultaneous continuum
refinement.

## A positive local detector event on growing default sources

For a degree-six interior site at c=1/8,U=1, the exact source moment bound
in `CHARGED_SPECTRAL_BANDS.md` supplies, for each charge separately,

    mu_q({35/256<=E^2<=2744/5}) >= p_*=25/100352.          (10)

The band lies in [1/4,47/2]. Choose n=4096, Delta=1/16 and the group
J={0,...,1024}. Its upper phase centre is 8pi>24, its period is 32pi>96,
and the band is at least 1/4 from every excluded centre on the circle.
Using pi<22/7 in (7) gives tail<=11/202. Therefore

    p_q(J) >= (191/202)p_*.                              (11)

Choose the cut radius by the exact rational search in `local_spectral.py`
so that (3/2)B_R(T)<=p_*/4. Since sqrt(2)<3/2, the cut instrument applied
to the same full ground has grouped probability at least

    [(191/202)-(1/4)]p_* > p_*/2.                        (12)

This is a finite-support, finite-time complete detector event with a
strictly positive original probability independent of the outer source
volume. Every outer source contains the cut cube, and the observed site
has all six source neighbours. Equation (12) is not a claim about the
cut cube's own ground. A supplied marginal or service error reduces the
lower bound by its value in (6). Both charge records remain in the same
complete instrument.

The ratio lambda=U/c is eight in this example. Homogeneity gives
K(c,lambda c)=c K(1,lambda); ordinary source times are unit-hopping times
divided by c, and all preparation weights are unchanged. At fixed lambda,
a fixed physical energy window under c=v/a corresponds to a shrinking
unit-hopping window. Equations (1)--(9) measure that window with explicit
resolution, locality and original-probability errors. They supply no
positive lower weight in it. That requires actual ground response
information on the scaling packets. The uniform band (10) lies at fixed
unit-hopping energies and is not a finite physical mass under that limit.

## A tapered complete clock with a positive inverse-window bound

Clock amplitudes can reduce alias leakage while preserving completeness.
Let a_l be the coefficients of (1+z+...+z^(m-1))^4, padded by zeros to
N=4m terms, and S=sum_l a_l^2. Their sum is m^4, so Cauchy--Schwarz gives
S>=m^7/4. Replace (1) by

    K_(q,j)=(NS)^(-1/2) sum_l a_l exp(2pi i j l/N)
                                      tau_(-l Delta)(A_q).             (14)

Parseval gives sum K*K=sum_l(a_l^2/S)tau_(-l Delta)(1)=1. In (5), the
uniform time average is replaced by the probability weights a_l^2/S;
the same sqrt(2)B_R error applies to the complete instrument.

The scalar kernel is |D_m(theta)|^8/(NS), where D_m is the length-m
Dirichlet sum. For wrapped distance d and a=N Delta d/(2pi)>=1,
|D_m(2pi u)|<=1/(2|u|) and the two grid tails give

    tail_4 <= N^7/(128 S) [a^(-8)+a^(-7)/7]
            <= (4096/7) a^(-7).                         (15)

The amplitude normalization S is computed exactly by the finite producer.
No phase label or charge direction is omitted.

At the default source, `UV_TAIL.md` proves the original total charge tail
above Ecut=20 is at most

    v=1647837513/156250000000000 < p_*/16.

The band (10), truncated at20, consequently has each-charge original
weight at least p_*-v. Take N=8192, m=2048, Delta=1/16, d=1/4 and
J={0,...,1664}. The last centre is (13/2)pi, separated from20 by more
than1/4. The lower excluded wrapped centres are also separated by1/4
from the known band. Equation (15) supplies

    t_4=(4096/7)(11/224)^7 < p_*/16.

Choose the cut so its complete error delta is at most p_*/16. Then the
local measured event has original probability at least

    q_low=(1-t_4)(p_*-v)-delta.

Applying (8) to that measured record, retaining its complete error and
the same ultraviolet budget, gives the true original spectral weight

    mu_q([0,(13/2)pi+1/4])
       >=q_low-delta-t_4/2-v >3p_*/4.                   (16)

The second delta is intentional: one comparison establishes the
measured lower bound; the inverse inference uses the measured value
with its own allowed ideal-versus-cut error. The finite producer retains
both charges and all16384 outcomes, and gives the exact rational radius
and lower bounds. Marginal and service errors require their additional
budgets from (6). This is a growing-source energy-window readout on the
specified ground, without conditional renormalization or a mass-shell
assignment.

## Exact source control

The finite control rebuilds the default full sixteen-state centered dimer,
its unique ground, every energy projector, and both local charge directions.
The source response has energies sqrt(17)/8 plus or minus 1/16, with
original weights 1/4 minus or plus 1/(4sqrt(17)). Exact source projector
contractions establish these values. Both the uniform and fourth-convolution thirty-two-step clocks with Delta=1
retain all sixty-four charge/phase records. Rational square-root and
trigonometric enclosures compute their probabilities. The operator
completeness proof (2) applies to the same actual source clock; separate
matrix tests compare (1) directly with (3) and retain excited inputs as a
control against a false ground-relative spectral interpretation.

The dimer control does not compute the infinite centered ground. No
physical particle, mass, scattering process or Lorentz symmetry is
assigned by these records.
