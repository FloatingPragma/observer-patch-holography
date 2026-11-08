# Prepared local access to original Lorentzian spectral weight

Use the complete centered CAR source, its chosen invariant ground and the
same source clock as `GROUND_MARGINALS.md`, `LOCAL_SPECTRAL_READOUT.md`
and `UV_TAIL.md`. Both original charge measures are retained. This
construction measures a Lorentzian energy score without a periodic
energy identification. It also replaces full-region state tomography by
a single pressure diagnostic when only that score is required.

## Complete interference records

For a unit finite-support one-particle packet g write A_+=a*(g), A_-=a(g)
and tau_t(A)=exp(itK) A exp(-itK). At each supplied source time t the four
Kraus maps are

    J_(q,b)(t) = [A_q + b tau_(-t)(A_q)]/2,
                     q in {+,-}, b in {+1,-1}.                 (1)

Summing both signs cancels their cross terms; summing both charges then
uses CAR to give sum J*J=I on every input. In particular, this is a
complete instrument with arbitrary untouched references. The signed
score b has expectation

    C_g(t)=sum_q Re omega(A_q* tau_(-t)(A_q))
          =integral cos(tE) dnu_g(E),                          (2)

where nu_g=mu_+ + mu_- is the positive ground-relative energy measure of
total mass one. Equality to the spectral expression uses the actual
ground state. On an excited input the same complete instrument measures
energy differences. No ground energy is supplied to its controller.

For eta>0 the desired original score is

    Q_eta(g)=integral eta^2/(eta^2+E^2) dnu_g(E)
            =eta integral_0^infinity exp(-eta t) C_g(t) dt.    (3)

It is also -eta times the imaginary part of the original signed retarded
resolvent at i eta. Addition/removal energies have opposite signs in that
representation, which does not alter this even function.

Use rational eta,Delta>0 and integer m>=1. Put r=1/(1+eta Delta),
w_j=(1-r)r^j for j=0,...,m-1. Choose time j Delta with that probability,
apply (1), and record time, charge and sign. The remaining probability
r^m selects a dummy branch with the two CAR maps A_q and score zero.
There are 4m+2 records; their sum is I. All source implementation failures
are additional retained zero-score records. A branch is never normalized
by its own success probability.

The infinite geometric sum has the scalar kernel

    Re eta/[eta+(1-exp(-i Delta E))/Delta].

Its difference from eta^2/(eta^2+E^2) is at most Delta |E|/2:
|1-exp(-ix)-ix|<=x^2/2 and the first denominator has real part at least
eta. If integral E^2 dnu_g<=v^2, the finite score s_m therefore satisfies

    |Q_eta(g)-s_m| <= Delta v/2 + r^m =: e_time.         (4)

This is a global source moment bound; it needs no ultraviolet alias cut.
The bound in `UV_TAIL.md` supplies v for the specified invariant packet
family. The exact local CAR identity supplies
v^2=(U^2+d c^2)/4 for a single mode at a degree-d site, even without
particle-hole symmetry. Its second moment is summed over both charges.
Separate charge weights need not equal one half for a general packet.

## One complete error budget

Let g have support on s sites in a radius-r0 cube, and cut the same source
bonds at radius R>r0. Write T=(m-1)Delta. The existing source locality bound
is

    B_R(T)=min(2,6cs(2R+1)^2 T 2^(-(R-r0))
                                             (exp(48cT)-1)).

For (1), only its time-evolved term changes. Summing the four records
before taking the norm gives

    sum_(q,b) (J-J_R)* (J-J_R)
       = (1/2) sum_q (tau_(-t)(A_q)-tau^R_(-t)(A_q))*
                         (tau_(-t)(A_q)-tau^R_(-t)(A_q))
       <= B_R(T)^2 I.                                  (5)

The probability-weighted time register leaves this bound unchanged; the
dummy maps are identical. Thus the complete output half-trace error is
at most B_R(T), with no number-of-times penalty. A signed score in [-1,1]
costs twice that error.

Suppose the supplied cut-region state, together with any compared
reference, is within half-trace distance delta_R of the true ground
marginal. Add complete preparation, source-service and accessible local
repair errors epsilon_p, epsilon_s, epsilon_r. A measured signed score
interval [l,u] then encloses Q_eta(f), for a unit target packet f with
||f-g||<=epsilon_f, after enlarging it by

    e = 2 epsilon_f + e_time
            + 2 min(1,B_R(T)+delta_R+epsilon_p+epsilon_s+epsilon_r).
                                                               (6)

The first term follows from the original retarded CAR isometry and the
positive contraction eta^2/(eta^2+L^2), rather than an identification of
finite cube coefficients with momentum. A genuine magnetic-momentum
packet may be truncated in spatial l2 norm and normalized with its
explicit epsilon_f. Both the Bloch chart and any unitary boundary
transport belong to that packet specification. A cube Walsh vector is
not promoted to a translation eigenvector.

A finite algebraic clock can make the local score observable itself
algebraic. Shift the cut Hamiltonian by a scalar to H with
||H||<=h=c times the number of retained bonds + U times the site count/4.
The exactly unitary Cayley step

    V=(I-i Delta H/2)(I+i Delta H/2)^(-1)

has ||V^j-exp(-ij Delta H)||<=j Delta^3 h^3/12. Splitting each time step
into n pieces divides this bound by n^2. The Heisenberg operator error
has an additional factor two. Hence a complete instrument service
allowance 2 T Delta^2 h^3/(12n^2) suffices. The same finite coherent
compiler can approximate these algebraic controls, adding its complete
error. This is a constructive program, with no practical cost guarantee
for a large collar and no claim that its native tape has been executed.

## A single-observable ground certificate

Choose the normalized finite packet g with algebraic coefficients and
a finite algebraic clock. Let S_R be the signed-score observable of
the complete local instrument. It is an even algebraic Hermitian matrix
with ||S_R||<=1. The dummy records have zero score. No expectation of an
isolated cut ground is inserted in (6).

Instead add h sum_x theta_(b x)(S_R) as a variational diagnostic only,
where b=1 in an ordinary translation chart or b=2 for the commuting
magnetic translations. The source preparation and clock are unchanged.
For a cube of side L divisible by b and support coordinate spans s_i<L,
let n_i=floor((L-1-s_i)/b)+1 and

    nu=1-[n_1 n_2 n_3/(L/b)^3].

The missing diagnostic density is nu/b^3. The finite minimum enclosures
[ell_0,u_0], [ell_+,u_+], [ell_-,u_-] give, for every ground invariant
under those translations,

    b^3(ell_+-u_0)/(h L^3)-3b^3 c/(hL)-nu
       <= omega(S_R)
       <= b^3(u_0-ell_-)/(h L^3)+3b^3 c/(hL)+nu.        (7)

Intersect with [-1,1]. This is the tiling and pressure-secant argument
of `GROUND_MARGINALS.md`, with the magnetic-cell density retained.
Algebraic finite matrices permit certified minimum-eigenvalue intervals.
Only one diagnostic observable is required, instead of 16^s coefficients
for a marginal on s sites. Equations (4)--(5) and the algebraic-clock
error then transfer (7) directly to Q_eta. This route certifies an ideal
source expectation. An independently measured finite prepared state
uses (6); equation (7) does not prepare the infinite ground or compare
its inaccessible exterior/reference state.

There is a rigorous stopping interpretation. The right derivative of the
concave pressure per site at h=0, multiplied by b^3, is the infimum of
omega(S_R) over all translation-invariant source grounds. To see the reverse inequality to
the variational upper bound, take invariant grounds of the perturbed
interaction along h down to zero. Compactness of local states gives a
limit; the local ground inequality passes to the unperturbed interaction.
The variational bound is e(h)>=e(0)+h omega_h(S_R)/b^3. Every limit
therefore has a score at
least the infimum over unperturbed grounds. Together with the upper
variational inequality this proves equality. Tiling enclosures at
rational h and sufficiently large L converge to each secant.

Consequently enumerating rational h down to zero, increasing L and
certified eigenvalue precision halts with a lower certificate above z
whenever every such ground has omega(S_R)>z with a uniform strict
margin. If distinct grounds straddle z, no such universal certificate is
asserted. No uniqueness, uniform spectral gap or effective selected-state
subsequence is assumed. The search is computationally finite on a
positive instance but can have prohibitive cost.

## Original spectral windows and a terminating fixed-source example

For W=N eta with N>0, t_N=1/(1+N^2), positivity and total mass one give

    nu_f([0,W]) >= max(0,(Q_eta(f)-t_N)/(1-t_N)).         (8)

Every score interval from (6) or (7) can be used before this conversion.
The result concerns the total original charge measure, including any
zero-energy component. An isolated pole, positive lower energy threshold
and equal separate charge weights do not follow from (8).

For the actual default degree-six local probe c=1/8,U=1, the exact
second moment is 35/128 on every source ground. Jensen gives

    Q_1 >= 1/(1+35/128)=128/163.

The same second moment directly gives the stronger analytic bound
nu([0,1])>=1-35/128=93/128 by Markov. The construction below certifies
a lower weight through a complete operational score and a terminating
state diagnostic; it does not improve that analytic moment inequality.

Take Delta=1/256, m=2048 and v=3/4. The exact rational time error in (4)
is less than 1/256. Choose the collar and Cayley subdivision so that each
complete half-trace error is at most 1/512. Then the total ideal-score
error is less than 3/256<1/64, and every invariant ground has
omega(S_R)>3/4. Equation (7)'s search therefore terminates with a strict
lower bound greater than 3/4. Subtracting 1/64 and using N=1 proves
original total weight greater than 15/32 in source energies [0,1].
This example proves termination of the operational certificate for a
source-defined broad window; the
finite search and its large native readout tape are not executed here.

The finite calculation executes the ideal complete algebraic instrument
on all sixteen rows of the actual default dimer and checks all time,
charge, sign and dummy probabilities. It compares its score with the
exact full-source Lorentzian measure, retaining every zero and nonzero
record. That finite witness does not represent a thermodynamic ground.

At fixed U/c=8 with c=v_physical/a, fixed physical eta corresponds to
unit-hopping eta proportional to a. A fixed physical momentum packet
also spreads over an increasing spatial region. Equations (4)--(8) give
an explicit experiment and a semidecision for a strict original-weight
floor. Their uniform Jensen floor is O(a^2). They supply no positive
scaling-uniform floor, no effective selected-ground marginal modulus,
and no interacting relativistic dispersion. Such a floor requires the
actual low-energy correlation or memory estimate in that same source.
