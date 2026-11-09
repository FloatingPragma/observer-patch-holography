# Capacity closure does not select the attenuation law

The retained pixel closure has a unique root in a certified interval after
its census, scale, response kernel and charged action law are specified.
That result is compatible with distinct attenuation laws having
self-consistent roots. The following perturbation preserves the complete
inner heat/scale system and instantiates this distinction with the actual
recorded bounds, without solving or fitting a replacement numerical root.

## The supplied law and its perturbation

Let alpha lie in one of the two certified intervals D in
`model/pixel_mass_closure/results.json`. The unchanged inner equations give
P=phi+sqrt(pi) alpha, u(alpha), mZ(alpha), v(alpha) and the running couplings.
Their verified endpoint signs and sign-definite implicit derivatives select
the same smooth branches on D. None of these inner equations uses the
charged attenuation chi. The outer map is g_0(alpha)=1/A_0(alpha), where
A_0 includes the full specified charged-fermion kernels.

For a supplied scalar s>=0 change only

    chi_s(u)=u^2/[512(1-u/77)]+s,
    rho_s=sqrt(2) exp(-chi_s),
    r_j(s)=1+rho_s cos(delta+2pi j/3),
    m_j(s)/v=r_j(s)^2 *
       [exp(3kappa)/(2*6^14*product_j r_j(s)^2)]^(1/3).       (1)

The phase delta and kappa retain their source laws. The determinant
normalization is unchanged. At s=0 this is precisely the owning action.
At s>0 it is a different specified interaction/controller, not a second
output of the same fixed program. No measured mass enters s. The exact
rational s=10^-14 is chosen from the certified analytic stability margin.
The original producer and retained receipts are untouched.

## Uniform perturbation bound on the complete closure

On each complete interval D, outward evaluation verifies u in (0,1),
all three r_j(0)>1/100, and A_0>101. For every s>=0,
r_j(s)=(1-exp(-s))+exp(-s)r_j(0), so r_j(s)>1/100.
Since rho_s<sqrt(2)<3/2,

    |d_s log r_j| < 150,
    |d_s log m_j| = |2d_s log r_j-(2/3)sum_k d_s log r_k| < 600. (2)

For the exact charged-lepton response kernel in the owning closure,

    K(m)=(2/pi) integral_0^1 x(1-x)
                  log[1+(mZ/m)^2 x(1-x)] dx,
    |d_log(m) K| <= 2/(3pi).

The quark and gauge terms do not change at fixed alpha. Summing the three
leptons gives |d_s A_s|<1200/pi<400. For 0<=s<=10^-4,
A_s>A_0-400s>100. Thus

    sup_(alpha in D) |g_s(alpha)-g_0(alpha)| <= s/25.        (3)

The complete response kernel is retained; no mass logarithm is substituted
for its finite-mass integral. The unchanged determinant fixes the product
of masses at a fixed alpha, so its leading summed mass logarithm cannot
select s by itself either. Equation (3) bounds the exact response.

## Two self-consistent responses from the certified margins

Let J contain g_0(D), as given by the owning centered interval image, and
put mu=min(J.lo-D.lo,D.hi-J.hi)>0. Both actual receipts have mu>9.9*10^-15.
For s=10^-14, equation (3) gives 4*10^-16<mu. Therefore g_s(D) lies
strictly inside D. Continuity and the intermediate value theorem prove
that g_s has at least one fixed point in D. This does not assert a numerical
location or uniqueness for that explicit perturbed value of s.

The original contraction constant L<0.00734 also gives, for any such root
alpha_s and the unique original root alpha_0,

    |alpha_s-alpha_0| <= [s/25]/(1-L).                     (4)

The whole-domain u enclosure gives a chi_0 interval of width below3*10^-18.
Adding s=10^-14 makes that entire chi_s interval strictly larger than the
original one. Hence every perturbed closed root has strictly smaller
attenuation exp(-chi_s) than the original root. Its positive-root Koide
ratio Q_s=(1+exp(-2chi_s))/3 is also strictly smaller; the executable
certifies disjoint outward Q intervals. Consequently the two closures have
different dimensionless mass-ratio shapes, regardless of their common
unit normalization. This is a distinction between two specified laws.

## Complete native realization throughout the perturbed root box

For e=600s<1, equation (2) gives the conservative rational mass enclosure

    m_j(s) in [m_j(0).lo*(1-e), m_j(0).hi/(1-e)].          (5)

The owning rational `pixel_mass_closure.memory.certify` is applied to these
enlarged whole-domain mass boxes. It verifies three separated cubic
projector-weight roots and two admissible quadratic controller branches,
with 0<p<1/12, 0<q<1/6 and total mass/v<1. The generic complete-record
identity therefore realizes the perturbed spectra with the original C/E
grammar, disabled records and original duty weights. The required p,q,g
are recomputed from the changed coefficient; old center controls are not
reused. No perturbed native tape or interacting spectral experiment is
claimed to have been executed.

## Parametric uniqueness and the selection boundary

There is also a qualitative parametric contraction statement. Positivity
of the masses and response denominators, the implicit derivative bounds,
and the analytic heat sums on these compact positive domains imply that
g_s and its alpha derivative vary continuously with s near zero. Choose
a sufficiently small neighborhood with
sup_D |d_alpha g_s-d_alpha g_0|<(1-L)/2 and value perturbation<mu.
Every map there is a strict self-map and a contraction, hence has its own
unique root in D. This argument supplies no explicit parameter radius for
that derivative estimate, so it is separate from the explicit-s existence
certificate above.

The full original law fixes chi uniquely; no ambiguity of that complete
specification is asserted. The claim excluded is that the pixel fixed-point
condition, shared geometry and complete-memory realizability select this
law over nearby controller laws. The conditional face-incidence producer
also treats its affine source/path weights as supplied architecture and
has multiplier countercontrols. Its fixed-point theorem proves uniqueness
after those weights are fixed, not their selection by capacity closure.

An extra equation or physical response derivation could constrain s. No
such additional selecting condition is used by the actual closure code.
The present construction does not prove that either conditional closure
describes Nature or that its scale design is physically selected.

## Owning evidence

- `model/pixel_mass_closure/{PROOF.md,closure.py,results.json}`: complete
  heat/scale branches, charged response, strict image and contraction bounds.
- `model/pixel_mass_closure/memory.py`: exact full-box controller inverse.
- `model/charged_lepton_response/p_dependence.py`: bounded outward shape
  evaluation; no replacement root solve.
- `code/particles/leptons/derive_charged_face_incidence_conditional_theorem.py`:
  supplied source weights and multiplier countercontrols, including the
  explicit bare-to-endpoint response premise.
