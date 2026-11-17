# Complete continuum boost-orbit apparatus and locality bound

## Construction, composition, and entropy

Let the continuous positive source and its unitary multiplier V be exactly those in `../spectral_flattening/`. For the complete ordinary Dirac Poincare representation U_D define U_F(g)=V U_D(g) V† on the full one-particle space. Fix the original finite-entropy regional projection chi_B and define

R_(g,B) = U_F(g) chi_B U_F(g)†.

This is an explicit linear CAR apparatus map. It has a second-quantized implementation in the source particle/hole Fock representation because U_F preserves its vacuum projector P. It is not the region V chi_B V† of the transported Dirac net. In particular, at g=identity the apparatus is chi_B itself.

Composition U_F(g)U_F(h)=U_F(gh), each R being an orthogonal projection, nesting for one fixed observer, and vacuum entropy invariance are exact. The compressed covariance on R_(g,B) is unitarily equivalent to chi_B P chi_B. Thus every such boosted proper ball has the same finite vacuum entropy and the same area coefficient if that asymptotic coefficient exists for the rest ball. No observer-dependent entropy subtraction is made. This transports a real region algebra with a real source unitary; it does not assert that a simple local sampling of the field on a tilted plane implements that apparatus. The latter implementation was tested separately and failed.

This is an observer architecture, not yet a native OPH record-generation law. The word "proper" here uses the declared background Minkowski geometry and actual source Poincare action; its geometric accuracy must be tested by localization, not merely assigned to the algebra.

## A full three-dimensional continuum bound

The same V is smooth at p=0: its scalar part is a smooth function of R^6 |p|^6, and its nonscalar part is beta(alpha dot p) R^3 |p|^2 times another such smooth function. At large |p| its angular limit is bounded and its derivative of order k is O(|p|^-k). Consequently for any unit vector n,

M6 = integral_R3 ||(n dot grad_p)^6 V_(R=1)(p)||_op d^3p < infinity.

Rotational covariance makes M6 independent of n. Both the regular origin and the integrable |p|^-6 ultraviolet tail matter. This is a defined finite constant of the source, not a freely adjustable geometric or gravitational coefficient.

For x != 0, integration by parts six times in direction n=x/|x| gives the distributional convolution kernel v_R(x), apart from any distribution supported at x=0, the bound

||v_R(x)||_op <= M6 R^3 / [(2 pi)^3 |x|^6].

The identity follows equally by the distributional derivative rule; it does not assume that V itself is integrable. All supported-at-zero terms vanish between separated spatial sets. The Schur/Young bound then gives, for sets A,B with Euclidean separation at least r,

||chi_A V_R chi_B|| <= integral_(|x|>=r) ||v_R(x)|| d^3x
                         <= [M6/(6 pi^2)] (R/r)^3 = delta_3(r).

The same holds for V†. This is a bound on the full operator including arbitrarily high momentum modes. It is not an average over low modes and does not invoke a finite Fourier band limit. The constant M6 is proved finite here; its numerical value has not been evaluated, so this theorem alone does not provide a numerical three-dimensional tolerance at a specified distance.

The ordinary complete Dirac boost has finite propagation along the characteristics of the Dirac equation. A ball of radius ell_B on its boosted rest plane reaches at most exp(|z|) ell_B in either null coordinate, where z is rapidity. For two apparatus balls whose geometric spacelike separation margin is

g = d - ell_A - exp(|z|) ell_B > 0,

expand the unboosted A by r_A=g/2 and B by r_B=g/[2 exp(|z|)]. Their ordinary Dirac apparatuses still do not overlap (one can take strictly smaller radii and pass to the limit). Insert these enlarged projections on either side of the Dirac boost in

chi_A V T(d) U_D(z) V† chi_B.

The term with both enlarged projections is zero by ordinary Dirac causal support. The remaining two terms are bounded by the quasi-local tails of V and V†. Translations commute with V. Therefore

||chi_A T(d) U_F(z) chi_B|| <= delta_3(g/2) + delta_3(g/[2 exp(|z|)]).

This gives controlled approximate cross-observer CAR locality as R/g tends to zero, with complete modes retained. It does not give strict microscopic locality, nor a unique observer-independent local algebra for every spacetime region. The bound scales with the geometric gap, not merely with the ball radius; making a ball large while holding two boundaries a distance R apart does not make this error small.

## An explicit one-dimensional constant

In the transverse-zero sector V(p)=exp[A theta(p)], A=beta alpha_x, A^2=-I, and theta(p)=arctan(p^3)/2 for R=1. Thus

V'''' = V { I[(theta')^4-3(theta'')^2-4 theta' theta''']
              + A[theta''''-6(theta')^2 theta''] }.

Using the absolute rational derivatives of theta gives the rigorous envelope

||V''''(p)|| <= [45 p^2 + 544.5 p^5 + 563.0625 p^8
                   + 1341 p^11 + 288 p^14 + 180 p^17]/(1+p^6)^4,

for p>=0. Twice its integral is exactly

M4_bound = 155 + (5739/256) pi.

The beta integral integral_0^infinity p^a/(1+p^6)^4 dp = B((a+1)/6,4-(a+1)/6)/6 evaluates each term. The convolution tail obeys

delta_1(r) <= M4_bound/(3 pi) (R/r)^3.

This analytic full-operator bound needs no numerical ultraviolet cutoff. The cross-observer estimate uses the same two enlargement radii as above. `BOUND.json` evaluates it at declared gaps/rapidities; these are upper bounds, not observed overlaps or fitted error models.

## Finite implementation and its limitation

`run.py` retains every mode of a finite Fourier CAR system and uses an explicitly self-adjoint, vacuum-preserving finite closure of the Dirac boost. The closure pinches the scalar {x,p}/2 into positive/negative momentum blocks and is conjugated by the exact V. Its algebra, group composition, complete frame, and entropy checks are exact finite-source tests. Held-out Gaussian wave packets test whether it approximates the physical continuum boost on their specified domain.

The full finite-cutoff spacelike operator norm also includes modes hitting the momentum boundary during a boost. The ordinary Dirac control has sizable nonzero overlaps from that closure, even though its continuum geometric algebra is exactly local. These finite norms therefore cannot be transferred as continuum locality measurements or quietly subtracted. They remain recorded failures of that approximation for a complete-mode locality test. The continuum bound above avoids this particular cutoff closure.

## Remaining gravity obligations

This branch repairs an observer-orbit entropy component and proves a qualitative full three-dimensional macroscopic localization bound. It does not prove exact geometric local nets, curved-spacetime observer transport, actual fixed-volume stationary entropy for matter plus geometry, a dynamical lapse equation, native apparatus/event implementation, interacting chiral or charge compatibility, or the measured gravitational coefficient in the same matter clock. A subsequent joint preparation law must use actual CAR entropy and test unused regions and observers. Neither Poisson nor n=-sigma has been installed here.
