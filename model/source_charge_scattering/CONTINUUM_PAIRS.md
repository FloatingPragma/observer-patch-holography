# Separated charge-pair preparations on the same field source

The [bounded source creator](../spin_charge_dispersion/PROOF.md) supplies
normal operators B_+=B and B_-=B* in the joint configuration-field algebra,
with charges±1 and norms at most1/2. They use one ordered source ground Ω,
positive energy L and commuting spatial translations U(x). At a momentum
scale r>0, its complete three-outcome instrument is

    B_+, B_-, B_0=sqrt(I-B_+*B_+-B_-*B_-).

The source order bound m>1/10 follows from the replayed infrared/biased
state certificate in `spin_charge_dispersion/order.py`. Its bounded
creator theorem gives, for each sign, an original probability at least

    kappa=m0²/17576=1/1757600                              (1)

in the specified positive joint energy-momentum band. In particular the
full unfiltered creation norm squared is at least kappa. This is an
input proved on the same source, without a fitted particle wavefunction.

## 1. A nonzero pair at arbitrarily large separation

Set A=B*B. The whole equal-time configuration algebra is commutative,
including spatial translates and their adjoints. Therefore

    w(x)=||alpha_x(B) B Ω||²
        =omega(A alpha_x(A))
        =<AΩ,U(x) AΩ> >=0.                                (2)

This uses the actual field translation group. It does not use a separate
coordinate system or assume that the state factors. Source charge
covariance makes the prepared vector in (2) an exact charge+2 vector.
The conjugate preparation has charge-2 and the same norm, because source
conjugation fixes Ω and commutes with the spatial translations.

Average (2) over the cube [-R,R]^3. The joint spectral calculus of the
three commuting translation generators gives the averaging multiplier

    product_(j=1)^3 sinc(R P_j).

Its absolute value is at most1, and it tends pointwise to the projection
onto joint momentum zero. Dominated convergence on each vector proves
strong convergence of these averages to the invariant projection P_inv.
Consequently

    lim_(R->infinity) average_(x in [-R,R]^3) w(x)
        =||P_inv AΩ||² >=|<Ω,AΩ>|² >=kappa².               (3)

The middle inequality uses the invariant unit vector Ω; it does not
assume P_inv is one-dimensional. This is a mean statement, not a
clustering limit. In particular, the limit need not equal omega(A)².

The contribution of any fixed bounded set of separations to the average
tends to zero: 0<=w(x)<=||A||²<=1/16 and its volume fraction tends to
zero. Thus beyond every prescribed separation there exists a translation
with w(x)>kappa²/2. Strong translation continuity makes w continuous.
A nearby rational translation, with slack in the distance, therefore has

    w(x)>=kappa²/4.                                      (4)

Equation (4) is an existence result at arbitrarily large separations.
It asserts neither that every large separation works nor an effective
sufficient separation. The packet includes exact periodic and
nonergodic algebraic controls that prevent either unjustified inference.
They are implication controls, not substitutes for the spin source.

Apply the complete three-branch preparation, then its spatial translate.
All nine Kraus operators are

    K_(s,t)=alpha_x(B_s) B_t,     s,t in{-1,0,+1}.

Their adjoint products sum to I by completeness, even without appealing
to commutativity. Retain all nine classical records and every native
implementation failure. The (+1,+1) and(-1,-1) records have the original
probability (4); the other seven records are not removed or used to
renormalize that bound. This constructs a separated pair of source
charge excitations. Two units of internal charge are not identified
with two particles or with a Fock occupation number.

## 2. An energy bound on the entire prepared pair

The parent creator has

    B=(1/2pi) integral exp(-i theta)
        integral_0^infinity exp(-u) sin(u X_theta) du dtheta,

where X_theta=sqrt(r)[Q3(F)cos(theta)+Q2(F)sin(theta)].
For ell=L/r the inherited all-input energy-form kick estimate gives

    sqrt(q_ell(exp(±iu X_theta)psi))
        <=sqrt(q_ell(psi))+sqrt(3)|u| ||psi||.

Taking the sine difference and using the energy-form triangle inequality
preserves that upper bound. The Laplace and angular integrals are
convergent in the closed form norm: integral exp(-u)du=1 and
integral u exp(-u)du=1. Thus

    sqrt(q_ell(B psi)) <=sqrt(q_ell(psi))+sqrt(3)||psi||.   (5)

The same estimate holds for B*, all translates, and every form-domain
input. In particular q_ell(BΩ)<=3 and ||BΩ||<=1/2. Applying (5) a second
time gives the original, unnormalized pair bound

    q_L(alpha_x(B)BΩ) <=(sqrt(3)+sqrt(3)/2)² r=27r/4.     (6)

No state-specific form assumption is inserted between the two creations.
The spatial translation commutes with L, so the bound is independent of
which large separation satisfies (4).

For the upper energy cut E_max=54r/kappa², spectral Markov gives
original weight above E_max at most kappa²/8. Combining this with (4)
proves original pair weight at least kappa²/8 in

    0<=L<=54r/kappa².                                    (7)

The lower endpoint is zero. Neither (6) nor (7) excludes a zero-energy
pair component or proves a positive two-particle threshold. The nonlinear
creators need not have total momentum confined to the original single
field annulus. Every energy and momentum complement is retained.

## 3. Operational attachment and scattering boundary

The parent bounded creators and their adjoints are common strong limits
of finite bounded source histories. Their translated products and the
continuous positive square-root complement therefore belong to that same
history construction. Every finite source approximation has the complete
nine-branch instrument. The existing native finite-instrument compiler
supplies its arbitrary-input/reference error bound with all first-failure
labels. This does not supply an effective infinite-ground preparation
modulus or a numerical value of a sufficient separation in (4).

For an implemented full experiment whose unhalved trace error is delta,
the original probability bound for a stated event is reduced by at most
delta/2 before any conditioning. The source-continuum approximation must
be controlled separately from the finite native pulse error. Existence
of the limiting nonzero event does not assert that its full tape or a
laboratory experiment has been executed.

The construction supplies nonzero separated source charge-pair records
with finite energy on the actual field clock. Scattering additionally
requires identified stable one-excitation channels, an asymptotic
identification into this same Hilbert space, quantitative incoming and
outgoing wave-operator tails, and the associated complete detectors.
No such particle assignment, asymptotic factorization or scattering
matrix follows from the mean existence argument. The exact finite
response and the conditional Cook attachment are given in [PROOF.md](PROOF.md).
