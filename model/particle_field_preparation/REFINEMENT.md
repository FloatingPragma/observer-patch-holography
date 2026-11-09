# Conserved-energy gap control and its refinement domain

The [finite polynomial preparation](PROOF.md) admits a stronger gap bound
than the short-time first-energy-moment estimate when the declared
plaquette coefficient is sufficiently large. This is a different fixed
member of the same coordinate-field Hamiltonian family. Its input
coefficient and preparation cost are explicit, with no particle-mass
calibration. A source with plaquette coefficient 1/100 is not silently
replaced by this member.

## The existing admissibility theorem on this source

The [retained overlap implementation](../common_source/overlap_regulator.py)
uses the bound of
[Hernández, Jansen and Lüscher, section 2.3 and Appendix C](https://arxiv.org/html/hep-lat/9808010v2):
for Wilson parameter s=0, unitary represented links and every represented
plaquette within epsilon of identity in operator norm,
`A†A>=1-30 epsilon`, where `A=I-aD_W`. Their finite lattice is periodic.

The source's three spatial directions are the constant fourth-direction
sector of a periodic four-dimensional Wilson lattice with fourth link
identity. The fourth covariant derivative and all temporal plaquette
defects vanish. The retained spatial beta and the Euclidean gamma5 are
different matrices, but both squared Hermitian kernels on this sector
equal `A†A`. The theorem therefore applies to the actual spatial Wilson
singular values. Its proof uses unitarity and represented plaquette norms,
so the unitary reducible M15 representation of `S(U(3)×U(2))` is admitted.
It does not require an SU(3)-fundamental color space.

Taking epsilon=1/40 gives the required squared gap 1/4. Define the actual
nonnegative magnetic potential

```
S(U)=sum_plaquettes [1-Re Tr_M15 U_p/15]
    =sum_plaquettes ||R(U_p)-I||_F²/30.
```

There are P=3V=E positive plaquettes on the periodic cubic lattice. If
`sigma_min(H_W)<1/2`, at least one represented plaquette has operator-norm
defect greater than 1/40. Thus, as multiplication operators,

```
1_Bad <= 48000 S.                                      (1)
```

This is a sufficient control criterion; a small Wilson gap does not by
itself imply loss of overlap locality. The same cited paper's section 2.4
gives a separate isolated-low-mode route under its projector-locality and
spectral-separation conditions. Those conditions are not supplied by the
present prepared state.

## A plaquette moment from the same finite word amplitude

For a threshold delta>0, let Large be the event that some represented
plaquette has operator-norm defect greater than delta. For every gauge
frame of a configuration in Large, one of its four link factors has
operator-norm distance greater than delta/4 from identity. This follows
by telescoping the four unitary factors, including repeated edges in a
one-cell torus. The normalized polynomial preparation in (1) of
[PROOF](PROOF.md) therefore obeys

```
Prob(Large) <= (60/59)^(2E) (1+8n)^(900E)
               exp[-n² delta²/480].                    (2)
```

On its complement, each plaquette contribution to S is at most delta²/2,
because `||R(U_p)-I||_F²<=15||R(U_p)-I||²`. On Large it is at most two.
Consequently

```
<S>_0 <= E delta²/2+2E Prob(Large).                      (3)
```

The integer planner chooses delta=2^(-r), a dyadic upper bound on
Prob(Large), and then L=n² using

```
Prob(Large) <= min(1,2^b_delta),
b_delta=2E+900E ceil(log2(1+8n))-floor(n² delta²/480).
```

For target eta, it makes each term in (3) at most `eta/(4*48000)`.
This gives `S0<=eta/(2*48000)`. The amplitude and deterministic finite
Householder preparation are the same construction as in PROOF; no
additional projection or repeated postselection is inserted.

## All-time bound from the actual conserved Hamiltonian

Keep the initial N=60 fermion preparation and no extra scalar interaction.
On this conserved occupation sector write

```
H_m = H_E+K_m-kappa sum_plaquettes ReTr_M15(U_p)/15,
H'_m = H_m+kappa P = H_E+K_m+kappa S,
||K_m|| <= B=N(2+epsilon_m)/a.
```

The scalar shift does not change evolution or detector probabilities.
For the same prepared state let `E0=(29/10) E L²` and use (3) for S0.
Positivity of electric energy and conservation of H'_m give

```
kappa <S>_t <= <H'_m>_t+B
             =<H'_m>_0+B <= E0+kappa S0+2B.
```

Together with (1), this yields the time-independent estimate

```
Prob(Bad at t) <= min(1,48000[S0+(E0+2B)/kappa])          (4)
```

for every real model time. The state has finite electric energy and lies
in the self-adjoint Hamiltonian domain, so energy conservation is applied
on its actual domain. This is not an external Gibbs law or an assumed
stationary reference.

An explicit admissible coefficient is

```
kappa = 2*48000*(E0+2B)/eta.                             (5)
```

The planner then proves (4)<=eta. This coefficient is fixed before the
observation begins. The bound requires no clock rescaling, cooling during
the experiment, stochastic gauge averaging or repeated gap measurements.
For kappa=1/100 the same conservative energy estimate is trivial. That
failure of the estimate is not a proof of physical gap loss in that
smaller-coupling program.

With eta=2^(-b), the exact-sign comparison for any chosen finite duration
T becomes

```
complete_final_record_TV
 <= min(1,(N T/a)[epsilon_m+(2+epsilon_m)2^(-floor(b/2))]). (6)
```

It compares the same prepared state, electric term, plaquette coefficient
and clock; only the sign polynomial changes. As in the owning proof, it
has one evolution interval and a complete final detector. Intermediate
interventions that inject energy require their explicit work allowance
in (4). A classical copy of a pointer is not a new physical gap bound.

The same all-time estimate holds directly on every finite source word
code `V_K` containing the prepared state. Compressing the positive
operator inequality (1) gives

```
P_K 1_Bad P_K <= 48000 P_K S P_K.
```

The compressed electric energy is positive, the kinetic norm stays at
most B, and the initial expectations are unchanged. Conservation of the
finite Hamiltonian `P_K H_m P_K` therefore gives (4) for the complete
two-outcome detector with effect `P_K 1_Bad P_K`. Its other effect is the
identity on the code minus that effect. This is not a statement about
the spectral gap of the truncated Wilson matrix `P_K H_W P_K`.
The certificate bounds this effect without materializing its matrix or
compiling its detector; a finite precision implementation needs its
operator-error allowance.

The executable certificate uses a 40-bit eta target, m=128 and T=1 as one
finite example. Its code dimension, electric-energy bound and selected
kappa are enormous. The analytic all-time gap statement does not imply
that a practical simulation or observation can attain that example.

## Conditional local roughness from energy density

There is a volume-independent version of the energy argument for local
plaquette roughness, under a different state hypothesis. Keep the same
periodic Hamiltonian, fixed spacing, polynomial order and coefficient
kappa>0. Let V be the number of cells, and suppose the initial joint
state has a source marginal invariant under every lattice translation.
Its electric energy is finite. Suppose also that the kinetic term obeys
`||K_m||<=b V` on an invariant sector containing the prepared state.
The full Fock carrier permits the uniform choice
`b=60(2+epsilon_m)/a`; a smaller conserved occupation density can improve
this bound. An untouched reference may be correlated with the source;
only its source marginal needs translation symmetry. Observation
evolution has no intervening drive or measurement.

Write the actual initial densities and individual plaquette costs as

```
e_E = <H_E>_0/V,     s_0 = <S>_0/V,
s_(x,mu,nu) = 1-Re Tr_M15(U_(x,mu,nu))/15 >=0,
S = sum_x sum_(mu<nu) s_(x,mu,nu).
```

Translation symmetry is preserved by the same Hamiltonian. Positivity
of H_E and conservation of H'_m give, for every model time,

```
<S>_t/V <= s_0+(e_E+2b)/kappa.                          (7)
```

For any fixed orientation mu<nu, translation symmetry makes all V
expectations of that orientation equal. The other two orientation sums
are nonnegative, so each individual plaquette satisfies

```
<s_(x,mu,nu)>_t
 = V^(-1) sum_y <s_(y,mu,nu)>_t
 <= <S>_t/V <= s_0+(e_E+2b)/kappa.                     (8)
```

There is no division by three in (8). Equality of the three orientation
means would require an additional rotation-symmetry hypothesis.
For the rough-plaquette effect
`R_p=1_{||R(U_p)-I||>1/40}`, the exact pointwise inequality
`R_p<=48000 s_p` follows from
`s_p=||R(U_p)-I||_F²/30`. Therefore

```
Prob_t(rough plaquette p)
 <= min(1,48000[s_0+(e_E+2b)/kappa]).                   (9)
```

If e_E, s_0 and b are bounded uniformly along a specified state family,
the right side is independent of volume at fixed kappa. A small bound
requires an actual small magnetic density and the stated electric and
kinetic allowances; the algebra does not supply those state estimates.

The same argument holds directly on a finite total-word code V_K. Its
projection P_K commutes with translations, so the compressed Hamiltonian
preserves the symmetry. Electric compression is positive and kinetic
compression preserves the norm bound. Use the costs `P_K s_p P_K` and
the complete two-outcome detector effect `P_K R_p P_K`, whose complement
is `I_(V_K)-P_K R_p P_K`. Compressing the positive inequality proves
(9), with the actual compressed initial energy densities. This effect
need not be a projection and is not a roughness test on a separately
truncated link or Wilson matrix. The proof retains the complete field
and reference output; it does not condition on the smooth outcome.

The specified filled-origin input breaks translation symmetry. Its
gauge amplitude is symmetric, but its interacting joint source marginal
is not, so it does not instantiate (8). Averaging its full density over
all translations defines a distinct, admitted mixed preparation to which
the finite lemma applies. That averaging preserves its total energy and
magnetic expectations; it does not establish volume-uniform density
bounds for the planner whose L itself grows with volume. The local
energy-density strategy in
[the native spin limit](../native_repair_spin_limit/PROOF.md#exact-conservation-and-the-energy-estimate)
and [its ground and thermal references](../spin_reference_limit/PROOF.md#charge-density-from-source-energy)
is a useful proof pattern. Those theorems use a two-flavor onsite penalty
and a controlled spin embedding. They do not instantiate the present
gauge-coordinate hypotheses.

### The missing local operator estimate

Roughness of one plaquette is not a local Wilson-gap event. Summing (9)
over every plaquette to guarantee a global gap restores a volume factor.
The global admissibility theorem alone therefore cannot convert (9)
into a volume-independent sign-response estimate.

The local spectral quantity required by an alternative route can be
stated exactly. At a coordinate configuration U, put
`Pi_low(U)=1_{|H_W(U)|<1/2}` and let P_X project onto the finite one-body
spin/internal modes in a spatial set X. The scalar spectral theorem gives

```
||(sign(H_W)-q_m(H_W)) P_X||_HS²
 <= epsilon_m² rank(P_X)
    +(2+epsilon_m)² Tr(P_X Pi_low P_X).                (10)
```

Equation (10) is an algebraic consequence of the owning polynomial
bound, at every finite coordinate configuration. It isolates local
low-mode weight rather than requiring that the entire lattice have no
low mode. No bound on `Tr(P_X Pi_low P_X)` follows here from (9).
Spatial decay of the low-mode projectors and control of their dependence
on the quantum links are also separate hypotheses; a pointwise matrix
kernel bound does not by itself control the electric dynamics.

A sufficient dynamical target is the following local Duhamel bound.
For `Delta H=H_sign-H_m`, an effect `0<=A_(X,R)<=I` on the observed
cells and an untouched reference R, let `rho_(m,V,R)(u)` be evolution
under H_m from the same prepared joint state. One needs a bound of the
form

```
|Tr rho_(m,V,R)(u) i[Delta H tensor I_R,
                    (alpha_sign^(T-u) tensor id_R)(A_(X,R))]|
 <= c_(X,T) [epsilon_m+f(eta_local)],    0<=u<=T,        (11)
```

uniformly in the exterior volume and reference, with `f(eta)->0` as
eta->0. Here eta_local must be derived from this state's local
roughness/low-mode estimates, and c_(X,T) must have a proved volume bound.
The same-state Duhamel identity would then bound the difference of
complete local final records by T times (11), before the separate volume,
cutoff and apparatus errors. Neither (9) nor (10) proves (11). A
uniform interacting relativistic refinement also needs control as
spacing and polynomial order change, and nonvanishing interacting
responses on that trajectory.

Finally, the Wilson gap is a regulator singular-value gap. It is not
a physical many-body excitation gap or a particle mass. The
[free flat source](../interacting_extension/chiral_lattice/PROOF.md#complete-free-spatial-overlap-control)
has `|H_W|>=1` and a zero-energy overlap mode at spatial momentum zero.

## The complete finite-code comparison

The Haar coordinate Hilbert space is the analytic comparison parent.
The finite prepared state belongs to the actual source word code, and
its finite evolution is the compression of the parent's polynomial.
For a local observation use the following order of comparisons:

1. Periodic exact-sign dynamics to the periodic polynomial, using (6).
2. Periodic polynomial to its interior open-box polynomial, using the
   [same-field locality bound](../particle_field_causality/PROOF.md).
3. Open-box polynomial to its finite representation compression, using
   the owning complete-history cutoff estimate.

The open-box reference is the entire box-plus-reference marginal of the
periodic preparation. It need not be a box-only Gauss singlet; boundary
charge correlations must be retained. The preparation has degree at most
L on each link, and taking a marginal preserves that finite local word
support. With all outgoing links retained, initial total degree at most
`3|box|L` is sufficient for the cutoff bound. Boundary link registers on
deleted edges are electric spectators. Periodic wraparound must stay
outside the selected collar.

An open Wilson box can have boundary modes, so the periodic gap theorem
is not applied to its truncated kernel. The separate comparisons above
avoid that error. Their complete-output bounds add, without conditioning
on a detector result. In particular, the large kappa in (5) must also be
used in the locality and cutoff constants; the constants of the 1/100
control cannot be imported unchanged.

## What the other finite theorems can and cannot supply

| Reused result | Valid input to this field route | Additional premise for a relativistic claim |
| --- | --- | --- |
| [Spatial free-overlap branch](../interacting_extension/chiral_lattice/PROOF.md#complete-free-spatial-overlap-control) | Exact free Fourier gap, full momentum census and low-momentum Dirac-symbol estimate | The flat coordinate fiber is not the normalizable interacting state; convergence of its complete quantum response is required |
| [Euclidean electric bundle](../interacting_extension/chiral_electric/PROOF.md) | Gauge-covariant derivative on a fixed-rank admitted chiral bundle | Common canonical many-body chiral Hamiltonian and its state; the bundle projector is not the spatial modified axial charge |
| [Fixed-spacing locality](../particle_field_causality/PROOF.md) | Full-Fock graded influence bound and complete local-record volume comparison | Constants controlled along the actual joint degree, spacing, volume and coupling family; no boost symmetry follows from a cone bound |
| [Free Weyl continuum](../continuum_free_limit/README.md) | A rigorous free continuum on its own centered source and scaling | It has a different action/state and does not supply this interacting gauge-coordinate limit |
| [Ordered canonical field](../spin_canonical_field/README.md) | Positive energy, translation and canonical limiting observables on its ordered spin source | A state/current/clock-preserving connection to this coordinate field, including its surviving interactions |
| [Finite causal observer nets](../../Lean/QFT/FiniteCausalObserverNet.lean) | Compatible local observable laws when its region and repair hypotheses hold | Its declared region system does not prove spacetime dimension, Lorentz transformations or continuum dynamics |

The public coordinate proof and the supplied particle-extension source
use the same separation between a free flat fiber and the fluctuating
gauge parent. Reusing that free calculation is appropriate; replacing
the interacting state by its flat fiber is not a source-preserving step.

The closest justified continuation is a common family of quantum local
correlations with finite approximation moduli at fixed physical scales.
Its dynamics must control canonical chiral sectors, surviving gauge and
matter interactions, wavefunction/parameter normalization and the actual
energy-momentum response. Only then can relativistic transformations and
particle thresholds or resonances be attached to these records.
Equation (5) supplies a sufficient all-time finite-gap family, but does
not prove those uniform estimates. In particular, increasing kappa may
freeze a sector rather than preserve the interactions needed in the
target limit. Architecture choices are permitted; their measured outputs
and any refinement trajectory require their own verification.
