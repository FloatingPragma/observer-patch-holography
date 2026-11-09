# Spatial coframes and the response of the interacting source

A declared spatial coframe specifies six strain responses of the actual
finite-range Wilson gauge action. The derivative acts on its Dirac,
Wilson, electric and magnetic terms, with the original quantum links
and every fermion occupation retained. A complete quantum-flux control
shares its Hamiltonian, Gibbs target and clock with the
[state preparation](../particle_field_state_response/PROOF.md).
This is an action-derived mechanical response to a specified geometry
input. Identifying that input with the observers' proper geometry is a
separate premise.

## 1. Source and constant spatial geometry

Use the [coordinate field](../particle_observer_outputs/FIELD_ROUTE.md)
with all sixty Dirac/internal modes per cubic cell. Its allowed
[finite-range Wilson choice](../particle_smooth_sign_field/PROOF.md)
has flat one-body action

```
h_I=(beta+H_W)/a
   =sum_i alpha_i D_i/a + beta sum_i L_ii/(2a),
alpha_i=i beta gamma_i,
D_i=(T_i-T_i^dagger)/(2i),
K_i=I-T_i,   L_ij=(K_i^dagger K_j+K_j^dagger K_i)/2.
```

T_i transports matter on direction i with the actual coordinate link.
Each T_i is unitary on the periodic parent carrier. A gauge change U
conjugates all T_i, so D_i and L_ij transform by the same conjugation.
L_ii=2I-T_i-T_i^dagger; this makes the displayed flat action exactly
the original Wilson action. The smooth-sign and retained Legendre
actions are distinct choices and are not relabeled as this action.

Declare a constant real coframe V with positive determinant. A coordinate
edge has physical coordinate vector a V e_i in this declared dictionary.
Set

```
G=V^T V,       v=det(V)>0,
h_V=sum_(A,i) alpha_A (V^(-T))_(Ai) D_i/a
       +beta sum_(i,j)(G^(-1))_(ij)L_ij/(2a).         (1)
```

The bare a in (1) is fixed during a strain experiment. In particular,
the Wilson coefficient has the metric dependence of a fixed regulator
length times a Laplacian. A simultaneous change of that regulator
length would be a different declared family.

All h_V are Hermitian. Their full second quantizations use the entire
Fock carrier. No momentum, occupation or mirror sector is discarded.
The source derivative on the Fock carrier is dGamma(partial h_V), by
linearity of finite second quantization.

## 2. Electric and magnetic terms on the same dictionary

Anchor the left electric generators E_i^r(x) of the three outgoing
links at each site x. Under a site gauge change these generators transform
in the same adjoint representation. Use its fixed invariant inner
product and define the symmetric electric Gram form

```
E_ij(x)=sum_r (E_i^r(x)E_j^r(x)+E_j^r(x)E_i^r(x))/2,
H_E(V)=(g_E/2)sum_(x,i,j) [G_ij/v] E_ij(x).          (2)
```

The coefficient g_E is the declared original electric normalization.
At V=I this is the sum of the original electric Casimirs. Cross terms
are between outgoing links at the same site; no parallel transport or
extra link multiplication is inserted into them. The finite carrier
uses the unchanged coordinate-word compression P_K, which preserves
gauge transformations and the electric generators. Positivity follows
by factoring the positive matrix G/v and writing the form as a sum of
squares of real linear combinations of the E_i^r.

On the uncut compact-group carrier, smooth functions are a common form
core and the positive quadratic form defines its Friedrichs extension.
For V in a compact set of positive coframes its form norms are equivalent
to the flat electric form norm. This is a form statement; no bounded
locality theorem is used to assert an operator-domain identity. The
finite matrix statements below require no unbounded differentiation.

For magnetic terms, let P_i(x) be the positive oriented plaquette in the
plane perpendicular to i, based at the same site x. Set C_i=I-R(P_i),
where R is the actual fifteen-row representation, and define the real
symmetric coordinate multiplication functions

```
Q_ij(x)=Re Tr[C_i(x)^dagger C_j(x)]/30 - delta_ij,
H_B(V)=kappa sum_(x,i,j) [G_ij/v] Q_ij(x).            (3)
```

Conjugation at the common base preserves every Q_ij. Unitarity gives
Q_ii=-Re Tr R(P_i)/15. Therefore (3) at V=I is exactly the original
magnetic potential, with its specified energy zero. The off-diagonal
Gram terms and metric-dependent subtraction are declared geometry
couplings. They are not consequences of knowing the flat potential.
The positive Gram part is bounded below and the subtraction is a
bounded finite-volume scalar.

The complete finite source is

```
H_K(V)=P_K[H_E(V)+dGamma(h_V)+H_B(V)]P_K.             (4)
```

Compression follows the complete parent expressions. In particular,
products of separately compressed plaquette links are not substituted
into C_i^dagger C_j. Such a substitution can destroy the flat identity
at cutoff endpoints. Equation(4) is a Hermitian, gauge-preserving
family on one finite carrier with one clock. At V=I it is precisely
the chosen finite-range Wilson source, not an inferred continuum field.

## 3. Six strain operators and the dilation identity

Let V(s)=I+s S with real symmetric S. At s=0,

```
dV^(-T)/ds=-S,      dG^(-1)/ds=-2S,
dv/ds=Tr S,         d(G/v)/ds=2S-(Tr S)I.            (5)
```

Differentiating the actual action gives K(S)=dH_K(V(s))/ds at zero:

```
K(S)=P_K{dGamma[-sum_(A,i) alpha_A S_Ai D_i/a
                         -beta sum_(i,j)S_ij L_ij/a]
 + (g_E/2)sum_(x,i,j)[2S_ij-(Tr S)delta_ij]E_ij(x)
 + kappa sum_(x,i,j)[2S_ij-(Tr S)delta_ij]Q_ij(x)}P_K. (6)
```

The three matrices S=e_i e_i^T and three matrices
S=(e_i e_j^T+e_j e_i^T)/2 give six independent symmetric strain
operators K_ii and K_ij. Their convention is
K(S)=sum_i S_ii K_ii+2sum_(i<j)S_ij K_ij. Each is Hermitian and
gauge invariant on the whole retained carrier. The code checks every
one of the sixty native spin/internal columns with noncommuting weak
links, including every shear. The coordinate configuration is an
algebraic control, not a normalizable quantum-link preparation.
Its displayed sixty-column total adds the magnetic scalar on the
one-particle sector. In the full field, only h_V is second-quantized;
the magnetic term multiplies the Fock identity and is not multiplied
by fermion number.

For isotropic V(s)=exp(s)I, write H_D for the first term of (1) and
H_W for its Wilson term. Equations(1)-(3) give exactly

```
H_K(s)=exp(-s)(H_E+H_D+H_B)+exp(-2s)H_W,
K_iso=sum_i K_ii=-H_E-H_D-H_B-2H_W.                 (7)
```

The Wilson regulator contribution is explicit. Replacing the right
side by -H without that extra H_W would give the wrong response for
this declared fixed-regulator family.

Constant physical frame rotations have a separate finite Ward identity.
For R in SO(3), G(RV)=G(V), and the source Dirac spin representation
S_R obeys S_R alpha_B S_R^dagger=sum_A R_AB alpha_A and
S_R beta S_R^dagger=beta. Hence

```
H_K(RV)=Gamma(S_R) H_K(V) Gamma(S_R)^dagger.          (8)
```

This is covariance of the declared frame coupling. Its spin lift and
gauge rotations act on different factors. It establishes neither
Lorentz boosts nor the observed spin of a particle. The code checks
an exact quarter-turn on a rational nonorthogonal coframe, independently
of its gauge-covariance check.

## 4. Same-state work, pressure and complete detector laws

For a differentiable controlled V(t), Schrödinger evolution under (4)
gives d Tr[rho(t)H_K(V(t))]/dt=Tr[rho(t)partial_t H_K]. At a fixed
prepared density, the initial work slope in direction S is Tr(rho K(S)).
The [same-action energy construction](../particle_gravity_response/PROOF.md)
and its finite [coarse-energy theorem](../../Lean/QFT/CoarseEnergyWard.lean)
apply to the full deformed action. Their sum identity does not require
reproving an energy law or identifying a number current with stress.

For the ideal Gibbs target rho_beta(V)=exp(-beta H_K(V))/Z(V), define
F(V)=-log Z(V)/beta. Finite trace differentiation gives

```
dF/ds=Tr(rho_beta K),
dS(rho_beta)/ds=-beta²[Tr(rho_beta H K)
                       -Tr(rho_beta H)Tr(rho_beta K)]. (9)
```

The entropy formula uses the same H, beta and K at the base coframe.
It holds with noncommuting H and K because H commutes with its Gibbs
density; it is not obtained by replacing a general modular operator
with beta H. These equilibrium derivatives compare separately prepared
states at the same beta. Isolated unitary driving does not in general
follow that Gibbs family; a bath, cooling law or quasistatic preparation
would require its own implementation and work accounting. The positive
polynomial preparation in the companion
package is rho_M, a controlled approximation to rho_beta. Its modular
Hamiltonian is -log rho_M, not beta H. Work and detector probabilities
computed on rho_M retain their certified target-state error.

The declared physical volume of one coordinate cell at bare a=1 is
v=det V. Along the isotropic curve dv/ds=3v. Thus the operational
pressure for this geometry coupling and Gibbs target is

```
p_beta=-(dF/ds)/(3v)=-Tr(rho_beta K_iso)/(3v).        (10)
```

The energy density is Tr(rho_beta H)/v in the same energy convention.
Both come from a specified state and action variation. A fixed run
with an arbitrarily chosen volume function supplies neither (9) nor
(10). This finite response is not a cosmological equation of state,
and v has not been identified with observer-derived proper volume.

## 5. The unchanged208-state quantum-flux control

Use exactly `field_route.exact_rotor(6)`: every one of the sixteen
four-Dirac occupations and thirteen fluxes -6,...,6 is retained.
The source coefficients and q6=6 charge are unchanged. Write

```
H=H_E+H_W+H_D,
H_E=I tensor E²/2,
H_W=onsite tensor I + ((hop+hop^dagger)/2) tensor(S6+S6^dagger),
H_D=((hop-hop^dagger)/2) tensor(S6-S6^dagger).         (11)
```

This is the one-direction Wilson expression inherited from (1),
with declared magnetic coefficient kappa=0, compressed only after its
parent shift is formed. Setting transverse
classical links to identity and removing transverse terms gives (11).
It is an algebraic directional reduction, not a claimed normalizable
reducing subspace of the full three-dimensional quantum-link system.
In particular the full six-component spatial proof and this finite
thermal control have the explicitly different carrier domains stated here.
At nonzero kappa, identity plaquettes would retain the magnetic constant
and its coframe response; they cannot be silently erased. The sixty-column
magnetic covariance control uses kappa=1/7 and is not the kappa=0 thermal
control.

Under V_parallel(s)=diag(exp(s),1,1), the retained electric coefficient
G_11/v is exp(s). Its actual action and slope are

```
H_parallel(s)=exp(s)H_E+exp(-s)H_D+exp(-2s)H_W,
K_parallel=H_E-H_D-2H_W.
```

The isotropic curve instead has K_iso=-H_E-H_D-2H_W, exactly as (7).
Both perturb the same original H at s=0. The companion target uses
beta=1/64 and the same unaltered source clock. Its source Gibbs density
contains all Gauss charge sectors, with a conjugate reference retained
for a joint invariant preparation. Every strain operator commutes with
the source Gauss generator, so its detector also preserves the joint
constraint without discarding any source charge sector.

The exact bounds ||H_E||<=18, ||H_W||<=4 and ||H_D||<=2 imply
||K_parallel||,||K_iso||<=28. Each strain therefore has a complete
two-outcome apparatus with effects (I+K/28)/2 and (I-K/28)/2 and
record labels +28,-28. Its mean is Tr(rho K), retaining all states,
outcomes and references. It is a mean-reading instrument, not a
projective claim that K has only two eigenvalues. The standard complete
pointer dilation realizes these effects as an admitted apparatus.
The longitudinal and isotropic menus are alternative experiments on
identically prepared inputs; their operators need not commute, and
their mean records are not assigned simultaneous classical values.

If the certified preparation trace distance to the Gibbs target is
epsilon, the complete detector total variation is at most epsilon,
the strain-mean error is at most56epsilon, and the pressure error at
v=1 is at most56epsilon/3. Actual preparation/compilation and dynamical
errors, if present, are added under the same contract. The code prints
exact rational work and pressure records of the prepared rho_M and
keeps the ideal Gibbs interpretation separate.

## 6. What the geometry consumer additionally requires

There is a precise limit to recovering geometry from a flat Hamiltonian.
For any real xi, H_tilde(V)=H_K(V)+xi log(det V)I has exactly the same
flat source H_K(I), states and all fixed-flat detector records. Its
strain shifts by xi Tr(S)I; its isotropic pressure shifts by -xi/v.
Thus flat records alone cannot determine a unique metric coupling or
gravitational vacuum stress. This is not a uniqueness requirement for
the construction: the coframe family (1)-(4) may be declared and its
consequences calculated. Claiming that native observers reconstruct
that family requires the corresponding operational identification.

The [ten-probe Einstein consumer](../matched_einstein_limit/PROOF.md)
is reused rather than re-proved. A spatial mechanical tensor from (6)
does not supply its Lorentz stress, curvature, proper clocks, modular
diamond kernel, fixed-proper-volume area law or shrinking-ball bounds
on one interacting state. The separate free-Weyl geometric state is
not substituted into these interacting formulas. No Einstein equation,
physical Newton constant, particle mass or species identity follows
from this finite coframe response alone.

## 7. Causal control on compact coframe families

The constant-coframe Wilson family also admits the
[refinement-uniform causal argument](../particle_physical_causality/PROOF.md).
Each lattice cell owns its three outgoing link factors. Consequently every
electric cross term in (2) is onsite in this tensor decomposition. Its
positive form defines a self-adjoint onsite Friedrichs generator. Add the
bounded onsite matter terms, including any self-loop bonds on collapsed
periodic axes, to that generator. The interaction picture preserves every
remaining support, norm and parity. Bounded spectral clippings and their
finite-volume strong-unitary limit give the same bound for the uncut
electric carrier. No electric coefficient or electric-energy expectation
enters the interaction rate.

Use the periodic or infinite parent with unitary shifts. An open-box
restriction retains/deletes terms of this parent interaction decomposition
with its stated boundary onsite terms. Replacing the parent products by
products of truncated nonunitary shifts would define another boundary
action and is not assumed here.

Write B=V^(-T), C=G^(-1)=B^T B and A=G/v. Combine the Dirac and diagonal
Wilson terms of each nearest bond. Apart from its unitary link, its
one-body block is

    [-beta C_ii - i alpha dot B_i]/(2a).

Its spin singular values are |C_ii plus/minus sqrt(C_ii)|/(2a), each
twice. All fifteen internal rows give full-Fock bond norm
30 max(C_ii,sqrt(C_ii))/a. A self-loop can be absorbed onsite; the same
upper bound also follows from zero spin trace and the trace-norm triangle
inequality. Two incident bonds per direction, each of support size two,
give weighted rate 120 lambda sum_i max(C_ii,sqrt(C_ii))/a.

For i<j, decompose the mixed Wilson term into root stars on
{x,x+e_i,x+e_j}. A star has coefficient beta C_ij/a times
(k_i^dagger k_j+k_j^dagger k_i)/2, where k_i is the row of I-T_i.
If these three cells are distinct, gauge changes on this tree reduce its
coordinate factor to the matrix with eigenvalues 3/2,-1/2,0. Its full-Fock
norm is 60|C_ij|/a. Collapsed periodic cells are covered by the universal
120|C_ij|/a bound: on the fifteen-row coordinate factor each k_i has
Hilbert-Schmidt norm squared at most60, and the symmetrized product has
trace norm at most60. The beta factor has two signs of multiplicity two.
There are at most three incident roots per pair and support size at most
three, so the universal weighted contribution is
1080 lambda^2 sum_(i<j)|C_ij|/a. The distinct-cell constant540 requires
that additional geometric premise and is not used by the code.

The diagonal magnetic interactions retain the original conservative
16 lambda^2 kappa sum_i A_ii bound. For i<j, ||Q_ij||<=2, so the combined
pair 2 kappa A_ij Q_ij has norm at most4 kappa |A_ij|. Both based
plaquettes together use only outgoing links owned by
{x,x+e_1,x+e_2,x+e_3}; this actual tail-cell support has size at most four
and diameter at most two. A cell belongs to at most four such rooted
stars. Their weighted rate is therefore at most
64 lambda^2 kappa sum_(i<j)|A_ij|.

For kappa(a)<=kappa0/a and lambda>1, a rational majorant is

    W(V)=120 lambda sum_i max(1,C_ii)
          +1080 lambda^2 sum_(i<j)|C_ij|
          +kappa0 lambda^2 [16 Tr(A)+64 sum_(i<j)|A_ij|].

Here sqrt(C_ii)<=max(1,C_ii). The complete weighted interaction norm is
at most W(V)/a, and W(I)=360 lambda+48 kappa0 lambda^2 exactly recovers
the existing flat Wilson rate. This is an extension of that action, not
a transfer of its sharp constants to the distinct smooth regulator.

The exported compact-family example uses kappa0=1/7, so its a=1 member
includes the actual sixty-column magnetic control at kappa=1/7.
The directional thermal experiment keeps its separate kappa=0 action.

For constant coframes with m^2 I<=G<=L^2 I, m>0, every entry of C has
absolute value at most m^(-2), while every entry of A has absolute value
at most L^2/m^3. Thus one bound uniform over the compact family is

    W_K=360 lambda max(1,m^(-2))+3240 lambda^2/m^2
          +240 kappa0 lambda^2 L^2/m^3.

The code checks membership through every nonempty principal minor of
G-m^2 I and L^2 I-G, fourteen minors in total. Leading minors alone
would be insufficient for the semidefinite inequalities.

Use Euclidean coframe distance on an open carrier and quotient Euclidean
distance, minimized over periodic lifts, on a torus. Plain unwrapped
coordinate distance across a periodic seam is not an allowed separation.
For graph distance D, both conventions obey d_V<=LaD. The same complete
local-history estimate therefore gives

    half_trace(output_U,output_I)
      <=min(1,2(2k+1)|X| exp[-(log(lambda)d_V/L-2W_K T)/a]).

The clock is the original Hamiltonian duration. A coframe-coordinate cube
of physical half-width R is contained in a coordinate cube of half-width
2R/m, so |X|<=(4R/m+3)^3/a^3 for a<=1. The same monotone-envelope proof
gives an all-smaller-spacing threshold with this prefactor and margin
[(lambda-1)/lambda]d_V/L-2W_K T. The conservative physical speed is at
most2L W_K/log(lambda); replacing log(lambda) by (lambda-1)/lambda gives
the exact rational upper bound in the certificate. Both compared finite
implementations add their own complete-experiment errors.

This result concerns constant spatial coframes and the same ordinary
Wilson interaction. It supplies uniform local-record control under
declared metric variation, without a convergent state family, proper
geometry, Lorentz boosts, a chiral limit or Einstein dynamics. The
208-state thermal control is not substituted for a spatial refinement
state by this bound.
