# From original energy transfers to excitation and mass readouts

## The fixed source and the input readout

Use the actual finite physical carrier of the
[primary Wilson contract](../particle_common_contract/PROOF.md), its
Hamiltonian H with unchanged energy zero, its stationary polynomial density
rho, and its full conjugate reference with idle Hamiltonian zero. Write E0
for the lowest eigenvalue on this physical carrier. The selected preparation
proves

```
[rho,H]=0,     Tr rho=1,     Tr rho(H-E0)<=epsilon_n=2^(-n-4).
```

The [common filtered algebra](../particle_primary_field_algebra/PROOF.md)
and [complete filtered apparatus](../particle_filtered_probe/PROOF.md)
supply the same bounded physical probe B and its original source-time
correlation `F(t)=Tr rho B* alpha_t(B)`. Its finite positive measure mu
records energy transfers, which can have either sign on a mixed input.
All filters, failure outcomes and untouched quantum references belong to
the stated experiment. The [positive-window theorem](../particle_positive_window/PROOF.md)
and [complete spectral readout](../spin_shell_readout/PROOF.md) give existing
weight/error methods on their stated sources; this proof does not repeat
their apparatus or transfer a spin-source observation to the primary state.
A supplied primary window weight must have its own complete source,
sampling, clock and apparatus allowances subtracted.

## A positive coupling with final excitation energy

Choose a common eigenbasis of H and rho, with eigenvalues E_i and p_i.
For degeneracies this is only a proof coordinate choice. Set

```
w_ij=p_i |B_ji|^2,
omega_ij=E_j-E_i,       z_ij=E_j-E0,       x_i=E_i-E0>=0.
```

The positive finite measure assigning weight w_ij to `(omega_ij,z_ij)`
has first marginal mu and second marginal

```
nu(S)=Tr[B rho B* 1_S(H-E0)].
```

Both original masses equal `s=Tr rho B*B<=C^2` when `||B||<=C`.
They are not divided by s. In particular s=0 gives two zero measures
and needs no separate conditional-state definition. Since `z-omega=x_i`,

```
integral |z-omega| dw
 =sum_i p_i x_i <i|B*B|i>
 <=C^2 Tr rho(H-E0) <=C^2 epsilon_n.                 (1)
```

This is a same-source transport bound. For every Lipschitz real or complex
function g for which the finite integrals are defined,

```
|integral g(z) dnu-integral g(omega) dmu|
 <=Lip(g) C^2 epsilon_n.                             (2)
```

The argument does not measure the initial energy, dephase the canonical
reference, prepare a ground vector, or counter-evolve the reference.
It is a positive mathematical coupling of two source-marginal spectral
laws. The complete apparatus remains the one used to establish the input
correlation. Source stationarity is required; stationarity of its retained
purification is neither used nor asserted.

If B is a triangularly filtered probe of half-width T with unsmeared norm
cap C, the reused bound `integral |omega| dmu<=2C^2/T` yields

```
integral z dnu <=2C^2/T+C^2 epsilon_n.                 (3)
```

Thus filtered final excitation measures have a uniform first moment for
fixed C,T. On the common subsequence where the original correlations
converge, (2) implies that mu and nu have the same weak limit: bounded
Lipschitz tests determine weak convergence of finite measures, their
masses agree, and the moment bound controls the tails. This identifies the
limiting positive-energy response with final ground-relative excitation
weight on that same source. It supplies no nonzero limiting weight.

## Window weights and the energy origin

If `mu([l,u])>=w`, `0<=l<=u`, and r>0, the original weight with x_i>r is
at most `C^2 epsilon_n/r`. All remaining pairs in the transfer window
have final excitation in `[l,u+r]`. Therefore

```
nu([l,u+r]) >=max(0,w-C^2 epsilon_n/r).                (4)
```

All weights are original. The primary choice `r_n=2^(-floor(n/2))`
makes both the energy enlargement and its weight allowance vanish;
the latter is `C^2 2^(-ceil(n/2)-4)`. A computed uniform positive w
would therefore survive this conversion. The package supplies the
actual epsilon and allowance, without inventing that primary observation.

The true E0 enters the proved cooling inequality and the definition of
nu. It is not numerically computed for the enormous primary matrix.
Equation (4) infers a ground-relative interval from transfer records and
the proved cooling certificate; it needs no numerical E0 estimate.
It is not a direct absolute-energy readout. A different experiment that
measures final H-energy in `[a,b]` and supplies a certified ground interval
`E0 in [g_l,g_u]` gives

```
E_final-E0 in [max(0,a-g_u), b-g_l].                  (5)
```

Without that ground interval, the direct-subtraction adapter returns
unknown. This separate requirement is not added to (4). Adding a scalar
to H and E0 leaves (1)-(4) unchanged, but the actual source energy zero
and controlled clock convention retain their contract values.

The [primary finite-survival theorem](../particle_primary_survival/PROOF.md)
supplies a proved positive original band `[delta_n,W_n]`, where W_n is
the full certified source width U_n-L_n. Since E0>=L_n and E_j<=U_n,
every final excitation obeys `0<=z<=W_n`. Also z>=omega. Consequently

```
nu([delta_n,W_n]) >= mu([delta_n,W_n]).                (5a)
```

This broad finite band needs no Markov subtraction. Its source-authenticated
consumer retains both the ideal triangular filter's lower weight and the
positive lower weight of its actual finite Cayley/quadrature blank map,
as separately proved by the survival theorem. Both use spectral projectors
of the unchanged H; the finite gates do not replace that generator.
The consumer keeps these original weights factored without constructing
their enormous denominators. Both the lower endpoint and lower weights approach
zero, and its upper endpoint diverges; (5a) supplies no fixed-band
noncollapse or rest mass. The tighter upper endpoint in (4) has the
explicit cooling cost needed for a fixed-energy inference.

## What turns excitation energy into rest energy

The finite primary box is open. Its coordinate or Fourier probe labels
are not conserved translation momentum. The
[periodic joint translation readout](../particle_joint_translation_readout/PROOF.md)
uses a different torus, state and joint field/reference translation.
Its canonical input's zero total crystal momentum does not assign zero
momentum to an individual excitation. The common primary time algebra
alone supplies no canonical spatial momentum operator or vacuum momentum.

When commuting energy and spatial translation generators and a density
invariant under all of them are supplied, a mixed-state joint correlation
measures `(E_final-E_initial,p_final-p_initial)`.
Cooling alone does not select the initial momentum. If a separate source
estimate gives initial momentum within radius p0 of the admitted vacuum
momentum, except with probability eta_p, then observed transfer magnitude
`|q| in [q_l,q_u]` gives

```
|p_final-p_vac| in [max(0,q_l-p0), q_u+p0],
additional original weight loss <=C^2 eta_p.         (6)
```

This uses the triangle inequality and `B*B<=C^2 I`. Combining (4) and
(6) uses one joint event and subtracts both tail allowances. Separate
marginal observations require their own intersection bound; they cannot
be multiplied as if independent. Fourier labels without a proved spatial
generator do not satisfy the premise of (6).

Suppose the retained same-branch event has final excitation energy
`E in [e_l,e_u]`, momentum magnitude `p in [p_l,p_u]`, speed
`v in [v_l,v_u]`, and a separately justified squared-shell error

```
|E^2-v^2 p^2-M^2|<=eta,       M>=0.
```

Then elementary monotonicity gives

```
max(0,e_l^2-v_u^2 p_u^2-eta) <= M^2
 <= e_u^2-v_l^2 p_l^2+eta.                            (7)
```

The code uses inherited outward rational square roots for M. A negative
upper bound or reversed interval rejects inconsistent premises. M is rest
energy; the model inertial mass is M/v^2. A broad band, even with positive
weight and positive (7), does not establish an isolated particle pole,
a common shell across different momenta, or a physical species. Those
attachments must be provided before calling (7) a particle mass result.
The matched quadratic gauge/matter characteristics alone do not provide
this interacting spectral shell.

The independent three-level control has H=diag(0,3,4),
rho=diag(15/16,1/16,0), and B=|2><1|. Its complete complement is
diag(1,0,1). The original response has weight 1/16 at transfer energy1,
while the final excitation is at4. Its transport cost3/16 saturates (1).
Declaring commuting control momenta p_initial=p_final=1 gives transfer
momentum0; the transfer invariant is1 while the final excitation invariant
is15 at v=1. These are separate finite controls, not primary momentum or
particle assignments. The exact source/reference instrument retains the
failure weight15/16 and all reference degrees of freedom.

## Units are explicit comparison inputs

The primary clock is `exp(-itH)` in source action units with hbar=1. The
unit proof and arithmetic are reused from
[source_mass_shell/CALIBRATION_PROOF.md](../source_mass_shell/CALIBRATION_PROOF.md)
and [units.py](../source_mass_shell/units.py). If one source time unit is T
seconds and the declared source signal speed is v, matching it to c gives

```
metres per source coordinate =c T/v,
GeV per source energy =hbar_GeV_s/T.                   (8)
```

An absolute source length or tick duration may be declared for comparison.
`declared_SI_comparison` accepts either, supplies the other by (8), and
returns GeV intervals immediately for any supplied source-energy interval.
For a valid rest energy, the corresponding mass unit is GeV/c^2, conventionally
reported in GeV when c=1. The executable one-metre example is an explicit
unit declaration used only in the separate comparison control, with no
physical default for the primary. It is not a prediction of a metre-sized
source feature or a derived clock calibration.

Specifying c alone fixes the length/time ratio. It leaves the simultaneous
rescaling `T->lambda T`, `length->lambda length`, `GeV/unit->GeV/unit/lambda`
free. This permits a declared comparison scale; it does not establish
agreement with observed masses. Matching the finite quadratic source speed
does not prove the interacting theory's Lorentz speed. The optional
same-generator cesium route
has the separate physical-process and clock premises proved in the reused
unit theorem. Neither route imports a fitted Yukawa coefficient as a
measured primary mass.

The primary report has a proved finite excitation band, without a numerical
full response, momentum-resolved rest interval or physical particle
assignment. Its mass table is empty because those mass-spectrum results
are absent, independently of the permitted unit declaration. Supplied
spectral provenance is recorded and
bound to the primary member; string identities do not prove the supplied
measurement. Architectural inputs and empirical calibration remain visible.
