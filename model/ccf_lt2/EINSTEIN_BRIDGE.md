# Einstein and Newton implications for the common source

The complete CCF-LT-2 Hamiltonian supplies native matter, field and independent per-vertex
energy response. Its Einstein attachment uses the established
[matched regional theorem](../matched_einstein_limit/PROOF.md) and
[quantitative Lean theorem](../../Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean).
The mathematical implications apply unchanged once their inputs are
constructed from this source. Finite lapse response alone does not
instantiate those inputs.

## Exact tensor theorem

Use one symmetric geometry tensor G, one symmetric stress tensor T and ten
regional entropy and area readouts in a common orthonormal frame. For the
ten fixed unit timelike probes u_r and C_ell=4 pi ell^4/15, require

```
|s_r - 2 pi C_ell T(u_r,u_r)| <= e_bulk,
|a_r + C_ell (G(u_r,u_r) - Lambda)| <= e_area,
|s_r + eta a_r| <= e_stationarity.
```

Here eta and ell are positive, the three errors are nonnegative, and the
same s_r and a_r enter all three estimates. The theorem gives

```
max_mu,nu |G_mu,nu + Lambda g_mu,nu - (2 pi/eta) T_mu,nu|
 <= 255 (e_bulk + eta e_area + e_stationarity)/(32 pi eta ell^4).
```

The exact decoder norm is 17/8. Ten probes determine all ten independent
symmetric components; deleting any probe loses that determination. This
route derives tensor balance without a separate null-balance assumption.
Its [executable tests](../matched_einstein_limit/test_frame.py) verify the
rational decoder, extremal norm, each deleted-probe control and the matched
residual identity. These algebraic tests are not a source refinement.

For one positive proper-radius sequence ell_n tending to zero, require
all three errors divided by ell_n^4 to vanish, eta_n to converge to a
common eta_*>0, and tensor and vacuum-reference readouts to converge in
common frames. Identify the geometry limit with G[g], for example through
local C2 convergence to a nondegenerate Lorentz metric. Then

```
G[g] + Lambda_* g = 8 pi G_N T,       G_N = 1/(4 eta_*).
```

Constant eta_*, stress Ward conservation, contracted Bianchi conservation,
metric compatibility and connectedness imply constant Lambda_*. The
vacuum reference fixes its value. The theorem does not select physical
units or a cosmological constant. First variations around one vacuum only
yield the linearized relation; nonlinear dynamics requires finite contrasts
or a controlled integration along a common state family.

The [composition theorem](../../Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean)
also provides a valid conditional interface. Its continuum premise takes
null balance as an argument, which the matched regional route can supply.
Its convention G=kappa T+referenceLambda g uses the negative of the
physical Lambda above.

## Source attachments

| Required input | What the finite construction establishes; required attachment |
| --- | --- |
| Common source and quotient | Full native polynomial and local-capacity columns preserve Gauss. Supply the finite presentation, quotient, boundary and any repair service with its exact descent and boundary properties. Free unitary evolution is not a repair iteration. |
| Common metric and readouts | Graph coordinates are supplied chart data. Attach geometry, event, stress, entropy, scale and modular readouts to one quotient and compatible refinement maps. The scalar smooth finite-menu theorem keeps its own source and observation domain. |
| Symmetric Lorentz stress | All matter and field energy shards and local geometry energy enter the total response. The constant-coframe family supplies six strain and three shift derivatives of H_S. Identify their physical Lorentz stress, spatially varying geometry, boundary and apparatus terms. |
| Regional modular entropy | The adopted cap6 Gauss-zero Gibbs family at beta=1/1024 has exact support ranks and positive floors on all geometry factors; supported regional logarithms exist. Establish the physical vacuum, geometric 2 pi kernel and central/edge normalization. |
| Generalized-entropy variations | Construct the allowed fixed-proper-volume entropy/area variations and stationarity error. The actual 9/32 off-diagonal entry gives a 9/16384 countervariation to unrestricted full-Gibbs/diagonal-area stationarity. A physically justified regional variation family remains required. |
| Common positive area coefficient | Derive the same eta_n from regional entropy/area response and a positive spacetime-independent limit. Supplied lapse stiffness or fitted temperature is not this derivation. |
| Ten-probe estimates | Calculate common s_r, a_r, G and T with independent bulk, area and stationarity errors. Defining a tensor to cancel the residual supplies no estimate. |
| One physical refinement | Specify angular/radial, field-cutoff, matter, reference and clock maps together and prove o(ell^4) errors. The angular/radial shell has shrinking cells and exact chain maps; state, action, cutoff, instrument and proper-metric compatibility remain unproved. |
| Conservation and scale | Establish Lorentz-stress Ward conservation and actual geometric Bianchi identities on the same limit. Calibrate a common clock/length scale independently of the tensor fit and retain source scalars and references. |

The local geometry coupling removes independent neutral-stress
factorization. A complete small diagnostic exhibits distant charged-clock
response to neutral energy. Literal full-native columns verify that the
same coupling enters the master source. Propagated master-source detector
records, regional entropy estimates and the physical refinement require
calculation. Neither the bounded scalar coordinates nor their absence as
microscopic tensor variables proves or refutes an emergent tensor limit.

The complete partition is K_v=A_v+M_v[g_v+N_v(h_v−e_v)+e_v−F]+F,
with one independent capacity factor at each vertex and one common numeric F.
The exact energy currents J_vw=i[K_v,K_w] satisfy dK_v/dtau=−sum_w J_vw.
Their regional and block-map identities instantiate the existing finite Ward
theorem. They are not covariant conservation of an identified Lorentz stress.
The [local bridge proofs and checks](bridges/README.md) give complete evidence
and the [stationarity countervariation](bridges/PROOF.md). Retaining a purifier
does not identify source-marginal, regional and joint pure-state entropies.

## Newton transfer from the same residual

The [gravity transfer](../openwave_gravity_capacity/PROOF.md) uses the same
Einstein component error epsilon_E. Supply a static spherical weak field,
pressureless matter and the same G_N, with the cosmological background
negligible or consistently subtracted. At first order,

```
Laplacian(Phi) = 4 pi G_N rho + (E00+E11+E22+E33)/2.
```

If every component is bounded by epsilon_E, the scalar source error is at
most 2 epsilon_E. If it vanishes outside radius R, regularity at the origin
and a vanishing potential at infinity give, for r>=R,

```
|delta Phi'(r)| <= 2 epsilon_E R^3/(3 r^2).
```

Finite-box boundaries, pressure, anisotropy and nonlinear remainders add
their own errors. The [Lean scalar bound](../openwave_gravity_capacity/Transfer.lean)
formalizes the four-component inequality; the analytic weak-field and
radial integration steps are in the cited proof and executable tests.
Particle acceleration, light paths and separated clock comparison require
the same physical propagation and instrument dictionary.

The continuous commuting-source square-completion comparison
`-g^2 T^T L_D^-1 T/(2 kappa)` does not exactly eliminate the bounded,
noncommuting quantum geometry. Its two-site inverse is not a
three-dimensional Newton kernel. The finite capacity potential selects no
physical cosmic expansion or vacuum amplitude.
