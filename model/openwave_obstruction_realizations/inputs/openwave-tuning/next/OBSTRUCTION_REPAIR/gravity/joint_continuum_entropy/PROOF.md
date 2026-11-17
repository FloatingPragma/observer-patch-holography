# Actual common-region entropy equations with independent coframe and lapse

The massless continuous source is D_sigma=e^(-sigma/2) D e^(-sigma/2), H_uv=D_sigma+beta R^3 D_sigma^4, and H_flat=H_uv/sqrt(1+R^6 D_sigma^6). D_sigma anticommutes with beta, so the positive spectral denominator commutes with H_uv and its projector is unchanged. The lapse is introduced independently by H_N={N,H_flat}/2, N=exp(n). This is a declared self-adjoint static source coupling. Its low-momentum part has the standard lapse/coframe Dirac form, while a full curved-spacetime covariance theorem for its nonlocal ultraviolet part remains unproved.

At the flat reference, write A=alpha dot (p+q), B=alpha dot p, E_uv(p)=sqrt(|p|^2+R^6 |p|^8), T(p)=H_uv(p)/E_uv(p). For a unit Fourier coframe mode,

delta D=-(A+B)/2,
delta(D^2)=-[(|p+q|^2+|p|^2)I+2AB]/2,
delta H_uv=delta D+beta R^3 (|p+q|^2+|p|^2) delta(D^2).

The complete negative-projector derivative is

delta P_sigma(p+q,p)=-[delta H_uv-T(p+q)delta H_uv T(p)]/[2(E_uv(p+q)+E_uv(p))].

For the independently varied lapse, the actual flattened energies E_flat=|p| enter:

delta P_n(p+q,p) = [E_flat(p+q)-E_flat(p)] [P(p+q)-P(p)]
                    /[2(E_flat(p+q)+E_flat(p))].

Uniform lapse gives zero exactly. The code integrates these vertices with the actual cell form factors w_h(p+q)w_h(p), retaining successive continuum alias shells. Uniform coframe is independently checked against the direct derivative of the source symbol. A small complete finite Galerkin source supplies an additional actual nonlinear entropy check of both vertices; its complement convention is stated separately from the production continuum approximation.

For each regional cell algebra, all G1/G2/H multiplicities give the actual covariance C_B and modular matrix K_B=log[(I-C_B)/C_B]. The derivative is Tr(K_B delta C_B). Reconstruction/aggregation of the complete modular matrix is checked against a direct full regional matrix before larger balls are used. No sector is omitted for geometric or lapse response.

The proper-volume constraint uses the genuine moving-cell pullback proved in the preceding packet. The coefficient is average_union sigma, computed from the integrated cell form factor and the cell centres. Thus the geometric row is the independently measured coframe response plus this geometric average times the independently measured unit inward shape response. The residual uniform row remains recorded and is never subtracted.

Matter is the exact neutral four-particle/four-hole Slater source from normalized continuum positive/negative projected radial seeds. Its full source Gram matrices and flattened-clock energies are measured. Projecting those real source modes into the observer cells gives the actual regional covariance update. All representation components are included for displaced centres. At the source centre, cubic symmetry also allows the actual finite Gaussian entropy difference to be evaluated and compared with its tangent. Changing packet width changes its actual state and energy; no arbitrary excitation-strength parameter is substituted.

The finite preparation ansatz uses three cubic coframe modes and three independent cubic lapse modes, each the average of cosines along the three axes. Every basis function, every ball size and every training/held-out designation is explicit. The linear preparation law solves the six actual entropy rows from two radii and three centres, then applies that same single field to an unused radius and other centres. It is an entropy-based preparation rule; neither a Poisson equation nor a relation between lapse and coframe is an input. A small training residual alone has no physical acceptance force.

The translation-averaged matter entropy row is measured directly from the same state. All selected geometric modes have zero translation mean, and a true uniform fixed-proper-volume coframe/lapse row vanishes in the infinite-volume coordinate limit. Thus a positive mean matter row remains an explicit obstruction to a globally static scalar equilibrium on the torus. It is not absorbed into a finite-size numerical K(0). An eventual open/infinite construction or homogeneous temporal source response must handle it.

The normalized cell algebra also traces subcell modes, so its global observer state can have a bulk entropy density. That density is measured and reported without subtraction. The parent fixed-region entropy convergence theorem applies to nested scalar cell refinements, but supplies no numerical convergence rate or derivative/large-ball interchange. These observer, finite-box, weak-source and shape errors must be controlled before the tangent equations can be interpreted as a gravitational equation.
