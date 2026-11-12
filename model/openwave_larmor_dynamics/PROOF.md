# Spin and dipole precession with a unique quasistatic self-field

## Construction and complete records

The matter operator is the complete positive kernel `K(A)=|D_W(A)|` of
[the stationary magnetic construction](../openwave_intrinsic_magnetic/PROOF.md).
Its exact unitary reduction retains all sixteen native Clifford rows as
four identical four-row blocks, both signed Dirac branches and every
folded valley label. The one-particle preparation occupies one block;
three spectator blocks and the empty Fock vacuum have zero probability.
Their states are present in the number-conserving Fock extension.

Matter has open Dirichlet links on a side-three cube of spacing one.
The real self-field is periodic on a side-nine cube, restricted to
mean-zero transverse vectors. The parameters are mass 8, scalar mass
well 0.01, Wilson coefficient 1/2 and field coupling 0.1. A fixed external
uniform field enters the same Peierls links. Its apparatus energy is not
part of the evolving matter/self-field system. The energy is

    E(psi,A) = <psi, |D_W(A+A_ext)| psi> + a^3 ||C A||^2/(2 kappa).

This is an explicitly quasistatic transverse extension: at each time A
minimizes this energy at fixed normalized psi, and

    i psi_dot = K(A[psi]+A_ext) psi.

It includes neither a propagating Maxwell kinetic field nor longitudinal
electric backreaction. The external scalar well supplies confinement.
It is not an implementation of the interacting Peter-Weyl Hamiltonian.
The state-dependent field makes this a nonlinear semiclassical mean-field
law. It is not a linear completely positive channel on arbitrary
entangled inputs. A prescribed frozen field has the linear second-
quantized kernel on the full Fock space; that reference-uniform channel
property does not transfer to the self-consistent state-dependent law.
Measurement completeness here concerns each stated trajectory and its
separate instruments, without an arbitrary-reference dynamics claim.
A computational one-particle rotating phase removes mass 8 during
propagation; multiplying by exp(-8 i t) restores the physical unitary.
Every reported measurement is insensitive to that common phase.

The complete separate instruments measure the three canonical spin
components, the instantaneous positive-kernel energy and the three
magnetic moment components. They are separate interventions, without an
assertion of a simultaneous spin/energy/moment measurement. All 108
signed spectral modes of the occupied block are included in each energy
instrument; the unprepared spectator blocks have zero probability. The
moment instrument diagonalizes the full derived observable, rather than
assigning a moment from a spin label.

## Current for arbitrary superpositions

The gap `|D_W| >= m` holds for every real field. Differentiating
`K^2=D_W^2` gives

    K X + X K = D_W delta D_W + delta D_W D_W.

In a complete Dirac eigenbasis with signed eigenvalues lambda_i,

    X_ij = (lambda_i+lambda_j)/(|lambda_i|+|lambda_j|) (delta D_W)_ij.

The denominator never vanishes. For density rho, let W be the same
self-adjoint divided-difference map applied to rho. Every actual bond
current is `J_e=-Tr(W partial D_W/partial A_e)`. This includes the Wilson
term, cross-sign coherences and degenerate subspaces. Replacing W by rho
is valid within a positive Dirac spectral subspace; it fails on the
arbitrary signed-state control.

Gauge covariance yields the exact local continuity law
`rho_dot + D*J = 0`, including the open matter boundary. Away from a
stationary state J need not be divergence free. The field equation is
therefore `C*C A=(kappa/a^3) P_T J`, with the mean-zero transverse
projection P_T. The harmonic and longitudinal currents are not relabelled
as transverse field sources. The moment is `mu=(1/2) sum_e r_e cross J_e`;
it equals the expectation of `-partial K/partial B` on the same state.

## A global finite contraction

All norms on the field coordinates use the counting edge norm. Uniformly
in A and external B, the matter matrix satisfies

    ||D_W|| <= M = m + eta max|x|^2 + 6r/a + 3/a,
    ||D_W'[h]||^2 <= 2(1+r^2)||h||^2,
    ||D_W''[h,h]|| <= a(1+r)||h||^2/2.

The first derivative bound follows from the four-row Frobenius norm of
each directed bond and its adjoint. The second follows from the block
row-sum bound. Direct sums of the four identical native blocks have the
same operator bounds. The positive Sylvester inverse has norm at most
1/(2m), so differentiating K twice gives

    ||K''[h,h]|| <= C ||h||^2,
    C = [(1+(M/m)^2)2(1+r^2) + M a(1+r)/2]/m.

Here `M=1403/100` and `C=5302169/2048000`. The nonzero transverse
Laplacian gap is `4 sin^2(pi/9)`. From `3<pi<22/7` and
`sin x >= x-x^3/6`, a rational lower bound is
`239557386916/562711519881`. Thus

    q = kappa C/(a^3 lambda_lower) < 0.608135 < 1.

The map `A -> (kappa/a^3)(C*C)^+ P_T J(psi,A)` is a global contraction
for every normalized psi. Its unique fixed point is the unique global
field minimizer: the field-energy Hessian dominates the matter Hessian.
For any trial A the exact fixed-point residual r gives
`||A-A_exact|| <= ||r||/(1-q)`. Displayed residuals are floating
calculations, without an outward-rounded numerical enclosure claim.

Every fixed-state current solve has this contraction. The implicit
exponential midpoint solve also changes its unitary half-step state with
the field; a separate global contraction bound for that coupled midpoint
map is not asserted. Its executed convergence and residuals are checked.
The exact continuous quasistatic dynamics conserves norm and total
energy by stationarity of the field energy.

## Time evolution and comparisons

Each numerical step uses its own self-consistent midpoint current and
an exact dense spectral exponential of the complete K. Step sizes 20,
10 and 5 cover a full cycle of duration 2200 at B=0.05. The half-field
experiment doubles the duration. Reversing B reverses the rotation; zero
B preserves the initial stationary spin. Spin and the independently
computed intrinsic dipole precess at matching frequencies. Complete
energy and moment records are printed at both endpoints.

The full-cycle frequency is compared with the independently reoptimized
stationary energy splitting at the same external field. The frozen
initial field is a distinct control. Its exact linear spectral evolution
is compared with a direct matrix exponential; its spectral residual gives
the usual finite-time unitary error formula `t ||H U-U Lambda||` plus
basis-unitarity error. Freezing the field leaves a nonzero self-field
residual during the rotation and does not solve the coupled construction.

A split doublet precesses at its energy difference, approximately
`2 |mu| B`. The projected canonical spin has length s below 1/2.
Consequently the static definition `g=2 E |mu|/s` differs from the
doublet-precession definition `g=2 E Omega/B` by the factor 2s in the
weak-field limit. Both conventions and both supplied-mass and measured
bound-energy references are retained. The finite lattice supplies no
continuous spatial rotation generator that removes this distinction.

Time-step differences, energy drift and residuals are finite numerical
controls. No useful interval bound for the complete long-time nonlinear
trajectory, continuum g=2 theorem, physical electron identification or
anomalous moment prediction is asserted.
