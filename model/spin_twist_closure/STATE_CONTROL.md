# Selected-ground regional energy and state information

Use the same source and the same volume-then-selection-field ground
prescription as `PROOF.md`. Its finite-ground origin gives a regional
second-moment bound beyond the bound obtained from local stability alone.

## The exact finite ground equation controls a regional second moment

Split a finite source containing a region Lambda as

\[
 H=H_\Lambda+H_{\Lambda^c}+V,
 \qquad -JbI\le V\le0,
\]

where b counts the crossing bonds. Onsite staggered-field terms belong
entirely to one regional Hamiltonian. Let E, E_Lambda and E_complement
be the three ground energies, and set
`Delta_Lambda=H_Lambda-E_Lambda`,
`Delta_complement=H_complement-E_complement`. These are commuting
positive operators. Variational bounds on the negative crossing
interaction give

\[
 -Jb\le c:=E-E_\Lambda-E_{\Lambda^c}\le0.
\]

For a true finite ground vector Omega, the eigen-equation gives

\[
 (\Delta_\Lambda+\Delta_{\Lambda^c})\Omega=(cI-V)\Omega,
 \quad \|(cI-V)\Omega\|\le Jb.
\]

Commuting positivity implies
`Delta_Lambda^2<=(Delta_Lambda+Delta_complement)^2`. Consequently

\[
       \omega(\Delta_\Lambda^2)\le J^2b^2.               \tag{1}
\]

The same assertion holds for ground mixtures by linearity. At fixed
Lambda and positive selection field this is a uniform bound in the
ambient volume on a fixed local operator. It therefore passes to the
selected infinite-volume state. Removing the field next is legitimate:
the regional Hamiltonian changes in norm by at most `B|Lambda|/2`, and
its minimum eigenvalue changes by no more than that norm. Thus (1)
holds in exactly the specified zero-field ordered state.

For an open cube of side L, b=6L^2 and |Lambda|=L^3, giving

\[
 \left\|{H_\Lambda-E_\Lambda\over|\Lambda|}\Omega\right\|^2
 \le {36J^2\over L^2}.                                   \tag{2}
\]

The local-reset bound in `PROOF.md` is also valid, so the producer
retains the minimum of `18J^2(L-1)/L^2` and (2). Combining (2) with
the same exact tiling bracket
`E_Lambda/|Lambda|-3J/L<=e0<=E_Lambda/|Lambda|` yields

\[
 \|(H_\Lambda/|\Lambda|-e_0)\Omega\|\le9J/L.             \tag{3}
\]

Including the forward crossing bonds costs at most another `3J/L`.
Thus the scalar total bond-energy coefficient has an explicit
microscopic ground-vector convergence rate. No rate for general
finite-volume marginal states is obtained from (3).

Equations (1) and the first-moment reset estimate also give a complete
spectral tail bound. For an excitation threshold u>0 above the finite
regional ground,

\[
 \omega(1_{[u,\infty)}(H_\Lambda-E_\Lambda))
 \le\min\{1,Jb/u,J^2b^2/u^2\}.                            \tag{4}
\]

The sharper local-reset second moment `J^2 M b` may replace `J^2b^2`
where smaller. The effect includes every regional energy outcome and
does not condition on low energy. These are regional energies; they
are not vacuum-relative particle masses.

## Why this does not select the full marginal

On the same zero-field source, rotate the maximally ordered ground by
pi about internal axis two. Every local bond Hamiltonian and every
product of them is unchanged by this symmetry, so the two ground
phases have identical complete local-energy statistics. Their one-site
states, however, are

\[
 \rho_\pm={I\pm2m\epsilon_x\sigma^1\over2},\qquad
 {1\over2}\|\rho_+-\rho_-\|_1=2m>{11\over50}.             \tag{5}
\]

The strict bound uses the same replayed source order certificate
`m>11/100`. A local spin pointer distinguishes these states. Only the
plus-field prescription is the chosen state; the rotated phase is a
countercontrol to deriving that prescription from energy data. In
particular, identical energy densities, tails and all energy moments
cannot by themselves give a trace-distance bound to one chosen phase.

A genuine gap g>0 above a finite Hamiltonian's ground subspace gives
`Tr((I-P_ground)rho)<=Tr(rho(H-E0))/g`. If that ground is unique, the
half trace distance to its pure state is at most the square root of
this ratio. The regional Hamiltonian in (4) is not the global
interacting reference Hamiltonian; its ground projector does not in
general coincide with the desired regional marginal. Inserting the
actual finite Doeblin lower gap into (4) gives only the bound one on
the demonstrated cubes. The code evaluates this with exact rational
source parameters. This diagnoses the bound, not the true gap.

Lieb--Robinson propagation can transport a known initial local-state
error through a finite experiment; it does not bound that initial
error. For n regional spins, a sufficient finite observable criterion
is explicit control of all Pauli expectation differences d_P. Pauli
orthogonality and `||X||_1<=sqrt(2^n)||X||_2` give

\[
 {1\over2}\|\rho-\sigma\|_1
   \le {1\over2}\sqrt{\sum_{P\ne I}|d_P|^2}.              \tag{6}
\]

The finite native compiler provides these complete local Pauli
instruments. A uniform bound on their difference from the selected
infinite state at the required growing supports is an additional
state estimate. Energy-only inequalities (1)--(4) do not supply it.

## The remaining force estimate uses a different source average

For the actual uniform field `P_p(G)=-p sum_x G(px)S_x^2`, the exact
current identity in `spin_shell_diagnostic/README.md` gives

\[
 i[L_p,P_p(G)]
   =\sum_{x,i}[G(px)-G(p(x+e_i))]C_{x,i}^{,2}
   =-p^{-2}\,p^3\sum_{x,i}D_i^pG(px)C_{x,i}^{,2},       \tag{7}
\]

where `D_i^pG=(G(p(x+e_i))-G(px))/p` and
`C_{x,i}^2=(S_x cross S_(x+e_i))^2`. Under the same beta translation,
the axis-two current has no staggered sign. The zero-integral weighted
mean-ergodic argument therefore makes the unamplified average in (7)
vanish on the ground vector, but it supplies no rate that survives the
factor `p^-2`. The rate (3) concerns a different, total-energy operator.
It cannot be substituted for a current correlation estimate.

Closing a field equation requires controlling (7) after subtracting
the proposed field force, on the same energy-filtered ground vectors
and histories. The corresponding first-order current residual in
`spin_shell_diagnostic` must have squared norm `o(J^2 p^3)` for its
specified packet normalization. The scalar energy density and its
rate prove neither that cancellation nor saturation of the canonical
covariance inequality. The field-generator identity is thus retained
with its actual current, without replacing it by a free-field force.
