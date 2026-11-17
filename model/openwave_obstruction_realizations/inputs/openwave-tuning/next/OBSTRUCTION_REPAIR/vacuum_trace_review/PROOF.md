# Independent full canonical trace and vacuum consumers

The consumed source is the frozen four-cell, two-wall compact U1 ring,
with all sixteen CAR modes and the infinite electric rotor. This review
imports no producer code. It independently rebuilds its bands, complete
canonical trace, rotor tail and quadrature calculation. Coulomb pointwise
ground bounds, other-number-sector bounds, preparation residuals and the
earlier charged-current residual bounds remain consumed inputs.

## A different complete trace calculation

At zero input mass, momenta at holonomy theta are
k_a=(theta+2 pi a)/4, a=0,1,2,3. Put W=1-cos(k), and

    s_small = 2 W^2 / (sqrt(1+4 W^2)+1),
    s_large = (sqrt(1+4 W^2)+1)/2.

These are the two singular values of the complete two-wall mass matrix
W I-J. The four bands per momentum are

    +/-8 sqrt(sin(k)^2+s_small^2),
    +/-8 sqrt(sin(k)^2+s_large^2).

No band is discarded. The stable expression for s_small also covers W=0.
`run.py` evaluates these formulas with 45-digit interval arithmetic at
every one of 1024 holonomies and compares with the producer's sixteen
stored eigenvalues. Its trace calculation instead uses the coefficient

    [z^8] product_(i=1)^16 (1+z exp[-beta(e_i-anchor/8)]).

Descending polynomial recursion retains every one of the 12,870 Fock
occupations, including degeneracies. It is an independent algorithm from
the producer's explicit subset sums. The eight-particle anchor cancels
exactly in the later eigen-count expression; no energy counterterm or
physical parameter is changed.

The positive Coulomb variance B implies lambda_j(H_free+B)>=lambda_j(H_free)
for all ordered Fock eigenvalues. The independently reviewed Coulomb input
strengthens the ground lower bound only. Therefore this calculation
subtracts one free ground Boltzmann contribution and replaces it with the
stronger ground bound. All excited levels remain. Interval arithmetic
controls the analytic bands, exponentials and polynomial sums; it does
not independently certify the consumed Coulomb lower bounds.

## Complete infinite-rotor trace and the quadrature

Write H=T+V(theta), with T=(E+qbar)^2 and bounded Hermitian V=H_matter+B.
The finite-dimensional [Golden–Thompson inequality](https://arxiv.org/html/1408.2008)
holds without assuming that the two terms commute. To extend it here,
compress to increasing Fourier subspaces of T. The variational eigenvalues
converge monotonically from above because their union is a form core and
V is bounded. The finite heat traces increase to the full heat trace.
On the right side, the T semigroups converge in trace norm and the
uniformly bounded V exponentials converge strongly. Their product trace
therefore converges, proving the required infinite-rotor inequality.

The diagonal kernel of each occupation-dependent shifted rotor is
(2 pi)^-1 sum_l exp[-beta(l+qbar)^2]. Gaussian Poisson summation expresses
it as a cosine series with positive Fourier coefficients. Its maximum
over all real shifts is at zero. Thus the interacting trace is bounded
by the unshifted rotor trace times the holonomy average of the full Fock
potential trace, even when the two Fock matrices do not commute.

For the rotor tail after |l|=16 this review uses a separate geometric
bound: for l=17+k, l^2>=289+35k. Consequently

    tail <= 2 exp(-289 beta)/(1-exp(-35 beta)).

The pointwise trace is integrated with a global eigenvalue Lipschitz
bound. The producer uses 16 sqrt(2). This review additionally establishes
the sharper constant16 on the actual source. The one-body boundary
derivative J=i(C-C dagger) obeys exactly

    J dagger=J,   J^3=64J,   tr J=0,   tr J^2=256.

Hence its spectrum is +8 twice, -8 twice and zero twelve times. The
complete N=8 second-quantized derivative has extrema +/-16, attained
by occupying the two like-signed modes and six zero modes. A cell phase
unitary gives the same derivative spectrum at every holonomy. B is
independent of holonomy. The actual interacting potential therefore has
derivative norm16, and the nearest-grid quadrature factor is

    exp(beta 16 pi/1024).

The code reports both the original and sharper factors. It bounds the
actual potential, rather than differentiating a possibly nonsmooth
pointwise comparison state.

## Eigen count and two preparations of the same vacuum

If exp(beta c) Tr exp(-beta H)<2, there is at most one eigenvalue at or
below c, including multiplicity. A variational state below c and the
other-number-sector lower bounds establish one unique full-source vacuum.
The independent initial run checks the producer's beta=5/4, offset1/8
witness for all four preparation cases.

`consumers.py` adds the rational verification offset1/6 using the sharper
quadrature bound. Increasing a cut only multiplies the already enclosed
count by exp(beta times the cut increase). All four new cuts must still
pass the below-two test and the other-sector exclusions before any
consumer is calculated. This offset is a proof parameter, not a source
coupling or a measured mass choice.

For each mass member, the K4 residual gives a ground interval [L,U] and
a shared excited-state cut c=U+1/6. Its neutral spectral gap is at least
1/6. The complete K2 and K4 preparations are two vectors of this same
infinite-rotor Hamiltonian, not different sources. For a preparation with
Rayleigh upper mu and full residual sigma, two valid weight bounds are

    p_excited <= sigma^2/(c-mu)^2,
    p_excited <= (mu-L)/(c-L).

For the second inequality, use the spectral energy expectation and the
fact that (mu-E0)/(c-E0) decreases as E0 increases when mu<c. It can be
much stronger for the coarser K2 preparation even though its direct
residual bound was zero as a fidelity claim. That older failed estimate
is retained in the producer's freeze.

The bounded CAR current, charged resolvent identity and trace distance
give the same state-transfer error 2 sqrt(p_excited)/g and energy-transfer
error (mu-L)/g^2 as the producer's proof. This review recomputes them with
outward interval arithmetic and both charged-sector lower bounds. The
recorded prepared-current error is retained as a consumed contribution.

This proves a finite-volume vacuum and response for the named one-spin,
one-dimensional control source within the explicit input arithmetic
admissions. It identifies neither a three-dimensional lepton pole nor
its physical mass ratio, continuum chiral protection or a shared gravity
clock. None of those missing transfers is inferred from the trace pass.
