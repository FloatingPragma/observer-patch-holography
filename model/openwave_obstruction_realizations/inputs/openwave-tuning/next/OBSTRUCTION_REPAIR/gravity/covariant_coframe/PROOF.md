# A positive connection completion and a geometric CAR pullback

## Fixed source and the missing connection

Let D=Σ alpha_j sin(k_j), L=2Σ(1−cos(k_j)), and beta anticommute with D.
The flat source used in the earlier packets is H=D+beta R³L². Define

    W=L−D²=Σ_j(1−cos(k_j))² ≥0,
    v=exp(−sigma),
    D_sigma={v,D}/2,
    Q_sigma=D_sigma²+vWv,
    H_sigma=D_sigma+beta R³Q_sigma².

These are actual finite lattice operators. Q_sigma is positive, local,
self-adjoint and commutes with beta. H_sigma is self-adjoint. At sigma=0,
Q=L exactly, so the flat covariance, modular Hamiltonians, entropy-area
coefficient and neutral preparation from the prior source are unchanged.

For a smooth continuum conformal spatial metric g=e^(2sigma)delta in three
dimensions, the Dirac operator transformed to the flat half-density Hilbert
space is e^(−sigma/2)D e^(−sigma/2)={v,D}/2. The equality uses the continuum
first-order product rule. Its square contains the spinor connection terms.
The lattice W is order h²p⁴ relative to the continuum Laplacian after restoring
h. Thus Q_sigma completes the covariant Dirac square at resolved momenta,
while retaining a positive lattice term that removes the flat high-momentum
doublers. Exact covariance at finite lattice spacing has not been proved.

For one Fourier coframe mode, with k'=k+q and A(k)=Σ alpha sin(k),

    delta D(k',k)=−[A(k')+A(k)]/2,
    delta Q(k',k)=−[L(k')+L(k)] I
                  +[A(k')²+A(k)²]/2−A(k')A(k),
    delta H=delta D+beta R³[L(k')+L(k)] delta Q.

The last two terms in delta Q are the omitted connection. In the continuum
they equal q² I/2 minus the spin matrix from the antisymmetric product of
q and k. They vanish identically for q=0; no uniform-coframe calibration or
Einstein coefficient is used. `dense_control.py` differentiates the complete
finite operators independently and checks finite entropy differences.

## Deterministic moving observer algebra

For a declared spatial vector field xi, let

    A_xi=(xi·partial+partial·xi)/2,  A_xi^dagger=−A_xi,
    U_epsilon=exp(epsilon A_xi).

On a finite lattice, symmetric derivative stencils define an exactly
anti-Hermitian matrix. Its second quantization preserves particle number.
The observer modes b=P_B U_epsilon c form an exact CAR factor of the global
finite Fock space, with covariance P_B U_epsilon C U_epsilon^dagger P_B.
No stochastic region selector, Shannon entropy, or per-ball entropy offset
is introduced.

For g=e^(2 epsilon sigma)delta, the first-order proper-volume constraint of
a continuum ball B fixes a uniform inward coordinate dilation

    xi(x)=−c_B x,  c_B=mean_B(sigma),
    delta V=3 integral_B sigma + integral_boundary xi·n=0.

xi may taper outside a stencil neighbourhood of B; the first-order regional
response depends only on the values within that neighbourhood. Each ball
has its own observer map applied to the SAME global source metric. For the
axial q mode, c_B=W_B(q)/V_B. The script computes the actual entropy derivative
from the commutator [A_xi,C], rather than substituting a fitted shape slope.

For a stencil partial f=Σ_m d_m[f(x+m)−f(x−m)] and a unit inward dilation,
the translation invariant covariance obeys

    delta C(x,y)=−Σ_(j,m) d_m[(x_j−y_j+m)C(x−y+m e_j)
                           −(x_j−y_j−m)C(x−y−m e_j)].

`run.py` evaluates this expression on the true regional covariance and
modular Hamiltonian. The dense control exponentiates the finite matrix and
checks the sign convention using H_coord=U^dagger H U: then the moving
observer covariance is exactly unchanged. This checks the implementation,
not equivalence of the source coframe deformation to a continuum coordinate
transformation.

## Explicit limits

A fixed-rank finite CAR factor does not by itself define literal geometric
volume. Its projector has constant trace even when its mode support moves.
The geometric volume in the construction is the ball boundary volume above;
identifying the moving mode factor with all local field observables in that
ball requires localization and refinement estimates. The finite unitary has
tails beyond the generator stencil at higher orders. Uniform-sigma
cancellation is therefore an essential measured diagnostic, not an axiom.

This branch has lapse N=1. It supplies neither a lapse equation, temporal or
boosted stationarity, nor a physical Newton coefficient. A positive-energy
source on a periodic static scalar background retains the mean-mode
compatibility problem. There is no autonomous event selection law here.

The beta R³Q² term breaks the naive axial symmetry that commutes with the
kinetic alpha matrices. A free low-momentum cone is not a proof of radiative
mass protection. No interacting chiral Ward identity, overlap/Ginsparg-Wilson
completion or counterterm ledger is established by this packet. The source
has exact total fermion number, but no joined fractional-colour gauge service
or transfer to the electron/clock unit map has been supplied.
