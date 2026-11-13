# Local gauge records on three full Haar refinements

## Source and state

For L=1,2,3 take the open cubic graph on {0,...,L}³. Its physical space is
the Gauss-invariant subspace of the product of L²(SU(2)) over every link.
All representations and intertwiners are retained. On the Casimir domain,

    H_L = alpha sum_e C_e² + beta sum_f(1-W_f),
    alpha=32/L, beta=L³/32, a=1/L, g²=64/L².

The magnetic part is bounded and nonnegative. The electric operator is
self-adjoint and positive; bounded perturbation gives the full source on
the same operator domain. These coefficients and the Haar Hamiltonian are
declared architecture. A native twelve-port update producing this source
is not established by this calculation.

Let Omega be the constant Haar vector. It is a supplied probe, not an
interacting vacuum. Its source energy is beta F, F=3L²(L+1); its preparation
energy above the ground lies in [0,beta F]. The fixed physical loop is the
unit-square outer-face perimeter, subdivided into ell=4L links. Define
A=chi_(1/2)(U_loop)/2 and B=chi_1(U_loop)/3=(4A²-1)/3. Both are Hermitian
gauge-invariant contractions. Haar orthogonality gives means zero and
second moments 1/4 and 1/9. Their ranges are [-1,1] and [-1/3,1].

## Local full-source propagation

Electric evolution preserves the support of any bounded loop operator:
the electric factors are on individual links. Multiplication by a face
Wilson function is bounded. The interaction-picture Duhamel identity
therefore gives, for either O=A or B,

    ||tau_H(t)(O)-tau_E(t)(O)||
      <= min(2, 2 beta N_inc |t|) = delta_L.

Only faces sharing a link with the loop enter the integrand; every other
face commutes with its electric evolution. The exact graph enumeration
gives N_inc=5 for L=1 and N_inc=8L-4 for L>=2. No global low-energy
projection or representation cutoff occurs. This is the local
interaction-picture argument used by the retained full-Haar propagation
theorem, with the free electric evolution as its comparison.

Each link traversed once by a spin-j Wilson loop has Casimir j(j+1).
Consequently A Omega and B Omega are exact electric eigenvectors with
energies alpha ell 3/4=96 and alpha ell 2=256, for every L. Free real
correlations divided by their initial second moment are cos(96t) and
cos(256t). Cauchy-Schwarz bounds the unnormalized interacting error by
delta_L/d, with d=2 or3. The normalized error is at most d delta_L.
Omega is not stationary under H: a normalized correlation need not lie
in [-1,1], and the receipt does not impose that restriction. Its elementary
absolute bound is d. Pairwise normalized refinement differences are at
most d(delta_L+delta_M). Rational Taylor remainders enclose the cosines.

The two common native times 1/8192 and1/4096 are supplied intervention
times. The spin-one observable is a held-out Wilson character using the
same source and times. It is not identified with physical stress. At the
later time, both normalized real correlations exceed0.97 on all three
refinements. The perturbative error grows as L⁴ at fixed time. A time
window of order L^-4 controls this estimate; no fixed-time continuum
convergence follows from it.

## Complete temporal acquisition and work

For either contraction O choose eta=1/2 and the three Kraus operators

    K_+=(I+eta O)/2, K_-=(I-eta O)/2,
    K_0=sqrt((I-eta² O²)/2).

They are smooth bounded multiplication operators with sum K_s* K_s=I.
They preserve Gauss invariance and the electric form domain. A supplied
blank three-label apparatus realizes their Stinespring isometry. The map
extends by the identity to every external quantum reference. All nine
sequential outcomes (s0,st) are retained, with probabilities

    p(s0,st)=Tr[K_st U(t) K_s0 rho K_s0* U(t)* K_st*].

Their weighted first instrument is
sum_s s K_s rho K_s = eta {O,rho}/2. The second weighted effect is eta O.
Thus, exactly for any state and the full interacting source,

    E[s0 st]=eta² Re Tr[rho O(t) O].

This relation includes the first acquisition's backaction. It requires
no postselection or infinitesimal probe limit. The executable receipt
evaluates all nine zero-time Haar probabilities and a complete
noncommuting two-state control with an entangled reference. At positive
times the full-Haar result is an analytic bound on the weighted moment;
individual infinite-carrier probabilities are defined, not claimed to
have been separately numerically evaluated. Different refinements,
times and characters require separately prepared probes.

For the electric carré du champ Gamma and any finite-form-energy input,

    sum_s Gamma(K_s) = eta² Gamma(O)/[2(1-eta² O²)].

Magnetic multiplication commutes with the Kraus operators. Summing the
complete instrument cancels the derivative cross terms, so its exact
mean source work is alpha times the expectation of this nonnegative
function. For eta=1/2 its coefficient is at most1/6. The inherited Haar
identity Gamma(A)=ell(1-A²)/4 gives work <=alpha ell/24=16/3.
The character chain rule gives Gamma(B)=16ell A²(1-A²)/9 <=4ell/9,
and work <=2alpha ell/27=256/27. Two readouts cost at most32/3 or512/27,
uniformly in the three refinements and in the input state. Free source
evolution contributes zero mean work. Probe preparation, initialized
apparatus records and intervention times are explicit resources; no
free reset or autonomous physical clock is asserted.

## Scope

This certificate gives full-Haar local memory, complete temporal quantum
records and refinement-independent readout work bounds. The correlation
error has an explicit nonuniform refinement dependence. It supplies no
physical time calibration, continuum Yang-Mills existence, Clay mass gap,
particle rest mass, Nature match or native-process admission certificate.
