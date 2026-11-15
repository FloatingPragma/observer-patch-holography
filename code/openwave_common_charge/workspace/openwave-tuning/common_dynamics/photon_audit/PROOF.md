# Independent quadratic photon and frozen-data audit

The source is the explicit noncompact quadratic amendment on the actual
carrier complex. The audited producer manifest is
`43fca237a676954e4a213b0888e2cdfee70bb9d8cefd8fbf113f5d2a60cbc81e`.
Its 12 files and 336 source pins are verified before and after execution.
No producer entry point that writes receipts is invoked. The audit owns
only this directory and preserves the separate frozen quantum-reference
audit.

The executable checks continuous-field identities, the exact rational
spectral certificate, all-direction frequency coefficients, the native
spectral experiment and explicit mutations. It does not certify a local
edge encoder, moving charged matter, native harmonic packets or an
official OpenWave row.

## Continuous source and common normalization

The integer matrices have dimensions B:96 by 288 and A:288 by 208,
with BA=0, rank B=95 and rank A=190. The physical cycle space has 193
dimensions. Its orthogonal complement to the 190 coexact modes consists
of three harmonic canonical pairs.

For a neutral real charge vector q=G-n, put
L=BB^T, G_L=(L+11^T/96)^(-1)-11^T/96 and P_q=B^T G_L q.
Direct matrix checks give BP_q=q and A^T P_q=0. The projections
P_L=B^T G_L B and P_T=I-P_L decompose each actual routed current J.
Thus P_(q+BJ)-P_q=P_L J and the transverse electric update is P_T J.
Every one of the 48 size-two routes is checked. Independent fields per
particle are outside this construction.

The Hamiltonian coefficients a=1/2 and b=3 are supplied architecture.
With rho=2/3, rho a b=1 and rho a=1/3. The independent cell calculation
uses the 12 by 36 zero-momentum incidence matrix. For the matrix C of
the six positive inter-cell edges and P_cycle=I-B_0^+ B_0, it obtains
C^T P_cycle C=(3/2)I. The harmonic flux frame
W=(2/3)P_harmonic C consequently satisfies C^T W=I and
W^T W=(2/3)I. This construction of W uses the cell potential corrector,
independently of its extraction from macroplaquette derivatives.

The conductivity 3/2 divided by the 12 vertices per cell gives acoustic
Laplacian coefficient 1/8. A specified root has squared projection 1/12
on the normalized constant mode. Its Green pole is therefore
(1/12)/(1/8)=2/3, and multiplication by a gives charge coefficient 1/3.
These are cell and leading-symbol identities. A quantitative real-space
remainder and a measured force require separate transfers.

All 190 positive eigenvalues enter the free Gaussian covariance. The
audit reconstructs its coexact factors and adds the three harmonic
packets with electric variance 1/2. Their product obeys
Cov(theta)Cov(E)=P_T/4. The harmonic coordinate variance evolves from
1/2 to 0.52 at time 0.4 under a=1/2; it cannot be assigned a stationary
normalizable ground state. The four reported charge energies and all
three shared-field route probabilities are independently reproduced.
Those route probabilities are logical continuous-field comparisons;
the eight-bit native photon program does not implement them.

## Exact frozen spectral certificate

The immutable data contain a 288 by 190 integer frame Q_num and 190
rational eigenvalue approximants, each with denominator 2^44. The
audit does not rerun an eigensolver to select these integers. Instead it
forms their exact Gram and curvature residuals. Three base-2^15 limbs
replace the producer's two base-2^22 limbs. Each int64 product has an
explicit absolute overflow bound, and arbitrary-precision integers
recombine the products. This establishes the exact squared Frobenius
norms of Q^TQ-I and KQ-Q Lambda for K=AA^T.

Let epsilon and r bound these operator norms through their Frobenius
norms, and m=min Lambda. With s=r/m, define

    Y=K Q Lambda^(-1),
    gamma=epsilon+2(1+epsilon)s+s^2,
    U=Y(Y^T Y)^(-1/2).

Then ||Y-Q||<=s, ||Y^T Y-I||<=gamma<1/2 and
||U-Y||<=gamma. The columns of Y are in im K. Their Gram lower bound
makes them independent, so U is an exact orthonormal frame of the full
190-dimensional coexact space. In particular BU=0 exactly. Rounded Q
itself is not asserted to satisfy exact Gauss support.

Using ||K||<=1504/17 and l=max Lambda,

    ||U^T K U-Lambda|| <= r+(1504/17+l)(s+gamma)=eta.

The audit reconstructs every rational coefficient in this expression.
The certificate gives eta<m/2. For positive matrices, the square-root
perturbation estimate follows from the Sylvester equation
sqrt(A)X+X sqrt(B)=A-B and its exponentially convergent integral.
It yields

    ||sqrt((3/2) U^T K U)-sqrt((3/2)Lambda)||
      <= (5/4) eta/sqrt_lower(m)=delta_omega.

The actual one-photon Hamiltonian in the exact frame U need not be
diagonal. The native diagonal approximant acts in that same specified
frame. The vacuum energy error is at most 190 delta_omega/2, and cap
one adds at most delta_omega. The full bound, including the vacuum
scalar, is 96 delta_omega, approximately 2.2783059109 times 10^-5.
The 65 padding states retain their declared zero Hamiltonian.

A coefficient mutation violates the frozen Gram bound. A different
control multiplies one selected eigenvector by -1: it preserves the
coexact projector exactly and changes the integer frame bytes. This
demonstrates why valid eigensolver output does not identify a unique
immutable encoding.

## Frozen replay interface

The producer's `fock.spectral_certificate()` selects rounded frame
integers through runtime `numpy.linalg.eigh`, and its native phase
generator calls that function. Degenerate eigenspaces and eigenvector
signs permit different valid selections across numerical environments.
This is an execution-data portability issue; it does not invalidate the
frozen integer certificate.

The additive `frozen.py` interface supplies:

- `spectral_certificate()`: the immutable rational frame, spectrum and
  bounds, with no eigensolver call.
- `phase_targets(time)`: outward phase intervals for that fixed spectrum.
  This mathematical target does not certify a native word at another time.
- `photon_phase(Fraction(1))`: the common eight-bit phase word extracted
  from both emitted nine-bit experiments, with a recalculated all-string
  Walsh certificate.
- `word(setting, Fraction(1))`: the emitted X or Y experiment, its exact
  stored bound and its fixed-spectrum phase target.

Native replay rejects a float time and exact times other than one. These
rejections are executed. Runtime recompilation or a different time needs
its own word and certificate. A consumer cannot transfer this default
word by changing its clock label. The adapter verifies all frozen input
hashes, returns copied mutable arrays, and writes no provider output.
`IMPORTS.json` gives its transitive source ancestry for portable export.
The separate facade control disables NumPy's eigenvalue and eigenvector
solvers during all replay APIs. It also mutates returned frame and word
objects and verifies that subsequent calls retain the frozen data.

## All-direction low-frequency theorem

The audit reconstructs the 36 by 26 local face symbol from the actual
size-three geometry and integer translation terms. A rank-factor
pseudoinverse, independent of the producer's `pinv` calculation,
gives rank A(0)=22, ||A(0)||_F^2=88 and
||A(0)^+||_F^2=3383/240<16. Thus the positive singular values at zero
momentum exceed 1/4 and ||A(0)||<10.

With P=I-A_0 A_0^+, Q=I-A_0^+ A_0 and linear term E_1, all six real
quadratic coefficients of (P E_1 Q)(P E_1 Q)^* are checked exactly.
They equal W(|k|^2 I-kk^T)W^T, and every linear term lands in the
three-dimensional harmonic space. Therefore the two positive leading
singular values are sqrt(2/3)|k| in every direction. Axis sampling is
not used as proof of rotational independence.

Integer sums of the actual translation lengths give, for r=|k|,
||E||<=11r=eta_k and ||A(k)-A_0-E_1||<=(17/2)r^2.
The skew-adjoint left and right block eliminators have norm at most
4 eta_k each. The first-order terms involving E cost at most
8 eta_k^2. The second derivative integral remainder costs at most
32 eta_k^2(10+eta_k). Hence

    D(r)=8 eta_k^2+32 eta_k^2(10+eta_k)+(17/2)r^2.

The high singular-value block is bounded below by 1/4-eta_k.
The explicit admissions eta_k<1/8,
r+D(r)<1/4-eta_k and D(r)<r/2 separate it from exactly two positive
low singular values. The incidence symbol has rank 12 for sufficiently
small nonzero k, since in-cell connectivity makes a possible left-null
vector constant and any nonzero inter-cell phase forces that constant
to vanish. BA=0 then bounds the face-symbol rank by 24. Together these
facts exclude additional positive low singular values.

Multiplication by sqrt(ab)<5/4 yields
|omega(k)-|k||<=(5/4)D(|k|). The audit reproduces the exact bounds at
h=2^-18, 2^-20 and 2^-22. These are Bloch-fibre estimates at k=hp,
not native executions of those spatial lattices. Allowed momenta of a
finite periodic box must be admitted separately.

The finite-photon phase comparison uses an explicit spectral dictionary
and the same source vacuum phase on both sides. It neither establishes
absolute continuum vacuum-energy convergence nor supplies a local-edge
field observable transfer. Nonzero creation-operator matrix element in
a normalized spectral mode does not by itself fix the normalization,
smearing, preparation or native pointer cost of an electric field.

## Native phase and complete records

The audit reconstructs the six-dimensional source generators directly
from their response and control receipts. It exponentiates these
matrices and applies each pulse with an independent integer bit-index
gather. It does not call `compiler.execute`, `compiler._native` or the
producer's experiment runner. The emitted diagonal word is tested on
vacuum, selected one-photon states and padding columns. Its full
256-state phase error is checked through every Walsh coefficient and
the pinned scalar control verifier.

The emitted 4,879 and 4,880 pulse words are executed on the nine-bit
prepared state. An independent formula retains the scalar vacuum phase
and gives photon amplitudes
(exp(-iE0)+exp(-i(E0+omega_1)))/2 and their difference divided by two.
The clock state is (cos t, i sin t), with t=-1/16 for X and
t=pi/4-1/16 for Y. No output phase alignment is performed.

Every fine C/E effect is compared with the live native source law.
Their two coarse effects are diagonal in the data basis. Independent
bitwise stochastic contractions therefore produce all 512 joint
outcomes, including the 65 spectral padding states in each clock
sector. The probability sum is compared with the unnormalized source
norm. Export complements and reset failures retain their declared
outcome space; there is no success division. Exact scalar/certificate
inputs and the producer's complete prefix-complement specification
support the failure accounting. The audit does not claim an independent
enumeration of every microscopic failure-prefix tape.

The audit reconstructs nine blank and nine pointer-ready reset failures,
the native phase/clock/pointer errors and the rational geometry error.
The reset minimum-eigenvalue argument applies to correlated initial
carriers through the pinned rank-one conclusive reset theorem. Every
sum equals the frozen complete-record bound.

Deleting the photon frequency changes a retained measured probability
by approximately 0.08277. Deleting the common vacuum scalar changes
the tested state vector by approximately 0.8894 while leaving the
tested vacuum/one-photon probabilities and independent clock unchanged.
The experiment measures an excitation frequency. A coherent reference
between vacuum and padding, or a controlled-source intervention, can
detect the retained vacuum phase; that is a different preparation or
instrument. The audit executes this distinction rather than claiming
that the measured Ramsey record determines absolute vacuum energy.
One time sample also does not identify an unaliased frequency: the same
cosine record admits integer-period frequency shifts. The positive
frequency interval is supplied by the declared Hamiltonian certificate.
A frequency-estimation experiment requires a specified time family and
photon quadrature settings; changing X to Y here changes only the clock
analysis.

## Exact cap boundary and a failing interaction transfer

For cap one, a_i=|vac><i| gives
[a_i,a_j^*]=delta_ij I-P_top(delta_ij I+a_j^*a_i).
The matrix-unit checks include diagonal and off-diagonal cases. A
single top state has canonical-commutator defect norm 2. Summing the
literal truncated quadratics gives diagonal energies E0 on vacuum and
omega_j/2 on a one-photon state. The declared top correction adds
(omega_j+sum omega)/2 and gives exactly E0+omega_j.

This free-Hamiltonian identity does not certify a Weyl displacement.
For a unit one-mode displacement, the cap-one vacuum survival is
cos(1)^2, while the infinite Gaussian survival is exp(-1). The executed
probability discrepancy exceeds 0.075. Intermediate occupation tails
are required even when the initial and final desired record domain
contains at most one photon. Harmonic translations require their own
native approximation as well.

## Smallest spatially distinct charged transfer

A size-three periodic cell lattice is the smallest member whose
positive and negative axis steps have distinct endpoints. It has 27
cells, 324 carrier ports and 972 edges. One concrete neutral experiment
holds a compensating charge in cell (0,0,0) and places the moving charge
in cell (1,1,0). Its two x-direction endpoints, (2,1,0) and (0,1,0),
are distinct from each other and from the fixed charge. The audit
reconstructs their actual root-port identifiers from the size-three
graph rather than reusing the size-two route table.

The separate spatial-domain control constructs all 702 size-three face
columns and computes rank 646 modulo the prime 101, which is a rational
rank lower bound. Connectivity gives incidence rank 323. The exact
incidence and three independent period identities give the matching
upper bound 972-323-3=646. Thus cap one has 647 physical states and
requires ten spectral bits, with 377 padding states. Adding the same
clock requires 2,048 complete outcomes per pointer law.

The actual diagonal compiler rejects a 1,024-entry phase target because
its declared limit is eight bits. This rejection is executed. The
existing Pauli-string primitive admits ten-bit labels, and an exact
template count gives (3b-5)2^b+6=25,606 pulses at b=10. An additive
Walsh compiler can therefore retain the same primitive ancestry, with
a fresh certificate for all 1,024 phases. A size-three spectral member
also needs its own rational frame bound, vacuum-plus-excitation factor
646/2+1=324, eleven blank resets and eleven pointer-ready resets.
These resource identities are not a size-three native implementation.

The corresponding transfer must parameterize the routed Dirac/Fock
source on this graph, use P_L J and P_T J from the same current, retain
fermionic signs and the harmonic packets, and implement a complete
occupation/current/field/clock instrument. A charge-energy finite
difference can be read by a specified coherent branch and inverse-route
analysis. An unperturbed force-displacement claim additionally requires
the actual evolved density/current response. These are distinct
measurements. Neither follows from a static minimum or a change of
clock label. Real-space asymptotics and uniform preparation/tail/native
budgets are required along an increasing-volume/refinement family.

Reproduce with the repository virtual environment and single-thread
BLAS: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
PYTHONDONTWRITEBYTECODE=1 reverse-engineering-reality/.venv/bin/python
openwave-tuning/common_dynamics/photon_audit/audit.py`. Run
`facade_control.py`, `spatial_domain.py`, and then `freeze.py` after the audit succeeds.
All official row-status outputs are empty.
