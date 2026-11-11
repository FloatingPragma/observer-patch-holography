# Finite charged response including virtual photons and pairs

The source is the unchanged orthogonal Wilson
cube at g=1/64, total photon cap3 and completion c=21/2, with its original
Psi/rho/Y frame, beta-sea scalar and cost4|R|. The charged-band package is
used only for exact matrix/control helpers. Its different completion23/2
is not substituted. Every low state below is in the reducing native code;
full native complements are separately excluded by energy bounds.

## Lowest charged rank-two cluster

Use E_q=E_free+4|q|+m lambda_q, with lambda_0=0 and
lambda_1=sqrt5/2. Divide excitation differences by m; the raw reference
energy is exactly4, not4m. Let lambda1=sqrt5/2,
lambda2=sqrt(17-4sqrt3)/2, v=1/16+1/(12sqrt6),
qstar²=(13+5sqrt6)/48 and C0=8,C1=67/8.
The exact fixed-sector perturbation bound is nu_q=24g qstar+Cq g².
On the free lowest two-dimensional charged eigenspace with photon vacuum,
the odd phase has zero compression. The literal finite-coordinate inequality
|1-<exp(igQ)>|<=g²v/2 gives

    ||P(H-H0)P|| <= d1=g²(12v+C1).

Each full-Fock bond has norm2 times its scalar phase coefficient. The
factor12 is therefore12edges times2 times the1/2 Taylor factor.
Min-max on the codimension-two code subspace, followed by the unchanged
native-complement bound, yields

    gap >= min(lambda2-nu1,21/2-Delta_upper)-lambda1-d1 > 0.

Only one nu1 occurs. The coarse2nu1 criterion did not establish this
result. The rank-two multiplicity refers to the fixed-charge active source,
not physical spin/species, arbitrary reference spectators or Lorentz mass.

## Actual perturbation and virtual spaces

On the physical code write H/m=H0+g A+g²B+remainder, where

    A = sum_e Q_e dGamma(J_e),
    B = -1/2 sum_e Q_e² dGamma(E_e) + 1/2 rho_ref^T G rho_ref.

Je and Ee are literal derivatives of the orthogonal Wilson bonds in the
fixed charge frame. The compensated density retains the root reference.
The frozen total-cap boundary completion belongs to H0. It is not removed.

A maps a neutral free sea to a particle/hole pair and one photon. Acting
on the lowest charged free two-space also moves the extra particle or
creates an additional particle/hole pair. The actual binary64 diagnostic
executes these full-CAR contractions with all five field modes and their
free energy denominators. It includes the two-photon-capacity boundary
through the literal K3 source, although this second-order virtual space
has only one photon. Photon-vacuum kinetic compression alone omits the
negative resolvent contribution and is not the response calculated here.

Let P denote the neutral free vacuum or the charged lowest free two-space.
The second-order matrix is

    M_q = P B P - P A (H0-Eq)^(-1)_(one-photon) A P.

The reference-corrected charged-minus-neutral transition has target
lambda1 I + g²(M1-M0 I). The diagnostic values are not certificates until
finite-precision and remainder enclosures are completed.

## Photon-parity Feshbach remainder

Let E/O denote even/odd photon parity, retaining all56 photon states.
The even perturbation has norm <=g²b_q, b_q=12qstar²+Cq.
The odd perturbation has norm <=g a, a=24qstar, and differs from gA by
at most gamma g³, gamma=4qstar³. These follow by diagonalizing each literal
finite Q_e separately; no finite CCR or commuting-coordinate identity is
used. On P, the direct even fourth-order remainder is <=3v²g⁴ because
<vac,Q_e⁴ vac>=3v² holds exactly at cap3.

Suppose certified compressions give ||AP||<=alpha, ||R0,odd AP||<=beta,
||BP||<=bP, ||PBP||<=D. Here R0,odd=(H0,odd-Eq)^(-1).
One can use bP<=12sqrt3 v+sqrt(Cq ||P C P||), where C is the positive
Coulomb coefficient. Define d=g²D+3v²g⁴. All low eigenvalues are <=Eq+d.
The free odd gap is at least2 in both charge sectors. The even complement
gap is delta0=sqrt5 (neutral) and delta1=lambda2-lambda1 (charged).

Set

    muO = 2-d-g²b_q,
    w = g²b_q+(ga)²/muO,
    muE = delta_q-d-w,
    s = g alpha+gamma g³,
    rho4 = qstar^4 g^4,
    t0 = g²bP+rho4+(ga)s/muO,
    r = d+s²/muO+t0²/muE.

The rho4 allowance bounds the full even Taylor remainder acting on P;
it is distinct from the sharper direct P-P cosine compression.
When both mu values are positive, two successive Schur complements show
|lambda-Eq|<=r. Eliminate the odd subspace first, then the excited even
subspace. Both resolvents are positive on the whole low interval.

The source-specific weighted-resolvent improvement is

    z = g beta[1+(r+g²b_q)/muO]+gamma g³/muO,
    t = g²bP+rho4+(ga)z.

A uniform low-space effective-matrix error relative to Eq+g²M_q is

    epsilon_q = 3v²g⁴
       +[2g alpha gamma g³+(gamma g³)²]/muO
       +g² beta²(r+g²b_q)[1+(r+g²b_q)/muO]
       +t²/muE.

The terms are respectively the direct cosine remainder, the odd sine
remainder, the shifted odd resolvent, and the excited-even elimination.
For the resolvent term use
Rodd-R0=R0[(lambda-Eq)-V_even,odd]Rodd and
||Rodd AP||<=beta[1+(r+g²b_q)/muO].
The second Schur complement uses the entire excited-even carrier,
including all extra fermion pairs and photon states. Its derivative is
negative semidefinite, so spectral inertia brackets the ordered exact low
energies by the ordered eigenvalues of g²M_q plus/minus epsilon_q.
Native complements remain reducing and are excluded separately.

No interacting infinite-volume/continuum result, empirical mass fit,
Standard Model assignment or source-derived GeV unit is asserted.


## Exact coefficient replay and its certificate boundary

`algebraic_intervals.py` encloses every admitted nonnegative rational power
by integer arithmetic. For x >= 0 and positive integer q, set S=10^60 and
k=floor(root_q(floor(x S^q))). Then k/S <= root_q(x) <= (k+1)/S;
if k^q=x S^q the interval is exact. Raising these positive endpoints to
p and inverting their order for p<0 encloses x^(p/q). Negative bases with
noninteger exponents and intervals containing a reciprocal pole are rejected.
This includes the coefficient 2^(3/4) that the restricted evaluator rejected.
The tests verify the enclosure by exact rational q-th powers, including
negative rational exponents, perfect roots, nested field variances and poles.

The canonical charged-band helper executes every frozen scientific hash
admission. Its retained source and reference receipts define the construction.

The binary64 candidate columns U are rational inputs, not exact eigenvectors.
Let eta bound the Hermitian row norm of U*U-I and r bound the row norm of
HU-UD. With eta<1/2, polar orthonormalization gives Q with ||U-Q||<=eta.
The dimension is32, hence ||HU-UD||<=sqrt32 r<6r. Therefore

    ||H-Q D Q*|| <= 6r + (||H||+||D||) eta.

The replay uses the exact Hermitian row bound for H and maximum absolute
diagonal entry for D. It admits this one-body error and computes the virtual
contractions in exact rational arithmetic with explicit dyadic term rounding.

The wide proposed many-body compressed-matrix transfer allowance is not a
proved enclosure. It additionally needs projection estimates for the exact
sea/charged eigenspaces, the induced Fock-space map and resolvent comparison.
The output types it as a proposed allowance and marks both that transfer and
the interacting spectral remainder as unproved. Successful arithmetic replay
therefore does not promote the diagnostic second-order eigenvalues to
certified interacting energies or masses.
