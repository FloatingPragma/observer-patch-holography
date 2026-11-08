# Conserved charge, local contrast and finite return survival

The source is the selected eight-site cube of `native_ordered_ground`:
12 antiferromagnetic bonds with coupling 1, selection field 1/4,
Hamiltonian H, unique normalized ground Ω, and its declared source clock.
No energy, state or time is replaced by a fitted particle model. The
`protected_defect_catalog.charge` construction pulls the full fermionic
spin charge through the actual spin inclusion and Marshall map:

    Q = (1/2) sum_x epsilon_x X_x,    [H,Q]=0,    Q Ω=0.       (1)

Each epsilon_x is the actual bipartite Marshall sign. Every complete C/E
repair branch, including its record, commutes with the corresponding full
source charge. Ground nondegeneracy and the exact charge-reversing source
automorphism prove QΩ=0. The charges here are internal source charges;
no electric charge or physical spin is identified.

## 1. A local complete instrument creates coherent charge channels

Set Z=Z_0 and A_q=Z(I-q epsilon_0 X_0)/2 for q=±1. Direct source algebra
proves

    [Q,A_q]=q A_q,    A_+ + A_-=Z,
    A_+* A_+ + A_-* A_-=I.                                  (2)

Measure only local X_0, retain both outcomes, and apply local Z_0. This
implements the two Kraus operators A_q. With

    w_q=||A_q Ω||²,    a_q=A_q Ω/sqrt(w_q),

both branches have w_q>4/25. Their probabilities are about .835475 for
q=-1 and .164525 for q=+1. Conditioning names a branch; it does not delete
its original probability. Fine measurement of all eight local X values
and grouping their charges is a different instrument: it generally
dephases inside a charge sector. That instrument is not used to prepare
a_q. The unconditioned local creation leaves every exterior marginal
unchanged. A conditioned branch can steer an entangled exterior marginal,
so exact exterior equality is not asserted for a_q.

Define the return survival on the same clock by

    c_q(t)=<a_q, exp(-itH) a_q>,    S_q(t)=|c_q(t)|².          (3)

The arbitrary common energy origin cancels in S. Each sector Q=q has
dimension 56. Conservation of Q does not imply a_q is an eigenvector of
H. The verified energy variances are

    q=-1: [0.338256843566, 0.338256843755],
    q=+1: [0.921555894395, 0.921555894820].                   (4)

These positive intervals prove transitions of the prepared vector into
its orthogonal complement within the same conserved charge sector.

## 2. Exact complete spectral and source-moment certificates

Each source bond is minus a unit projector. Their sum lies between -12I
and 0. The selection field has norm at most 1. Thus

    -13 I <= H <= I,    spectral width <=14,
    ||H+6I|| <=7.                                           (5)

The retained 256-column integer eigensystem proposal is accepted by the
independent residual and Gram certificate of `native_ordered_ground`.
It certifies every level, the ground gap and the ground-vector error e_g.
For its integer ground vector g, b_q=(2A_q)g and n_q=b_q*b_q, the normalized
branch vector error is bounded by

    e_q = 4 e_g sqrt((g*g)/n_q).                              (6)

This follows from ||A_q||<=1 and normalization of the nonzero branch.
Every integer matrix power in

    m_n = b_q* (H+6I)^n b_q / n_q                           (7)

is evaluated exactly. The degree-N Taylor scalar p_N(t) approximates the
proposal's unitary amplitude with error at most (7t)^(N+1)/(N+1)!: the
integral remainder contains a unitary factor of norm 1. The true amplitude
therefore lies within that radius plus 2e_q. Squaring its modulus with
outward rational square-root bounds gives the survival interval.
No eigenvalue is selected, dropped, rounded into a mass, or assigned an
empirical label.

A separate computation retains all 256 approximate spectral weights

    v_j = |V_j* b_q|²/n_q.

Their sum differs from 1 by at most the certified Gram defect d. Let O be
the polar orthogonal factor of V. The same certificate proves
||H-O diag(E_j)O*||<=r and ||V-O||<=d. Duhamel's identity bounds the unitary
operator difference by tr. Consequently the full spectral scalar
sum_j v_j exp(-it(E_j+6)) differs from the true amplitude by at most

    tr + d(2+d) + 2e_q.                                     (8)

An independent scalar Taylor remainder is added. Every tabulated short
time has overlapping source-moment and full-spectrum enclosures. The
spectral proposal is evidence subjected to exact inequalities, not an
unverified diagonalization.

The variance from (7) has error at most 6*7² e_q. The first term changes by
at most 2*7²e_q, and the square of the first moment by at most 4*7²e_q.
No domain beyond this finite matrix is involved.

## 3. A certified first retention time

For the exact probability measure ν_q of H in a_q,

    S_q(t)=integral integral cos((E-E')t) dν_q(E)dν_q(E'),
    1-S_q(t) <= Var_q(H) t².                                (9)

On 0<t<3/14 every nonzero energy difference satisfies
0<|(E-E')t|<3<pi. Differentiation of the finite sum gives S'_q(t)<0
because the variance is positive. This proves strict monotonicity only
on that stated interval. Exact endpoint enclosures and rational bisection
therefore certify the unique first time at which survival equals .99:

| Source charge | First .99-retention time, source clock |
| --- | --- |
| -1 | [0.1726210202, 0.1726210331] |
| +1 | [0.1044566248, 0.1044566376] |

This is an observable overlap-retention threshold, not a decay lifetime.
For q=-1 the same complete calculation gives

    S_-(2) in [0.301809595920, 0.301809595922],
    S_-(5) in [0.780935158862, 0.780935158865].                (10)

The rebound is certified from degree-160 integer source moments. It rules
out treating the whole response as monotone exponential decay.

## 4. Complete echo and local-contrast readouts

After creation outcome q and time evolution, apply Z, then inverse ground
preparation and a complete computational-basis readout. Designate the
prepared-input basis word as return; retain all other words and every
native first-failure label. Charge orthogonality gives

    <Ω,Z exp(-itH) A_q Ω> = w_q c_q(t).

Thus the original joint probability of charge q and return is

    R_q(t)=w_q² S_q(t),
    Pr(q,nonreturn)=w_q-R_q(t).                              (11)

The four coarse probabilities sum to one. Return at t=0 is w_q²,
not w_q: the reference-ground echo tests only the part that can be
recombined with the ground. This factor is retained in every receipt.
The ideal ground-preparation/inverse unitary is supplied by the finite
native preparation theorem. The local X pointer and Z word are native
controls. Complete H evolution and inverse compilation use the declared
coherent controller service. This package executes exact matrix and
interval certificates, not the full echo pulse tape.

If complete implementation has total-variation event error eta and the
weight estimate has error delta_w, with true and estimated w at least b,
the reconstructed ratio has conservative error

    |Rhat/what² - R/w²| <= eta/b² + 2 delta_w/b³.             (12)

This follows by separating numerator and inverse-denominator errors,
using R<=1 and the derivative bound |(x^-2)'|<=2/b³. Complete source-history
comparison and control errors enter eta before division. No failure is
conditioned away.

A conserved accessible contrast survives the return loss. For every t,

    sum_x epsilon_x [Pr_q(X_x=+) - Pr_Ω(X_x=+)] = q.         (13)

Some site therefore has probability contrast at least 1/8 in magnitude;
the site need not be fixed. With original preparation weights,

    max_x |Pr(q,X_x=+) - w_q Pr_Ω(X_x=+)| >= w_q/8 >1/50.   (14)

There is a fixed complete detector that reads this distinction without
locating the maximizing site. Average equally over the eight declared
local-X settings, retain x, and score the observed eigenvalue s by
epsilon_x s. Equivalently use a uniform three-bit setting register and
retain its outcome. The conditioned score mean is q/4, and

    Pr(score=+ | q)=1/2+q/8,
    Pr(q,score=+)=w_q(1/2+q/8).                              (15)

The ground reference gives 1/2. All 32 charge/setting/binary records are
retained; the four charge/score coarse outcomes are displayed. The
original-weight binary contrast is exactly w_q/8 for every source time.
The bound persists under complete repairs preserving (1). Arbitrary
controls that break Q are outside that assertion. If the joint detector,
reference detector and weight estimate have errors eta_joint, eta_ref,
delta_w, the contrast error is at most their sum. This makes charge
persistence measurable even when S is small.

Equations (13)-(15) establish finite internal-charge contrast. They do not
identify a fixed defect location, a propagating relativistic wavepacket,
a particle species, or an infinite-volume superselection charge. Local
contrast, return survival and conserved charge answer different questions.

## 5. What a continuum decay calculation must establish

For a finite charge sector let p_j be its weights grouped by distinct
energies. Direct integration of the finite frequency sum gives

    lim_(T->infinity) (1/T) integral_0^T S_q(t) dt
      = sum_j p_j² >=1/56.                                  (16)

Simultaneous recurrence on the finite phase torus makes S arbitrarily
close to 1 at arbitrarily large times. Both statements prohibit an
irreversible finite-source decay law tending to zero. They do not prohibit
loss of a local signal or finite-time relaxation.

There is a precise source-limit route. Suppose the normalized prepared
channel spectral measures ν_n converge weakly to a probability ν on the
same physical energy convention. Then c_n(t)->c(t) for every fixed t
because exp(-itE) is bounded and continuous. Tightness gives convergence
uniform on compact time intervals: truncate to a compact energy set,
use equicontinuity there, and bound both tails. If ν has an integrable
density, approximate it in L1 by compact step functions; each step's
Fourier transform tends to zero, so c(t)->0. This proves return dephasing
in that limit. The infinite-volume limit precedes large time; (16) cannot
be bypassed by a large finite simulation.

Even this limit does not identify particle decay. A decay claim needs
source-derived asymptotic incoming/outgoing particle observables, an
unstable sector with appropriate thresholds, and the transition records
into its decay products. A resonance-width interpretation additionally
requires a controlled continuation or equivalent scattering statement.
Conserved total charge can pass into several products while (13) holds.
No such particle or resonance assignment is supplied by the finite
retention times above. In particular, the two charge branches are not
identified as a measured particle and antiparticle.
