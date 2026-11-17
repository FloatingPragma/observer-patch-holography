# What an observer criterion does and does not select

Platform: row 11 untested. This is a class-T theorem about the declared coherent
generation-register construction (class F). It supplies a conditional analyzer
selection theorem and identifies the remaining operational selection law. It
does not identify three physical charged particles.

## 1. Exact source response before choosing an analyzer

Let U be the actual 432-state tick and p the existing two-object preparation.
Keep the family register F=C^3. Normalize the source's three orthogonal flag
vectors to H,+,-. Define A_j=(<j| tensor I)U(|p> tensor I). Direct exact
application of `register.full_tick`, itself using `Observer.repair`, gives

    A_H = I/3,       A_+ = R/(3 sqrt(2)),       A_- = R^-1/(3 sqrt(2)),
    R^3=I,          Tr(R)=Tr(R^2)=0.

The source clock multiplies the last two amplitudes by z and z^-1, |z|=1.
The complete prepared output is an isometry. The receipt retains all 144 raw
object-coordinate Kraus maps K_y=(A_y+sqrt(2)B_y)/24 as integer tables, including
zero maps, and checks sum_y K_y^dagger K_y=I exactly. It also checks U on every
one of its 432 input columns.

For the three flags,

    sum_j A_j^dagger A_j = (2/9)I;
    G_ij = Tr(A_i^dagger A_j)/3 = diag(1/9,1/18,1/18)_ij.

The orthogonal archive complement has effect (7/9)I. Replacing the flag basis by
any orthonormal basis preserves this full effect. Thus selecting a particular
analyzer output always comes with the other flag outputs and the entire
remaining archive; neither their weight nor their reference correlations are
discarded by the theorem. For a two-outcome read, the complement effect is
I-M_c^dagger M_c.

## 2. A continuum of covariant analyzer responses

For any unit vector c=(c_H,c_+,c_-) in flag space, its selected Kraus map is

    M_c(theta) = conjugate(c_H) I/3
              + conjugate(c_+) z R/(3 sqrt(2))
              + conjugate(c_-) z^-1 R^-1/(3 sqrt(2)).                 (1)

Every such M commutes with R. Family C3 covariance therefore imposes no
additional restriction on c. In the real parity-even subfamily

    c=(sqrt(t),sqrt((1-t)/2),sqrt((1-t)/2)),  0<=t<=1,

write M=aI+b(zR+z^-1R^-1). Then

    a=sqrt(t)/3,    b=sqrt(1-t)/6,
    b^2/a^2=(1-t)/(4t)  (t>0).                                  (2)

The source's balanced ratio b^2/a^2=1/2 holds exactly at t=1/3. Its signed
circulant quotient is (1+t)/(6t); it is the mass-style square-root quotient
only where the three eigenvalues are nonnegative and the squared-eigenvalue
mass interpretation is supplied. No mass interpretation is used here.

An exchange of H and + is not a symmetry of this prepared response implemented
by unitary transformations on the input/output family: A_H^dagger A_H=I/9
while A_+^dagger A_+=I/18. Unit conjugation cannot change these scalar effects.
This excludes that proposed flag-permutation symmetry, not every possible
larger observer symmetry or a different preparation. The family third-turn
and interchange of + and - must not be treated as an S3 permutation of flags.

## 3. A positive conditional selection theorem

Suppose the specified observer task is an **unbiased flag eraser**: maximize the
worst probability that a selected analyzer outcome accepts any one of the three
normalized orthogonal flag states. This task's objective is

    J(c) = min_j |c_j|^2,       sum_j |c_j|^2=1.

**Theorem.** max J=1/3, and equality holds exactly when every |c_j|^2=1/3.

**Proof.** The least of three nonnegative numbers is at most their average.
If the least is 1/3 and their sum is one, all three must be 1/3. Conversely
that choice attains the bound. Every relative phase is unrestricted. QED.

This is a useful conditional derivation of the equal analyzer magnitudes.
It becomes a source-selection theorem only after the source specifies this
actual task, normalized-flag comparison, preparation and implementation. The
flagship and generation-register source explicitly retain the Fourier analyzer
as a construction input; the audited sources do not supply this additional law.

Even requiring a Hermitian family compression leaves a phase family. Set

    c=(1,exp(-i phi),exp(i phi))/sqrt(3).

Equation (1) is the balanced Hermitian circulant with effective phase theta+phi.
Thus unbiased erasure plus family covariance and Hermiticity do not fix the
phase. An operationally required phase convention or coherent prediction task
could fix it; it would need to be specified and realized on the same source.

## 4. Two other exact observer objectives

The mean selected probability on the uniform family reference is

    p(c)=Tr(M_c^dagger M_c)/3
        = |c_H|^2/9 + (|c_+|^2+|c_-|^2)/18
        = (1+|c_H|^2)/18.                                    (3)

Cross terms vanish by Tr R=Tr R^2=0. Maximizing this success probability gives
|c_H|^2=1, p=1/9; the compression is a scalar. The equal analyzer gives p=2/27.
Equation (3) is independent of phase. If every trial has an equal declared
cost, maximizing the retained outcome per trial is this objective. We do not
assume arbitrary physical analyzer implementations have equal work costs.

A different task maximizes the least **original branch weight** reaching the
selected output after commitment. Let w=(1/9,1/18,1/18) and p_j=|c_j|^2. If
min_j w_j p_j >= m, then p_j>=m/w_j and hence

    m <= 1/(sum_j 1/w_j)=1/45.

Equality is attained uniquely in magnitudes at p=(1/5,2/5,2/5). Thus normalized
flag erasure, original-weight minimax and mean throughput are three different
well-defined questions with different optima. This calculation does not refute
other observer criteria. It makes the objective that would justify the desired
one explicit, without choosing it by the known lepton masses.

## 5. Committed-record predictive objectives cannot select a lift phase

Consider any complete instrument {K_r}, including all failure and timeout
labels. Replace K_r by exp(i phi_r)K_r. Each completely positive branch map

    rho -> (K_r tensor I_ref) rho (K_r^dagger tensor I_ref)

is unchanged, for every input and untouched reference. The same statement
holds if a record has several Kraus operators all multiplied by its common
phase.

**Adaptive-objective theorem.** Compare two models with these identical
committed instruments at every available action. Fix the same initial state,
interventions and controller. The controller may choose its next action from
the entire past record and source-independent randomization, and may retain
quantum workspace. Then every complete finite transcript and its terminal
quantum/reference state have identical unnormalized outputs. Consequently any
record-measurable predictive loss, consensus score, entropy objective or
cost function has the same distribution in both models, provided the counted
costs for the compared actions also agree.

**Proof.** At zero steps the states agree. Condition on a retained history h.
The identical controller applies the same instrument or mixture of instruments
to the same conditional state. Each next committed branch map is identical,
so every extended history h,r has the same unnormalized state. This proves the
claim inductively without conditioning away any zero, failure or timeout
branch. Sum the identical histories for any coarsening, score or randomized
stopping policy. A finite horizon includes its timeout as an outcome. The
countably infinite version follows from equality of cylinder measures when
the stopping/output variables are measurable. Tensoring the argument with the
retained workspace proves the reference-sensitive statement. QED.

Taking the infimum over the same admissible controller class preserves this
equality. Least-resource prediction therefore cannot break this equivalence
using only the committed task and equal resource counts. This conclusion has
two necessary boundaries:

1. A coherent read before commitment can distinguish the phases. Such a read
   changes the accessible task/instrument and is outside the theorem.
2. A phase-dependent physical cost may distinguish implementations. Equal pulse
   counts alone do not prove equal work, storage or initialization costs.

The existing `model/observer/source_selection.py` executes precisely such an
equal-count example on the source C/E instruments: its two scheduled phase
words have equal declared pulse/service counts, committed CP maps coincide,
and its coherent read distinguishes 3/4 from 1/4; commit-first gives 1/2 under
both. Those are existing source results read here, not a fresh execution claim.
The source resource function explicitly excludes physical initialization,
storage and readout costs. The analogous 432-state phase control below is
exactly executed algebraically; physical implementation cost equality for that
larger register is not claimed.

The equal counted schedules are not a theorem of equal *minimal* compilation
cost. In particular, the existing witness deliberately retains zero-angle
phase services in the untwisted schedule. A compiler permitted to eliminate
those services changes that cost question. A global least-resource selection
claim requires the source's admissible implementations and cost functional,
then an optimization over them. This packet neither assumes that unrestricted
optimum ties nor claims that the present per-record phases minimize it.

## 6. The ambiguity acts on the actual lepton observable

On every output of the actual tick apply the phase

    D=(-1)^(C record bit).

This multiplies H by +1 and both repair flags by -1. For every raw object
coordinate y and all 432 input columns, its committed branch is only multiplied
by this sign. Every committed CP map is therefore unchanged. Coarsening the
records preserves the equality, and the native C record alone already resolves
the sign. The clock commutes with D, so the statement holds at every theta.

For the equal analyzer write

    M_+ = aI+b(zR+z^-1R^-1),
    M_- = aI-b(zR+z^-1R^-1),
    a=1/(3sqrt(3)), b=1/(3sqrt(6)), E_+/-=M_+/-^dagger M_+/-.

Both selected effects have trace 2/9. They need not have the same spectrum:

    Tr(E_+^2)-Tr(E_-^2) = 48 a b^3 cos(3theta)
                       = 4sqrt(2)cos(3theta)/243.              (4)

To check (4), expand (aI+bX)^4-(aI-bX)^4 with
X=zR+z^-1R^-1. Tr X=0 and Tr X^3=3(z^3+z^-3). At the inherited theta=2/9,
cos(2/3)>=1-(2/3)^2/2=7/9, so the difference is at least 28sqrt(2)/2187>0.
This is a strict exact spectral witness with no phase scan or mass scoring.

If the flag record is committed first, the selected effect is instead
sum_j A_j^dagger A_j/3=2I/27 for both signs. The distinguishability therefore
depends on access to coherent archive combinations, rather than a relabeling
of the committed records. No physical mass is assigned to either effect.

## 7. A3 entropy selection does not silently add a Fourier analyzer

The basis adopts minimum relative entropy to a **specified** reference under
**specified** observer-visible constraints. Let the feasible set be invariant
under a unitary D and let D commute with the reference tau. Functional calculus
and trace cyclicity give

    Drel(D rho D^dagger || tau)=Drel(rho || tau).

Hence the entire orbit of any minimizer is a minimizer. If the minimizer is
unique on a convex feasible set, it must commute with D. For independent flag
phases this implies a diagonal state. It does not uniquely select a coherent
rank-one equal analyzer. If feasibility is restricted to pure analyzers, the
problem is no longer that convex state-selection problem; with reference I/3
all pure states have relative entropy log(3), irrespective of their phases or
magnitudes. Maximizing Shannon entropy of the *coefficient probabilities*
does select equal magnitudes, but is a different explicitly supplied objective.

Coherent constraints or a reference that carries a phase can break this
invariance. A3 then makes a conditional prediction relative to those inputs.
Their source realization, rather than an assumption of an external spacetime,
is the remaining question for the observer-first derivation.

## 8. The next constructive obligation

An adequate next law can specify a reproducible coherent observer task for
which (i) normalized flag erasure is required, (ii) the relative phases are
operationally fixed and survive the complete record comparisons, (iii) the
preparation/readout/reset cost and failures are accounted for, and (iv) the
resulting observable agrees with the common rest/spatial-energy dictionary.
No uniqueness theorem for every imaginable architecture is required. The
specified law must reproduce the intended records and physical interpretation
on one candidate. The present result derives the balance conditional on one
precise task and proves exactly what committed prediction cannot select.

## 9. Post-run constructive extension: a calibrated unbiased analyzer

Added analytically at 2026-11-14T16:25:57Z, after the preregistered certificate
and replay. This theorem is not presented as an additional preregistered run,
a new numerical gate or a native implementation. It uses the exact response
maps already certified above.

There is a concrete target-free task that fixes the remaining analyzer phases
within this declared source. Use its **unclocked** tick U0 and prepare the
normalized invariant family vector

    v0=(1,1,1)/sqrt(3),     R v0=v0.

The projection of the complete output on the three-flag subspace is

    [ |H>/3 + |+>/(3sqrt(2)) + |->/(3sqrt(2)) ] tensor v0.     (5)

Fix the unbiased condition |c_j|=1/sqrt(3) from section 3. Among these analyzers
maximize the unconditioned selected-output probability on preparation (5).
Write c_j=exp(i alpha_j)/sqrt(3). The selected amplitude is

    [exp(-i alpha_H)/3
      + exp(-i alpha_+)/(3sqrt(2))
      + exp(-i alpha_-)/(3sqrt(2))]/sqrt(3).

The triangle inequality gives maximum amplitude (1+sqrt(2))/(3sqrt(3)), hence

    maximum probability = (3+2sqrt(2))/27.                    (6)

Every coefficient in this three-term sum is strictly positive in magnitude.
Equality in the triangle inequality therefore requires all three phases to
agree. Conversely equal phases attain (6). The analyzer projector is uniquely

    |Phi><Phi|,    |Phi>=(|H>+|+>+|->)/sqrt(3),

up to its physically irrelevant global phase. This is a conditional phase
selection theorem derived without a lepton mass or target phase.

All outputs remain. The complementary outcome has probability
1-(3+2sqrt(2))/27=(24-2sqrt(2))/27; no successful branch is renormalized. Keeping
the other two orthogonal flag outcomes and the original archive complement is
a complete implementation of the same measurement, if those read services
are supplied. On arbitrary inputs the complement effect remains I-M^dagger M.

The task's timing is essential. **Calibrate at U0, then freeze the analyzer
before applying K_theta**: its phases remain equal and the later response has
the actual theta. If one instead reoptimizes the same objective after applying
K_theta, equation (5)'s flag amplitudes become 1/3, z/(3sqrt(2)), z^-1/(3sqrt(2)).
The maximizing analyzer is then (1,z,z^-1)/sqrt(3); equation (1) loses both
clock phases. Retuning would cancel theta and restore the phase-zero response.

Thus this construction replaces a free analyzer phase by an explicit
preparation, calibration and freezing protocol. It does not establish that
the native process supplies that protocol or requires the normalized-flag
fairness objective. It also does not select K_theta, convert the exclusive
count into a phase, identify squared amplitudes with physical rest masses,
or authenticate preparation/control costs. Those remain separate source and
physical-identification obligations. The calibrated eraser is a concrete
candidate observer task for a subsequent source realization.
