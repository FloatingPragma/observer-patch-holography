# Charged records about the complete scalar ground

This calculation uses the unchanged Hamiltonian

\[
 H=W(\operatorname{diag}(1,2,3))+W^\dagger+N_H
 \tag{1}
\]

on all `3*2^18=786432` scalar and canonical-anticommutation-relation (CAR)
occupation states. The scalar has a vacuum and two single-occupation
components. Its pair cost is one; the common fermion offset is zero.
The reference is the verified unique full-carrier ground |G>, of charge
`Q*=N_F-2N_H=6`, energy about −7.94630530665, and zero gauge/spin quantum
numbers. The [owning proof](../particle_scalar_ground/PROOF.md) certifies
that ground and its preparation. Its rational singlet trial has trace
distance below `5.1*10^-6`; this error is retained in every spectral weight
reported here. The [307-state prepared-sector response](../particle_observer_outputs/CHARGED_HIGGS_RECORDS.md)
is a different experiment and is not used as this vacuum spectrum.

## One complete compensated probe

The existing full-CAR current theorem gives the transformation of every
creation and annihilation operator under the actual gauge generators.
For the eighteen modes i define

\[
 K_+=\frac1{\sqrt{18}}\sum_i c_i^\dagger\otimes|i^*\rangle,
 \qquad
 K_-=\frac1{\sqrt{18}}\sum_i c_i\otimes|i\rangle.
 \tag{2}
\]

The reference carries respectively the dual or direct gauge representation.
It therefore compensates the field probe on all input occupations, rather
than treating a bare charged operator as a gauge-invariant local pointer.
The complete effects are

\[
 K_+^\dagger K_+=I-N_F/18,\qquad
 K_-^\dagger K_-=N_F/18,
 \qquad\sum_{\pm}K_\pm^\dagger K_\pm=I.
 \tag{3}
\]

These are the existing universal CAR identities, applied to the complete
carrier and its retained quantum references. No channel is accepted or
discarded. A finite unitary dilation realizes this declared instrument;
the source energy supplied during a charged intervention belongs to its
apparatus/controller account. A zero reference Hamiltonian is the declared
energy convention here. Giving that reference a different energy requires
including it in an apparatus energy-transfer calculation.

The apparatus reads addition/removal, family, and the whole left or
singlet gauge multiplet. Those family/multiplet projectors commute with
the gauge action. It does not assume a preferred weak-gauge component
can be read independently. The subsequent energy apparatus is the spectral
measure of the same full H, with the same `exp(-itH)` clock. All unrecorded
probe/reference degrees of freedom remain in the instrument's output.

The exact ground is a gauge and spin singlet, and H commutes with those
actions. The four left modes and two singlet modes of each family thus
have identical diagonal spectral laws within each multiplet. The code
calculates one representative and restores multiplicities four and two.
The original joint weight of a spectral group is its single-mode weight
multiplied by `multiplicity/18`; this factor is never divided out of the
complete histogram. The representative moment equalities are also tested
directly against the literal source.

## Exact finite response domains

The full ground proof writes the empty-scalar part as u and its occupied
image as Au. In its six-cell basis the types are
`dual,dual,fundamental,fundamental,fundamental,fundamental`.
A single CAR operator changes only one cell, either to its other active
type or to an inactive singleton. Addition reaches Q*=7 and removal
Q*=5. The code takes the literal Hamiltonian closure of every original
occupation in `c_i^(dagger)u` and `c_i^(dagger)Au`, including every
nonzero source image. The exact ground has positive coefficients in the
owning cell frame, so this support is the complete exact-ground probe
support, not a numerical-threshold choice.

The representative reducing blocks have dimensions:

| Family | Left removal | Left addition | Singlet removal | Singlet addition |
|---:|---:|---:|---:|---:|
| 1 | 165 | 80 | 53 | 240 |
| 2 | 60 | 212 | 144 | 93 |
| 3 | 60 | 212 | 144 | 93 |

Every block is built from the actual full-CAR vertex and scalar-number
term. Its complement is a reducing subspace of the same H and has zero
weight for this particular input. It is retained for other inputs. The
block projector is a computational reduction, not an extra proposed gauge
measurement. Physical energy projectors are the full spectral projectors
of H; all energy outcomes outside the reported blocks remain in that
complete apparatus. In particular, directions of energy one and groups
whose measured weight interval includes zero are not removed.

## Exact admission of the complete spectra

The matrices in these blocks are integer and symmetric. The full ground
comparison implies `||H||<9` on the entire carrier. A floating eigensolver
proposes a complete real frame V and a sorted diagonal energy matrix Λ;
both are quantized to dyadic rationals before admission. Exact integer
products calculate

\[
 \epsilon\geq\|V^TV-I\|_2,\qquad R=HV-V\Lambda.
 \tag{4}
\]

Every fixed-width product has an explicit integer-overflow bound checked
before multiplication; residual sums and squares use Python integers.
The symmetric maximum absolute row sum bounds the first norm. The product
of maximum absolute row and column sums bounds the squared residual norm.
No floating eigenvalue, residual or success flag supplies a premise.

For `epsilon<1`, the polar frame `U=V(V^TV)^(-1/2)` is orthogonal and

\[
 \|U-V\|\leq d={\epsilon\over2-\epsilon},\qquad
 \|HU-U\Lambda\|\leq\delta
    =\|R\|+(B+\|\Lambda\|)d.                    \tag{5}
\]

The verifier takes B from the exact maximum absolute row sum of its
admitted matrix, so (5) does not depend on an unchecked norm assumption.
Indeed the singular values of V lie between `sqrt(1-epsilon)` and
`sqrt(1+epsilon)`, and `1-sqrt(1-epsilon)<=epsilon/(2-epsilon)`.
Conjugating H by U and applying the Hermitian eigenvalue bound gives
every sorted exact eigenvalue in `[Lambda_i-delta,Lambda_i+delta]`.
Overlapping intervals are retained as one complete spectral group.

Let J index one such group and P its exact projector. Its rank is |J|,
because its interval is separated from the remaining intervals. Put
`g=min_{i in J,j outside J}|Lambda_i-Lambda_j|-delta>0`. Expanding the
column residual in an exact H eigenbasis gives

\[
 \|P-U_JU_J^T\|\leq
 {\|R_J\|_F+(B+\|\Lambda\|)\sqrt{|J|}\,d\over g}.
 \tag{6}
\]

The spectral subspaces have equal rank, so the projector difference is
bounded by their maximum sine angle. Adding `d(2+epsilon)` bounds
`||P-V_JV_J^T||`. All square roots in these bounds have rational
enclosures. This controls complete projectors, including degeneracies.

The prepared trial has form `sqrt(alpha)u/sqrt(n) -
sqrt(beta)Au/sqrt(nr)`, where `n=||u||^2`, `r=||Au||^2/n`,
`delta_r=sqrt(1+4r)`, `alpha=(1+1/delta_r)/2` and `beta=1-alpha`.
For x and y obtained by applying the same charged probe to u and Au,
respectively, every retained quadratic form is

\[
 {\alpha\over n}x^TMx+
 {\beta\over nr}y^TMy-
 {2\over n\delta_r}x^TMy.                         \tag{7}
\]

The integer/rational terms in (7) are calculated exactly. Its square-root
coefficient is interval enclosed. With `M=V_JV_J^T`, (6) supplies the
projector error; adding the ground-preparation trace distance supplies
the exact-ground weight interval. Multiplication by `multiplicity/18`
then gives its original complete-instrument probability. Every group
is exported, even when its lower probability bound is zero.

For example, the first-family left removal has a certified excitation
interval within `[0.0710017,0.0710521]` and single-mode spectral weight
within `[0.6615519,0.6616432]`. First-family left addition has a line within
`[0.4872697,0.4872936]` with weight within `[0.2473969,0.2474326]`.
These are finite source-energy differences with operator residues, not
one-particle mass shells. The complete output also contains many other
lines, and a largest residue is not a species assignment.

## Moments, time evolution and the original probability sum

Replacing M in (7) by `H^k` gives exact unnormalized spectral moments,
enclosed through order four. The true-ground error is at most D for k=0
and `2D*9^k` otherwise, where D is the retained preparation trace distance.
This follows from the complete bounded observable on the original input,
without dividing by a probe-success probability. The spectral groups also
define the same-H correlation law

\[
 C_i^{\pm}(t)=\langle G|c_i^{\mp}e^{-it(H-E_0)}c_i^{\pm}|G\rangle
             =\int e^{-it\omega}\,d\mu_i^{\pm}(\omega).
 \tag{8}
\]

The weights and frequencies in (8) use one ground, Hamiltonian and clock.
Addition and removal are distinct measures. Their total weights obey the
exact CAR identity `mu_i^+(R)+mu_i^-(R)=1`. Summing the complete
family/multiplet instrument gives one, including all dark and zero-weight
outcomes. Pointer compilation or imperfect preparation adds its complete
instrument error budget; no large compiled primitive tape is asserted.

## Why an energy offset does not make these conjugate particle laws

The full finite singlet reference is not a proven Lorentz vacuum. Its
addition and removal spectra need not be a physical CPT pair. This is
visible in quantities that a number-energy shift cannot repair.
For `H_epsilon=H+epsilon Q*`, the same reference's addition excitation
energies shift by `+epsilon` and removal energies by `-epsilon`. Every
original spectral weight and each conditional variance is unchanged.
The original pair cost and Y are held fixed in this comparison.

All six multiplets have certified unequal addition/removal total weights.
Even after normalization solely for a shape diagnostic, their conditional
variances differ. For first-family left probes,

\[
 \operatorname{Var}(\mu_+)-\operatorname{Var}(\mu_-)
       \in[-0.12243,-0.08993].                    \tag{9}
\]

Thus no offset can make even those two normalized spectral shapes equal.
The mean-aligning offsets also disagree between families: approximately
−0.171 for first-family left and +0.673 for second-family left. These
diagnostics do not replace the original unnormalized records in the
reported instrument.

The Q*=6 reference stays the unique global ground for
`|epsilon|<ground_gap/12`, since `Q*` ranges from −2 to 18 and commutes
with H. Beyond that sufficient window the global ground must be
recomputed. The no-translation statement concerns the same reference;
it is not a theorem excluding every offset-dependent vacuum. It is also
not a claim of physical CPT violation. Identifying observed conjugate
particles requires the spatial dynamics, relativistic transformations
and a physical charge-conjugation map for these same probes.

The specified rational Y is a declared fixture. Family singular-value
transport gives the same computational route for the coherent source
matrix, but its different coefficients require their own ground and
charged comparison. No result here fits a lepton ratio or supplies a GeV
conversion. Bare Y controls the full interaction; a charged readout is
the spectral response of the prepared interacting state, not a bare
coefficient renamed as a mass.
