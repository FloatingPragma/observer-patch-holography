# Exact flavour freedom and a complete current-loop readout

This construction uses the actual fifteen-component native current tensors,
their three charged Higgs lines and the optional same-Higgs LH invariant.
The conditional action supplies three family copies, canonical positive
kinetics and physical CAR statistics. Those architectural choices do not
derive a physical family assignment or a chiral quantum measure.

## Family-basis equivalence on the complete finite action

The five inequivalent internal multiplets Q,u^c,d^c,L,e^c each have a
three-dimensional family factor. Every native current acts as G⊗I3.
Block-diagonal family unitaries commute with all twelve full 45-coordinate
currents and preserve the identity kinetic metric. The executable checks
all 540 commutators of the 45 anti-Hermitian family generators with those
currents. Their full CAR lift is a unitary canonical basis change on every
occupation sector, not a truncation to one-particle states. Preparations,
interventions, quantum references and detectors must be transported with
the action. A fixed detector's entries need not be invariant.

In four-component family notation these coordinate changes act as

    Yu -> UQ† Yu Uu,  Yd -> UQ† Yd Ud,  Ye -> UL† Ye Ue.

The all-left conjugate singlets use the complex-conjugate right matrices.
The symmetric neutral coefficient transforms as Mnu -> UL^T Mnu UL.
For anti-Hermitian infinitesimal generators A, the exact tangent equations
are

    dYu=-AQ Yu+Yu Au, dYd=-AQ Yd+Yd Ad, dYe=-AL Ye+Ye Ae,
    dMnu=AL^T Mnu+Mnu AL.

These are proved basis equivalences of the finite canonical action. They
do not discard a physical strong-CP parameter or establish invariance of
an unconstructed chiral functional measure. Gauge theta/global-group
choices are separate unspecified inputs outside this flavour count.

## Exact quotient dimensions

Three complex charged matrices contain 54 real coordinates. The family
unitary group has dimension 45. For generic nonsingular distinct masses,
the quark stabilizer is the common baryon U1. The charged-lepton matrix
has three singular-vector U1 stabilizers. Thus the orbit dimension is
45-4=41 and the quotient dimension is 13: nine singular values, three
quark mixing angles and one quark CP phase. The executable constructs the
complete 54-by-45 real tangent matrix at an exact rational complex point,
proves rank 41 and supplies its four independent null directions. The
explicit stabilizers give the upper bound; the attained rank proves the
generic count, without claiming that degenerate strata have that dimension.

A complex symmetric three-by-three Weinberg coefficient adds twelve real
coordinates. Generic charged and neutral matrices remove the three lepton
phase stabilizers, leaving baryon U1. The complete 66-by-45 tangent matrix
has exact rank 44. Hence the generic quotient has 22 coordinates: the
charged thirteen, three neutral singular values, three lepton mixing angles,
one lepton Dirac phase and two Majorana phases.

Rank two of a generic complex symmetric matrix is complex codimension one,
or real codimension two. The executable also verifies orbit rank 44 at an
exact rank-two neutral point. Its quotient therefore has 64-44=20 real
coordinates. The zero singular value removes one mass and one Majorana
phase, not only one mass. With a massless first column, multiplying both
massive PMNS columns by the same unit phase p multiplies Mnu by p*².
The common charged-lepton basis change UL=Ue=pI restores Mnu and Ye.
Only the relative phase of the two massive columns survives. This is
checked exactly with p=(3+4i)/5.

The actual normal-ordering calibration packet therefore consumes nineteen
independent empirical flavour coordinates: nine charged singular values,
four CKM parameters and six oscillation inputs. Its rank-two architecture
also chooses one independent unmeasured Majorana phase. Two printed zero
Majorana phase entries include a redundant common phase at rank two.
The zero lightest mass and ordering are architectural choices; their source
derivation, scale running and matching error are unproved. The Higgs norm
and GeV unit are additional inputs outside the quotient count.

## What the native identities constrain

For each charged channel, let C0,C1 be the literal oriented fifteen-row
Higgs tensors. Their actual internal current equations are verified before
lifting. For every complex family matrix Y, the full coefficient is

    K_a = C_a⊗Y* + (C_a⊗Y*)^T.

Since the currents act as G⊗I3, the internal Ward identity factors for
arbitrary Y. No family entry is restricted by this equation. The complete
27-column coefficient embedding has a positive diagonal Hilbert-Schmidt
Gram, so all 27 complex charged directions are independent. Each complete
45-row neutral-background coefficient obeys transpose symmetry and electric
neutrality; epsilon_spin⊗K is antisymmetric on all 90 CAR modes. The
three-channel complex counterexample is also checked directly against all
72 full-width current/Higgs identities. This is a full-operator statement,
not an equality of compressed energy tables.

The actual LH invariant obeys all twelve native current identities. Its
square tensors with any complex symmetric family coefficient, leaving six
complex neutral directions. The native charge, Higgs and canonical kinetic
identities select allowed operator types; they impose zero additional linear
constraints on those family coordinates. A specified source dynamics can
impose further equations, which must be evaluated separately.

## A source equation, its boundary and an exact counterexample

The existing source trial literally returns Yu=a F3, Yd=b F3 and Ye=J F3.
Executing those producer functions proves

    [Yu Yu†, Yd Yd†]=0.

With nondegenerate singular values this gives aligned family frames and
zero quark Jarlskog invariant. It is a consequence of that shared-feature
trial. It does not follow from the native gauge/Higgs identities alone.

An exact admitted counterexample takes Yu=diag(1,2,4),
Yd=V diag(2,3,5), Ye=diag(1,3,7), with all three mixing sines
3/5,5/13,8/17 and CP phase pi/2. These values are algebraic test settings,
not empirical mass or angle choices. The full native tensors obey every
required identity while [Yu Yu†,Yd Yd†] and

    J=Im(V11 V22 conjugate(V12 V21))

are nonzero. The exact determinant and J are retained. The complete native
45-coordinate mass-frame transformation gives V in each of the three
colour copies of the actual charged weak current. Thus the counterexample
and readout amplitudes are attached to the source current, not inserted
as an unrelated mixing table.

## One extra observable for one CP direction

Let a=V11 V22 and b=V12 V21. A separate declared two-path filter service has
successful map diag(a,b) and failures diag(sqrt(1-|a|²),0) and
diag(0,sqrt(1-|b|²)). Their squared effects sum to I. Prepare the path
qubit in its plus state and retain both Pauli-Y outcomes after success.
The complete four probabilities are

    P(Y+)=(|a|²+|b|²)/4-J/2,
    P(Y-)=(|a|²+|b|²)/4+J/2,
    P(fail0)=(1-|a|²)/2,  P(fail1)=(1-|b|²)/2.

They are nonnegative and sum to one, so J=P(Y-)-P(Y+) uses all attempted
trials without postselection. Independent row and column rephasings multiply
a and b by the same phase, leaving J unchanged. Complex conjugation preserves
every singular value and |Vij|², but reverses J and exchanges the two Y
probabilities. Incoherent current counts alone cannot distinguish this pair.
A source-derived value of this one signed contrast would constrain one
remaining real CP direction when the three nonzero mixing angles are fixed.
The chosen common-feature trial predicts zero, whereas the exact admitted
counterexample predicts a nonzero value.

The literal current readout is a separate complete instrument. Use two distinct
colour-register copies of the native charged-current block. Path zero
prepares input families (0,1) and path one prepares (1,0). Apply V⊗V and
record every output family pair (r,s), followed by both path-Y outcomes.
Writing a_rs=Vr0 Vs1 and b_rs=Vr1 Vs0, its eighteen Kraus rows are

    K_(rs,+)=(a_rs,-i b_rs)/sqrt(2),
    K_(rs,-)=(a_rs,+i b_rs)/sqrt(2).

Unitarity gives sum K†K=I2, checked exactly on every path coordinate;
therefore arbitrary untouched references are retained. For the plus path
preparation, probabilities are |a_rs∓i b_rs|²/4. Their sum is one, and
J=P((0,1),Y-)-P((0,1),Y+) without conditioning on that output pair.
All sixteen nontarget records are kept, including any zero values. The
colour labels distinguish the two occupied mode registers, so no identical
single-mode occupation or extra particle family is introduced. Both paths
make the same total charge transfer; any compensating apparatus reference
is shared and is not a flavour-dependent calibration.

The four-record filter and eighteen-record current instrument share the
specified target signed contrast. Neither equality of their complete
path/reference channels nor a classical coarse-graining relation between
them is asserted.

These are exact finite instruments with declared coherent services. This construction
does not claim native pulse compilation, laboratory weak-interaction
accuracy, physical mixing parameters or source-derived empirical Yukawas.
The analogous neutral phase distinction is explicit: Majorana column
phases cancel from U diag(m²) U†, while the remaining relative phase changes
|Mnu_ee|² in the charged-lepton mass basis. Oscillation records alone do not
constrain that lepton-number-violating direction.
