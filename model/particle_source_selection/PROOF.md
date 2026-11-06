# Classical records and coherent family phases

The source is the complete record instrument in
`model/common_source/native_family_record_memory.py`. It retains 128 raw
labels, twelve old-memory inputs and their work-state maps. Its fair-coin
event rates are (p=1/12), (q=1/6), with joint rate (pq=1/72). Therefore

\[
2\sqrt{p/q}=\sqrt2,\qquad p+q-2pq=2/9.
\]

These identities are exact. Their interpretation as a coherent modulus and
an oriented phase requires additional structure. The following two arguments
settle that requirement for functions of the complete classical record
instrument and for cyclic averaging of its real family kernel. They do not
exclude a source law with an independently specified coherent reference.

## Phase nonidentifiability from complete classical records

Let a finite instrument have maps

\[
\mathcal I_r(\rho)=\sum_a K_{ra}\rho K_{ra}^{\dagger},\qquad
\sum_{ra}K_{ra}^{\dagger}K_{ra}=I.
\]

For arbitrary real numbers \(\theta_{ra}\), replace
\(K_{ra}\) by \(e^{i\theta_{ra}}K_{ra}\). Every map \(\mathcal I_r\)
is unchanged, as is the complete classical-quantum output
\(\sum_r |r\rangle\langle r|\otimes\mathcal I_r(\rho)\).
This is equality on the full matrix space, so it holds with arbitrary
entangled references. Repeated composition and any feedback depending on
retained classical records preserve this equality outcome by outcome.

Consequently, a well-defined function of these instrument maps cannot depend
on the free Kraus phases. In particular, a coherent sum
\(\sum_{ra}c_{ra}K_{ra}\) does not follow from the instrument alone when its
value changes under such a phase transformation. This is an insufficiency
theorem for the specified input data, rather than a claim that the physical
source cannot carry coherent phase information.

A source-selected Stinespring isometry together with an accessible coherent
environment reference can contain that information. Measurements of
off-diagonal environment observables can distinguish the phased isometries.
Those observables are absent from a committed classical record. In the
native family instrument, old memory is explicitly discarded and record
coherences are not retained. Recovering their phases from event probabilities
alone is therefore not a valid inference.

A convention can choose one representative Kraus family or set a numerical
phase equal to a function of the rates. The theorem permits such a choice;
it proves that the classical records do not establish its physical
identification or distinguish it from the other coherent lifts. In
particular, setting a probability equal to an angle in radians is an
additional rule. If a construction is equivariant under complex conjugation
and its input matrices are real, its output must be real. Real principal
matrix square roots, real matrix products and Reynolds averages respect this
condition. A nonreal chiral response therefore requires an input or a
selection rule that supplies the missing oriented complex structure.

## Exact counterfamily with the native event rates

The following auxiliary construction uses the same numerical rates \(p,q\),
with a regular three-cycle \(R\), an identity branch and coherent record erasure.
It is a mathematical counterfamily, not an extension or intertwining of the
native 128-record instrument and not an implemented native controller.
Set

\[
K_0=\sqrt q I,\quad K_1=\sqrt p\,zR,\quad
K_2=\sqrt p\,\overline zR^2,\quad
K_3=\sqrt{1-q-2p}\,I,\qquad |z|=1.
\]

The fourth branch retains the probability omitted by the first three, and
\(\sum_jK_j^\dagger K_j=I\). All four retained classical outcome maps are
independent of \(z\). A unitary Fourier measurement on the first three
coherent environment coordinates has one outcome with Kraus operator

\[
L_0=C(z)/\sqrt3,\qquad
C(z)=\sqrt q I+\sqrt p(zR+\overline zR^2).
\]

The other two Fourier outcomes and the fourth original outcome are retained;
the unread channel is unchanged. The selected outcome is a filter, and its
success probability is \(\operatorname{tr}(C^2\rho)/3\). Its probabilities
and spectrum depend on \(z\). Assigning its eigenvalues to square-root masses
would be an additional physical interpretation.

For both \(z=1\) and \(z=(99+20i)/101\), the certificate proves strict
positivity of \(C\) by exact leading principal minors. The first spectrum has
a repeated eigenvalue and the second has three distinct eigenvalues. Both
satisfy

\[
\frac{\operatorname{tr}C^2}{(\operatorname{tr}C)^2}
=\frac{q+2p}{3q}=\frac23.
\]

The complete classical instrument, exact event rates and this balanced
invariant therefore fail to select the phase or spectrum, even after
granting the cyclic carrier. The rational unit phase is an algebraic control;
it was not chosen against particle data. This counterfamily proves no
selection rule from the classical instrument can distinguish these lifts.
It does not equate the coherently measured instruments, which are different.

## Real cyclic averaging has a repeated eigenvalue

Let \(R\) be a real orthogonal regular order-three action on a real
three-dimensional carrier. Its invariant projector is
\(P=(I+R+R^2)/3\), of rank one. The complement is a real two-dimensional
rotation plane. For every real symmetric \(X\), its Reynolds average is

\[
\overline X=\frac13\sum_{j=0}^{2}R^jX(R^T)^j
=aP+b(I-P),\qquad
a=\operatorname{tr}(PX),\quad b=(\operatorname{tr}X-a)/2.
\]

To see completeness, a symmetric matrix commuting with a nontrivial planar
rotation is scalar on that plane, and the singlet-to-plane blocks vanish.
Thus the characteristic polynomial is \((x-a)(x-b)^2\). This holds for
every real symmetric input, rather than a finite menu of examples.
Functional calculus of this averaged operator preserves the repeated pair.

The native stationary projector average \(Y\) is real symmetric. Exact
replay proves that its discriminant is nonzero and that it commutes with none
of the twenty inherited order-three rotations. Averaging \(Y\) over any one
of these cyclic actions produces a repeated pair, so averaging cannot give a
three-distinct-eigenvalue charged-family operator. Complex Hermitian
circulants can have three distinct eigenvalues; the second counterfamily
member verifies this distinction explicitly.

## Scope and finite-stack reuse

The exact positive-circulant invariant is the algebra used in
`Lean/ObserverPatchHolography/KoideCirculant.lean`. It fixes a modulus under
its positive-chamber and mass-interpretation premises, without selecting a
phase. The source rotation and frame come directly from
`native_gauge_object.geometry()`; no numerical fitted basis is introduced.
The stationary weights are reconstructed from all native record weights,
rather than read from a stored result.

The phase argument applies to all finite retained classical instruments and
classical-record feedback. It does not exclude coherent source extensions,
phase-sensitive initial preparation, a selected environment reference,
different matter carriers, or a different source-to-mass map. The twirl
argument applies to real symmetric kernels on a real regular three-cycle;
it does not apply to a complex Hermitian response with an imaginary
antisymmetric component. No measured spectrum, continuous source phase,
absolute mass unit or cosmological population distribution follows here.
