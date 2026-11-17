# Native-current three-family oscillation records

## Declared coefficient and complete native source

Let Y be the positive real3×3 matrix returned by
`common_source.native_repair_mass.architecture()`. Its entries derive from
the declared five-state record-memory feedback, with native event probability
1/12 and stationary weights proportional to (1,1/12,...,1/12⁴). It is a
stationary mean-coefficient construction, with the inherited controlled
small-step interpretation. No equality to a finite noisy history is assumed.
Let P₂ be the actual native rank-one port projector and declare one coherent
pulse V=I+(−i−1)P₂. The neutral coefficient is

    Mν = V* Y V†.

This is a symmetric, explicitly chosen member of the allowed family-tensor
space. The same native Y supplies the charged Dirac coefficient. The pulse,
unit Higgs component h=(0,1), and equality of charged/neutral singular values
are architecture choices made without a neutrino mass, angle, or mass-splitting
input. This source is distinct from the analyzer-effect charged-lepton mass
architecture used for the observed charged-lepton ratio comparison.

Let D be the actual oriented native charged-lepton Higgs tensor and
ell=LH·h the actual native neutral invariant. The full45×45 internal-family
kernel is

    K = D⊗Y + (D⊗Y)^T + ell ell^T⊗Mν.

The 90-mode spin pairing is epsilon⊗K*, which is antisymmetric. All12
identities G^T LH + LH g_H=0 and the electric Ward identity Q^T K+KQ=0 are
checked exactly. The pairing and kinetic source are the free same-Higgs
quadratic construction; all quark, singlet, charge-conjugate and opposite
helicity sectors are retained. No dynamical Higgs or interacting vacuum is
substituted for this declared background.

## Deriving the current frame

Choose a real ordered orthonormal eigenframe F of Y with distinct positive
roots λ_i. The charged left and right frames are F, including the actual
internal tensor sign. The neutral Takagi frame is VF, since

    (VF)^T Mν (VF) = diag(λ_i).

Transforming the full native weak raising current in this same45-dimensional
basis gives the charged-to-neutral matrix, up to its common internal sign,

    C = F^T V F.

Thus no imported Pontecorvo--Maki--Nakagawa--Sakata matrix enters the
construction. The numerical eigensolver displays this exact spectral
definition and verifies the full-frame intertwiner. Eigenvector signs
change row and column phases only. The matrix predicts the declared model's
own mixing values, not the measured neutrino parameters.

## Complete finite observer instrument

For a positive momentum p and the negative-helicity branch, the free Dirac
or Majorana dispersion gives E_i=sqrt(p²+λ_i²). Its actual left-current
positive-band weight is w_i=(E_i+p)/(2E_i). The source and return analyzer
therefore use

    Kp = diag(sqrt(w_i)) C†,
    F0 = sqrt(I−Kp†Kp), F1 = sqrt(I−Kp Kp†).

The first instrument has outcomes Kp and F0. On success, the branch evolves
with U=diag(exp(−itE_i)); the return instrument has Kp† and F1, followed by
the three orthogonal weak-family record projectors. Both failures retain
their full quantum state. For initial family a the original-input law is

    P_b = |[Kp† U Kp]_(ba)|², b=0,1,2,
    P_fail0 = 1−||Kp e_a||²,
    P_fail1 = ||Kp e_a||²−||Kp† U Kp e_a||².

These five probabilities sum exactly to one. The apparatus may be implemented
by the finite contraction dilation, with its environment and references
retained. It is a declared spectral/current measurement service; no local
hardware implementation theorem is assumed. Spectator sectors are identity
factors in the clock and belong to the complement of the prepared branch.

Writing A=C U C† and r=max_i(1−w_i), commutation of U with the diagonal
weights gives ||Kp†UKp−A||<=r. The successful probability vectors differ in
l1 norm by at most2r. Their lost mass is at most2r. Hence the complete law is
within total variation2r of the ideal three-family oscillation law with
zero failures. The approximation has an explicit positive-momentum domain;
the replay uses p=1,3,10 and four source times for all three initial families.

Conjugating C gives the opposite weak-helicity current law. For the Majorana
coefficient these are helicity channels, not distinct conserved particle
species. A nonzero probability contrast is a source CP diagnostic. Column
Majorana phases cancel from all five records. A zero port-phase control has
C=I and exactly zero cross-family appearance. Floating-point values are
diagnostics of the exact formulas, not interval certificates for a fit.

The construction supplies neither an absolute neutrino GeV scale nor
measured mixing angles or mass splittings. Its parameter-free family
geometry produces a concrete comparison target conditional on the declared
pulse and neutral coefficient.
