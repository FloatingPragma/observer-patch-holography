# Static hypercharge response of the interacting source

The complete interacting source has an exact obstruction to transferring the
free folded-matrix magnetic response to its prepared charged lines. Any
central hypercharge Peierls background in its hopping terms is absorbed by
independent unitaries on the retained quantum links. Every energy eigenvalue
is unchanged, including every degeneracy. This is a result for the declared
finite source and its hypercharge current. It does not identify a physical
electromagnetic field or imply that every fixed-preparation response vanishes.

## Source and background

The Hamiltonian is the same
\[
 H=H_0+\lambda T+4N_{\rm ref},\qquad \lambda=2^{-34},
\]
used by the prepared `ec` and `L` records. It contains 240 fermionic modes,
eight five-state quantum scalars and 24 complete 119-state Peter-Weyl links.
All three local Yukawa coefficients and the scalar and link excitation
penalties equal one. The reference is the retained static fermionic system
at the preparation cell. The underlying exact action is
[`peter_weyl_links/model.py`](../coherent_defect_wave/chiral_joint/peter_weyl_links/model.py);
the coefficients and full prepared-state bounds are inherited from
[`interacting_charged_response`](../coherent_defect_wave/chiral_joint/interacting_charged_response/PROOF.md).

On each oriented link let
\[
 K_e=-G^Y_{e,\mathrm{left}}.
\]
For every native matrix coordinate of that link,
\[
 [K_e,U^{(r)}_{ab}]=Y_r U^{(r)}_{ab}.
\]
Here \(Y_r\) is the actual native hypercharge, without fitting or changing
its normalization. `ec` has \(Y=1\), and `L` has \(Y=-1/2\). Removing a
source fermion and retaining its reference partner gives the opposite
source charge. The identity applies to occupied-to-occupied link
transitions as well as transitions from the blank state.

Declare a real external link angle \(A_e\) by multiplying every actual
forward hopping path in representation \(r\) by
\(\exp(iY_r\sum_{e\in\mathrm{path}}A_e)\), with its conjugate on the
adjoint term. No hopping, reference, scalar or link state is removed.
Set
\[
 V(A)=\exp\!\left(i\sum_e A_eK_e\right).
\]
Different link generators commute. The coordinate identity gives the
finite operator identity
\[
 H(A)=V(A)H(0)V(A)^\dagger. \tag{1}
\]
Each \(K_e\) commutes with the link excitation penalty. The other terms of
\(H_0\), the source spin, and the reference act on different factors.
Hypercharge is central in the native gauge algebra, so \(K_e\) also
commutes with every endpoint Gauss generator. Consequently (1) holds on
the complete physical source-reference sector, with all scalar and link
occupations retained. A fixed weak-isospin direction is not central and
does not satisfy this preservation argument.

Equation (1) permits arbitrary real angles on individual links. It is
stronger than covariance under a gradient of site phases. The executed
profile has angle \(3\pi\) on one edge of an elementary square and zero
elsewhere. Its `ec` Wilson phase is \(-1\). No site gauge gradient has that
plaquette phase.

## Current and static response

For a real profile \(f_e\), put \(A_e=\theta f_e\) and
\(K_f=\sum_e f_eK_e\). Differentiating the actual Peierls Hamiltonian gives
\[
 J_f=\left.\partial_\theta H\right|_0=i[K_f,H],\qquad
 D_f=\left.\partial_\theta^2H\right|_0=-[K_f,[K_f,H]]. \tag{2}
\]
The first derivative convention here defines the current; an opposite
convention changes its sign without changing the conclusion. For every
complete energy projector \(P_E\), including a degenerate projector,
\[
 P_EJ_fP_E=0. \tag{3}
\]
Let \(Q_E=1-P_E\), and let
\(R_E=Q_E(E-H)^{-1}Q_E\) be the reduced resolvent on this finite spectrum.
Expanding both commutators between energy projectors proves
\[
 P_ED_fP_E+2P_EJ_fR_EJ_fP_E=0. \tag{4}
\]
For an intermediate energy \(E'\ne E\), the first term contributes
\(-2(E-E')P_EK_fP_{E'}K_fP_E\); the second contributes its negative.
Thus the static second-order response cancels exactly. Equation (1)
proves the stronger all-orders statement: every source energy and every
spectral multiplicity is independent of the background. No spectral
approximation or restriction to the charged energy bin is used.

The current itself is nonzero. Literal occupied-column controls evaluate
the first and second background derivatives independently from the
diagonal link conjugation and check (2). Omitting either derivative
fails these controls. The cancellation in (4) includes the diamagnetic
term and the complete complementary spectrum.

## Complete prepared records

Conjugate the complete inherited preparation instrument by \(V(A)\) on
each retained source copy, retaining every failure output and apparatus
reference. The native removal pump acts on source matter and the static
reference, so it commutes with this link-only unitary. Functional calculus
and (1) give
\[
 P_I(H(A))=V(A)P_I(H(0))V(A)^\dagger,\qquad
 e^{-itH(A)}=V(A)e^{-itH(0)}V(A)^\dagger.
\]
The interval \(I\), reference cost, source clock parameter and all energy
bin complements are preserved. Source spin and the existing location
effects commute with \(V(A)\). The paired experiment therefore has the
same complete joint record law for both prepared charged families. This
statement includes preparation failure, pump failure, energy-bin
complement, nonorigin/cloud records and the retained quantum systems.

For the atomic projector \(P_0\), every link is blank and \(V(A)P_0=P_0\).
Thus any fixed input supported in \(P_0\), followed by the native pump,
the background Hamiltonian and the existing complete arrival effects,
has background-independent arrival probabilities. This stronger
fixed-input statement uses exact blank-link support. It is not assigned
to the inherited energy-filtered interacting line.

A fixed preparation need not be its paired preparation. The code gives
an exact physical control state formed from the native local lepton-pair
pulse plus its image under one actual hopping edge. Its Gauss constraints
are checked directly. Its unnormalized squared norm is \(9/8\). Holding
the normalized coherent state fixed changes its energy mean by
\(-\lambda/3=-1/51539607552\) at the test background; transforming the
state by \(V(A)\) restores the original energy mean exactly. This control is not substituted for either
inherited charged line. Equation (1) alone supplies no assertion about
all fixed-state transients or all magnetic observables.

## Link action and the magnetic-response boundary

The declared \(H_0\) contains a link excitation penalty and no bare link
plaquette term. For the closed `ec` loop \(W\) around the test square,
\[
 V(A)WV(A)^\dagger=-W.
\]
The code checks all affected Gauss generators on
\((W+W^\dagger)|\mathrm{blank}\rangle\). This vector has squared norm two;
the difference under conjugation has squared norm eight. Adding a fixed
nonzero coefficient of \(W+W^\dagger\) to the source would therefore
invalidate the step establishing (1). This comparison operator is absent
from the executed source. A consequence derived through virtual hopping
from the unchanged source inherits (1); naming it an effective plaquette
term does not remove the conjugation identity.

The explicitly different longitudinal Coulomb construction uses
\[
 H'=H_{\mathrm{local}}+\lambda T+4N_{\rm ref}
       +\frac{\kappa}{2}\rho^*G\rho,\qquad \kappa>0,
\]
where the density \(\rho\) includes the actual matter, scalar and
reference hypercharges and \(G\) is the finite graph pseudoinverse.
Its implementation is
[`openwave_current_coulomb/source.py`](../openwave_current_coulomb/source.py).
Every \(\rho_x\) acts on factors distinct from the quantum links, so
\([K_e,\rho_x]=0\). The density-quadratic energy therefore commutes
with \(V(A)\), and
\[
 H'(A)=V(A)H'(0)V(A)^\dagger.
\]
The static current identities (2)–(4) hold with \(H'\) as well.
Replacing the unit-link binding cost by this longitudinal energy does
not remove the static hypercharge spectral obstruction. The executor
checks the actual alternative column action on 21 complete source and
reference inputs, including vacuum, single-particle and filled-reference
occupations. It also exhausts all 119 link labels at fixed density.
The operator proof holds for every positive \(\kappa\), independently
of these fixtures. It transfers no ground, physical preparation,
charged energy bin or transport probability from \(H\) to \(H'\).

The folded 16-component free-envelope calculation of a response tending
to \(g=2\) uses a prescribed classical connection and a Clifford mass
operator. It has no transfer theorem to this interacting state with its
retained quantum links. The exact static hypercharge splitting here is
zero. Calling that result an electromagnetic \(g\) would require a photon
direction, a physical field intervention and an energy-to-mass relation
on this same source. Those identifications are not made by this theorem.

The separate
[`intrinsic magnetic construction`](../openwave_intrinsic_magnetic/PROOF.md)
uses the positive kernel \(|D_W(A)|\), a confining scalar mass well,
a Wilson regulator and a semiclassical transverse Maxwell field with
energy proportional to \(\|CA\|^2\). Its stationary numerical current
and self-consistent field belong to that stated construction. Its
classical connection and Maxwell energy do not satisfy the independent
PW119 link-conjugation premises used here. Neither its preparation
nor its intrinsic dipole is transferred to the interacting charged line
or to the longitudinal Coulomb alternative.

`magnetic.py` executes all 59 full coordinate matrices, their 1,156
nonzero entries, central endpoint-current commutators, occupied 240-mode
Hamiltonian columns and the nonzero-flux controls. The tests also check
complete axial paths at side four. The analytic operator identities carry
the all-state statement; the basis controls verify the literal source
implementation without materializing its full matrix.
