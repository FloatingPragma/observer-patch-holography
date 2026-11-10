# Exact charged energy denominator

## Source and comparison domain

The source is the complete 30-mode fermion Fock space, five-state quantum
scalar and 119-state compressed Peter–Weyl links at the unit local and link
coefficients. The controlled comparison is
\(H_\lambda=H_0+\lambda T\), where \(T\) is its complete native hopping.
The [atomic control](../atomic_reference_control/PROOF.md) proves the
physical atomic ground projector \(P_0\), the positive gap \(g\),
\(P_0TP_0=0\), and the full second-order operator

\[
 K_2=-P_0TQ[Q(H_0-VE_0)Q]^{-1}QTP_0,
 \qquad Q=I_{\rm physical}-P_0.
\]

Every distinct edge contribution is retained. On the two-cell-per-axis
source the 24 kinetic paths contain one distinct link each. After one hop,
their excited-link sectors are mutually orthogonal and are preserved by
\(H_0\), so all cross-edge terms in this second-order expression vanish.
The complete contribution of an edge with axis \(a\) is

\[
 K_{2,e}=-A_eI+J_e\sum_i R_{a,ii}S_{x,i}S_{y,i},
 \qquad
 J_e=\frac14\sum_f d_f\int_0^\infty e^{-t}b_f(t)^2\,dt.
\]

Here \(R_a\) is the spin rotation by \(\pi\) about axis \(a\). The
factor \(1/4\) is \(4|c_e|^2\), with the actual \(|c_e|^2=1/16\).
The internal sum contains the blank-link normalization \(1/d_f\) from
the full Peter–Weyl coordinate. Summing its two internal indices gives
\(d_f^2/d_f=d_f\). Both flux orientations and their fermion signs are
included in the inherited exact contraction. All quark channels remain in
the scalar coefficient \(A_e\). No Hamiltonian coefficient changes in
this calculation.

The local semigroup is the complete unrestricted operator
\(E(t)=e^{-t(h-E_0)}\), including charged local sectors. On the physical
ground spin doublet,

\[
 P_g c_{s,a}^\dagger E(t)c_{u,b}P_g
 =\delta_{ab}\left[a_f(t)\delta_{su}I+
 b_f(t)\sum_i\sigma_{i,us}S_i\right].
\]

Hermiticity makes \(a_f,b_f\) real. Actual source charge conjugation
gives the addition response with the corresponding transposed Pauli
indices and minus sign, as specified and tested by the atomic control.
The pair denominator has the strict lower bound
\(1+2h_{\min}-2E_0\geq1-2\delta>0\). Thus each Laplace integral and
every pair resolvent below converges on its complete charged domain.

## Exact quark response cancellation

Let \(S_\ell\) and \(S_q\) be the native total lepton and quark spin
currents. Every local Yukawa term pairs two fermions into a spin singlet
within one of these two sectors. `model.separate_spin` verifies all six
full normal-ordered identities \([h,S_{\ell,i}]=[h,S_{q,i}]=0\), with
all 30 fermion modes and all scalar states. It also verifies the 90
opposite-sector field commutators, including
\([c_{q,s,a},S_{\ell,i}]=0\).

The exact local ground belongs to the unique maximal quark module of
Casimir 54. Its multiplicity is one and the commuting quark-spin
representation is trivial. Consequently \(S_\ell=S\) on the ground
doublet. Every quark response
\(P_gc_{q,s,a}^\dagger E(t)c_{q,u,b}P_g\) commutes with this irreducible
spin one-half representation. Schur's lemma makes it scalar on the ground
doublet for every \(t\geq0\), hence

\[
 b_{u^c}(t)=b_Q(t)=b_{d^c}(t)=0.
\]

This identity makes only the spin-dependent quark response vanish. Quark
occupations, all quark Yukawa couplings and the quark contribution to the
scalar energy are retained. The scalar and quark dynamics also remain in
the lepton charged resolvent.

## Complete charged lepton blocks

The particle-hole transformation, actual maximal quark module and local
ground are the frozen [local reference](../local_reference/PROOF.md).
Its carrier is \(140\times64\times5=44800\) states. An elementary
lepton field acts trivially on the quark factor and therefore preserves
this module exactly. The full source ground has zero amplitude in every
other quark module, and local evolution preserves these modules.

The source highest-spin block has dimension 16; its ground projector is
rank one. A native \(e^c\) annihilator reaches transformed lepton number
four, hypercharge \(-1\), weak spin zero. A native lower-component
\(L\) annihilator reaches transformed number two, hypercharge \(+1/2\),
and highest weak weight \(+1/2\). Starting from source spin \(+1/2\),
native up and down removal reaches magnetic spin weights zero and one.

`model.target` retains every coordinate of each charge and magnetic weight,
then computes the exact kernel of the weak raising current. Its dimension
is independently checked by adjacent weight multiplicities. The dimensions
are 17 and 7 for \(e^c\), and 33 and 12 for \(L\). Complete unreduced
weight spaces have dimensions 135, 48, 134 and 45, respectively. The actual
full CAR map, including particle-hole signs and scalar identity, is
evaluated before projection. Exact reconstruction of its whole image
proves zero transfer into every excluded weak representation. The other
weak component follows from the exact native weak covariance.

Quarks precede leptons in the reduced tensor convention. Their transformed
occupation is 12, so the odd lepton Jordan–Wigner factor is \(+1\).
All scalar occupations and all six-mode lepton occupations in each target
weight are retained. No vacuum approximation is made. For each target,
the full Hamiltonian intertwiner, its physical Gram metric, and the
metric adjoint of the elementary transfer are checked exactly.

Let \(D_s\) be these complete transfer matrices, \(h_s\) their target
generators, and \(h_g\) the 16-dimensional source generator. The signed
response is

\[
 b_f(t)=e^{E_0t}\operatorname{Tr}\left[
 P_g\{D_\uparrow^\ddagger e^{-th_\uparrow}D_\uparrow
       -D_\downarrow^\ddagger e^{-th_\downarrow}D_\downarrow\}\right],
\]

where \(\ddagger\) denotes the actual Gram adjoint. This follows by
taking the two native spin expectations in the ground-up state; their
difference is \(b_f(t)\), with no factor of two. The zero-time moments
are independently checked against \(\beta_{e^c}\) and
\((1-\beta_{e^c})/2\) from the actual ground spin current.

## Algebraic ground projector and pair resolvent

Write \(e=E_0\). Its minimal polynomial is

\[
 f(e)=[e(e-1)]^3-313[e(e-1)]^2+31708e(e-1)-1041468,
\]

and its selected embedding is the smallest real root. All calculations
take place in the exact number field \(\mathbb Q[e]/(f)\). For the
squarefree source characteristic polynomial \(p\),

\[
 P_g=\frac{p(h_g)/(h_g-e)}{p'(e)}.
\]

Polynomial division defines the numerator without a singular matrix
inverse. Its trace is checked to be one. The complete target squarefree
annihilators \(q_f\) have degrees 9 and 17. Their action on every target
state is checked directly. They include dark states and spectral
degeneracies; no per-vector eigensolver output is used.

The signed moments
\(m_j=\operatorname{Tr}P_g(D_\uparrow^\ddagger h_\uparrow^jD_\uparrow
-D_\downarrow^\ddagger h_\downarrow^jD_\downarrow)\),
for \(0\leq j<\deg q_f\), determine the exact rational response

\[
 R_f(z)=\frac{N_f(z)}{q_f(z)}
       =\sum_r\frac{w_r}{z-\epsilon_r},\qquad
 b_f(t)=e^{et}\sum_rw_re^{-\epsilon_rt}.
\]

The sum is over spectral projectors, including degeneracy, and its weights
may have either sign. It introduces no basis-dependent amplitudes.
To evaluate the pair denominator, define

\[
 g_f(z)=-R_f(2e-1-z).
\]

The strict pair gap implies \(q_f(2e-1-z)\) is invertible modulo
\(q_f(z)\). `pair_functional` computes this inverse in the exact number
field and verifies both its multiplicative identity and

\[
 q_f(2e-1-z)g_f(z)+N_f(2e-1-z)=0\pmod {q_f(z)}.
\]

For the resulting polynomial representative \(g_f(z)=\sum_jg_jz^j\),
applying the signed moment functional gives

\[
 I_f=\sum_jg_jm_j
 =\sum_{r,s}\frac{w_rw_s}{1+\epsilon_r+\epsilon_s-2e}
 =\int_0^\infty e^{-t}b_f(t)^2\,dt,
 \qquad J_e=\frac{I_{e^c}+2I_L}{4}.
\]

This evaluates the actual denominator; it does not replace it by an
average excitation energy. The receipt records its exact degree-at-most
five algebraic coefficient vector and an outward-rounded rational
interval obtained from the certified ground-root isolation. A separate
two-pole exact test checks the reflection sign and both spectral sums;
singular and dimension-inconsistent inputs fail closed.

## Scope and replay

The scalar shift \(-A_eI\), every edge, and every complementary state
are retained by the inherited full effective operator. The canonical
effective-Hamiltonian error remains
\(6(|\lambda|\tau)^3/g^2\) on the stated atomic-control interval.
This packet sharpens \(J_e\) within that same operator and does not
extend the admitted hopping range by itself. It proves no unit-hopping
ground, thermodynamic phase, relativistic mass or empirical species
identity.

The mutually exclusive write and check modes bind the complete imported
Python closure, the package's own code and proof, and validated immutable
provider receipts. Check mode recomputes the exact identities and values
and compares the whole receipt without writing. No numerical eigensolver
or floating diagnostic defines a certified result.
