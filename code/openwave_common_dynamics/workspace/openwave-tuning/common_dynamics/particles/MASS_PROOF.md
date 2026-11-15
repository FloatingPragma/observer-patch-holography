# A source-local massive amendment

The frozen source is W=SR, where S is the native holding/transport involution
and R=2AA*−I is the reflection about the uniform twelve-event state at each
port. This construction proposes one additional parameter m and the word

\[
W_\mu=SR\exp(-i\mu R),\qquad \mu=hm,\qquad |\mu|\le1/2.
\]

Here h is the spacing used by the common clock t=nτ₀h. The mass coefficient
m is an architecture input. This proposal does not derive an observed
particle mass, choose a physical vacuum, or change predecessor receipts.

## Exact finite spectrum and all sectors

Let f be an eigenvector of T=A*SA with eigenvalue t=1−λ/12. The span of
Af and SAf is invariant. In that ordered, generally nonorthogonal basis,

\[
W_\mu=\begin{pmatrix}0&-e^{i\mu}\\e^{-i\mu}&2t\cos\mu\end{pmatrix},
\qquad G=\begin{pmatrix}1&t\\t&1\end{pmatrix}.
\]

Its determinant is one and its eigenvalues are e^{±iω}, with
cosω=t cosμ. For t=1 the two-coordinate basis degenerates: its actual
one-dimensional uniform state has phase e^{-iμ}. It cannot be counted twice.

On the finite 96-port source, the physical 1152-mode spectrum is:

| Subspace | Phase | Multiplicity |
|---|---|---:|
| Constant Af | e^{-iμ} | 1 |
| Divergence-free antisymmetric moving cycles | e^{iμ} | 193 |
| Remaining R−, S+ subspace | −e^{iμ} | 768 |
| Nonconstant transfer modes | e^{±i arccos(t cosμ)} | 95 pairs |

There are 288 moving edges. The antisymmetric cycle dimension is
288−96+1=193. The S+ dimension is 576 holds+288 symmetric edges=864;
its R− intersection has dimension 864−96=768. These counts sum to1152.
No stationary or zone-edge modes are discarded.

The encoded 2048-dimensional source also retains896 padding dimensions.
The32 invalid addresses have uniform-event phase e^{-iμ}; their remaining
480 directions and the384 valid-address event-padding directions have
phase −e^{iμ}. The executable checks all invalid-address uniform columns
and event-padding representatives in addition to the complete16-event coin.
Padding is a decoder outcome, not a physical species.

Quasienergy always means −Arg(Wμ)/(τ₀h), with the common energy contract's
branch convention. The phases at −1 are retained with their principal-zone
alias. For μ>0 the cycle energy is −m/τ₀, while −e^{iμ} has energy
(π−μ)/(τ₀h). A positive rest frequency is not a proof that the empty Fock
state is the energy ground state.

## Native word and charge conservation

Let U prepare the source's uniform twelve-event vector from event0. The
entire modified coin is compiled as the single source word

\[
U\operatorname{diag}(e^{-i\mu},-e^{i\mu},\ldots,-e^{i\mu})U^*.
\]

The emitted diagonal angles are [−μ,π+μ,…,π+μ] on all16 encoded events.
The common swap word is unchanged. Keeping its scalar phase is necessary
when controlling the whole word. At m=0 this recovers the full base coin,
including padding and its global phase.

`mass_coin.native_error` checks the emitted Pauli-angle templates against
exact rational μ and certified π intervals. It adds the native control
errors, two source-preparation rotation errors, and the frozen swap error.
The numerical column residual is a control, not the error certificate.

The Fock lift Γ(exp(−iμR)) is constructed from the same onsite Givens word
and16 occupation phases [−μ,+μ,…,+μ]. The vacuum phase is one. It preserves
the total occupation at every port, so it commutes with both modular and
integer Gauss operators. Spectator parity strings remain in the Givens
implementation. `massive_coupled.py` emits this word for all128 addresses,
including padding, and composes its full-sector operator error with the
actual autonomous field source. It does not replace the source with a
separately written Hamiltonian.

## Complete continuum generator

Use the frozen continuum projector P=ker(W(0)−I), of rank26, and its three
source-derived matrices Hₐ=gₐe*+egₐ*. They obey gₐ*gᵦ=δₐᵦ/48. The integer
numerator checks in `mass_coin.py` establish

\[
RP=PR=2ee^*-P,\qquad RH_a+H_aR=0.
\]

The first derivative of the amended word gives, on all of P,

\[
K_m(p)=i\sum_a p_aH_a-imRP.
\]

Thus its signed energy matrix is [−ΣpₐHₐ+mRP]/τ₀. On the scalar and
longitudinal vector directions it has energies

\[
E_\pm=\pm\sqrt{m^2+|p|^2/48}/\tau_0.
\]

The other24 directions have energy −m/τ₀: two transverse directions at
nonzero momentum and22 additional stationary directions. At p=0 there is
one +m branch and25 −m branches. This is not a two-component Dirac field.
The inherited scalar/vector rotation representation, investigated by DYN10,
does not become spin one half because a mass term was added.

The coefficient of the native displacement projection is still PHₐP via
the predecessor identity PJₐP=−Hₐ. The mass term changes the state dynamics;
it introduces no new current readout. The resulting scalar amplitude obeys
the Klein–Gordon dispersion on the admitted propagating subspace. All other
components and their records remain part of the source.

## Uniform complete-state limit

Let |p|∞≤K, hK≤1/2, and |hm|≤1/2. Initial states lie in P, possibly
entangled with an arbitrary unchanged reference. Let X be the same optical
corrector used in the frozen massless proof. Because [R,P]=0, the new
mass term has no Q-to-P first-order block; the defining equation for X
continues to hold with K replaced by K_m. Put

\[
x=5K/4,\quad L=K+|m|,\quad b=K/4+|m|.
\]

The derivative of W(hp)e^{-ihmR} has norm at most L; its quadratic remainder
has norm at most h²L²/2. Also ||X||≤x and ||K_m||≤b. For J_h=P+hX,

\[
\|W_\mu(hp)J_h-J_he^{hK_m}\|
\le h^2\{x(L+b)+(L^2/2)(1+hx)+b^2/2\}
\le h^2 C_m,
\]

where

\[
C_m=x(L+b)+13L^2/16+b^2/2.
\]

The last step uses hx≤5/8. Unitary telescoping and both endpoint correctors
give the complete-state bound

\[
\boxed{\|W_\mu(hp)^nP-Pe^{nhK_m(p)}\|
\le(5/2)hK+nh^2 C_m.}
\]

It controls every terminal instrument, including source/event outputs and
an arbitrary reference, in total variation after capping at one. There
is no success normalization or optical projection of the actual output.
At m=0, C_m=77K²/32 exactly, recovering the frozen massless theorem.
For fixed t=nτ₀h the bound tends to zero linearly in h. Numerical controls
compare complete144-component matrices, including p=0, with this inequality.

For the acoustic branch λ(hp), the inherited bound gives
λ≤15h²K²/32 and |λ−h²|p|²/8|≤3h⁴K⁴/8. Elementary cosine remainders,
using cosω=(1−λ/12)cos(hm), yield

\[
\left|(\omega/h)^2-m^2-|p|^2/48\right|
\le h^2\left\{\frac34(m^2+5K^2/64)^2+\frac{m^4}{12}
+\frac{5m^2K^2}{128}+\frac{K^4}{16}\right\}.
\]

For clarity, 2(1−cosω) differs from μ²+λ/6 by at most μ⁴/12+λμ²/12;
the difference between ω² and2(1−cosω) is at most ω⁴/12. The stated
window gives ω<π/2 and ω²≤3(μ²+λ/6), which yields the displayed bound.
Physical squared frequencies divide this expression by τ₀².

The continuum theorem concerns the declared ideal family. Uniform native
compilation and clock/error resources for every refined member remain an
implementation obligation. A sharp intermediate position measurement can
leave the bandlimited P domain; exact finite continuation is supplied
separately rather than silently reusing this theorem after that operation.

## Coupled scope

The concrete combined proposal is
V_m=U_E U_B Γ(S)_F Γ(R exp(−ihmR)). Its emitted native program, field
state, current, preparation and detector are on the same source. The
three-step development experiment uses h=1/8, m=0 or1, α=.17 and one face
with β=.41. It retains all288 electric registers and all position/event
branches. The source-derived mass response is measurable in position.

This finite result does not establish an interacting infinite-volume pole,
a ground-state vacuum, a homogeneous Maxwell limit, pair creation, a
two-charge force law, or observed species. Those claims cannot inherit the
free dispersion merely by naming the same mass parameter.
