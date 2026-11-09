# Momentum on a retained-origin observation domain

## 1. Source and comparison algebras

The source is the exact open Wilson/Gauss member of the
[common contract](../particle_common_contract/PROOF.md), with side
\(L_n=2^{2n+1}\), spacing \(a_n=2^{-n}\), physical density \(\rho_n\),
source Hamiltonian \(H_n\), and ground excess
\(\epsilon_n=2^{-n-4}\). Its full canonical conjugate reference has zero
Hamiltonian. Source stationarity does not mean that this purification is
stationary under source-only evolution.

The [spatial algebra construction](../particle_primary_spatial_algebra/PROOF.md)
supplies, at each regulator, an infinite physical comparison algebra
\(\mathcal G_n\), integer translations \(\beta_m\), a commuting infinite
time action \(\alpha_t^\infty\), and the extension \(\widetilde\omega_n\)
of the actual finite state by invariant exterior blanks. The physical
quotient imposes the local Haar constraints. No isometric embedding of the
finite physical matrix algebra in that quotient is assumed. The finite
source time is compared with the infinite action on interior observables
by the owning bulk bound.

The averaging argument for a fixed regulator appears in
[the spatial ground proof, section 6](../particle_spatial_ground_state/PROOF.md).
Here the regulator changes with \(n\), the averaging is an explicit retained
observer register, and physical spatial and time smears supply continuous
actions. The synthetic two-point translation instance in
`Lean/Variational/TranslationInvariantComposedInstance.lean` is not a
primary-source momentum theorem.

## 2. One origin throughout a complete history

Put
\[
 b_n=\min\{\lfloor(L_n-1)/2\rfloor,2^{\lfloor3n/2\rfloor}\},\qquad
 F_n=\{b_n,\ldots,L_n-b_n-1\}^3,\quad M_n=|F_n|.
 \tag{1}
\]
An independent classical register \(X\) is uniform on \(F_n\). Its
Hamiltonian is zero. The exact input is
\[
 \Omega_n={1\over M_n}\sum_{x\in F_n}|x\rangle\langle x|
                  \otimes|\Psi_n\rangle\langle\Psi_n|.
 \tag{2}
\]
The reduced source input is exactly \(\rho_n\). No source preparation,
reference dynamics, clock, or source coupling is replaced. The origin is
an additional declared observation service. Its finite label allocation
is \(M_n\); this is a specification, not execution of the large carrier.

On the infinite comparison algebra define
\[
 j_n(A)=\bigoplus_{x\in F_n}\beta_x(A).
 \tag{3}
\]
This is an isometric unital star homomorphism: every component is an
automorphism, and each product uses the same component. For a contained
local seed, the corresponding finite block is its actual translated
ambient operator restricted to the finite physical space. It is physical
because the seed is even and gauge invariant. Equation (3) is not an
assertion that every infinite physical operator has a finite norm
embedding. Time-evolved and quasi-local representatives use the comparison
in section 4.

For a complete local instrument \(\{K_r\}\), each translated family obeys
\(\sum_r K_{r,x}^\dagger K_{r,x}=I\). A history with fixed source time gaps
has finite branch operator
\[
 L_{\mathbf r,x}=K_{r_k,x}e^{-it_{k-1}H_n}\cdots
                       K_{r_2,x}e^{-it_1H_n}K_{r_1,x}.
 \tag{4}
\]
The same \(x\) controls every intervention. The complete output retains
\(X\), all outcomes including failures, the final source, and the unchanged
reference carrier. Its diagonal record block is
\(M_n^{-1}(L_{\mathbf r,x}\otimes I)|\Psi_n\rangle\langle\Psi_n|
(L_{\mathbf r,x}^\dagger\otimes I)\). Coherent fixed pointer entries use
\(\widetilde\omega_n(L_{\mathbf s,x}^\dagger L_{\mathbf r,x})/M_n\).
Summing over \(x\) is the pointer marginal of this complete output. No
success probability divides any weight.

Conditioning only on the input label \(X=x\) gives precisely the complete
fixed-origin translated experiment; its joint weight is \(1/M_n\) times
that law. This is the comparison map to fixed-origin records. It does not
identify their individual states with their average. A fresh origin at
each step is a different experiment. For the state \(|10\rangle\) and
one-site number readouts, a common uniform origin gives only records
\((1,1),(0,0)\), each with probability \(1/2\). Fresh origins give all four
pairs with probability \(1/4\). Also
\(\langle j(N^2)\rangle=1/2\), whereas the arithmetic operator average
\(\bar N=(N_0+N_1)/2\) has \(\langle\bar N^2\rangle=1/4\).

## 3. A quantitative translation-invariant observation state

The relative-observable functional is
\[
 \eta_n(A)={1\over M_n}\sum_{x\in F_n}\widetilde\omega_n(\beta_x A).
 \tag{5}
\]
For an integer displacement \(m\), write \(q_n=L_n-2b_n\). Reindexing the
sum, without an invariance assumption on \(\widetilde\omega_n\), gives
\[
 |\eta_n(\beta_m A)-\eta_n(A)|\leq\|A\|d_n(m),\quad
 d_n(m)=2\left[1-\prod_{i=1}^3
        \left(1-{|m_i|\over q_n}\right)_+\right]
 \leq 2\sum_i\min\{1,|m_i|/q_n\}.
 \tag{6}
\]
The final upper bound may be capped by 2. For fixed physical displacement
\(y\), take \(m_i=\lfloor y_i/a_n\rfloor\). Then \(d_n(m)\to0\).
Also \(b_n/L_n\to0\) and \(a_nb_n\to\infty\). Thus the retained origins
occupy an asymptotically full fraction while their physical collar diverges.

Use the regular algebra construction from the spatial proof: compact
spatial smears have cell-integral weights and the bound
\(\|A_n(f)-A_n(g)\|\leq\|f-g\|_1\); fixed triangular time smears supply
uniform time continuity. Include the physical one-cell collective seed
\[
 C_x={1\over\sqrt{42}}\sum_{\alpha=1}^2\sum_{i=1}^6
       c^\dagger_{+,\alpha i,x}c_{-,\alpha i,x}
 \tag{7}
\]
from [the refinement response](../particle_primary_refinement_response/PROOF.md),
its adjoint, and the charge and plaquette seeds. The dimension-six native
block has integer native charge 1; the bilinear has gauge charge zero.
Its complete Fock norm is at most 1. Adding this seed is an explicit
extension of the menu; membership in the unsmeared point algebra is not
assumed. Norm-continuity follows by the same bounds with one owner cell.

Take a countable dense menu of rational compact spatial smears, positive
rational time widths, rational translates and times, and all star words.
In the bounded product modulo the norm-null ideal their closure is a
separable regular algebra \(\mathcal A_{\rm rel}\). Rounded translations
form an exact point-norm continuous \(\mathbb R^3\) action \(\beta_y\)
on this algebra; they do not form an exact group on the entire product.
The source time action \(\alpha_t\) commutes with it. Boundedness and a
common diagonal subsequence give a state \(\eta\) on the algebra. Equation
(6) and spatial smear continuity imply \(\eta\circ\beta_y=\eta\).
This constructs a state on the declared relative-observation domain. It
does not prove invariance of the fixed-origin limiting state.

## 4. Same clock, boundary error, and positive time energy

For support radius \(R\), finite source-time extent \(T\), and an elementary
seed occupying \(s\) owners, the least origin collar has the conservative
graph-distance bound
\[
 r_n=b_n-\lceil R/a_n\rceil-2.
\]
The actual crossing interaction norm is bounded by
\(C_n=(1272/7)L_n^2/a_n\). With \(W_*=5232/7\), the inherited Duhamel/Lieb-Robinson
comparison gives, uniformly over the same retained origins,
\[
 e_n(R,T,s)=\min\{2,2sT C_n
             \exp[-r_n\log2+2W_*T/a_n]\}.
 \tag{8}
\]
The formula is applied once \(r_n\geq1\); outside that domain the norm cap
2 applies. At zero time the contained local representatives coincide.
For fixed \(R,T,s\), equation (8) tends to zero since
\(b_n\asymp2^{3n/2}\), \(T/a_n\asymp2^n\), and
\(\log C_n=O(n)\). Normalized spatial and time weights add no count of
origins or cells. For a fixed word, telescoping bounds its difference by
the sum of its factor errors times the other factor norm caps.

For positive contraction effects with norm difference at most \(e_n\),
the two square-root columns differ by at most \(\sqrt{2e_n}\). Their
complete instrument outputs, including any quantum reference, differ in
half trace norm by at most \(\min\{1,\sqrt{2e_n}\}\). A fixed history
telescopes these bounds. The output is block diagonal in \(X\), so its
half trace distance is the average of the origin-conditioned distances,
bounded by their maximum. This statement retains the origin-conditioned
canonical reference; it is not a statement about source marginals alone.

The finite source marginal is stationary. Equation (8) passes that
stationarity to every limiting relative word and gives
\(\eta\circ\alpha_t=\eta\). For a uniformly bounded finite block observable
\(B_n=\bigoplus_xB_{n,x}\), the source-only energy is
\(I_X\otimes H_n\), its ground floor is the same \(E_{0,n}\), and its
ground excess in (2) is \(\epsilon_n\). The original downward-transfer
mass below \(-\delta\) is therefore at most
\[
 \|B_n\|^2\epsilon_n/\delta.                 \tag{9}
\]
There is no \(M_n\) factor. Applying the
[positive-energy proof](../particle_joint_reconstruction/POSITIVE_ENERGY.md)
to each regular word, using the uniform comparison (8), gives nonnegative
time spectrum in the canonical Gelfand-Naimark-Segal representation of
\((\mathcal A_{\rm rel},\eta)\). The quantum reference is not part of this
energy generator. Exact finite observer outputs retain it, but no limit
of the growing full reference or origin-register density is asserted.

On the dense cyclic domain the canonical actions are
\[
 U(t)\pi(A)\Omega=\pi(\alpha_t A)\Omega,\qquad
 V(y)\pi(A)\Omega=\pi(\beta_y A)\Omega.
 \tag{10}
\]
Invariance makes both well-defined unitaries. Point-norm continuity makes
them strongly continuous. The algebra actions commute, hence so do the
unitaries. Their joint spectral theorem gives
\(U(t)=e^{itH_{\rm rel}}\), \(V(y)=e^{iy\cdot P}\), with
\(H_{\rm rel}\geq0\), \(H_{\rm rel}\Omega=0\), \(P\Omega=0\), and strong
commutation of all four generators. Fixing \(\Omega\) fixes their additive
origins. The spatial representation can be trivial, and the entire probe
response can collapse. These statements supply no nonzero momentum
weight, dispersion relation, Lorentz covariance, particle pole, or mass.

## 5. Original response weights and the collar

For any original positive success spectral measure \(\mu_x\) with total
mass at most 1, put \(\theta_n=1-M_n/L_n^3\). The actual full-box and
retained-origin measures obey exactly
\[
 \mu_{\rm box}=(1-\theta_n)\mu_F+\theta_n\mu_{\rm outside},\qquad
 \mu_F(S)\geq\max\left\{0,{w-\theta_n\over1-\theta_n}\right\}
 \quad\text{if }\mu_{\rm box}(S)\geq w.       \tag{11}
\]
This comparison uses the actual finite \(H_n,\rho_n\) on both sides.
It does not use an infinite comparison spectrum. The success marginal of
uniform origin with the two-outcome local instrument
\(C_x,\sqrt{I-C_x^\dagger C_x}\) equals the averaged success measure in
the refinement-response theorem. Its complete failure instrument is
specified separately and need not equal the single aggregated failure
map of that theorem.

Consequently its original moving band
\([1/(100\beta_n),1446/a_n]\), with full-box weight \(3/100\), has the
lower bound (11). At level 13, \(q_n/L_n=127/128\), and the bound is
\(353639/51209575>0.0069\). At level 16 it exceeds 0.0185; its limit is
0.03. At level 0 all eight origins are retained. The finite registered
event floor \(1/50\) transfers by the same formula for its stated clock
resolution. Each bound keeps the original probability normalization.
No fixed spatial or time smear survival, or fixed finite-energy limit,
follows from these unsmeared moving bands. Fixed pointer-history limits follow from the regular
algebra and the common subsequence, independently of that survival issue.

## 6. Executed controls and theorem boundary

The code regenerates the source member and the owning collective response,
computes exact cube overlap, collar, original-weight loss and a dyadic
bulk-error bound without allocating large boxes. Small quantum controls
retain an entangled reference and repeated origin-conditioned histories;
they verify completeness, reference preservation under the unread complete
channel, and the difference from fresh origins. A constant-probe control
has the same invariant cyclic construction and zero momentum, establishing
that momentum existence does not imply momentum response.

The limiting state is selected by a common subsequence, with no effective
subsequence, uniqueness or nonzero field asserted. The source Hamiltonian
has no finite periodic translation symmetry. The origin service supplies
an explicit relative-observation domain whose continuum translation group
has canonical commuting generators; it supplies no isomorphism to the
fixed-origin state representation. Finite Fourier apparatus labels and
these limiting generators require a separate convergence/response map
before a mass-shell interpretation.
