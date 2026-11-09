# Selected phase work and its three-dimensional gradient limit

The unchanged primary Wilson source gives an exact spatial comparison of
complete local phase experiments across capacity. Its unknown selected-state
hopping coefficient cancels from their work ratios. The limit is a positive
three-dimensional Dirichlet work form, with a controlled lattice error and
a smaller phase-nonlinearity error. This is a limit of energy moments. It
does not establish nonzero probability in a fixed positive-energy band,
a relativistic dispersion relation or a particle mass.

## 1. Actual state, clock and complete profile apparatus

Use the [primary contract](../particle_common_contract/PROOF.md):

\[
 a=2^{-n},\quad R=1/a,\quad L=2R^2,\quad V=L^3,\qquad
 H_n=\widehat H_{L,n+1}/a.
 \tag{1}
\]

The full sixty-mode matter carrier, quantum links, open Hamiltonian bonds,
product word cutoff, boundary Gauss constraints, selected polynomial density
\(\rho_n\), idle conjugate reference and source clock are retained. Put
\(n_y=N_y-30I\), where \(N_y\) counts all sixty native modes. The
[selected energy inequality](../particle_primary_energy_tightness/PROOF.md#6-selected-hopping-energy-and-a-complete-density-sum-rule)
supplies a coefficient floor with a sharper uniform form on the actual
capacity sequence:

\[
 h_n=-\frac aV\operatorname{Tr}(\rho_n H_{\rm hop}),\qquad
 h_n\ge h_* >0.0424,\qquad h_n\le30B/V<90,
 \tag{2}
\]

with \(B=3L^2(L-1)\). No exact value of \(h_n\) is asserted. The
variational vector enters the energy inequality only; the expectation in
(2) is the actual preparation. Proper cube symmetry gives equal expected
total hopping energies in the three directions on this open source.

For the coefficient improvement, let \(g=(u_0+720)/8\) be the
authenticated tiled-trial gain per cell. The exact lower floor is
\(-g-P/(7V)-a\epsilon/V\), with \(\epsilon=a/16\). The selected
sequence has \(P/V=3(1-a^2/2)^2\) and \(V=8/a^6\). Thus

\[
 h_n\ge h_*+\frac{3a^2}{7}-\frac{3a^4}{28}-\frac{a^8}{128},
 \quad h_*=-g-\frac37,
 \quad \frac{3a^2}{7}-\frac{3a^4}{28}-\frac{a^8}{128}
             \ge\frac{281a^2}{896}>0.
 \tag{2a}
\]

This retains the relation between finite boundary size and cooling rather
than subtracting their separate uniform worst cases.

Choose fixed positive integer physical half-widths
\(b=(b_1,b_2,b_3)\), a fixed nonzero rational amplitude \(s\), and a
member with \(R\ge\max b_i\). Write \(m_i=b_iR\). On the finite cyclic
address list define

\[
 t_m(k)=\max(1-d_L(k,0)/m,0),\quad
 f_b(k)=\prod_{i=1}^3t_{m_i}(k_i),\qquad
 U_{X,\sigma}^{b,s}
 =\exp\left(i\sigma\pi s a\sum_y f_b(y-X)n_y\right).
 \tag{3}
\]

The source Hamiltonian has no wrapping bonds. Only the declared apparatus
coefficient list wraps. The complete instrument has Kraus operators
\(U_{X,\sigma}^{b,s}/\sqrt{2V}\), for every address and both signs.
Its original spectral mass is one. The address, sign, quantum reference
and both complete energy records are retained. Initial energy recording
preserves the stationary source marginal and generally changes its
correlations with the idle reference. No successful probability is divided
out. For an origin separated from the boundary by the profile widths,
(3) is supported in an ordinary spatial box of those physical widths.

Every coefficient in the phase, divided by \(\pi\), is rational with
denominator dividing \(\operatorname{den}(s)b_1b_2b_3R^4\). Thus (3)
is a finite number-phase program, using roots of unity of order dividing
twice that denominator. It commutes with the link and Gauss projectors.
It is not assumed to preserve its input state.

## 2. Exact selected work and finite spatial errors

For a fixed direction \(i\), let

\[
 k_i(b,s)=\sum_{z\in\mathbb Z_L^3}
       [1-\cos(\pi s a(f_b(z+e_i)-f_b(z)))],\qquad
 \overline k=\tfrac13\sum_i k_i.
 \tag{4}
\]

Only hopping changes under (3). Averaging the retained signs removes its
sine term. Averaging the retained origins gives the same coefficient
\(k_i/V\) on every actual bond in direction \(i\), including bonds
near a face. This translation is of apparatus settings; the state need
not be translation invariant. Equal directional hopping expectations then
give the exact original signed energy-transfer moment

\[
 W_n(b,s)=-\frac{\overline k}{V}\langle H_{\rm hop}\rangle_{\rho_n}
          =h_n\frac{\overline k}{a}.
 \tag{5}
\]

The elementary finite tent sums are

\[
 S_2(m)=\sum_k t_m(k)^2=\frac{2m^2+1}{3m},\quad
 S_4(m)=\sum_k t_m(k)^4=\frac{6m^4+10m^2-1}{15m^3},
\]
\[
 \sum_k(\Delta t_m)^2=2/m,\qquad
 \sum_k(\Delta t_m)^4=2/m^3.
 \tag{6}
\]

These identities include the case \(L=2m\). For
\(D_p=\sum_{i,z}|\Delta_i f_b(z)|^p\), equation (6) gives exactly

\[
 D_2=\sum_i\frac2{m_i}S_2(m_j)S_2(m_k),\qquad
 D_4=\sum_i\frac2{m_i^3}S_4(m_j)S_4(m_k),
 \quad \{j,k\}=\{1,2,3\}\setminus\{i\}.
 \tag{7}
\]

The all-real Taylor inequality
\(0\le u^2/2-(1-\cos u)\le u^4/24\) proves

\[
 0\le \frac{\pi^2s^2aD_2}{6}
          -\frac{W_n(b,s)}{h_n}
 \le\frac{\pi^4s^4a^3D_4}{72}
 <\frac{25s^4a^3D_4}{18}.
 \tag{8}
\]

This error concerns an expectation in the original selected state. It
contains no trial substitution, spatial covariance approximation or
spectral-gap premise. For fixed \(b,s\), \(D_4=O(a)\), so the right
side is \(O(a^4)\).

Define the ordinary compactly supported physical tent
\(f_b^{\rm cont}(x)=\prod_i\max(1-|x_i|/b_i,0)\). Direct integration
and the exact lattice sums give

\[
 \int_{\mathbb R^3}|\nabla f_b^{\rm cont}|^2\,dx
   =\frac89\sum_i\frac{b_jb_k}{b_i}=:I_b,
\]
\[
 aD_2-I_b
 =\frac{4a^2}{9}\sum_i\frac1{b_i}
                  \left(\frac{b_j}{b_k}+\frac{b_k}{b_j}\right)
       +\frac{2a^4}{3b_1b_2b_3}.
 \tag{9}
\]

Consequently the exact output obeys the explicit estimate

\[
 \left|\frac{W_n(b,s)}{h_n}
       -\frac{\pi^2s^2}{6} I_b\right|
 \le \frac{5s^2}{3}(aD_2-I_b)
        +\frac{25s^4a^3D_4}{18}.
 \tag{10}
\]

Its lattice error is \(O(a^2)\), and its phase error is \(O(a^4)\),
with every coefficient in the right side rational and computed. For the
unit-width tent the limiting coefficient is \(4\pi^2s^2/9\).
Equation (2) supplies a positive lower coefficient and a finite upper
coefficient. For example, whenever \(|s|a^2/\min_i b_i\le1\), use
\(1-\cos(\pi u)\ge2u^2\), \(|u|\le1\), to obtain

\[
 W_n(b,s)\ge\frac{16s^2h_*}{27}\sum_i\frac{b_jb_k}{b_i}>0.
 \tag{11}
\]

The restriction is necessary for this finite-member lower bound: at
member zero and amplitude two the phases in (3) are exactly trivial.
It is satisfied eventually for every fixed admitted profile.

Every sequence of members has a subsequence along which \(h_n\) tends
to some \(h\in[h_*,90]\). On the same subsequence, the countable family
of fixed integer-width, rational-amplitude work measurements converges to

\[
 W(b,s)=\frac{\pi^2h}{6}
             \int|\nabla(s f_b^{\rm cont})|^2\,dx.
 \tag{12}
\]

This is a positive three-dimensional Dirichlet work form evaluated on
the stated profiles. Its coefficient need not have a unique subsequential
limit. The formula does not assert a continuum limit of all source
correlations or identify the work form with an effective field action.

The work limit can be recorded with finite source clocks. Let
\(\mathcal R_n\ge\|H_n\|\) be the authenticated source norm bound and
use the complete alias-free clock with tick \(1/(4\mathcal R_n)\)
from [joint spectroscopy](../particle_primary_joint_spectroscopy/PROOF.md).
Its signed energy labels have absolute value below \(14\mathcal R_n\).
If each clock has resolution \(\eta\), tail probability at most \(b\)
and the implemented two-clock law has half-trace distance at most \(e\)
from the ideal finite-clock law, the reported signed mean obeys

\[
 |W_n^{\rm reported}-W_n|
       \le2\eta+56\mathcal R_n(2b+e).
 \tag{12a}
\]

On the good event each transfer error is at most \(2\eta\). The bad
event has probability at most \(2b\), and its reported-versus-true
transfer discrepancy is below \(30\mathcal R_n\). The implemented
signed statistic has range at most \(56\mathcal R_n\), giving (12a).
These are original complete laws with the reference retained. Set
\(u_n=a^2(16h_*/9)/128\), \(\eta=u_n/16\), and give each clock
tail and implementation allowance \(u_n/(1024\mathcal R_n)\). The
actual finite Cayley/Fejer planner then gives error at most
\(11u_n/32=O(a^2)\). The receipt binds its clock sizes, substeps and
all \(2VJ^2\) fine records. Thus the finite reported means have the
same limit (12). No full histogram or finite-shot estimate is asserted.

## 3. Spatial comparisons that remove the unknown coefficient

For any two admitted profiles, their exact work ratio in (5) cancels
\(h_n\). Compare with the unit-width, amplitude-one profile. Let
\(Q=\pi^2s^2aD_2/6\) and
\(\eta=(5/6)s^2a^2D_4/D_2\); use subscripts zero for the reference.
Equation (8) gives, when \(\eta_0<1\),

\[
 \frac{Q}{Q_0}\max(0,1-\eta)
 \le\frac{W_n(b,s)}{W_n(1,1)}
 \le\frac{Q/Q_0}{1-\eta_0}.
 \tag{13}
\]

All endpoints are rational. Their common limit is
\(s^2\sum_i b_jb_k/(3b_i)\), independently of the value or convergence
of \(h_n\). Doubling the amplitude gives ratio four. Doubling all three
widths gives ratio two. Doubling one width gives ratio three halves.
Permuting the widths leaves the exact finite work unchanged. These are
comparisons of spatial apparatus work, not particle energy-momentum laws.

A second exact comparison concerns the complete volume-normalized real
density apparatus. Use cosine and sine coefficients at spatial phase step
\(aq\), common normalization \(30V\), and their complete failure
effect. For positive fixed \(q=\pi r\), rational \(r\), the source
Ward identity gives

\[
 M_{1,n}(q)=-\frac{1-\cos(aq)}{3(30V)^2}
                   \langle H_{\rm hop}\rangle_{\rho_n},\qquad
 \frac{M_{1,n}(\pi r)}{M_{1,n}(\pi)}
       =\frac{1-\cos(\pi ar)}{1-\cos(\pi a)}\longrightarrow r^2.
 \tag{14}
\]

Put \(\zeta_0=5a^2/6\), \(\zeta_r=r^2\zeta_0\). Rational lower
and upper bounds for this ratio are
\(\max(0,r^2(1-\zeta_r))\) and \(r^2/(1-\zeta_0)\).
At coarse members a frequency may alias to zero. The original mean of
this volume-normalized apparatus is \(O(a/V)\); its fixed positive-energy
tail vanishes under the selected cooling estimate. Ratios of its vanishing
moments do not supply a particle residue. The frequency settings in
(14) are not asserted to be canonical momentum outcomes.

## 4. Energy probability and the actual criticality boundary

The source [energy-form estimate](../particle_primary_current_tail/PROOF.md)
applies to every profile in (3). Indeed its generator satisfies

\[
 \|[F,[H,F]]\|\le30\pi^2s^2aD_2<300s^2aD_2.
 \tag{15}
\]

For each fixed profile this is uniform in capacity. It gives a common
energy-form domain and bounded first absolute energy moments for finite
histories. In particular downward transfer costs at most the selected
ground excess \(\epsilon=a/16\), so the original absolute transfer
moment of the complete experiment is at most \(W_n+2\epsilon\).
None of (8)--(15) bounds the selected energy-weighted ultraviolet tail.
The probability measures may converge to a purely elastic measure while
their positive work is carried by a vanishing probability at diverging
energy. Thus (12) cannot be passed through the limiting probability
measure without an additional uniform-integrability estimate.

The fixed dimensionless electric coefficient105 and onsite coefficient3
give no theorem here that the interacting dimensionless gap is positive
uniformly. The exact physical gap6 of the decoupled product anchor is
proved in [the source spectral calculation](../particle_primary_spectral_scale/PROOF.md).
The actual endpoint contains hopping that creates orthogonal vacuum
packets. For one dimensionless bond \(T_b\),

\[
 \|T_b\Omega_0\|^2=15,
 \tag{16}
\]

at every admitted cutoff. A star of \(z\) bonds has corresponding
squared norm \(15z\), because different-link packets are orthogonal.
Suppose this fixed-anchor perturbation is decomposed as
\(T_b=R_b+B_b\), with \(R_b\) purely form bounded by a positive
decoupled local operator \(G\) annihilating \(\Omega_0\):
\(|\langle v,R_bv\rangle|\le\alpha\langle v,Gv\rangle\).
Applying the inequality to \(\Omega_0+t v\), with real and imaginary
\(t\), shows \(R_b\Omega_0=0\). Therefore

\[
                   \|B_b\|\ge\|B_b\Omega_0\|=\sqrt{15}.
 \tag{17}
\]

The star version has lower bound \(\sqrt{15z}\). This prevents that
particular decomposition from hiding the vacuum-creating term in a
purely relative perturbation with arbitrarily small bounded remainder.
It does not exclude a dressed anchor, a convergent connected expansion,
or every spectral-gap stability theorem.

Theorems1--2 of
[Yarotsky, *Ground states in relatively bounded quantum perturbations of classical lattice systems*](https://arxiv.org/abs/math-ph/0412040)
require small perturbation constants, or a purely relative part together
with a sufficiently small bounded part. The source bounds supply no
verified threshold for this endpoint. The quasiparticle theorem in
[Yarotsky, *Quasi-particles in weak perturbations of non-interacting quantum lattice systems*](https://arxiv.org/abs/math-ph/0411042)
also requires an isolated particle subspace at the reference point.
Neither theorem is imported as an endpoint or continuum conclusion.
A proof of a uniform dimensionless gap would exclude fixed physical
positive-energy survival on the unchanged \(H_n=\widehat H_n/a_n\)
trajectory. A viable route would then require an explicitly different
critical source trajectory. No such gap or trajectory is derived from
the work form (12).
