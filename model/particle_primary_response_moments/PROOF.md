# A populated finite excitation band

The selected primary state at member zero has original neutral-pair
transfer weight at least \(1/2\) in \([2,10]\), and original final
excitation weight at least \(3/8\) in \([2,21/2]\), in its declared source
clock. A complete finite two-energy apparatus reports the enlarged transfer
interval \([15/8,81/8]\) with probability at least \(15/32\). These are
finite interacting-source bounds, with retained failures and no success
normalization. They supply neither a particle pole nor a rest mass.

## 1. The actual source and neutral pair

Use [the primary contract](../particle_common_contract/PROOF.md), including
the full 60-mode matter carrier, product link cutoff, Gauss restriction,
selected stationary \(\rho=Q^2/\operatorname{Tr}Q^2\), idle canonical
reference and source-only clock \(e^{-itH}\). The source has spacing
\(a\), side \(L\), volume \(V=L^3\), onsite term
\(3d\Gamma(\beta_D)/a\), and the specified Wilson hopping, electric and
magnetic terms at coupling one. No source coefficient is varied.

The native \(Q_6=6\) block has dimension one. In the actual Dirac chart,
choose the orthonormal positive-beta vectors
\(p_1=(1,0,1,0)/\sqrt2\), \(p_2=(0,1,0,1)/\sqrt2\), and
negative-beta vectors \(m_1=(-1,0,1,0)/\sqrt2\),
\(m_2=(0,-1,0,1)/\sqrt2\). At each cell define

\[
 B_{x\alpha}=c^\dagger(p_{x\alpha})c(m_{x\alpha}),\qquad
 P_{x\alpha}=B_{x\alpha}^\dagger B_{x\alpha}
             =n_{m_{x\alpha}}(1-n_{p_{x\alpha}}).
 \tag{1}
\]

Each probe is a gauge-invariant even partial isometry of norm one. It
preserves every native-block total occupation and acts within the actual
physical carrier. These selection rules do not assign a continuum spin or
species. The matching beta basis is a declared detector chart; an individual
probe is not asserted to be a cubic scalar.

The electric and magnetic terms commute with (1). The onsite commutator is
\(6B/a\). Let \(h\) denote the dimensionless parent one-body hopping
operator in this one-dimensional native block. Multiplication by its link
characters is unitary pointwise before product compression. For fixed
\(x,\alpha\), put \(u=hp\), \(v=hm\). Both vectors are supported on
neighboring cells and are orthogonal to \(p,m\). With internal degree
\(d_x\), the actual Wilson matrices give

\[
 \|u\|^2=\|v\|^2=h_x=d_x/2,\qquad c_x=\langle u,v\rangle.
 \tag{2}
\]

The equality uses \(S_j^\dagger S_j=(I-\gamma_j)/2\) and
\(S_j S_j^\dagger=(I+\gamma_j)/2\), with their assignment determined by
the oriented bond. Every beta eigenvector has zero gamma expectation.
Open boundaries can give \(c_x\ne0\); orthogonality of \(u,v\) is not
assumed. The parent hopping commutator is the second quantization of
\(M=|u\rangle\langle m|-|p\rangle\langle v|\).

## 2. An all-occupation commutator norm

The full Fock norm, rather than its one-particle restriction, satisfies

\[
 \|d\Gamma(M)\|^2
   =h_x+\sqrt{h_x^2-|c_x|^2}\le d_x.             \tag{3}
\]

For completeness, rescale to \(h_x=1\), change phases so
\(c=|c_x|/h_x\in[0,1]\), and choose four orthogonal modes
\(p,m,e_2,e_3\) with \(u=e_2\),
\(v=ce_2+\sqrt{1-c^2}e_3\). Then
\(d\Gamma(M)=c_2^\dagger c_m-c_p^\dagger
(c c_2+\sqrt{1-c^2}c_3)\). On occupation sectors zero and four its squared
singular values vanish. Sectors one and three each have eigenvalues
\(0,0,1,1\); sector two has
\(0,0,1+s,1+s,1-s,1-s\), \(s=\sqrt{1-c^2}\).
Thus its characteristic polynomial on the 16-dimensional Fock space is
\(\lambda^8(\lambda-1)^4(\lambda^2-2\lambda+c^2)^2\).
The degenerate endpoints follow by continuity; extra spectator modes do
not change the norm. Restoring the common squared norm proves (3).

The matter probe commutes with the tensor link projection. Therefore
compression of the parent commutator equals the commutator of the
compressed source, without a discarded cross term. Compression and Gauss
restriction contract the norm. Consequently the exact physical residual

\[
 R_{x\alpha}=[H,B_{x\alpha}]-\frac6a B_{x\alpha}
 \quad\hbox{satisfies}\quad
 \|R_{x\alpha}\|^2\le\frac{d_x}{a^2}.          \tag{4}
\]

At member zero every cell has degree three. The explicit chart gives
\(|c_x|^2=1/4\) for both pairs at all eight corners, so (3) gives the
sharper cap \(3/2+\sqrt2<3\). The public bounds below use the rational
cap three. For arbitrary side \(L\), the averaged cap is
\(6(1-1/L)/a^2\). No selected-state weight bound at other members is
inferred from this norm estimate.

## 3. Original selected-state spectral weight

Diagonalize the actual stationary state and Hamiltonian together,
\(\rho=\sum_i p_i|i\rangle\langle i|\). Define the finite positive
measure without dividing by its total mass:

\[
 \mu=\frac1{2V}\sum_{x,\alpha,i,j}
       p_i|\langle j|B_{x\alpha}|i\rangle|^2
       \delta_{E_j-E_i}.                       \tag{5}
\]

Its total mass is the original success probability for a uniformly chosen
probe address. [The occupation bound](../particle_primary_occupation_bound/PROOF.md)
authenticates the same member and selected preparation, and proves
\(\operatorname{Tr}\rho D<5\) at member zero. The pair inequality
\(\sum_{x,\alpha}P_{x\alpha}\ge2V-D\), on the full Fock carrier,
therefore gives

\[
 \mu(\mathbb R)\ge\frac{11}{16},\qquad
 \int(\omega-6)^2\,d\mu(\omega)
   =\frac1{16}\sum_{x,\alpha}\operatorname{Tr}\rho R_{x\alpha}^\dagger
      R_{x\alpha}\le3.                        \tag{6}
\]

Markov's inequality outside radius four gives the concrete source result

\[
                    \mu([2,10])\ge\frac12.     \tag{7}
\]

The receipt also retains the stronger exact rational bound obtained from
the authenticated defect expectation, without decimal rounding. The
product vector used inside the occupation operator proof is not substituted
for \(\rho\) in (5) or (6).

The [original-weight excitation transfer](../particle_primary_mass_readout/PROOF.md)
couples each transfer \(\omega=E_j-E_i\) to final excitation
\(z=E_j-E_0\). The average initial excitation carried by (5) is at most
\(\epsilon=1/16\), since the aggregate success effect is at most the
identity. Discarding initial excitation above \(r=1/2\) costs at most
\(\epsilon/r=1/8\). For the surviving events \(z\in[2,21/2]\), so

\[
           \nu([2,21/2])\ge\frac38.             \tag{8}
\]

This inference requires no numerical value of the actual ground energy.
It uses its defining lower floor and the proved cooling inequality.

## 4. Complete finite records

One exact declared probe service has seventeen Kraus outcomes at member
zero:

\[
 K_{x\alpha}=B_{x\alpha}/4,\qquad
 K_f=\sqrt{I-\tfrac1{16}\sum_{x,\alpha}P_{x\alpha}}.
 \tag{9}
\]

The pair effects commute, so the complementary root is a finite diagonal
functional calculus of occupations. Equation (9) is complete on every
source input and arbitrary reference. Its success spectral measure is
exactly (5). The failure instrument is specified, retained and not removed
by normalization. Its success measure equals a uniform choice of probe;
its single failure map differs from retaining sixteen separate failure
maps and no equality of those complete instruments is claimed.

Apply the complete energy instruments from
[joint spectroscopy](../particle_primary_joint_spectroscopy/PROOF.md)
before and after (9). The fine-history Kraus operators are
\(F_c(H)K_bF_a(H)\), with all \(17J^2\) labels retained and the canonical
reference untouched. Their source probability law is

\[
 \sum_{i,j}p_i|(K_b)_{ji}|^2|F_a(E_i)|^2|F_c(E_j)|^2. \tag{10}
\]

The initial energy readout preserves the stationary source marginal. It
changes the source-reference joint state, and its record is retained.
No joint-purification stationarity is invoked.

Reuse that package's source-only `clock_plan` at member zero with energy
resolution \(\eta=1/16\), each complete energy-tail allowance at most
\(1/128\), and each complete clock implementation allowance at most
\(1/128\). Its alias-free signed chart and Cayley bounds act on the actual
source Hamiltonian with no energy-origin shift. Two tails and two clock
errors cost at most \(1/32\). On successful energy-resolution events the
reported difference differs from the true transfer by at most \(2\eta\).
Consequently the complete finite apparatus has

\[
 \Pr\{b\ne f,\ \widehat E_c-\widehat E_a\in[15/8,81/8]\}
       \ge\frac{15}{32}.                       \tag{11}
\]

An additional certified complete-instrument error or statistical event
error subtracts from (11) separately. The exact probe root service is an
architectural declaration; no efficient compilation or execution of the
large primary histogram is asserted. The certificate gives a sufficient
finite clock plan and an analytic probability bound. It is not a list of
sampled outcomes. All-reference errors refer to the complete instrument
with the same finite reference; no limiting reference is assumed.

## 5. What the band establishes

Equations (7), (8) and (11) concern a neutral, number-preserving excitation
channel on the actual finite interacting member zero. They do not isolate
an eigenvalue or exclude several transitions and continua within the band.
The bounds do not assign a spatial momentum, rest frame, Lorentz spin or
particle species. A declared absolute unit conversion can express these
source energies in comparison units, without converting an unresolved
transition band into a measured mass. A finite particle comparison requires
a resolved branch with controlled finite-domain dispersion and errors; a
continuum claim additionally requires the appropriate uniform survival and
field reconstruction. The present member-zero estimate supplies neither.
