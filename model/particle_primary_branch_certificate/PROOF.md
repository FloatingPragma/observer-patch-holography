# An admitted virtual-pair bound for the interacting odd sector

The full weak-pair/singlet occupation sector of the primary eight-cell
source has a bottom energy bounded above by
\(-715.100189896245\ldots\). The certificate uses 43 real trial
coordinates, including two successive hops on distinct links. Every
vector is inside the actual link cutoff and physical Gauss space. The
selected polynomial density gives a strictly positive original
sector-and-final-energy event in the certified interval. No isolated
particle pole or physical mass is inferred from this variational bound.

## 1. The source, sector and analysis vector

Use member zero of the [common contract](../particle_common_contract/PROOF.md):
all 480 matter modes, all quantum links, the degree-one product word
cutoff, physical Gauss projection, selected polynomial density, idle
reference, energy origin and source clock. The interaction coefficient
is one. The reducing native occupation increment is

\[
 \Delta=(0,2,0,0,1),\qquad (Q_6)_r=(-4,-3,1,2,6).
\]

Let \(\phi\) be the normalized physical vector with the two weak
positive-beta modes and the charge-six positive-beta mode occupied in
one spin row at a corner. Every other matter mode is the filled
negative-beta anchor; all links are constant. It is a gauge singlet.
The [sector calculation](../particle_primary_sector_threshold/PROOF.md)
gives \(H_0\phi=b\phi\), \(b=-711\), and a lower bound
\(L_\Delta=-722.2210758583\ldots\) on the **entire** sector.
Write \(H=H_0+T+B\), where T is all native hopping and B is signed
magnetic multiplication compressed to the actual cutoff. The source
bound is \(\|B\|\le6/7\).

The analysis vector is not the original input density. Variational
bounds concern the same source operator. Section 4 derives an original
probability directly from the selected density and its spectral weight.

## 2. Distinct-link second actions inside the source

For each native representation r and extra beta-defect count
\(k\in\{0,2\}\), let \(u_{r,k}\) be the corresponding projection of
\(T\phi\). Seven are nonzero. Their norm squares and relative
onsite/electric energies are

\[
 n_{r,0}=3|\Delta_r|/4,\quad n_{r,2}=12d_r-n_{r,0},
 \qquad d_{r,k}=105C_r/2+3k.                         \tag{1}
\]

The code separately enumerates the actual 480-bit CAR basis and all
three doubled Wilson matrices, whose nonzero entries are Gaussian units.
It recovers every norm in (1). Native link matrix coefficients have
Haar inner product \(1/d_r\); all ten native or dual representations
have distinct signed central charges.

For each of the 66 unordered pairs of distinct cube edges define

\[
 v_{bc}=(T_bT_c+T_cT_b)\phi.                         \tag{2}
\]

Each link in (2) carries at most one native word. The source cutoff
therefore admits both intermediate and final actions exactly. Operators
T preserve Gauss, as do the projections to native electric eigenspaces
and onsite defect count. No compression invariance is assumed.

Different edge pairs are Haar orthogonal: a link present in only one
pair has a nontrivial one-word Haar integral. Project (2) to an unordered
native pair \(r\le s\) and extra defect count \(k=0,2,4\), and sum
over edge pairs. This gives 30 nonzero mutually orthogonal packets
\(w_g\), with norms \(N_g\) and relative energies

\[
 D_g=105(C_r+C_s)/2+3k.                              \tag{3}
\]

The exact CAR/Haar enumeration gives

\[
 \|v_{bc}\|^2=\begin{cases}1785/2&b,c\text{ share a vertex},\\
 900&b,c\text{ are disjoint},\end{cases}\quad
 \sum_gN_g=59220,\quad \sum_gN_gD_g=11350341/2.      \tag{4}
\]

There are 24 adjacent and 42 disjoint pairs. A distinct pair of hops
has a nonzero net site-number divergence: the graph has no two-edge
cycle. Consequently no resulting matter profile equals that of phi.
The code checks this on every retained coordinate.

The hopping matrix element between \(u_i\) and \(w_g\) is
\(c_{ig}=\langle u_i,Tw_g\rangle\). Its exact computation takes the
inner product of (2) with the two-hop image of the specified first
packet. Terms with any other second edge are Haar orthogonal. Every
aggregate c is real, and \(\sum_{i,g}c_{ig}=59220\), an independent
adjoint/norm check. Bounded Gaussian-integer arithmetic is used before
the rational Haar factors; the certificate rejects any coefficient
outside that exact domain.

Include the five constant-matter magnetic-loop packets from the
[sector trial](../particle_primary_sector_threshold/PROOF.md). Their
common norm is \(g_0=6/22050\), and their relative matrix is
\(A_{rs}=210C_r\delta_{rs}-M_{rs}/420\), where M is the authenticated
native character-product matrix. Bipartite matter parity makes T
vanish within the one-hop block and within the two-hop block. It also
eliminates the one-hop/two-hop matrix elements of H0 and B. A loop
packet has anchor matter profile and hence no H0 or B overlap with a
two-hop packet; T between these even matter parities vanishes.
These statements retain the complete operator outside the trial space.

Only the two-hop B block is bounded rather than evaluated: its quadratic
form is at most \((6/7)\sum_gN_g z_g^2\). Thus the following 43-coordinate
quadratic form is an **upper bound** on the expectation of the actual H
for each specified real coefficient vector. In the unnormalized
orthogonal basis its metric is \(\operatorname{diag}(1,n_i,N_g,g_0I_5)\).
The relative form has phi/one-hop entries \(n_i\), phi/loop entries
\(g_0\), one-hop diagonal \(n_id_i\), one-hop/two-hop entries c,
two-hop diagonal \(N_g(D_g+6/7)\), and loop block \(g_0A\).

For \(z\ge0\) eliminate the two-hop block and define

\[
 S_{ij}(z)=n_i(d_i+z)\delta_{ij}
       -\sum_g\frac{c_{ig}c_{jg}}{N_g(D_g+6/7+z)},\quad
 h=S(z)^{-1}n,\quad \ell=(A+zI)^{-1}{\bf1},
\]
\[
 a_g=-\frac{\sum_i c_{ig}h_i}{N_g(D_g+6/7+z)},\qquad
 q(z)=\sum_i n_ih_i+g_0\sum_r\ell_r.                 \tag{5}
\]

Exact Schur pivots show the excited form is positive at zero. Hence
q is positive and decreasing. Rational bisection gives a width
\(2^{-24}\) bracket for \(q(z)=z\). At its lower endpoint choose
\(\psi=\phi-\sum_i h_i u_i-\sum_g a_gw_g-\sum_r\ell_r b_r\).
For \(s=\|\psi\|^2-1\), the actual Rayleigh quotient obeys

\[
 \frac{\langle\psi,H\psi\rangle}{\|\psi\|^2}
 \le b+\frac{-q(z)-zs}{1+s}
 =:U_\Delta=-715.100189896245\ldots.                 \tag{6}
\]

All coefficients and the upper bound in the receipt are rational.
The norm is \(1.1107699157673947\ldots\). The 13-coordinate source
trial has upper \(-714.7851\ldots\), so (6) is a strict improvement
on that complete-sector upper bound. It gives neither an eigenvector
nor an isolated eigenvalue. The known scalar retained-penalty family
has comparison floors strictly below \(-715.4069\ldots\); (6) does
not bridge that particular isolation test. Those comparison ceilings
are not upper bounds on the actual complementary spectrum.

## 3. Compulsory degeneracy and its exact symmetry reduction

The filled-beta anchor has 240 matter particles. Every vector in this
native sector has odd total CAR parity. The central element of the
actual binary cube acts as that parity. The possible proper irreps
are therefore G1, G2 and H, with dimensions 2, 2 and 4. Since the full
source H commutes with these rotations, every energy eigenspace in
this sector has dimension at least two.

Let E_delta be the sector ground energy. For **any** rank-one P in
the full native sector, its orthogonal complement Q contains a
normalized sector ground vector. Therefore

\[
 \inf\sigma(QH_\Delta Q)=E_\Delta.                  \tag{7}
\]

A rank-one complementary floor above the sector ground is impossible
on this full sector. Equation (7) does not exclude an isolated multiplet
or a pole in a scalar cyclic spectral measure.

The compulsory partners can be removed exactly. The actual six-state
local kernel supplies G1 and H matrices by restricting its 48 literal
exterior rotation matrices to the spin-coefficient Casimir subspaces.
G2 is the A2 character times G1. The code checks all 576 matrix-coefficient
inner products and the corresponding group sums on the actual local
kernel. On the full source define

\[
 E^{\lambda}_{ab}=\frac{d_\lambda}{48}
     \sum_{g}\overline{D^{\lambda}_{ab}(g)}U_g.       \tag{8}
\]

Finite-group orthogonality gives
\((E_{ab}^\lambda)^*=E_{ba}^\lambda\) and
\(E_{ab}^\lambda E_{cd}^\mu
 =\delta_{\lambda\mu}\delta_{bc}E_{ad}^\lambda\).
Every matrix unit commutes with H. The eight diagonal units form a
complete orthogonal projector menu on the odd native sector, with all
other native sectors retained as complementary outcomes. The same
energy clock commutes with this menu. Thus
\(H|_\lambda=I_{d_\lambda}\otimes h_\lambda\).

For a proposed normalized vector in \(\operatorname{ran}E_{aa}^\lambda\),
take its rank-one P and \(Q=E_{aa}^\lambda-P\). The exact Feshbach
operator is

\[
 PHP-zP-PHQ(QHQ-zQ)^{-1}QHP,                         \tag{9}
\]

where the inverse exists. The forced irrep partners do not occur in
this Q. Equation (9) is an instantiated symmetry reduction of the
actual source. A lower bound on this reduced QHQ, or a controlled
resolvent through the trial window, is not supplied by (6). The trial
at one corner is not assigned to a particular total cube irrep.

## 4. Original density weight in the lower-energy sector interval

Let \(L_0\) be the authenticated lower bound on the **global** physical
ground and let \(D\le2^{b_D}\) bound physical dimension. The actual
selected density is \(\rho=p(H)\). Its
[polynomial modular derivative](../particle_primary_gravity_premises/PROOF.md)
satisfies \(0<-\partial_E\log p(E)\le\beta\). For every source
energy \(E\le U_\Delta\), comparison to a global ground eigenvalue
E0 gives

\[
 p(E)\ge\frac{e^{-\beta(E-E_0)}}D
 \ge2^{-b_D-\lceil2\beta(U_\Delta-L_0)\rceil}
 =:c_{\rm low}.                                    \tag{10}
\]

The first inequality uses \(p(E_0)\ge1/D\); the second uses
\(E_0\ge L_0\) and \(e<4\). This is a source-marginal operator bound
on \(1_{H\le U_\Delta}\), not on the retained pure source-reference
density. The certificate regenerates every term from the same member.
The energy-local dyadic exponent is substantially smaller than the
full-spectrum exponent 582407133; no floating exponential is evaluated.

There is a sector ground eigenspace of multiplicity at least two in
\([L_\Delta,U_\Delta]\). The exact source charge-conjugation unitary
preserves H, rho and the clock, and exchanges the two disjoint signed
native sectors. Consequently the ideal original probabilities are at
least \(2c_{\rm low}\) per sign and \(4c_{\rm low}\) for their union.

Apply the complete native occupation measurement and one complete
source-energy clock, with resolution eta=1/64 and the authenticated
tail plus implementation loss \(\delta<1/512\). Every energy
eigenvector in this interval is reported in
\([L_\Delta-\eta,U_\Delta+\eta]\) with probability at least
\(1-\delta\). Positivity of that actual finite-clock effect therefore
gives original probability

\[
 p_{\rm reported,sign}\ge2c_{\rm low}(1-\delta),
 \qquad p_{\rm reported,either}\ge4c_{\rm low}(1-\delta). \tag{11}
\]

The complete-clock loss is inside the eigenvector effect before the
tiny density factor is multiplied. No fixed error is subtracted from
the tiny final probability. All native sectors, clock outcomes and
reference outputs are retained. There is no initial energy measurement,
creation probe, conditional success normalization or transfer-law
interpretation in (11). These finite weights have no claimed uniform
physical refinement limit. The source-energy interval is separated
from the global ground but is not an identified GeV rest mass.
