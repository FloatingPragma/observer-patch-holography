# Neutral stationary band and complete native records

## Declared source and exact energy data

The carrier is the complete periodic two-cell-per-axis lattice with thirty
fermion modes and five scalar states per cell and 119 states per oriented
link. Its dimension is
\(D=(5\,2^{30})^8 119^{24}\). The physical space is the joint kernel of the
actual twelve Gauss generators at each cell. Declare

\[
 H=H_{\rm loc}+2N_{\rm nonblank}+\lambda T_{\rm per},
 \qquad\lambda=2^{-34},\quad \|T_{\rm per}\|\leq360.
 \tag{1}
\]

The local quantum Yukawa, scalar and representation operators are the native
ones with coefficients one. `source.full_column` acts on full carrier keys.
Periodic fermion boundary conditions change twelve seam signs. The link
penalty is two; this choice is independent of a measured mass. The link
projector commutes with every left/right native current, translations and
the bipartite parity. Source dimension, charged intermediate states and
nonblank-link complements are retained. The coefficient-one source is a
distinct member; its ground state or preparation is not transferred.

Let \(u=x(x-1)\). The complete local physical ground energy \(e_0\) is the
smallest real root of \(u^3-313u^2+31708u-1041468\). Its eigenspace is a
raw-lepton-number-one spin doublet. The first local excited energy \(e_\nu\)
is the smallest root of \(u^3-313u^2+31738u-1041088\); it has one raw-number-zero
and one raw-number-two spin-zero eigenline. The complete local spectrum,
including the excluded collective modules, is consumed from
`local_reference`. Every other local physical energy exceeds \(-11.13\).
The unrestricted local lower bound is \(l_0=-11.71\). `isolation.py` isolates
the algebraic roots and rounds outwards to exact dyadic endpoints, so all
thresholds below are rational certified bounds.

Put \(\Delta=e_\nu-e_0\). A support with \(m\) nonblank links touches at most
\(2m\) cells. Untouched physical cells have energy at least \(e_0\), while
touched cells have energy at least \(l_0\). Thus its energy is bounded below by

\[
 8e_0+m[2-2(e_0-l_0)].                                    \tag{2}
\]

This exceeds \(8e_0+1.04m\). At link penalty one the same bound is only
\(8e_0+0.0411705884m\); its one- and two-link lower bounds overlap the
neutral excitation window. That loss of a sufficient isolation bound proves
no actual spectral collision in the coefficient-one member.

## Complete band isolation and the reference ground

Restrict to conserved raw lepton number seven and baryon number zero. At
\(\lambda=0\), the lowest band consists of one local raw-number-zero neutral
eigenline, seven true ground doublets and blank links. Its energy and rank are

\[
 a=7e_0+e_\nu,\qquad r=8\,2^7=1024.
\]

The other local level bound, two-neutral threshold and (2) give a complete
orthogonal-complement separation \(d>0.0110299637715905\). Write its projector
as \(P\), \(Q=I-P\) in this physical charge sector, and
\(V=\lambda T_{\rm per}\), \(v=360\lambda\). Every native hopping word
contains a blank-link matrix element zero; hence \(PVP=0\). For any vector
\(p+q\), completing a square gives

\[
 \langle H-a\rangle\geq(d-v)\|q\|^2-2v\|p\|\|q\|
 \geq-\epsilon\|p\|^2,\qquad
 \epsilon=\frac{v^2}{d-v}.
\]

The min-max principle on \(P\) gives the matching upper bound \(a\), while
the complementary min-max bound is \(a+d-v\). There are exactly 1024 full-H
eigenvalues in \([a-\epsilon,a]\), with
\(\epsilon<3.9810\times10^{-14}\). No enormous matrix diagonalization is
used in this conclusion.

For the spectral projector \(E\) of that known interval, the exact Sylvester
equation is

\[
 (QH_0Q)(QE)-(QE)(H|_E)=-QVE.                              \tag{3}
\]

The left operator has spectrum at least \(a+d\), while the right operator
has spectrum at most \(a\). The convergent semigroup integral solving (3)
therefore gives \(\|QE\|\leq v/d=\eta\), independent of the dimensions.
The projectors have equal finite rank, so
\(\|E-P\|=\|QE\|\leq\eta<1.900\times10^{-6}\). This argument uses the
actual separated spectra and supplies the stated denominator \(d\).

For the full physical ground, the atomic projector has rank 256 and gap
\(\Delta\); (2) excludes all nonblank sectors. Its hopping compression also
vanishes. The same square-completion argument gives

\[
 8e_0-\epsilon_g\leq E_g\leq8e_0,
 \qquad\epsilon_g=\frac{v^2}{\Delta-v}.
\]

Uniqueness and reference momentum follow from a separate finite theorem
shortcut. The actual complete charged semigroup and native current tensors
give the second-order spin operator

\[
 K_2=-24A_2I+2J_2H_{\rm cube},\qquad
 J_2=\frac14\sum_f d_f\int_0^\infty e^{-2t}b_f(t)^2\,dt>0. \tag{4}
\]

This retains the scalar shift and all charged intermediate states. The
`atomic_reference_control` response bound has
\(b_L(t)\geq(1-\beta_{ec})/2-2/15=:b_* >0\) for
\(0\leq t\leq w=1/1104\). Only the link weight changes in (4).
Consequently \(J_2\geq w(1-2w)b_*^2/2>8.2569\times10^{-5}\).
The complete 256-spin cube spectrum has a unique ground and gap
\(\delta_{\rm cube}>4/5\). Bipartite parity commutes with the doubled link
penalty, so the full-complement canonical graph theorem gives
\(R_4\leq25v^4/\Delta^3\). The exact rational check proves

\[
 2\lambda^2J_{2,\min}\delta_{\rm cube}-2R_4
 >4.4542\times10^{-25}>0.                                 \tag{5}
\]

This inequality holds throughout the positive interval up to the declared
\(\lambda\). The cube ground has translation characters \((1,1,1)\).
Eight odd local ground factors acquire an even fermionic permutation sign
under each axis translation; the four wrapped factors also have even parity.
Thus periodic boundary conditions preserve the native ground character.
The isolated continuation preserves it, proving that the actual physical
ground of (1) is unique and has momentum zero.

Combining the neutral-band and ground bounds gives

\[
 E_\nu-E_g\in
 [0.0895553304760346,\,0.0895553304760794].                  \tag{6}
\]

## Explicit local profile and all momentum sectors

`local_profile.py` builds the actual sixteen-dimensional raw-number-zero
physical block with its native Gram matrix and full fermion/scalar embedding.
Let \(q\) be the square-free characteristic polynomial. The selected local
projector is the explicit polynomial

\[
 P_{\nu,\mathrm{loc}}=
 \frac{\sum_{i=0}^{\deg q-1}H_{\rm loc}^i
                  \sum_{j=i+1}^{\deg q}q_j e_\nu^{j-i-1}}
 {q'(e_\nu)}.                                           \tag{7}
\]

Polynomial remainders prove its idempotence, eigenvalue equation and rank
one. The spin-one block has no common selected factor. The actual physical
embedding annihilates all twelve local Gauss currents. Global raw lepton
number seven differs from the ground's eight; the exact source conserves it.
No external charged reference, assigned sterile coordinate or repair-based
definition of electric neutrality is introduced.

The complete local occupation detector resolves native \(N_{ec}=0,1,2\)
and scalar blank/excited. Raw lepton number zero gives \(N_L=N_{ec}\).
Their six probabilities are approximately

\[
 (0.204735733,0.185796878,0.256485576,
  0.233033289,0.0602544543,0.0596940694).
\]

Their exact algebraic numerators sum to the denominator. This neutral
profile is not an empty Fock state. Occupation measurement is performed
after the profile projection and can change the energy; no stationary
postmeasurement state is inferred.

The complete atomic band is indexed by a neutral cell and seven ground-spin
bits. Translating the seven odd factors carries their actual fermionic
permutation sign. Exact enumeration of all 1024 labels proves commuting
order-two translations with character \((1024,0,0,0,0,0,0,0)\).
Thus every momentum has rank 128. Native full-occupation translation controls
check the grading. The spectral gap and covariance preserve all eight ranks
at nonzero \(\lambda\). The antiperiodic member, whose odd translation
squares to minus identity, is not identified with this periodic member.

## A finite polynomial preparation and its cost

Let \([a_-,b]\) enclose the interacting band, \(c>b\) bound its complete
complement from below, and \(U\) bound \(\|H\|\). These are explicit
rationals in the receipt. Inside the joint physical, raw-number-seven,
baryon-zero, momentum-zero projector \(\Pi\), define

\[
 x(E)=\frac{2E-U-c}{U-c},\quad
 z_a=-x(a_-),\quad z_b=-x(b)>1,\quad
 f(E)=\frac{(-1)^mT_m(x(E))}{T_m(z_a)},\quad
 M(E)=1-(1-f(E))^K.                                      \tag{8}
\]

Outside \(\Pi\), set \(M=0\). Both \(\Pi\) and \(H\) commute. Choose an
exact dyadic \(s>0\) with \(s^2\leq z_b^2-1<4s^2\). The inequalities
\(T_m(z_a)\geq(1+s)^m/2\) and
\(d\log T_m(z)/dz\leq m/s\) give a stop amplitude
\(\delta\leq2/(1+s)^m\) and pass loss
\(p\leq m(z_a-z_b)/s\). On the band, \(0\leq1-M\leq p^K\).
On the complement,
\(|M|\leq K\delta/(1-K\delta)\), by the binomial/geometric majorant.
The exact integer search gives

\[
 m=20951,\quad K=8,\quad\deg M=167608,\quad
 \|M-E_{k=0}\|\leq\varepsilon=10^{-69}.                   \tag{9}
\]

This also proves \(\|M\|\leq1\). The Kraus pair
\((M,\sqrt{I-M^2})\) defines a complete finite instrument. With an arbitrary
reference system its unhalved trace error relative to the exact binary
spectral instrument is at most
\(2(\varepsilon+\sqrt{2\varepsilon})\), from the operator Hölder bound on
the square roots and the two Stinespring isometries.

For input \(I/D\), the accepted density is exactly

\[
 \rho_M=\frac{\Pi M(H)^2}{\operatorname{Tr}(\Pi M(H)^2)}.
                                                               \tag{10}
\]

It commutes with the full interacting Hamiltonian. The failure density also
commutes with it. Writing \(r_0=128\), the original-input success satisfies

\[
 \frac{r_0(1-\varepsilon)^2}{D}\leq q\leq
 \frac{r_0}{D}+\varepsilon^2,
 \qquad r_0/D\simeq2.85183176571\times10^{-126}.
\]

The unhalved conditional trace error from \(E_{k=0}/r_0\) is bounded by
\([4\varepsilon+2(D/r_0)\varepsilon^2]/(1-\varepsilon)^2
<7.02\times10^{-13}\). The conditional energy-complement weight is at most
\((D/r_0)\varepsilon^2/(1-\varepsilon)^2\). An explicit polynomial replaces
any unknown-lowest-projector prescription. Its degree counts source-H
applications in a vector recurrence; it is not a native primitive-gate
compiler or a feasible preparation claim. The complete carrier is not
materialized. `replay.py` evaluates (8) on all sixteen levels of the actual
local block with seven spectator grounds at \(\lambda=0\); it is a numerical
diagnostic. The nonzero-\(\lambda\) full-carrier bounds are (2)-(10).

## Complete records, localization and calibration

The final energy detector resolves the explicit band interval and its
complete complement. Inside the band, the position effects \(P_y\) specify
the neutral local eigenline at \(y\), true ground doublets elsewhere and
blank links; \(\sum_yP_y=P\). An additional effect \(I-P\) retains every
virtual nonblank and other profile. Translation invariance of (10) gives a
common probability \(w\) at each cell, with
\((1-\eta^2)/8\leq w\leq1/8\), conditional on the band. The native
occupation detector after \(P_y\) has the exact probabilities \(p_\alpha\)
from (7). If \(h\) is the conditional band weight, the complete 51-record law
is

\[
 (1-q),\quad q(1-h),\quad qh(1-8w),\quad
 \{qhw p_\alpha:y=0,\ldots,7,\ \alpha=1,\ldots,6\}.
                                                               \tag{11}
\]

Its sum is one. All outcomes are source-time independent before the final
measurements. The high-probability atomic profile has every local matter-plus-
scalar gauge current zero; virtual charge/flux outcomes reside in the retained
complement. Exact global neutrality and conserved raw number hold throughout.
This is closed-source stationarity, not a lifetime bound under an external
interaction.

A localized input \(P_y/128\), followed by the band filter retaining all
eight momenta, has ideal herald probability at least \(1-\eta^2\). Its
conditional finite-filter unhalved trace error is at most
\(4\varepsilon/(1-\eta^2)\). Every initial vector in \(P_y\) differs from
its normalized band projection in the transverse direction by at most
\(\eta/\sqrt{1-\eta^2}\). Band evolution differs from a global phase by at
most \(\epsilon|t|\). Therefore the combined probability of every other
location and the profile complement is bounded by

\[
 \left(\frac{\eta}{\sqrt{1-\eta^2}}+
                 \epsilon|t|\right)^2
 +\frac{2\varepsilon}{1-\eta^2},                         \tag{12}
\]

capped at one. Every momentum outcome is retained. On the native translation
orbit \(P_yP_kP_y=P_y/8\), and the squared Hilbert-Schmidt norm of their
commutator is \(7/32\); simultaneous sharp position and momentum are not
asserted. For source spacing \(\pi\), the centered position observable has
norm \(\pi/2\). In-band stationarity bounds its mean velocity by
\((\pi/2)\epsilon\). All momentum energies lie in the same narrow interval;
individual dispersion eigenvalues and a positive motion lower bound are
unresolved. At zero hopping the band is exactly flat.

Use the supplied empirical electron energy 0.00051099895069 GeV to calibrate
the actual local charged gap, whose selected energy is the smallest root of
\(x^4-2x^3-212x^2+213x+11114\). Dividing that empirical input by the charged
local gap gives one common energy conversion. Applying it to (6) gives

\[
 [6.616386592241406,6.616386592244709]\times10^{-5}
 \ \mathrm{GeV}.                                        \tag{13}
\]

This is a calibrated excitation energy at momentum zero relative to the
actual ground. No measured neutral mass is fitted. If this entire band were
assigned the relativistic law \(E(k)^2=m^2+c^2|k|^2\) at the source wave
number one, its width forces
\(c^2\leq2m_{\max}\epsilon+\epsilon^2\), incompatible with unit model speed.
A relativistic physical mass shell, physical spin, cosmological abundance,
neutrino or dark-matter identification, and an independent time/length
calibration are not consequences of this finite member.
