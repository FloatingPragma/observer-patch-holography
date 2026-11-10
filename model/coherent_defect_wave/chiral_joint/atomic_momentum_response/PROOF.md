# Source momentum and retained references

## Source, reference and clock

The source is the complete `peter_weyl_links.model.column` on an even periodic
cube of side \(N\), spacing \(\pi\), thirty native fermion modes and the
five-state quantum scalar at each site, and 119 states on every directed
nearest-neighbor link. All three native Yukawa coefficients, scalar and link
penalties, and the hopping coefficient are one. Every axial hopping word in
`periodic_family.kinetic_terms` is retained, including paths longer than one
edge when \(N>2\). The clock is \(e^{-itH_{\rm PW}}\). No scalar or link is
replaced by an expectation value.

The consumed `local_reference` theorem supplies the global zero-boundary-charge
physical ground projector \(P_g\) of the **single-cell** onsite Hamiltonian.
It has rank two, spin one half, raw lepton label one and odd fermion parity.
Its spin-invariant density is \(P_g/2\). The declared lattice preparation is
\[
 \rho_0=\bigotimes_x(P_g/2)_x\otimes|0_{\rm links}\rangle
 \langle0_{\rm links}|.
\]
This is a ground density of the decoupled onsite Hamiltonian with positive link
penalty. It is used here as a finite full-source quench preparation. No
stationarity or ground property under unit hopping is assumed.

The reference is physically retained. At each site introduce a static,
gauge-neutral fermionic spin doublet with precisely one reference fermion.
Choose the actual ground doublet basis so that its raising and lowering
operators have the usual spin-one-half matrices. The local preparation is
\[
 |s_x\rangle=(|g_\uparrow\rangle|\downarrow\rangle_R
             -|g_\downarrow\rangle|\uparrow\rangle_R)/\sqrt2.
\]
It is normalized, Gauss invariant, rotation invariant and even in the combined
fermions. A fixed reference raw lepton label minus one compensates the raw
source label; this is a reference assignment, not an identification with a
measured particle number. The product of these local singlets purifies
\(\rho_0\). All reference clocks are identity. A finite preparation isometry
and the readouts below are declared apparatus services. This package does not
supply an efficient local gate implementation of those services.

## Actual translations and their restriction

Use the actual antiperiodic (AP) source translations of the PW model. A mode
crossing a seam acquires a minus sign, and reordering occupied modes supplies
the full canonical anticommutation relation sign. Scalar and link coordinates
are permuted. These commuting unitaries obey
\[
 T_i^N=(-1)^{F_{\rm source}},\qquad T_i^{2N}=I.
\]
Consequently the complete source phase readout has \((2N)^3\) outcomes,
\(\exp(i\pi\ell_i/N)\). On an even source sector every \(\ell_i\) is even; on
an odd source sector every \(\ell_i\) is odd. Mixed-parity phase triples have
zero projectors. Zero outcomes are retained.

Each cell in the atomic ground subspace has odd parity. On each length-\(N\)
translation cycle, its graded block permutation contributes
\((-1)^{N-1}\) and its seam crossing contributes \(-1\). For even \(N\) their
product is one. Thus the restriction of every source translation to the
atomic subspace is the ordinary permutation of its \(N^3\) ground spin labels.
This argument applies to every odd local Fock occupation in each ground
wavefunction, irrespective of its number of fermions. It does not replace the
ground wavefunction by one bare fermion. The static site-reference array uses
the same AP convention and has the same spin-permutation restriction.

`model.shift_spins(mask, side, displacement)` implements this restriction:
site order is lexicographic \((x,y,z)\), and bit zero/one labels the ground
up/down state. `site_permutation` returns the corresponding site permutation.
The checks execute all 256 spin labels and all three actual source
translations at \(N=2\), and arbitrary odd local occupations at \(N=2,4\).
They also execute the central relation after a local parity change.

The full source commutes with these translations. The onsite source is an
identical even operator at every site. Under a one-step translation, the two
fermionic endpoints of every hopping word supply precisely the AP seam factor
of its translated kinetic coefficient, and its PW path is permuted to the
translated path. `native_translation_join` verifies every kinetic word at
\(N=2,4\), and full occupied-carrier columns with nonblank quantum links.
Translations permute the twelve Gauss generators among sites, preserving
their common kernel. They do not commute with each individual site generator.
The exact Fourier source readout is therefore a physical measurement on that
kernel, and all its probabilities are preserved by the full clock for every
real \(t\).

At side two, `rotation_translation_join` also executes every proper cubic
rotation on the actual eight-dimensional AP spatial one-fermion carrier.
For all 24 rotations and all three axes it checks
\(R T_iR^\dagger=T_j^{\epsilon}\), where
\(R e_i=\epsilon e_j\). The spin rotation commutes with spatial translation,
and the native internal index is unchanged. Second quantization lifts these
72 exact identities to every Fock sector. The checks retain the 72 scalar
permutation identities and all 1,728 oriented-edge identities. The latter use
the same full-state link orientation involution on both sides, supplying the
full 119-state bosonic lift. On the even side-two source sector,
\(T_j^{-1}=T_j\), so these rotations permute the three X phase sectors.
This is an operator symmetry statement; an interacting spectral gap or
isolated excitation cluster requires a separate theorem.

## Size-dependent source records

For a displacement \(a\in\mathbb Z_N^3\), its site permutation has order
\(d=N/\gcd(N,a_1,a_2,a_3)\) and \(N^3/d\) cycles. A spin basis vector is fixed
precisely when its label is constant on every cycle. The normalized character
and its complete Fourier readout are
\[
 \chi_N(a)=2^{N^2\gcd(N,a_1,a_2,a_3)-N^3},\qquad
 p_N(k)=N^{-3}\sum_a e^{-2\pi i k\cdot a/N}\chi_N(a).
\]
Grouping translations by their order and using elementary divisor inversion
eliminates complex arithmetic:
\[
 m_N(k)=\frac1{N^3}\sum_{d\mid N}2^{N^3/d}
 \sum_{\substack{e\mid d\\e\mid k_1,k_2,k_3}}
       \mu(d/e)e^3,\qquad p_N(k)=\frac{m_N(k)}{2^{N^3}}.
\]
The inner sum follows by imposing \(\gcd(d,b_1,b_2,b_3)=1\) on translations
\(a=(N/d)b\). The numbers \(m_N(k)\) are nonnegative integer multiplicities
and sum to \(2^{N^3}\). Exact cycle counts and all multiplicities are executed
at \(N=2,4,6\). For \(N=2\),
\[
 p_2(0,0,0)=23/128,\qquad p_2(k)=15/128\quad(k\ne0).
\]
The other 56 outcomes on the full \(\mathbb Z_4^3\) source phase grid have
probability zero. They are included in the receipt.

The pure singlet purification is invariant under **joint** source and
site-reference translation. Its complete joint records satisfy
\[
 \Pr(k_S,k_R)=\delta_{k_R,-k_S}\,p_N(k_S).
\]
Source momentum is mixed, while the sum is sharply zero. The explicit
\(N=2\) execution retains all 65,536 source/reference spin basis states, tests
the eight nonzero orthogonal branches, their sum, and all phase eigenvalues.
Changing the reference translation by a scalar \(N\)-periodic character
changes the joint label without changing any source probability. A control
with reference character \((1,0,0)\) executes this distinction. A reference
calibration cannot turn this mixed source momentum into a sharp source
momentum. These source and reference records remain distinct under the full
source clock.

## A normalized neutral Fourier pump

The actual native neutral composite is
\[
 J_{x,s}=\phi_{x,0}c_{x,s,11}-\phi_{x,1}c_{x,s,10}.
\]
It is an odd Gauss singlet, with its operator-valued scalar retained. On the
complete active fixed-spin carrier of two fermion modes times scalar dimension
five, direct exact matrices give
\[
 \operatorname{spec}\{J,J^\dagger\}=\{0^{(10)},1^{(4)},2^{(6)}\},\qquad
 \operatorname{spec}(J^\dagger J)=\{0^{(14)},1^{(4)},2^{(2)}\}.
\]
The other 28 local fermion modes are spectators. The graded embedding is
unitarily equivalent to this local matrix tensored with their identity.
In particular \(\{J,J^\dagger\}\le2I\), and \(J^2\ne0\). No canonical
fermion relation for this composite is used.

For each AP transfer phase \(\theta_i=\pi(2m_i+1)/N\), set
\[
 J_{\theta,s}^\dagger=N^{-3/2}\sum_x
     e^{-i\theta\cdot x}J_{x,s}^\dagger.
\]
Actual operator-column checks include both seam and interior sites, both
spins, every active occupation, every scalar state and both adjoints.
They give \(T_iJ_{\theta,s}^\dagger T_i^\dagger
=e^{i\theta_i}J_{\theta,s}^\dagger\).
Offsite odd composites graded-anticommute, so
\[
 \{J_{\theta,s},J_{\theta,s}^\dagger\}
  =N^{-3}\sum_x\{J_{x,s},J_{x,s}^\dagger\}\le2I,
 \qquad \|J_{\theta,s}\|\le\sqrt2.
\]
This full anticommutator inequality supplies the uniform-in-volume bound;
the local norm bound alone would not supply it.

Retain a separate initially empty neutral fermionic spin reference \(r\),
with raw lepton charge minus one for an occupied mode. The rectangular
success map from its vacuum to its one-fermion sector is
\[
 B_\theta=(J_{\theta,\uparrow}^\dagger r_\downarrow^\dagger
          -J_{\theta,\downarrow}^\dagger r_\uparrow^\dagger)/(2\sqrt2).
\]
It preserves Gauss, total spin and total parity. With all source modes ordered
before this pump reference, its matrix includes \(J^\dagger(-1)^{F_S}\)
before reference creation. The pump reference is ordered before the static
site-reference array. The initial source has even parity. The exact 160-by-80
active matrix checks the three spin intertwiners, parity and
\[
 B_\theta^\dagger B_\theta
  =\tfrac18\sum_sJ_{\theta,s}J_{\theta,s}^\dagger\le I/2.
\]
The effect formula is on the full **source** carrier with pump-reference
vacuum input. The same creation expression on arbitrary reference occupations
has additional occupation projectors and is not substituted for this
rectangular map. Those other reference input sectors are retained by the
unitary completion below.
The uniform bound follows from the Fourier anticommutators above. An exact
finite unitary implements success and failure:
\[
 \begin{pmatrix}\sqrt{I-B^\dagger B}&-B^\dagger\\
 B&\sqrt{I-BB^\dagger}\end{pmatrix}.
\]
Identity completes unused pump-reference sectors. The recorded failure map is
\(C=\sqrt{I-B^\dagger B}\). Functional calculus gives
\(CB^\dagger=B^\dagger\sqrt{I-BB^\dagger}\), so this block matrix is unitary.
The two output sectors are orthogonal, and \(C^\dagger C+B^\dagger B=I\).

## Actual ground weight and complete pump records

Let \(P_{\rm at}=\bigotimes_xP_g\otimes P_{\rm blank\ links}\). Cross-site
terms of \(B^\dagger B\) vanish after restriction by \(P_{\rm at}\), since they
change odd local parity at both sites. At the same site,
\(\sum_sJ_sJ_s^\dagger\) is an exact spin scalar. Schur's lemma on the unique
spin-one-half ground doublet therefore gives
\[
 P_{\rm at}B^\dagger BP_{\rm at}=wP_{\rm at},\qquad
 w=\tfrac18\operatorname{Tr}(P_gJ_\uparrow J_\uparrow^\dagger).
\]
This is a compression identity; invariance of the atomic subspace under the
full effect is not assumed.

The consumed `ground_composite_response` receipt computes this trace with
actual full-family ground projectors. Its individual setting has probability
\(1/4\) and success operator \(J_\uparrow^\dagger/2\), hence its joint weight is
\(\operatorname{Tr}(P_gJ_\uparrow J_\uparrow^\dagger)/32\). **The coherent
spin-singlet pump weight is four times that per-setting joint weight.** The
code reconstructs the exact algebraic weight from its complete cross-moments
and spectral projector polynomials, rather than substituting a decimal.
For the least real root \(e\) of
\[
 e^6-3e^5-310e^4+625e^3+31395e^2-31708e-1041468=0,
\]
the receipt records the rational function \(w(e)\) and the certified interval
\[
 \frac{11577757079357}{70368744177664}\le w\le
 \frac{5788878539679}{35184372088832}<\tfrac12.
\]
Both consumed receipts have their internal input hashes validated. Their
actual onsite coefficients and full native carrier are checked against this
source.

Write the even initial source phase as \(2k\) and the odd pump transfer as
\(\ell\), with all arithmetic modulo \(2N\). The complete nonzero records are
\[
\begin{array}{c|c|c|c|c}
\text{record}&\text{source}&\text{site reference}&\text{pump reference}
 &\text{probability}\\\hline
\text{success}&2k+\ell&-2k&-\ell&w p_N(k)\\
\text{failure}&2k&-2k&0&(1-w)p_N(k).
\end{array}
\]
The pump reference carries the declared inverse transfer character on its
occupied sector and phase zero on vacuum. Its \(N\)th translation power is
minus one when occupied, consistent with its odd parity. All other outcomes
of the success/failure flag and the three full phase triples have probability
zero and are retained. The receipt evaluates the complete records at
\(N=2,4,6\), for \(\ell=(1,1,1)\). The same formulas hold for every odd transfer
triple. Their probabilities sum to one because \(\sum_kp_N(k)=1\).

After success the **source** is odd and obeys \(T_i^N=-I\), so its momentum grid
is shifted by half a step in every direction. Failure has the even grid.
Source plus site-reference transfer is sharply \(\ell\); including the
compensating pump reference gives total phase zero. This calibration changes
neither the odd source grid nor its mixed momentum distribution. The full
clock commutes with source translations and acts trivially on the retained
references, so every displayed unconditional probability is preserved at
all real times. The pump and the phase projectors may be read in either
order after propagating their labels by these exact covariance relations.

These are finite operational momentum and transfer records around a
specified atomic preparation. No energy-momentum curve, interacting unit-hopping
vacuum, particle rest mass, continuum limit or empirical family identification
is inferred.
