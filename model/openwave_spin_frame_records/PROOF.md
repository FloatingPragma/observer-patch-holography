# Native field, frame and spin records

## The source and its actual rotations

The source is the periodic side-four kinetic/pairing construction of
`openwave_kinetic_phase_cycle`, with all 90 physical fermion modes per site,
all 3600 analyzer flags per physical mode, all three characters and every
massless spectator. The exact native frames in `openwave_lepton_mass_shell`
intertwine its charged branches with

\[
h_k(p)=-\tau_3\otimes\boldsymbol\sigma\cdot\boldsymbol\pi(p)
       +e_k\tau_1\otimes1,\qquad \pi_i(p)=\sin p_i.
\]

Their native charge is minus one, and the native spin currents restrict to
\(1\otimes\sigma_i/2\). The same measured coefficients \(e_k\), retrospective
phase choice and electron scale calibration enter this construction.

Let \(P_z(x,y,z)=(-y,x,z)\) modulo four and
\(r_z=\exp(-i\pi\sigma_z/4)\). The one-particle transformation
\(U_z=r_z\otimes P_z\), applied to every native internal and analyzer
coordinate, is the literal spatial quarter-turn and spin rotation.
The real-space finite-difference matrices satisfy
\(U_zKU_z^*=K\), and \(r_z\epsilon r_z^T=\epsilon\) preserves every onsite
pair term. The spin-independent analyzer, leg parity and spatial inversion
controls commute with this proper cube rotation. Thus every actual compiler
pulse and its complete output instrument are covariant. The corresponding
half-turn \(U_x=(-i\sigma_x)\otimes P_x\),
\(P_x(x,y,z)=(x,-y,-z)\), is also a symmetry.

On the full fermion algebra,

\[
\Gamma(U_z)^4=(-1)^N,\qquad \Gamma(U_z)^8=1.
\]

This follows from \(U_z^4=-1\) and the exterior-power representation; it
includes all occupation sectors and ancillary flags. The code checks the
complete 128-row spatial/spin matrices and a 16-row fermion lift with the
actual antisymmetric pair generator. Rotating only the two charged branch
coordinates does not establish the full-source identity. Continuous spatial
rotations are not symmetries of this lattice.

Choose the exact quadratic ground reference by filling all negative modes
and leaving the massless zero modes vacant. The occupied subspace is
invariant under proper cube rotations. Its determinant phase is one:
time reversal on each canonical Dirac or Weyl block has square minus one,
commutes with proper rotations, and preserves the negative subspace. The
restriction of a proper rotation is therefore a symplectic unitary, with
determinant one. Empty zero modes contribute one. The native particle-hole
change of coordinates has a filled two-spin reference on each right leg;
its determinant under a proper spin/spatial rotation is also one. Numerical
occupied-subspace determinants provide independent diagnostics on all three
native masses and the massless block. This reference has even parity; one
additional charged quasiparticle has odd parity.

## One preparation, autonomous evolution, complete records

Restrict the preparation, rather than the Hamiltonian, to the four momenta
\(p_x=p_y=0\). This gives a plane-uniform packet with a nontrivial longitudinal
profile. All other source coordinates remain present in their chosen ground
reference. A charged reservoir and an opposite idle reference supply a
charge-balanced creation instrument. For
\(f=|\beta+\rangle\otimes|+x\rangle\), prepare
\(\phi_p=P_+(p)f/2\), retaining the creation failure, and evolve with this same
\(h_k(p)\). This uses a declared positive-band preparation/filter service.
It does not identify its effects with unfiltered local particle positions.
Writing \(m=e_k\), \(E=\sqrt{1+m^2}\), \(u=m/E\), the original success
probability and the conditional plane-zero amplitude are

\[
w=(3+u)/4,\qquad
\psi_0(t)=\frac{2e^{-imt}+(1+u)e^{-iEt}}{4\sqrt w}
             |\beta+\rangle\otimes|+x\rangle .
\]

The two opposite nonzero sine momenta cancel the spin/chirality-odd component
at this plane. Its probability \(p=\|\psi_0\|^2\) is positive at every time,
with \(p\ge(1-u)^2/(16w)\). The profile oscillates with \(E-m>0\); the code
retains all four planes and four spinor outcomes at times zero, one and
three. The time-zero tail caused by the positive-band filter is included.
A proper quarter-turn leaves the axial spatial packet unchanged and rotates
its native plane spin vector. After \(j\) such rotations,

\[
s_i=\tfrac p2 n_i,\qquad
n=(\cos(j\pi/2),\sin(j\pi/2),0).
\]

These are matrix elements of the actual native spin currents through the
proved source frame. The field is measured on a spatial plane, not on a
point-localized three-dimensional particle.

Let \(P_{\rm exc}\) select this family's one-positive-quasiparticle sector on
the full chosen ground reference. For each detector axis \(a\), define the
positive effects

\[
E^a_\pm=P_{\rm exc}\bigl[|z=0\rangle\langle z=0|\otimes
 1_\tau\otimes(1\pm\sigma_a)/2\bigr]P_{\rm exc},\qquad
E^a_0=1-E^a_+-E^a_-.
\]

The complement includes every other sector; all analyzer outputs are retained.
The axis is chosen uniformly. Their conditional probabilities are
\(p^a_\pm=p(1\pm n_a)/2\), \(p^a_0=1-p\). Multiplication by \(w/3\), plus
the original failure \(1-w\), gives ten complete frame records.

Prepare two independent copies with this same instrument. There is no
copying of an unknown state. The observer's apolar tensor field is
\(Q_{ab}=s_as_b-\delta_{ab}|s|^2/3\), realized by products of native currents
on the two copies. For each uniformly selected axis pair \((a,b)\), measure
only the product outcome \(z\in\{-1,0,+1\}\), with effects
\(F^{ab}_z=\sum_{uv=z}E^a_u\otimes E^b_v\). The unmeasured degeneracies are
retained as quantum outputs. Its conditional law is

\[
q^{ab}_\pm=\tfrac{p^2}{2}(1\pm n_an_b),\qquad q^{ab}_0=1-p^2.
\]

The 27 success records weighted by \(w^2/9\), and the three explicit
preparation-failure classes, sum to one. This complete classical tensor law
has period two quarter-turns, while the oriented frame law has period four.
The full microscopic quantum output need not have the tensor's shorter
period. The tensor is a declared two-copy field readout; no fundamental
nematic field is inferred. Nontriviality is certified by

\[
\operatorname{TV}(q_0,q_{\pi/2})=w^2p^2/9>0,\qquad
\operatorname{TV}(p_0,p_\pi)=wp/3>0.
\]

Because \(u<m<1\), conservative all-time lower bounds are
\((1-m)^4/4096\) and \((1-m)^2/64\), respectively. The exact rational upper
endpoints of the native mass intervals give the receipt's strict bounds.

A retained pointer coherently controls four or eight source quarter-turns.
An odd success state gives the pointer phase minus one at \(2\pi\) and plus
one at \(4\pi\), with the same preparation failure retained. Rotating the
opposite odd reference too cancels the minus sign; dephasing the pointer
before readout gives equal probabilities. Both are executable controls.
These phases follow from the actual spatial/spin transformation above, not
from relabeling an abstract half-angle coordinate.

## Finite compiler transfer and spin integral

All pulse generators are even and covariant under both cube transformations.
The full-source parity response is therefore exact even with ancillary
excitation, provided all retained source modes are rotated. For the target
quadratic dynamics, the kinetic compiler bounds the unhalved complete-state
trace error by \(\varepsilon\le10^{-12}\) for all displayed times, including
references and nonvacuum analyzer outputs. One-copy record total variation
is at most \(\varepsilon/2\), and two-copy total variation is at most
\(\varepsilon\). Differences between two orientations lose at most
\(\varepsilon\) and \(2\varepsilon\), respectively. The positive lower
bounds above are much larger. No finite full-volume simulation is represented
as an execution of the enormous cycle count; the full-Fock theorem provides
the transfer. Declared ideal preparation, rotation and detector services have
no derived primitive repair implementation here.

For the preparation \(|\beta+\rangle\otimes|+z\rangle\), axial evolution
preserves the spin eigenvalue. The filtered native excess density obeys

\[
s_z(z,t)=\tfrac12\|\psi_z(t)\|^2,\qquad
\sum_z s_z(z,t)=\tfrac12.
\]

All position and energy outcomes, together with their original failure,
are retained. This is an excess canonical-spin integral. The unfiltered
source spin has ground fluctuations from the occupied sea. For one filled
negative charged block at momentum \(p\), its spin outcomes \(-1,0,+1\)
have probabilities \(a/2,1-a,a/2\), where
\(a=(\pi_x^2+\pi_y^2)/(2(|\pi|^2+m^2))\). Each occupied massless Weyl row has
outcomes \(\pm1/2\) with probabilities \((1\pm\pi_z/|\pi|)/2\). Vacant zero
modes contribute zero. Convolving every factor gives the receipt's complete
unfiltered spin law. Its ground mean is zero and its variance is positive.
Adding the spin-pure axial quasiparticle shifts this entire law by one half:
its orbital is a spin eigenvector orthogonal to the occupied subspace, so
the Slater-determinant characteristic function gains exactly
\(e^{i\theta/2}\). The resulting full state is not an eigenstate of the
unfiltered total spin. The lattice supplies no continuous spatial rotation
Noether generator, and no equality of this canonical spin integral with
orbital-plus-field angular momentum is proved.
