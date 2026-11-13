# Prepared charged-line localization under the actual clock

The [sequential motion theorem](MOTION.md) uses the preparation and error
bound below to certify atom-resolved arrivals after a retained position
intervention. The [static hypercharge theorem](MAGNETIC.md) derives the
background conjugation identity for the same complete Hamiltonian, its
currents and its paired observer records. Neither calculation changes the
source reference or imports a free-particle clock.

## Source and preparation

Fix the complete interacting side-two source of
[`interacting_charged_response`](../coherent_defect_wave/chiral_joint/interacting_charged_response/PROOF.md).
It has 240 native fermion modes, eight five-state scalar factors and
24 complete 119-state link factors. Its Hamiltonian, reference and clock are

\[
H=H_0+\lambda T+4N_{\rm ref},\qquad
\lambda=2^{-34},\quad \|\lambda T\|\le v=360\,2^{-34},
\qquad U(t)=e^{-itH}.
\]

The static full CAR reference is at cell zero. The actual elementary
removal pump \(B_r\), for \(r=e^c,L\), transfers a native fermion into
that reference. It is a contraction with a retained failure output.
The finite physical ground-preparation service retains its second source,
all environments and its failure record. Its conditional source fidelity
with the exact interacting ground is at least \(1-d^2\). A purification
can therefore be chosen within norm \(2d\) of that ground and an idle
environment vector. All inequalities below include that environment.

Let \(u\) be the inherited atomic comparison vector. It lies in the
tensor product of the eight actual local ground doublets, with every link
blank. It is an eigenvector of \(H_0\), although it need not be a product
of the eight spins. Decompose \(B_r u\) into all its exact atomic energy
components, and write \(\phi\) for the lowest visible component. Its
positive weight \(w=\|\phi\|^2\) is calculated from the complete native
charged resolvent. At every cell other than zero the matter is a gauge
singlet and all links are blank, so the complete single-charge location
projector \(\Pi_0\) satisfies \(\Pi_0\phi=\phi\).

The energy-bin projector \(P_I(H)\) is the actual full-H projector with
radius \(1/1024\). It is not a projection of a truncated charged model.
For the actual prepared purification \(\Psi\), define the unnormalized
successful line \(\chi=P_I(H)B_r\Psi\). The source's spectral-theorem
comparison gives

\[
\|\chi-\phi\|\le\epsilon,
\qquad
\epsilon=\eta_g+2d+v\sum_j d_j^{-1},                 \tag{1}
\]

where \(\eta_g\) is the complete interacting-to-atomic ground-vector
bound and \(d_j\) are the retained distances of every atomic response
energy to the appropriate bin edge. In (1), \(\phi\) includes the idle
environment vector. Every atom and both spin contributions enter the
sum. The code rounds \(\epsilon\) upward and the atomic weight interval
outward using exact rational arithmetic.

For \(w\in[w_-,w_+]\), the actual line probability \(p=\|\chi\|^2\)
satisfies

\[
p_-=(\sqrt{w_-}-\epsilon)^2\le p\le
\min\{1,(\sqrt{w_+}+\epsilon)^2\},\quad
\sqrt{w_-}>\epsilon.                                \tag{2}
\]

These probabilities are conditional only on the separately recorded
physical preparation. Multiplication by its success probability gives
the original-input probability. No successful branch is silently
renormalized in the complete output table.

## The prepared line includes a bound on screening

The all-state bin theorem in
[`interacting_charge_localization`](../coherent_defect_wave/chiral_joint/interacting_charge_localization/PROOF.md)
bounds each distant *single* charge, while retaining an unbounded cloud
outcome. Equation (1) yields a different, preparation-specific statement:

\[
\Pr(I-\Pi_0\mid\text{line},t=0)
\le\epsilon^2/p_-.                                  \tag{3}
\]

The effect \(I-\Pi_0\) includes every distant single charge and every
screened cloud. Thus no union bound over only the single-charge outcomes
is used. This conclusion depends on the actual preparation in (1).
An arbitrary vector of the same bin can lie entirely in a cloud sector;
a degenerate three-state countercontrol demonstrates that distinction.

Let \(a\) be the exact atomic energy of \(\phi\), including the static
reference cost. Since \((H_0+4N_{\rm ref})\phi=a\phi\), Duhamel's
identity on the full carrier gives

\[
\|U(t)\phi-e^{-ita}\phi\|\le |t|v\sqrt{w_+}.
\]

Unitary evolution of (1), followed by \(I-\Pi_0\), therefore proves

\[
\boxed{\quad
\Pr(\Pi_0\mid\text{line},t)\ge
1-\min\left\{1,
\frac{(\epsilon+|t|v\sqrt{w_+})^2}{p_-}\right\}.
\quad}                                               \tag{4}
\]

This inequality holds for the actual source clock and the retained
mixed-state preparation, using its purification. It does not require
\([\Pi_0,H]=0\), an isolated full-H charged eigenvalue, purity of the
source marginal or exclusion of dynamical screening. The fixed
compensator and the small declared hopping coefficient remain part of
the construction. The bound is finite in time; it supplies neither a
thermodynamic localization theorem nor an asymptotically free electron.

For a normalized prepared-state disturbance of norm at most \(\delta\)
and a Hermitian generator change of norm at most \(\eta_H\), another
Duhamel estimate replaces the square root of the leakage bound \(L(t)\)
by \(\sqrt{L(t)}+\delta+|t|\eta_H\). The receipt evaluates
\(\delta=10^{-3}\) and \(\eta_H=v/100\). These are disturbances of the
declared construction, explicitly distinguished from its original clock.
A large-drive countercontrol makes the guaranteed localization trivial.

## State persistence, the clock and one energy calibration

For \(\psi=\chi/\|\chi\|\), the two preparation errors at the endpoints
give

\[
\|U(t)\psi-e^{-ita}\psi\|\le
D(t)=\min\left\{2,
\frac{2\epsilon+|t|v\sqrt{w_+}}{\sqrt{p_-}}\right\}.
\]

Consequently

\[
|\langle\psi,U(t)\psi\rangle|^2
\ge\max\{0,1-D(t)^2/2\}^2.                           \tag{5}
\]

The same norm estimate bounds the characteristic function's difference
from \(e^{-ita}\). A controlled-source Ramsey apparatus measures its
real and imaginary parts. That coherent apparatus is a declared finite
service, and is not a self-starting, frequency-rigid oscillator. The
source clock itself is the unchanged \(H\) evolution. Removing the known
reference phase \(4t\) and the initial ground phase is an explicit
comparison; neither energy is dropped from the Hamiltonian.

The lowest visible atomic gaps are approximately \(0.691656680938395\)
for \(e^c\) removal and \(0.096888497892281\) for \(L\) removal. The
\(e^c\) hole is a colour and weak singlet with native hypercharge \(-1\).
One optional unit convention maps the certified rational centre of its
atomic gap to \(0.00051099895069\) GeV. The same conversion is applied to
the complete energy bin and the \(L\) channel. The latter is a native
doublet transition; it is not assigned a neutrino or muon mass. These
GeV values are calibrated source transition energies. A relativistic
mass shell and a measured-particle dictionary are additional assertions.

For every time setting the complete coarse records are

\[
1-q;\quad q p_{\rm failure};\quad q p_{\rm outside};\quad
q p_{\rm line}w_0(t);\quad q p_{\rm line}[1-w_0(t)],
\]

where \(p_{\rm failure}+p_{\rm outside}+p_{\rm line}=1\).
Every source copy, quantum reference and preparation output is retained.
Energy and position are sequential measurements; their projectors are
not assumed to commute.

## Native full-turn phase on the actual prepared state

The native spin action on each cell's thirty CAR modes has generator
\(S_z=(N_\uparrow-N_\downarrow)/2\). It commutes with all twelve local
gauge generators. On the complete 240-mode Fock carrier,

\[
e^{-2\pi i S_z}=(-1)^{N_\uparrow-N_\downarrow}
=(-1)^{N_\uparrow+N_\downarrow}=(-1)^N.              \tag{6}
\]

The scalar and gauge-link factors have integer spin and contribute no
central sign. The atomic cube ground is even, since each local
spin-half ground doublet is odd and there are eight cells. All terms of
\(H\) and the physical preparation filter preserve source parity. Each
successful elementary removal is odd; its full-H energy filter and all
subsequent clock evolution preserve that odd parity. Equation (6) thus
acts as \(-I\) on the *actual* prepared line, with \(+I\) after \(4\pi\).
This is an action on the same native source, rather than an assignment
of a spinor label to a separate matrix model.

A retained pointer begins in \(|+\rangle\). Controlled source rotation,
followed by its second Hadamard, gives deterministic pointer one after
\(2\pi\) and pointer zero after \(4\pi\), conditioned on the line.
The original-input probabilities are respectively \(q p_{\rm line}\);
every preparation, pump and energy failure remains a separate output.
The reference has odd parity as well. Rotating source *and* reference
by \(2\pi\) produces \(+I\), which gives the opposite pointer result.
The code checks this control and the loss of the deterministic result
under pointer dephasing.

The ray period is \(2\pi\), and the amplitude period relative to the
unrotated pointer is \(4\pi\). In a half-angle coordinate
\(\varphi=\theta/2\), these are \(\pi\) and \(2\pi\); this coordinate
change does not derive an apolar local tensor field. Continuous spin
rotation at intermediate angles is an allowed intervention and need
not be a symmetry of the finite hopping. No measured electron angular
momentum, spin-statistics continuum theorem or passive repair-only
controller is inferred from (6).
