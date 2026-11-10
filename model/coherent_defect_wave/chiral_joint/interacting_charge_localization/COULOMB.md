# Charged propagation and a Coulomb comparison

## What the native source proves

The fixed PW119 source gives a positive elementary odd response, actual
charge transfer to neighboring cells with a retained flux link, and the
single-unscreened-charge bounds in [PROOF.md](PROOF.md). Its unit link
excitation cost and `lambda=2^-34` make the selected low energy bins
strongly localized in that precise detector sense. The cloud outcome
retains all dynamical screening. Neither a free charged dispersion nor
a thermodynamic confinement statement follows from these finite bounds.

The charge used in the proof is the actual native hypercharge. The
electromagnetic generator `Q=T3+Y`, its normalization, an unbroken
abelian sector and its coupling to the surviving charged modes require
a declared electroweak field construction. The labels e^c and L alone
do not supply that construction or a physical charged-lepton dictionary.

## A finite spatial charged-band test

For a finite regulator with a translation-compatible physical charge
sector, choose gauge-invariant dressed operators D_x. Each D_x includes
the actual string or distributed flux and the retained compensator.
Use the actual interacting ground g and clock H of that regulator.
The matrix-valued positive spectral measure

`M_xy(J)=<g|D_x^dagger P_J(H-E_g) D_y|g>`

and its time transform give a direct test of spatial propagation. A
finite charged band would require energy intervals and nonzero residues
in the actual allowed momentum sectors, together with propagation of
wave packets built from those same intervals. The bin complement and
every pump failure remain part of the complete detector instrument.
Large time or large volume claims require corresponding error bounds.

A static compensator fixed at one cell changes the physical translation
problem: translating the source alone moves its joint Gauss constraint
away from that reference. A Fourier label on bare source matter is
therefore insufficient. Consistent alternatives include an explicit
boundary flux sector with its boundary moved outside the observation
region in a controlled limit, or a translation-covariant quantum
reference whose position, recoil and energy are retained. Translating
both the source and reference gives total momentum; identifying that
with the source excitation momentum requires a further decoupling
argument. Replacing a static reference by mobile compensating matter
changes the source and can instead produce a neutral composite band.

## The continuum charged threshold has an infrared boundary

Buchholz considers a local relativistic theory with its field-domain
assumptions, translation spectrum and Gauss law. For compactly supported
real test functions phi in the spacelike complement of zero, set
`F(phi_R)=integral d^4x R^-2 phi(x/R) F(x)`. The charged vector must
have convergent expectations of these operators and bounded asymptotic
mean-square fluctuations. Nonzero charge is determined by the limiting
Gauss flux. Under these hypotheses it cannot be an eigenvector of the
mass operator. See equations (1), (2a), (2b) and the argument in
[Buchholz, Physics Letters B 174 (1986), 331–334](https://doi.org/10.1016/0370-2693(86)91110-X),
also available in the
[author's paper](https://www.researchgate.net/profile/Detlev-Buchholz/publication/256576300_Gauss%27_law_and_the_infraparticle_problem/links/59e07400aca272386b73338b/Gauss-law-and-the-infraparticle-problem.pdf).

The theorem is conditional on a continuum charged sector and those
asymptotic assumptions. It is not a no-go theorem for finite matrix
eigenvalues. In a Coulomb field limit satisfying its hypotheses, the
appropriate target is a dressed or inclusive energy-momentum threshold,
with controlled detector resolution and charge transport. Requiring a
nonzero isolated charged mass-shell residue in that limit would impose
an incompatible criterion. A threshold does not by itself establish
propagation, Lorentz covariance or scattering.

## An explicit quadratic-Gauss comparison

A transparent comparison uses an incidence matrix D from edges to
vertices, a neutral charge vector rho and a real edge electric field E.
Declare `D E=rho` and quadratic electric energy `a_E ||E||^2/2`,
with `a_E>0`. On a connected finite graph write `G=(D D^T)^+` and
`E_L=D^T G rho`. Any admissible field is `E=E_L+E_T` with `D E_T=0`.
Orthogonality gives the exact variational identity

`min_(D E=rho) a_E ||E||^2/2 = a_E rho^T G rho/2`.

This replaces the linear cost of counting excited links by the energy
of flux distributed over all edges. On a translation-invariant cubic
torus, a pair `rho=q(delta_0-delta_y)` has energy
`a_E q^2(G(0)-G(y))`; the finite Green kernel is the Fourier sum with
denominator `4 sum_j sin^2(pi k_j/N)` and its zero mode omitted.
An inverse-distance interpretation requires a separate large-volume,
small-spacing estimate. The exact finite Green identity alone proves
neither a physical Coulomb law nor a quantum photon.

Adding canonical transverse coordinates A and the magnetic term
`a_B ||C A||^2/2`, with curl matrix C and `C D^T=0`, specifies a
Gaussian field comparison. After gauge and harmonic-zero-mode choices,
its transverse frequencies are the square roots of the positive
electric/magnetic quadratic-form eigenvalues. Both terms belong to
one declared H and clock. The current can be retained algebraically
through charge-covariant link coordinates and native matter hopping,
but the fully coupled ground and charged response must be derived for
that changed H. A quadratic field sector does not determine them.

A finite graph is not a finite-dimensional quantum carrier. Real
canonical link coordinates or oscillator modes have unbounded local
Hilbert spaces. A finite matrix cannot satisfy `[A,E]=iI` exactly, as
the trace of a commutator is zero. An occupation cutoff P introduces
boundary terms. For a declared coupled Hamiltonian H and a state psi,
one sufficient finite-time comparison, for t nonnegative, is

`||exp(-itH)psi-exp(-itPHP)Ppsi||`
` <= ||(I-P)psi|| + integral_0^|t| ||(I-P)HP exp(-isPHP)Ppsi|| ds`.

This states a required tail and boundary-amplitude control; it does
not assert those quantities are small. All cutoff losses must be
retained in the observer law. The Gaussian comparison here is analytic,
and is not offered as a finite executable replacement for PW119.

The repository supplies reusable, separately scoped mathematics:
[discrete_coulomb_green.py](../../../../code/electromagnetism/discrete_coulomb_green.py)
and its independent verifier check exact Green and Thomson identities
on a twelve-port, thirty-seam table. The
[Whitney quantum coefficient evaluator](../../../../code/electromagnetism/whitney_interacting_quantum.py)
uses a different finite geometry and continuous configuration variables;
it does not discretize the 56-dimensional wavefunction. The
[Whitney Maxwell proof](../../../../paper/tex_fragments/WHITNEY_MAXWELL_DYNAMICS.tex)
quantizes finitely many modes on an infinite-dimensional Fock space.
Those identities and their domain controls retain their own source,
state and clock. No charged-source transfer is established by citing them.

## A compact gauge architecture is another declared choice

An alternative declares electric Casimir energies and independent
magnetic plaquette couplings, with a representation cutoff and its
boundary controls, while retaining the exact endpoint current Ward
identities. This electric/magnetic Hamiltonian structure is the
[Kogut–Susskind construction](https://doi.org/10.1103/PhysRevD.11.395).
Changing the PW119 unit penalty, adding plaquette terms or changing
the link representation menu changes H. The actual vacuum, current
response, charged dressing and finite-clock records require a proof
for that member. A weak-coupling label is not such a proof.

Finite local gauge menus are not excluded from every emergent massless
phase. Fröhlich and Spencer prove a massless phase for their four-
dimensional U(1) lattice model and an intermediate massless phase for
their sufficiently large Z_N models. Their hypotheses concern those
specific spatial statistical models and couplings.
[Fröhlich–Spencer, Communications in Mathematical Physics 83 (1982), 411–454](https://doi.org/10.1007/BF01213610).
The PW119 menu and the one-edge compiled service are not identified
with either model. Importing that theorem requires an explicit action
or transfer-matrix equivalence and control of matter and regulator
changes. It does not follow from finite gauge invariance alone.

## Why matter loops do not provide a phase shortcut

The [native matter-induced square](../matter_induced_plaquette/PROOF.md)
uses unit hopping, a supplied one-fermion state and a retained dual
reference. It proves nonzero square-loop amplitudes and complete
finite-time records on that declared source. It does not prepare the
small-lambda interacting ground used here or isolate an autonomous
low energy gauge sector. Its reference orientation and reference
energy also differ from the elementary-removal setting.

In the small-lambda atomic split, P0 has every link blank. The
controlled low energy effective Hamiltonian acts on the 256 local
ground spins. A four-hop term in that space does not by itself create
a retained magnetic gauge coordinate. The frozen remainder control
`R4<=25(lambda*360)^4/g^3` bounds a full effective correction under
its finite-volume hypotheses. It is not an electric/magnetic coupling
ratio or a photon theorem. A loop amplitude alone cannot establish
a dense fluctuating-flux vacuum, Coulomb energy versus separation,
or a transverse massless spectral response.

The finite source route therefore has concrete additional observables:
reference-energy-controlled static charge energy differences, a
translation-compatible dressed charged spectral measure, and field
strength correlations with their transverse spectral support. Any
architecture choice used for these observables must be fixed before
computing its vacuum and complete detector law. Source current identity
and one Hamiltonian clock can be preserved by construction; their
preservation does not transfer a vacuum or a continuum theorem.
