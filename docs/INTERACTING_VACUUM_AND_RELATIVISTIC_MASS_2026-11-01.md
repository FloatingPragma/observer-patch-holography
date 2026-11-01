# An interacting reference and a controlled relativistic mass construction

The particle work now separates three calculations that have different
answers: the ground of the original three-dimensional cube, the charge gap
of its one-dimensional restriction, and a relativistic limit along an
explicitly constructed coupling trajectory. The first two use the original
couplings. The third specifies how those couplings change under refinement.
That additional choice is recorded rather than treated as a prediction of
the fixed source.

This continues the [resolved spectrum and massless-sector work](RESOLVED_SPECTRUM_AND_RELATIVISTIC_SECTOR_2026-11-01.md)
and the five-issue campaign. It advances the interacting-field and particle
tasks, [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142)
and [#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143).
The earlier results for [#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138),
[#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139)
and [#1145](https://github.com/FloatingPragma/observer-patch-holography/issues/1145)
retain their stated domains. No remote issue is closed by these local changes.

## The actual interacting cube ground

The [global ground certificate](../model/interacting_vacuum/README.md) covers
all 65,536 Fock states of the unchanged open `2×2×2` cube, at `c=1/8,U=1`.
Its ground energy is exactly

\[
E_G=-\frac{\sqrt3}{4}=-0.433012701892\ldots .
\]

The ground space has dimension five, all at particle number four, with
hidden internal spin two. This internal symmetry has not been identified
with physical Lorentz spin. Exact negative-band and contact-kernel
calculations determine the five states. In the other potentially competing
number sectors, integer Gram certificates prove lower bounds; floating
Cholesky proposes those certificates but does not supply the proof.

For a specified pure ground, local up-component addition and removal each
have weight `1/2` at positive excitation energy
`λ=√3/16=0.108253175473…`. Down-component addition has the exact orthogonal
decomposition

\[
\mu=\frac1{10}\delta_\lambda+
\frac1{10}\delta_{1+\lambda}+\frac45\mu_r.
\]

The residual is a normalized positive measure with mean
`1/2−√3/64` and variance `1069/4096`. All its energies exceed `7/1024`.
The displayed atom weights are exact symmetry contributions; coincident
residual weight has not been excluded. The residual spectrum is not fully
resolved. These are excitation energies above an actual finite ground,
which removes the earlier excited-reference negative-energy obstruction
for this calculation. A two-site-wide box does not determine a relativistic
dispersion or rest mass.

The declared repair completion acts as identity on the entire
`N=3,4,5` sum, including coherences and external references, because its
active support lies in `N=6`. This is passive preservation by that specified
repair law. Native ground preparation and an active correction mechanism
have not been constructed.

## What the original one-dimensional law generates

An exact onsite change of frame maps the periodic `(L,1,1)` restriction,
`L` a multiple of four, to the repulsive Hubbard chain with `t=c/2=1/16`
and `U=1`. The [charge-gap calculation](../model/mass_generation_audit/README.md)
uses the explicitly imported thermodynamic Hubbard result and certifies its
value with rational series and tail bounds:

\[
\Delta_c=0.7715517157380375\ldots .
\]

For the half-filled canonical reference, the original-H addition and
removal costs are respectively `0.8857758578690188…` and
`−0.1142241421309812…`. Consequently that reference is not the original
zero-chemical-potential ground across all number sectors. The centered
generator `K=H−(U/2)N` has equal charge thresholds
`Δc/2=0.3857758578690188…`. This reference convention must remain part of
any subsequent charged-particle interpretation.

Under refinement with fixed speed `v`, `t=v/(2a)`, and fixed `U`, the
charge gap tends to zero. Its positive lattice value therefore cannot be
carried unchanged into a relativistic mass table.

## A specified trajectory with relativistic rest masses

The [Hubbard relativistic construction](../model/hubbard_relativistic_sector/README.md)
uses exact thermodynamic quasiparticle energy and momentum formulas,
with the repulsive/attractive transformation and normalization stated
explicitly. The source is [Melzer's Hubbard scaling analysis](https://arxiv.org/html/cond-mat/9410043),
not an assumption that a nonzero gap is already relativistic.

Write `A=2πt/U`, let `M>0` denote a limiting rest **energy**, and choose

\[
a(A)=\frac{4v e^{-A}}{M\sqrt{2\pi A}},\qquad
t(A)=\frac{v}{2a(A)},\qquad U(A)=\frac{\pi v}{a(A)A}.
\]

The exact Bessel-series energy and momentum converge uniformly on bounded
rapidity intervals to

\[
E(\theta)=M\cosh\theta,\qquad
p(\theta)=\frac{M}{v}\sinh\theta.
\]

Thus this constructed limit supplies the missing mass-shell calculation:

\[
E^2-v^2p^2=M^2,
\qquad m_{\rm rest}=m_{\rm inertial}=\frac{M}{v^2}.
\]

Both members of the charge quasiparticle doublet have this same mass.
The estimate controls finite parameters: for all `|θ|≤1`, the absolute
error in `(E²−v²p²)/M²` is below `0.034` at `A=64` and below `0.017`
at `A=128`, and tends to zero. Rational Bessel enclosures and an explicit
series-tail bound certify these errors.
The theorem concerns thermodynamic quasiparticle branches of the centered
half-filled theory. The finite-ring sector selection rules still apply;
an isolated native one-particle detector pole is not inferred. The
one-dimensional field theory is not a three-dimensional particle sector.

There is a concrete trajectory through the original couplings. Set

\[
v=\frac18,\quad A_0=\frac\pi8,\quad
M_0=\frac{e^{-\pi/8}}\pi.
\]

Substitution gives `a(A0)=1`, `t(A0)=1/16`, `U(A0)=1` exactly. Along
this chosen continuation, the limiting values are

| Quantity | Value in the declared model units |
| --- | ---: |
| Charge-quasiparticle rest energy `M0` | 0.2149329913552644… |
| Charge-quasiparticle rest mass `M0/v²` | 13.75571144673692… |
| Opposite-charge partner rest mass | 13.75571144673692… |

This is a computed result for an anchored construction, **not a uniquely
predicted mass of the original fixed model**. Other coupling trajectories
through the same starting point can have different limits. In particular,
the fixed-`U` continuation has vanishing charge gap. The original finite
half-gap and this limiting rest energy are different quantities, and no
GeV conversion or measured-species fit has been supplied.

The nonuniqueness has an explicit witness: for any desired `M>0`, multiply
the displayed `a_M(A)` by `1+(M/M0−1)exp[−(A−A0)]` and retain the same
relations for `t` and `U`. The factor is positive, enforces the same starting
point, and tends to one. The resulting relativistic rest energy is `M`.
Thus computing the fixed bare point more precisely cannot determine the
missing refinement law.

## What has closed, and what has not

The finite cube now has an exact global interacting ground and positive
charged spectral data. The source-bound chain now has a certified
interaction-generated charge gap. A declared running-coupling continuation
has a controlled relativistic mass shell, so its rest mass is calculated
rather than guessed from a finite gap.

The remaining physical work is to derive or specify the actual
three-dimensional refinement law and its stable reference, preserve the
interactions and protected observables on that same carrier, construct the
measured charge/spin and detector correspondence, and fix the common units.
No electron, muon, tau, photon, or neutrino identification follows from
the equal-mass one-dimensional charge doublet. The required interacting
3+1-dimensional field theory and Einstein coarse-graining remain open.
All imported mathematical results, reference choices and refinement inputs
remain visible in the claim and premise records.

## Validation

The three new packages pass 83 focused tests and fresh receipt replays.
Integration regression passes 466 tests and 22 subtests; 23 program-class
checks also pass. The final map contains 62 realized scoped attachments
among 67 reviewed domains. All 84 canonical premise dispositions and 346
claim classifications retain their prior status. The source-hash check,
claim registry, public quantitative-surface policy and generated premise
views pass. All eight paper PDFs build with no overfull boxes, unexplained
underfull boxes, or reference, citation, glyph or font errors; their new
ground and mass sections were rendered and visually inspected.

These checks establish the recorded computations and their integration.
They do not establish the remaining physical hypotheses or close the five
GitHub issues.
