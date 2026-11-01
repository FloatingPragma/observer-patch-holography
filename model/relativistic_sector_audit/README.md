# A controlled massless relativistic sector of the original action

The original Weyl/contact family has a controlled **eight-valley massless
Weyl limit on every fixed finite-particle, bounded-momentum domain**. This
package constructs its isometry and controls the *full unprojected interacting
evolution*. It outputs model energies, helicities and zero continuum mass.
The contact interaction vanishes in this limit. A positive-energy interacting
vacuum, massive relativistic sector, physical species and SI units are not
thereby constructed.

```sh
.venv/bin/python -W error -m pytest -q model/relativistic_sector_audit
.venv/bin/python -W error -m model.relativistic_sector_audit --write
.venv/bin/python -W error -m model.relativistic_sector_audit --check model/relativistic_sector_audit/results.json
```

## Source and what is new

The action is exactly the supplied [Weyl/contact action](../interacting_field/DERIVATION.md):

\[
H_a=d\Gamma(h_a)+U(a)D,\qquad
h_a(k)=-\frac va\sum_{j=1}^3\sigma_j\sin k_j,\qquad
D=\sum_{x\in\mathbb Z^3}n_{x\uparrow}n_{x\downarrow}.
\]

The scaling choice is `c=v/a`, with positive fixed model speed `v`; each
`U(a)` is finite and nonnegative, matching the supplied source constructors.
There is **no scalar subtraction, chemical-potential
shift, mass term or time reversal**. In particular, this is a different
comparison from the reflected, critically tuned [moving pair](../moving_bound_pair/README.md).
The local action and its native implementation remain declared simulator
inputs, not inferred physical constants.
This is a growing-lattice scaling branch, distinct from the fixed eight-site
protected-cube experiment. The same local formula does not prove that these
wave packets survive that cube's five repair seams; no such protection or
identification with its eight measured response lines is assumed.

The earlier derivation already counted eight zeros and bounded the one-body
sine remainder. The new attachment is the simultaneous eight-valley
continuum-to-lattice Fock isometry, a countable-lattice contact estimate,
and an operator-norm dynamical comparison for every fixed `Nmax`, including
all contact escape and arbitrary untouched references. It also gives
controlled spectral projectors, explicit valley spin frames, a Lorentz
Clifford check and exact energy intervals.

`valley_certificate()` checks every one of the 128 Fourier columns of the
actual periodic `4×4×4` hopping matrix using Gaussian rational arithmetic.
No doubler or negative branch is removed. Its finite torus is a source
identity check, not the infinite-volume limit used below. The separate
scalar regulator has a supplied massive Klein–Gordon comparison; its
different action and mass input are not transferred to this source.

## Explicit continuum embedding and relativistic algebra

Let `C=[-Lambda,Lambda]^3`, `a Lambda < pi/2`, and
`nu=pi b`, `b in {0,1}^3`. The eight patches `nu+a C` are disjoint on the
momentum torus. Use measures `d^3p/(2pi)^3` and `d^3k/(2pi)^3` and set

\[
\mathcal K=\bigoplus_{b\in\{0,1\}^3}L^2(C,\mathbb C^2),\qquad
(W_af)(\nu+ap)=a^{-3/2}f_b(p),
\]

with zero output outside these patches. This is an exact isometry by
change of variables. There is no extra `1/sqrt(8)` factor. Its fermionic
extension `W_a^{<=Nmax}=direct_sum_{n=0}^{Nmax} wedge^n W_a` is isometric,
including superpositions of particle numbers; it does not perform a
projection during time evolution.

Write `s_j=(-1)^{b_j}` and `chi_b=-s_1 s_2 s_3`. The limiting blocks are

\[
h_b(p)=-v\sum_j s_j\sigma_jp_j,\qquad
E_\pm(p)=\pm v|p|,\qquad E^2-v^2|p|^2=0.
\]

Thus every block has **continuum mass zero**. There are four `chi=+1`
and four `chi=-1` blocks. The proper diagonal rotation
`R_b=-chi_b diag(s)` has an explicit constant SU(2) lift: the identity or
`i sigma_j` about its positive axis. Define spin matrices
`S_{b,j}=R_{b,jj} sigma_j/2`; they obey the usual rotation algebra and
`h_b=2 chi_b v S_b.p`. Positive and negative energy have mathematical
helicities `chi_b/2` and `-chi_b/2`, respectively. At `p=0` both spin states
have zero energy and helicity is undefined.

After these constant spin rotations an opposite-chirality pair has kinetic
matrices `alpha_j=diag(-sigma_j,sigma_j)`. With off-diagonal identity `beta`,
`gamma=(beta,beta alpha_1,beta alpha_2,beta alpha_3)` obeys all sixteen
Clifford identities for signature `(1,-1,-1,-1)`. This is the usual massless
Dirac/Weyl equation in time coordinate `v t`; it is more than a scalar
dispersion analogy. Pairing valleys does not select four physical Dirac
species. The standard Weyl/helicity interpretation is reviewed in
[Tong's field theory notes, §§4.4 and 4.7](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S4.html).
The regulator and its fixed momentum cutoff are not Lorentz invariant;
the limiting differential expression has this covariance.

## Uniform bounds including the actual interaction

On `|p_j|<=Lambda`, the sine remainder and the Pauli norm give

\[
\|h_a(\nu+ap)-h_b(p)\|
\le\frac{va^2}{6}\sum_j|p_j|^3
\le\frac{va^2\Lambda^3}{2}.
\]

For `a Lambda<=1` and `p!=0`, the positive/negative spectral projectors obey
`||P_{a,b,±}(p)-P_{b,±}(p)||<=a^2 Lambda^2/6`. Indeed, the vector sine
remainder is at most `(a^2 Lambda^2/6)|p|`, and half the distance between
the two unit vectors is at most their original relative error. No spectral
projector is assigned at the exact node. The energy squared also obeys

\[
-\frac{v^2a^2}{3}\sum_jp_j^4
\le E_{a,\pm}(p)^2-v^2|p|^2\le0.
\]

Here is the full infinite-lattice contact argument. Let `P_a=W_a W_a*`
and `Omega_a` be its Fock range projection through `Nmax`. The spin-complete
band has site-diagonal projection

\[
f_a=\langle\delta_x,P_a\delta_x\rangle
       =8(a\Lambda/\pi)^3,
\qquad d=\lfloor N_{\max}/2\rfloor.
\]

On the `N<=Nmax` sector, the positive countable partial sums defining `D`
are bounded by `d`; they converge strongly and `D^2<=dD`. For a band state,
annihilation is `a(P_a delta_x)` of norm `sqrt(f_a)` and leaves the remaining
particles in the band. Consequently

\[
\langle\psi,D\psi\rangle
=\sum_x\|a_{x\downarrow}a_{x\uparrow}\psi\|^2
\le f_a\langle\psi,N_\downarrow\psi\rangle,
\]

and likewise with up/down exchanged. Decompose by the two number operators
(which commute with `Omega_a,D`) to obtain `Omega_a D Omega_a<=f_a d`.
This handles odd `Nmax` correctly. Hence

\[
\|U(a)D\Omega_a\|\le |U(a)|d\sqrt{f_a}.
\]

This is the countable-volume version of the [finite dilute-band estimate](../dilute_momentum_ward/README.md),
with a proof rather than a limit assumed from finite tori. Its discrete
counterpart uses the exact discrete fraction, not the continuum volume.
Combining it with the exterior-power one-body bound gives

\[
\begin{split}
\|e^{-itH_a}W_a^{\le N_{\max}}
 -W_a^{\le N_{\max}}e^{-itd\Gamma(h_{\rm W})}\|
&\le\min\{2,|t|\epsilon_a\},\\
\epsilon_a&=N_{\max}va^2\Lambda^3/2
 +|U(a)|d\sqrt{8(a\Lambda/\pi)^3}.
\end{split}
\]

For each fixed `a` both generators are bounded on these domains. Duhamel's
formula integrates `(H_a W_a-W_a dGamma(h_W))`; the *left* evolution is the
full `H_a`. The band need not be invariant under contact. In particular,
the estimate includes every outside-band scattering state. It remains
valid after tensoring any untouched reference. Half trace distance of
corresponding complete normalized outputs is at most the vector bound,
and so is each probability error, capped by one. This does not automatically
compare arbitrary extra instruments which change the preparation domain.

At fixed `Nmax,Lambda,v,t`, `U(a)=o(a^{-3/2})` makes the bound vanish.
For `Nmax=1` the contact is exactly zero for *every* `U(a)`. At fixed
`U=1,Nmax=2,Lambda=v=t=1,a=1/64` the bound is less than **0.001237**;
the one-particle projector error is at most `1/24576`.
`continuum_budget()` computes outward rational bounds using a Machin-series
enclosure of pi. The critical pair scaling `U~a^{-2}` is not controlled by
this estimate. No finite-density or uniform full-Fock limit is claimed.

The finite negative control uses all eight zero nodes on a `4^3` torus.
An exact normalized opposite-spin pair has `||D psi||^2=1/8` and
`||(1-Omega)D psi||^2=7/64`, proving contact escape is nonzero. At time
`1/8`, direct propagation of all 8,128 two-particle states gives vector
error about `0.0441641`, below the rigorous `sqrt(1/8)/8`. This is a
floating matrix diagnostic beside the exact identities, not a new native
word or large-volume preparation execution.

## Resulting energies, masses and the remaining physical boundary

At `v=1,a=1/64`, all eight valleys have the following model energies.
Both signs are retained; the table displays their common positive magnitude.
The receipt contains exact rational enclosures behind every rounded entry.

| Momentum `p` | Lattice `|E_a|` | Continuum `|E|` | Continuum mass |
| --- | ---: | ---: | ---: |
| `(0,0,0)` | 0 | 0 | 0 |
| `(1,0,0)` | 0.999959310393 | 1 | 0 |
| `(1,1,0)` | 1.414156018578 | sqrt(2) | 0 |
| `(1,1,1)` | 1.731980331101 | sqrt(3) | 0 |

These are generator energies relative to the original empty-sector zero;
they are not fitted masses. At finite spacing, the nonzero mass-shell
defect is a regulator error, not a positive or tachyonic particle mass. A constant
Hermitian `2×2` matrix anticommuting with all three kinetic Pauli matrices
must vanish, checked as a rank-four exact linear system. This excludes
that single-valley mass insertion only; valley mixing, other interactions
and genuine many-body mass generation need their own construction.

The empty Fock state has energy zero but is not a ground state: on the
finite source a one-particle mode at `(pi/2,0,0)` has energy `-c`, and its
contact term vanishes at every `U`. This is a comparison across number
sectors, not number-preserving work extraction from the empty sector.
In the uncut continuum expression the negative branch is unbounded below.
A filled-vacuum/antiparticle representation and its interacting stability
have not been provided by this empty-Fock comparison. Nor has the source
number current been identified with physical electric charge, or the
valley spin with the observed particle spin.

Thus the positive result is a **controlled massless relativistic dynamical
sector** with explicit energy, spinor and error data. It is free in the
proved limit, includes four valleys of each chirality, and supplies no
massive physical species. Physical species assignments remain zero. The
earlier calibrated dimer masses, supplied scalar mass, and nonrelativistic
pair mass cannot fill those missing attachments.
