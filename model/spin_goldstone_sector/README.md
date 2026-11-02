# Ordered and gapless internal-spin sector

The fixed-spacing, three-dimensional antiferromagnetic spin limit constructed in
[FA90](../protected_spin_limit/README.md) has an ordered ground sector and genuine
positive excitation spectrum arbitrarily close to zero. This package gives an
explicit infrared integral certificate, a ground-state selection construction,
and a positive-energy transverse spectral-window bound. The spin is the source's
internal SU(2) spin. It is not an identification with spatial or Lorentz spin.

## Domain and imported mathematical results

The Hamiltonian is

\[
 H_J=J\sum_{\langle xy\rangle}(\mathbf S_x\cdot\mathbf S_y-1/4),
 \qquad J>0,\quad S=1/2,\quad x\in\mathbb Z^3.
\]

Each unoriented nearest-neighbor bond occurs once. The subtraction is a scalar
and does not change states, correlations, commutators or excitation energies.
Finite approximants below are even periodic cubes, with side at least four.
The boundary extension is explicit; its infinite local interaction is exactly
FA90's cubic interaction. The lattice spacing stays fixed. No limit in which
that spacing vanishes is taken.

We import the following precise results, rather than an assumption of magnetic
order or of a spin-wave approximation.

* Dyson, Lieb and Simon, *J. Stat. Phys.* **18**, 335–383 (1978),
  [original author-hosted paper](https://math.caltech.edu/SimonPapers/68.pdf),
  equations (16),(21),(42),(43), Theorems 4.1–4.2 and the antiferromagnetic
  Lemma 6.1/proof of Theorem 6.1: finite even cubic tori have the normalized
  Duhamel bound `D_beta(S_-k^a,S_k^a)<=1/(2 beta J E(k+Q))`.
  These infrared bounds hold for every spin magnitude, including 1/2;
  the larger-spin restriction in their eventual ordering theorem is not a
  restriction on this input. We derive the needed zero-temperature bound
  below, retaining the factor two in the Duhamel spectral formula.
  Kennedy, Lieb and Shastry's energy-sum-rule argument
  (*J. Stat. Phys.* **53**, 1019–1030, 1988,
  [DOI](https://doi.org/10.1007/BF01023854), Section 2, equations (1)–(4))
  motivates the improved order estimate. We reproduce its finite bound
  from the accessible DLS input and replace the rounded numerical integral
  by an exact rational upper bound.
* Koma, [arXiv:1712.09018v4](https://arxiv.org/abs/1712.09018v4),
  Theorems 2.1 and 2.2: for the reflection-positive isotropic hypercubic
  antiferromagnet, the limit of an infinitesimal staggered field realizes the
  maximum spontaneous staggered magnetization; if this is positive, an
  infinite-volume **pure** ground with that magnetization and gapless
  positive excitation spectrum exists. This includes S=1/2,d=3. No uniqueness
  of all infinite-volume ground states is required.
* Koma, [arXiv:2105.04970v2](https://arxiv.org/abs/2105.04970v2),
  equations (3.4), (4.1), (5.3)–(5.12): for the same field-selected ground,
  transverse susceptibilities obey `chi_a(k)<=1/(2J E(k+Q))`; smooth
  momentum packets and smooth strictly positive-energy filters pass to the
  infinite-volume ground with their two-point norms and first energy moments.
  His Theorem 2.1 is consistent with the explicit constants derived below.
  Our bounds concern packet spectral weight and energy, not a sharp pole.

The hypotheses are matched as follows. The local algebra is the full spin-1/2
matrix algebra; the graph is bipartite cubic; exchange is isotropic and positive;
there are no further-neighbor, frustration, anisotropy or magnetic terms in the
limiting dynamics. For reflection positivity, rotate the odd sublattice by pi about internal
axis 2: the real S1,S3 components then have negative exchange and the
imaginary S2 component positive exchange, exactly as in DLS Lemma 6.1.
For the Perron–Frobenius field argument, the rotation about axis 3 instead
makes the flip-flop exchange negative and the axis-1 staggered field a
uniform negative transverse field. `algebra.py` checks both distinct
frames and their field identities. In particular a positive field connects the entire
finite spin basis by single flips and has a unique finite ground. A field is
used only to select states; the limiting dynamics has field zero.

## Certified nonzero order

First we check the finite normalization. DLS defines
`D_beta(A*,A)=Z^-1 integral_0^1 Tr(exp(-s beta H) A* exp(-(1-s)beta H) A) ds`.
The real/imaginary reflection lemma and Gaussian bound give the stated
`1/(2 beta J E(k+Q))` bound. In the finite spectral expansion (their
Eq.21), multiplication by beta and zero-temperature passage leaves
`chi(k)+chi(-k)=2chi(k)`, by inversion symmetry. Ground-to-ground
matrix elements must vanish for k different from Q, since otherwise
`beta D_beta` would diverge. Consequently

\[
 \chi(k)\le{1\over4J E(k+Q)}.
\]

In the invariant zero-field ground projector,
`<S_x^a S_y^a>=-e_L/9` on a bond; the exact double commutator below is
`4J e_L E(k)/9`, and its two equal first-moment contributions give
`M1(k)=2J e_L E(k)/9`. Cauchy–Schwarz between the positive first and
inverse moments now gives, for every finite k different from Q,

\[
 g_L(k)\le\sqrt{{e_L\over18}}\sqrt{{E(k)\over E(k+Q)}}.
\]

This factor is not inferred from the weaker positive-field bound used
later; the factor two in the Duhamel limit is essential. The finite
spectral expansion works equally for a degenerate invariant ground projector.

Normalize momentum integration by `(2pi)^-3`. For a rotation-invariant
zero-field ground sequence, write `g` for one component's structure measure,
`b=g({Q})`, and `h(k)=-(cos k_1+cos k_2+cos k_3)/3`. The energy sum rule is
`integral h dg=e/9`; also `h(Q)=1`. Applying the imported infrared bound to the
positive part of h gives

\[
 b\ge e/9-\sqrt{e/18}\,I,\qquad
 I={1\over\pi^3}\int_{[0,\pi]^3}
 \left({C(x)\over3}\right)_+
 \sqrt{{3+C(x)\over3-C(x)}}\,d^3x,
 \quad C(x)=\sum_i\cos x_i.
\]

The product Néel trial vector gives `e>=3/4`. This is a variational energy
bound, not a claim that the product vector is a quantum ground. The displayed
lower bound is increasing for `e>=3/4` when `I<2/5`. Hence

\[
 b\ge {1\over12}-{I\over\sqrt{24}}
 >{1\over12}-{2\over5\sqrt{24}}>{1\over1000}.
\]

`infrared.py` proves `I<0.382612<2/5` with integer/rational arithmetic. Divide
the positive octant into 64 cubical cells per direction. The integrand is
increasing in C, so the lower-coordinate corner bounds every nonsingular cell
from above. Cosines at those corners are bounded by their degree-40 Taylor
polynomials plus the absolute degree-42 remainder, using rational bounds on
pi obtained from Machin's identity and alternating arctangent series. Each
square root is rounded upward by an exact integer-square-root inequality.
Permutation multiplicities enumerate all 64^3-1 regular cells.

For the single origin cell of width `h=pi/64`, use
`1-cos x>=2x^2/pi^2`, so the integrand is at most `pi sqrt(3)/|x|`.
Enclosing the cell in a ball of radius `sqrt(3)h` gives normalized integral
at most `6pi sqrt(3)/64^2<33/64^2`. This deliberately loose bound is included,
not dropped or sampled. The receipt retains the exact upper fraction.
The 32-division control fails to certify `I<2/5`; the checker does not promote
an insufficient mesh.

Here the order bound also holds for the full finite-torus staggered averages;
this does not follow merely from an atom in a weak limit. Explicitly the
finite-volume version of the imported infrared bound, after taking zero
temperature at fixed even side length L, gives
`b_L>=e_L/9-sqrt(e_L/18)*I_L`, where `I_L` is the normalized momentum sum
for the displayed integrand with the single singular point Q omitted. The
bound is valid for the rotation/translation-invariant ground projector.
We have `3/4<=e_L<=9/4`. Away from Q the Riemann sums converge normally.
Near Q the integrand is at most `pi sqrt(3)/|k-Q|`; the integer shell
`||n||_infinity=r` contains `24r^2+2<=26r^2` points. Its contribution
is at most `13 sqrt(3) r/L^2`. Summing to `r<=delta L/(2pi)` proves a
uniform O(delta^2)+O(delta/L) bound. Therefore `I_L` converges to I, and
the preceding monotone bound gives `liminf b_L>1/1000` for full-box averages.

Let `M=N^-1 sum_x (-1)^x S_x^1` and `b_N=<M^2>`, with `N=L^3`.
Take G to be a purification of the invariant finite ground projector; all
physical operators act on its spin factor. Thus no finite-ground uniqueness
is required. Symmetry gives `<M>=<M^3>=0`, `||M||<=1/2`.
The normalized biased vector `(1+M)G/sqrt(1+b_N)` then has
magnetization `2b_N/(1+b_N)` and excitation-energy expectation at most `3J/N`.
The latter follows from the exact double-commutator identity below.
Project this vector to energies `[0,J/sqrt(N)]` and normalize. The removed
norm squared is at most `3/sqrt(N)`, and the vector-distance error tends to
zero. Its excitation-energy **operator** norm is at most `J/sqrt(N)`;
therefore every unit-norm local A satisfies the ground inequality with defect
at most `J/sqrt(N)`. Taking a weak-* limit constructs an ordered ground.
The biased state and the projected state are invariant under even
translations and under one-step translation combined with a pi internal
spin rotation reversing axis 1; therefore the two sublattice means are
opposite, and their averaged magnetization is their local staggered
magnetization. Thus the maximum magnetization obeys

\[
 m_s\ge {2b\over1+b}>{1\over1000}=:m_0.
\]

We use the conservative exact lower bound `m_0`; no numerical value of the
actual magnetization is inferred. The pure-state conclusion follows from
Koma's stated theorem; the quantitative packet bounds below apply directly
to his field-selected state, which need not itself have been shown pure.

## Currents, moments and a strictly positive spectral window

For any internal axis a the directed bond current is
`j_xy^a=J(S_x cross S_y)^a`. Exact continuity and the first-moment identity are

\[
 {dS_x^a\over dt}=-\sum_{y\sim x}j_{xy}^a,\qquad
 [A_f^*,[H_J,A_f]]=-J\sum_{\langle xy\rangle}
 |f_x-f_y|^2\sum_{b\ne a}S_x^bS_y^b.
\]

The transverse bond has norm 1/2, so the expectation of the double
commutator is at most `(J/2) sum_edges |f_x-f_y|^2`. For a unit-normalized
plane wave this is `J E(k)`. For the axis-1 staggered field, transverse
first moments acquire at most `B/2`, which vanishes in the ordered limit.
All signs, all three axes, arbitrary complex two-site coefficients, and
current norm `J/2` are checked exactly in `algebra.py`.

Write `mu_k` for the axis-2 positive-energy structure measure and
`nu_(k+Q)` for axis 3. The Ward identity is
`m_s=-i omega([S_-k^2,S_(k+Q)^3])`. The reflection-positive susceptibility
bound and the first moments give

\[
 \int_{(0,\infty)}{d\mu_k(E)\over E}\le{1\over2J E(k+Q)},
 \quad \int E\,d\mu_k(E)\le J E(k),
 \quad \nu_{k+Q}((0,\infty))\le\sqrt{{E(k+Q)\over2E(k)}}.
\]

These statements are used first on the finite positive-field ground, where
the transverse means vanish and the ground is unique, then on smooth
packets in the cited order of limits. They do not posit a pointwise
infinite-volume plane-wave vector. Reflection gives equal bounds for both
terms in the Ward commutator.

Put `x=E(k),y=E(k+Q)`. Splitting each Ward term at positive energies a,b,
Cauchy–Schwarz bounds its low-energy part by `a/(2J sqrt(xy))`.
Its remainder is bounded by
`[mu_k([a,b])+Jx/b]^(1/2) [y/(2x)]^(1/4)`.
Consequently

\[
 \mu_k([a,b])\ge
 \left({m_s\over2}-{a\over2J\sqrt{xy}}\right)^2
 \sqrt{{2x\over y}}-{Jx\over b},
\]

whenever the parenthesis is positive. This estimate retains both Ward terms,
the high-energy tail, and the susceptibility contribution at low energy.

Take a smooth packet with momentum support `p/2<=|k|<=p<=1` and
`integral |fhat|^2 dk/(2pi)^3=1`. Then
`p^2/20<=x<=p^2/2`, `5<=y<=6`, and `p/2<=sqrt(xy)<=2p`.
For any known `0<m<=m_s`, choose

\[
 a=mJp/8,\qquad b=32Jp/m^2.
\]

The last inequality gives at least `(m^2/8) sqrt(x/y)>=m^2p/96`.
For example, `sqrt(2)>7/5` leaves a strictly positive rational margin in
this comparison. Let `g_hat` be smooth, equal to one on `[a,b]`, bounded
between zero and one, and supported in `(a/2,2b)`. If
`A_f=sum_x f_x S_x^2` and `A_(f,g)=integral g(t) alpha_t(A_f)dt`, then

\[
 \omega_0(A_{f,g}^*A_{f,g})\ge {m^2p\over96},\qquad
 {mJp\over16}\le
 {\omega_0(A_{f,g}^*[H_J,A_{f,g}])\over
   \omega_0(A_{f,g}^*A_{f,g})}
 \le {48Jp\over m^2}.
\]

The numerator is bounded by `Jp^2/2`; the filter also confines the entire
excitation spectrum to `(mJp/16,64Jp/m^2)`. This vector is nonzero and
orthogonal to the **entire** zero-energy subspace. Letting p decrease proves
that the field-selected ground's positive spectrum accumulates at zero.
Using `m=m_0=1/1000` removes an unverified order parameter from this result.
The constants are conservative bounds, not a measured propagation speed.
A similar mode near Q exists by the cited theorem; this package's explicit
constants and direct witness concern the packets near zero only.

## Transfer to the actual source and recorded repairs

For E bonds keep the FA90 trajectory `c=J E^4,U=J E^8` and its complete
C/E field-extension instrument. Initial spin density matrices are arbitrary
in FA90's finite-history comparison, including the selected finite states
above and the finite positive-field ground states. We inject them through
P (one fermion per site). No magnetic field is added to the source dynamics.
The source is not redefined, and its original bare one-bit repair is not
identified with the declared two-flavour extension.

For a finite positive-field ground, a unit-norm observable on s sites has
zero-field local ground defect at most `B s`; under the dressed low-band
Hamiltonian the extra defect is at most `2 eta=6J/E`. Hence first E tends
to infinity and then B to zero, or an appropriate diagonal sequence, gives
the selected zero-field ground on the same limiting spin algebra. The
undressed/dressed expectation error is at most `4/E^3`.
FA90's all-record trace error is
`8(1+Gamma E T)/E^3+6JT/E` for the specified independent Poisson edge clocks.
It tends to zero for every fixed local experiment and finite time interval.

Smooth packets are quasi-local because f decays rapidly; their energy filters
are quasi-local by the already proved Lieb–Robinson dynamics. First truncate
f and g(t), perform the source comparison at those fixed truncations, then
remove the truncations in norm. This constructs the same filtered correlation
and nonzero spectral weight as a limit of actual source histories. It does
not assert a uniform finite-box resolution of all arbitrarily small energies.
The selected initial states have been constructed mathematically; their
large-volume native pulse preparation is not claimed executed.
In particular this proves existence of an accessible ordered ground sector,
not order of every exact finite-source ground subsequence.

## Verification and limits

`infrared.py` is a rigorous integral upper bound, `algebra.py` checks the exact
current/rotation/Ward identities, and `bounds.py` checks all rational margins
and source-selection defects. Finite diagonalization controls independently
check moments, susceptibilities, and biased-ground energy. They are tests of
the formulas and normalization; finite matrices do not prove the infinite
order theorem. The receipt separates these new checks from inherited source
custody. Old experiments and receipts are not silently recomputed or repaired.

This closes nonzero internal Néel order and a genuine gapless spectral sector
for this fixed-spacing interacting branch. It supplies neither a sharp
one-particle pole nor Lorentz invariance, scattering completeness, a physical
photon, a measured particle species, a relativistic rest mass, or Einstein
dynamics. The earlier bounded-U continuum branch remains distinct.

From the repository root:

```sh
python -W error -m pytest model/spin_goldstone_sector -q
python -m model.spin_goldstone_sector --write
python -m model.spin_goldstone_sector --check
```
