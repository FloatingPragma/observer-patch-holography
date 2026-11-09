# Koide balance and the same flavor action

The Koide theorem supplies an exact relation between three positive masses
under a specified square-root-mass attachment. The finite connected-register
construction derives its balance within a declared response model. The
charged-lepton coefficients of the native flavor action do not satisfy that
relation. These statements concern different maps, which can be tested
without an empirical mass input.

## The exact identity and the finite balance theorem

For a regular three-cycle `R`, let

\[
 C=aI+bR+\overline b R^2,\qquad
 \lambda_k=a+2|b|\cos(\delta+2\pi k/3).
\]

The formal result in
[`KoideCirculant.lean`](../../Lean/ObserverPatchHolography/KoideCirculant.lean),
`koide_quotient_phase_independent`, proves

\[
 \sum_k\lambda_k=3a,\qquad
 \sum_k\lambda_k^2=3a^2+6|b|^2,\qquad
 Q={\sum_k\lambda_k^2\over(\sum_k\lambda_k)^2}
   ={1\over3}+{2\over3}{|b|^2\over a^2}.              \tag{1}
\]

Thus `Q=2/3` is equivalent to `|b|/a=1/sqrt(2)` for `a>0`.
This is a physical square-root-mass quotient when `lambda_k>=0` and
`m_k=s lambda_k^2`, with one common `s>0`. For the balanced carrier the
positive chamber is `|delta|<=pi/12` modulo `2pi/3`. Outside it, replacing
signed roots by their absolute values changes the denominator in (1).
The phase cancels in (1); the theorem does not select it.

The executable
[`derive_charged_koide_orientation_isometry.py`](../../code/particles/leptons/derive_charged_koide_orientation_isometry.py)
contains a stronger finite normalization result. Its source is

\[
 \mathcal A=B((1\oplus\chi\oplus\bar\chi)\otimes\mathbb C^2),
 \quad \rho=I_6/6,\quad
 Z_0=P_0\otimes I_2,\quad Z_c=P_c\otimes q_+ .       \tag{2}
\]

The orthogonal events `Z0,Zc` both have rank two. Conditioning on their
sum gives `p0=pc=1/2`; the excluded event has probability `1/3` in the
unconditioned six-state trace. In the normalized tracial GNS space the
three operators `I,R,R^2` are orthonormal. The specified square-root lift
therefore gives `a=sqrt(p0)`, `|b|=sqrt(pc/2)`, and hence the balanced
modulus in (1). This is a finite algebraic derivation for the declared
register, event and lift. Its hypotheses do not identify the output with
an arbitrary Yukawa coefficient.

The physical version of this construction requires a state-preserving
recoverable channel that intertwines the singlet and charged records,
acts isometrically on their GNS norms, and identifies the image with the
positive square-root-mass response after kinetic whitening. If the
canonical action coefficient is `Ye`, that response is, up to common
normalization,

\[
 C_{\rm mass}=(Y_e^\dagger Y_e)^{1/4}.                \tag{3}
\]

The fourth root in (3) is essential: the Dirac masses are proportional
to the singular values of `Ye`. They are not the squared singular values.
No running-to-pole identification follows from the finite GNS construction.

## A direct test on the declared native flavor action

The source in
[`native_flavor_orbit`](../native_flavor_orbit/PROOF.md) supplies the full
current-gated memory law, canonical kinetic form and four duty policies.
Its angle is `theta=pi/5`, reset port is `2`, and all record outcomes are
retained. For each policy its charged-lepton coefficient `Z_D` is positive
semidefinite. The interval source in
[`spectrum.py`](../native_flavor_orbit/spectrum.py) certifies all three
singular values, including the exact zero for rank-two duties. With
`m_i=singular_value_i(Z_D)`, the direct test is

\[
 Q_{\rm action}={\sum_i m_i\over(\sum_i\sqrt{m_i})^2}. \tag{4}
\]

The pure projection `source_mass_readout.summarize` in
[`source_mass_readout.py`](source_mass_readout.py) applies rational square-root
bounds to those source intervals. Its `certificate()` authenticates the
owning flavor receipt; `cold_replay=True` repeats its scientific producer.
The resulting coarse bounds below each strictly
exclude `2/3`.

| Complete duty policy | Certified enclosure for `Q_action` |
|---|---:|
| always | `[7675,7677]/10000` |
| C | `[5797,5799]/10000` |
| E | `[9502,9503]/10000` |
| CE | `[7754,7755]/10000` |

These are action-coefficient diagnostics, with no identification of their
entries as measured lepton masses. The exact current-preserving family
embeddings act by independent left/right unitaries on each Yukawa channel,
so they preserve its singular values. A common Higgs scale or clock
conversion multiplies all three masses equally and cancels from (4).
Neither operation converts these intervals into `2/3`.

There is a quantitative formulation of the failed balance premise. For
any positive `C`, put

\[
 E_0={(\operatorname{Tr} C)^2\over3},\qquad
 E_c=\operatorname{Tr}(C^2)-E_0.
\]

Then `Ec/E0=3Q-1`, so `Ec=E0` is equivalent to `Q=2/3`.
The coarse intervals above imply respectively

\[
 |E_c/E_0-1|>3/10,\quad >26/100,\quad >85/100,
 \quad >326/1000.                                   \tag{5}
\]

Consequently no attachment satisfying all the GNS record-isometry and
mass-identification conditions above can map (2) to these same four
coefficients. This conclusion does not exclude a different source
interaction or a different physical mass response. Such a construction
must exhibit its action and its map; replacing `Z_D` by a chosen balanced
circulant changes the specified action.

## Phase information and the conditional tau calculation

The actual geometric connection in
[`native_mass_holonomy.py`](../common_source/native_mass_holonomy.py)
uses spin-coherent projectors on the twelve native port directions. Every
oriented triangular face has exact spin-half Bargmann phase `pi/10`.
An integral face-boundary calculation exhausts all closed edge chains,
whose phases are integer multiples of `pi/10` modulo `2pi`.
Under the separate balanced circulant and `Y=C^2` attachment, equal-link
phase is one third of the loop phase. Positivity leaves three unordered
classes `delta=0,pi/30,pi/15`. Their source-certified mass ratios are,
respectively, approximately

\[
 (1,67.9411),\quad(6.33689,203.37238),\quad
 (109.88531,1969.28464).                              \tag{6}
\]

The specified connection therefore does not supply the proposed angle
`2/9`. A raw-event probability equal to `2/9` is not an equality to a
transport angle. Complete classical record channels are invariant under
individual Kraus-branch phases, as proved in
[`particle_source_selection`](../particle_source_selection/PROOF.md).
An accessible coherent reference can distinguish such phases; its
preparation and coupling belong to the source specification.

A phase law is unnecessary for predicting the third mass *conditional on*
physical Koide and two mass inputs. With `x=sqrt(me)`, `y=sqrt(mmu)`, (1)
gives

\[
 \sqrt{m_\tau}=2(x+y)+\sqrt{3(x^2+4xy+y^2)}.          \tag{7}
\]

The heavy-mass ordering selects the plus branch. The executable
[`koide_balance_comparison_certificate.py`](../../code/particles/leptons/koide_balance_comparison_certificate.py)
propagates its declared electron/muon input boxes through this formula
and checks that changing the tau reference does not change the output.
Its empirical ancestry is explicit: the balance is the published Koide
relation, and the comparison is conditional. Equation (7) does not derive
the two supplied masses or prove that (3) satisfies Koide.

There are two distinct sufficient routes. A source-to-mass attachment
proving physical Koide permits the conditional test (7) without a separate
phase theorem. A source prediction of both independent lepton ratios
requires an additional ratio or phase determination on the same physical
action. The finite GNS balance, the phase-independent Lean identity and the
native action diagnostics establish the precise premises for these routes;
they do not assign a physical mass to any of the fifteen internal states.
