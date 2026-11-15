# A common normalization for native matter and field records

The length scale `ell0`, time scale `tau0`, action unit `hbar` and charge
unit `q0` are shared positive conversion inputs. A cell has physical side
`ell0 h`, a completed matter macro has declared duration `tau0 h`, and the
attached native clock has relative phase `h/8` per macro. Hence

\[
x_{\rm phys}=\ell_0 h x,\qquad
t=n\tau_0h=8\tau_0\varphi_{\rm clock},\qquad
p_{\rm phys}=\hbar p/\ell_0.
\]

The factor eight follows from the actual clock word
`exp(-i h Z/16)` and its two eigenvalues. It is distinct from the number of
native pulses, addressed attempts, control strength, or elapsed processor
time. A service-duration model can assign costs to the word; it does not
select an absolute physical second. Clock resolution also depends on the
complete finite-shot experiment. A systematic mean-law error is not a
single-shot clock resolution.

If the source eigenvalue is `exp(i phi)`, a declared phase lift gives
`E=-hbar phi/(tau0 h)`. Changing the lift by `2pi` changes this energy by
`2pi hbar/(tau0 h)` while preserving the source eigenvalue. The branch,
reference state and finite phase-estimation resolution must accompany an
energy claim. The conserved diagnostic `1-Re V` is not this generator.
Its small-phase behavior is quadratic, so multiplying it by a fixed unit
does not reproduce a linear acoustic energy.

## Mass, speed and calibration

The scalar mass coin supplies

\[
E_\pm=\pm\frac{\hbar}{\tau_0}
 \sqrt{m^2+|p|^2/48},
\]

on its propagating subspace. Its other24 finite-energy components and all
other native outcomes are retained by the producer. The Dirac replacement
instead supplies the free limit with `|p|^2` coefficient one. Write the
corresponding dimensionless speed as `c_d`. Both cases use

\[
c_{\rm phys}=\frac{\ell_0}{\tau_0}c_d,\qquad
E_{\rm rest}=\frac{\hbar|m|}{\tau_0},\qquad
M_{\rm rest}=\frac{\hbar\tau_0|m|}{\ell_0^2c_d^2}.
\]

Thus the coin coefficient is a dimensionless rest-frequency parameter,
not a mass in kilograms. The scalar and Dirac proposals are distinct
sources; their speeds cannot be combined by applying a private time unit
to one row.

For nonzero `m`, choosing
`tau0=hbar|m|/E_target` and
`ell0=c_target tau0/c_d` calibrates a target rest energy and a target speed.
Agreement with those targets is then an input to the unit map. All mass
ratios, dimensionless charge couplings and native outcome probabilities
are unchanged by such conversion. The source mass parameter is declared;
no electron species or interacting mass pole is identified here.

A controlled single eigenphase measures the relative rest phase at
`E_rest/hbar`. A coherent superposition of positive and negative rest
branches has population terms `cos(2 E_rest t/hbar)`. The executable
uses the source's exact rest eigenphases and resolves Fourier bin one for
the relative phase and bin two for population beating. Calling the latter
the frequency `mc^2/hbar` loses a factor of two. A global phase alone is
not an observable internal oscillation; an actual reference and instrument
must supply the comparison.

## Canonical connection and unit electric flux

Let `A0` and `B0` be the actual zero-momentum face and incidence symbols.
The exact harmonic edge space is
`ker B0 intersect ker A0^T`. Let `M` sum the two positive external
channels in each coarse direction. Its restriction to this three-dimensional
space is invertible. The unique harmonic basis `W` with `MW=I` obeys

\[
B_0W=0,\qquad A_0^TW=0,\qquad W^TW=\rho I,
\qquad \rho=2/3.
\]

The executable reconstructs this basis from the actual36edge symbols
using exact rational arithmetic. It does not select an arbitrary
orthonormal basis and label its amplitudes as unit charge. The latter
would have crossing matrix `sqrt(3/2) I` and would change the coupling.

Assume the canonical unwrapped field algebra on this source. For edge
connection `theta` and electric field `E`, define the coarse connection
and flux by

\[
a=W^T\theta,\qquad p=W^TE/\rho.
\]

On the harmonic space,
`theta=W a/rho` and `E=W p`. Therefore
`E^T dtheta=p^T da` and `[a_i,p_j]=i delta_ij`.
The coordinate `p` is actual coarse electric flux since `ME=p`.
The coordinate `a` is the line integral on a unit routed coarse path,
after removal of exact face cycles. Unit matter transport couples through
`exp(i a)`. These statements fix the charge normalization.
With the producer's incidence orientation, `BF=G-n`, so the physical
electric flux for the charge convention `q=n-G` has the opposite sign.
That common orientation change leaves the canonical form and squared
coupling unchanged; it is not an additional fitted sign.

The electric form is `E^TE=rho p^Tp`. The exact quadratic magnetic
calculation on the same face complex yields

\[
\theta^TK(k)\theta
 =a^T(|k|^2I-kk^T)a.
\]

Seven rational momentum controls and polarization verify every quadratic
coefficient. This is a conditional harmonic field reduction, with its
state, cutoff and continuum hypotheses supplied separately. Finite Z16
registers do not possess the exact canonical commutator assumed here.

## Speed and charge coupling from the same field coefficients

For `N` field substeps per matter macro, with the same homogeneous
electric coefficient `alpha` and magnetic coefficient `beta`, define

\[
A_e=N\rho\alpha,\qquad B_m=N\beta.
\]

The long-wavelength harmonic Hamiltonian has the forms
`(A_e p^2+B_m(curl a)^2)/2` in the common scaled time. Its characteristic
speed and canonical charge coupling satisfy

\[
c_d^2=A_eB_m=N^2\rho\alpha\beta,\qquad
g_{\rm charge}^2=\frac{A_e}{c_d}
 =\sqrt{\frac{\rho\alpha}{\beta}}.
\]

Indeed, rescaling `a=g_charge A_c` and `p=Pi_c/g_charge` gives the
canonical Maxwell coefficients `c_d/2` for both terms. The unit routed
phase becomes `exp(i g_charge A_c)`. The symbol `g_charge` denotes gauge
coupling and is distinct from the magnetic Landé factor.

The scalar matched setting `N=1, alpha=1/8, beta=1/4` gives
`c_d^2=1/48` and `g_charge^2=1/sqrt(3)`. The Dirac matched setting
`N=8, alpha=1/16, beta=3/8` gives
`A_e=1/3`, `B_m=3`, `c_d=1` and `g_charge^2=1/3`.
The actual homogeneous curvature bound158 gives substep stability
`alpha beta 158=237/64<4` for the latter setting. It is a harmonic
stability statement, not a nonlinear or quantum vacuum theorem.
For the scalar setting the same conservative upper estimate is `79/16`,
which exceeds four. That bound does not certify stability of every
ultraviolet scalar-field mode; it also does not prove an unstable mode
exists. The small-momentum speed and full-spectrum stability are distinct.

The executed one-face member has `alpha=.17, beta_face0=.41` and207zero
face coefficients. It is not admitted to this homogeneous normalization.
Even a distinct homogeneous member with those two coefficients would
have speed squared `697/15000`, different from the scalar value `1/48`.
The executable rejects both substitutions.

Keeping `alpha beta` fixed preserves the speed while changing
`alpha/beta` changes the coupling. An exact control replaces the Dirac
substep values with `alpha=1/32, beta=3/4`: speed one is preserved but
`g_charge^2` is halved. The speed/stability choice therefore does not
derive or calibrate an observed charge coupling.

## Conditional Coulomb coefficient

An independent exact Schur calculation of the actual incidence symbol
gives its acoustic Laplacian eigenvalue `|k|^2/8+O(|k|^4)`. A cell-total
charge has normalized constant-vector coefficient `q/sqrt(12)`, so the
real static Green kernel has leading coefficient `(2/3)/|k|^2`. This
agrees with the unit-flux electric form above. Real minimization of
`||E||^2/2` under the Gauss constraint gives this quadratic kernel.
An integer instantaneous field does not automatically realize that
minimum, and a finite quantum state does not inherit it merely by having
the same Gauss label.

If this real static response is established on the same physical source,
with controlled spatial and thermodynamic limits, the pair potential is

\[
V(r)=\frac{\hbar\ell_0}{\tau_0}
       \frac{A_e q_1q_2}{4\pi r}.
\]

Here `q_i` are integer source charges; the torus has its compensating
reference and zero mode, so the `1/r` statement is an infinite-volume or
controlled intermediate-distance limit. If `Q_i=q0 q_i` denotes a
physical charge unit, matching the conventional expression gives

\[
\epsilon_0=\frac{q_0^2\tau_0}{\hbar\ell_0 A_e},\qquad
\alpha_{\rm EM}=\frac{q_0^2}{4\pi\epsilon_0\hbar c_{\rm phys}}
 =\frac{g_{\rm charge}^2}{4\pi}.
\]

Consequently the stated Dirac coefficients imply the conditional value
`1/(12pi)`, while the scalar coefficients imply `1/(4pi sqrt(3))`.
These are architectural consequences under the stated limit premises.
They are not measured constants, a calculation of a renormalized
interacting coupling, or an executed two-charge force. Rescaling `q0`
rescales the inferred permittivity by its square and leaves
`alpha_EM` invariant. A private charge unit cannot repair a dimensionless
comparison.

The actual finite source, state/reference and clock contracts govern
which of these formulas can be used. The separate Wilson-family
theorem stack supplies a useful contract pattern, but its vacuum,
representation and field results are not imported into this source.
Canonical field-state, interacting-limit, species and force claims are
represented by explicit premises in UNITS.json.
