# Positive spin gradients and complete nuclear records

## Declared action and input ancestry

The carrier, native colour, spin and weak currents, hopping and all Fock
complements are those of `openwave_axial_state`. Add the positive interaction

\[
 H_g=\lambda\sum_{a=1}^3(S_0^a-S_1^a)^2,
 \qquad \lambda=5/2.
\]

The proton, neutron, charged-pion and axial values, and the deuteron binding
energy, are five empirical calibration inputs. The electron energy supplies
the declared GeV unit. The rational gradient coefficient is a retrospective
construction choice motivated by the known odd-channel problem. The finite
menu 0,1,2,9/4,5/2,3,4 is retained, with five targets recalibrated for each
member. Agreement with the five inputs is not an independent prediction,
and the absence of a bound odd channel is not a blind validation.

All other source coefficients are positive. A positive global colour square
has coefficient four. It vanishes on every global singlet and is selected
by a sufficient finite-copy fragmentation bound, with coefficients one and
two retained as failed lower-bound controls. Failure of a lower bound is
not proof of a physical decay.

## Exact native operator and parameter reduction

Writing R for the signed native spin-raising current on the two cells and
Z for twice the signed spin-z current gives the integer identity

\[
4\sum_a(S_0^a-S_1^a)^2
=(R+R^\dagger)^2-(R-R^\dagger)^2+Z^2
=8\sum_x S_x^2-4S_{\rm total}^2.
\]

Each term in the first sum is the square of a Hermitian native current.
The code checks the identity on the full 2,024-state three-quark and
576-state quark-antiquark sectors before restricting to singlets. The
full source has the same positive-square definition on all occupations;
this is not a projector onto a desired composite state. Global colour,
spin and isospin, and exact cell reflection, are preserved.

Within the baryon singlet, the orbital coefficient is
`D=8U/3-4A/9`. At fixed total spin and isospin, the global spin and weak
terms are scalars. The actual native axial matrix element therefore
selects D for each declared gradient. A scalar pion equation then fixes A
and the bare mass, while the neutron-proton difference fixes delta. The
binding equation selects eta. This reduction avoids an unconstrained
multi-parameter fit and never rescales the native weak current.

The selected coefficients are rationalized at 44 binary places. All
reported certified intervals concern that exact rational action.
The complete baryon spin/isospin multiplets retain 120 states; complete
meson charge blocks retain 64. Integer projectors, positive exact Gram
factorizations and the finite spectral enclosure theorem certify every
level. The actual axial trial has integer native occupation amplitudes,
an exact residual and a separated ground interval. Its current enclosure
is obtained from the residual/gap bound, rather than floating agreement.

## Complete nuclear sectors and thresholds

Use the source-independent exact colour, spin, weak and cell-parity
projectors in `openwave_nuclear_distance`. Their ranks sum to 240 at
zero spin and weak projections, and to 2,520 after all magnetic
components are included. The gradient restricts in each block to
`2*S4-4*S*(S+1)`. Exact Gram reductions and outward spectral bounds
certify every retained level. Isospin components have the known
`-delta*Iz` shift. No unwanted spin or parity channel is removed.

The breakup comparison minimizes over complete three-quark multiplets
whose spins and isospins can couple to the specified six-quark channel.
It uses two independent copies of this identical action. In particular,
a spin-three state is compared with its compatible two-Delta channel,
rather than being declared unbound solely because it exceeds the
proton-neutron threshold. The selected source places the odd spin-zero,
isospin-zero and odd spin-one, isospin-one channels above their two-nucleon
thresholds. The even spin-one, isospin-zero binding is calibrated. A
spin-three, isospin-zero channel is also bound relative to its own
threshold and is explicitly unassigned.

Every other quark-number partition of six has a source lower bound.
For a local occupation n, native colour triality gives a Casimir floor
zero for n divisible by three and 4/3 otherwise. Minimizing

\[
 {A\over9}\bigl(x^2+(n-x)^2\bigr)
 +U\bigl(c(x)+c(n-x)\bigr)+\kappa c(n)
\]

over its two-cell occupations gives each non-three-quark fragment's
internal floor. Bare flavour shifts cancel at total isospin projection
zero. Positive hopping, spin, gradient and weak terms may be dropped in
this lower bound. All partitions are enumerated. A non-singlet
three-quark fragment costs at least the smaller of `A+3U` and
`5A/9+8U/3`, plus the global-colour adjoint cost `3*kappa` and its bare
mass. Actual singlet multiplet spectra supply the other three-quark
thresholds. These bounds establish the selected state's finite-copy
binding against every quark-fragment partition, without identifying
asymptotic scattering channels or proving spatial confinement.

## Complete preparation and readout

A normalized integer trial in the even spin-one, isospin-zero sector
defines a rank-one preparation filter on the full uniform input sector.
Its success is exactly the reciprocal of the full occupation dimension.
For each native occupation i, the success and failure probabilities are
`v_i^2/(norm*N)` and `1/N-v_i^2/(norm*N)`. The receipt lists every
nonzero numerator and explicitly assigns zero to every other success
branch, retaining all complementary failure probabilities. The full-Fock
complement has an identity outcome. An idle opposite-charge reference is
unchanged. Unital preparation and occupation readout have zero total
mean work on the stated uniform input.

Exact trial residual and spectral gap bound its distance from the
stationary ground line and therefore every complete subsequent
instrument, including an arbitrary external reference. The exact
isospin-zero quadrupole identity is `<n0*n1>/18`; its unit is graph
separation squared. No physical length, continuum orbital parity,
relativistic pole or complete observed nuclear spectrum is inferred.
