# A positive native-current hadron mass construction

## Source, calibration and complement

The carrier has the fifteen native internal coordinates, two spins and two
graph cells, hence sixty fermionic modes. Its colourful coordinates are the
actual Q doublet and the two native conjugate colour triplets. On those modes
use positive bare number masses `m_u=mbar-delta/2`, `m_d=mbar+delta/2`,
with the same coefficients on their conjugates. Uncoloured charged modes
have number energy one and neutral modes have zero energy. Their complete
Fock factors are retained. The coloured source is

`H=sum m_f N_f+t sum_edges,a(c_xa-c_ya)^*(c_xa-c_ya)`
`  +U sum_x C_x^2+Js sum_x S_x^2+Jw sum_x I_x^2`
`  +A sum_x B_x^2+eta sum_edges(I_x+I_y)^2`.

Here C is the actual native colour current with its orthonormal Cartan metric,
S is colourful spin, I is the native weak current on Q, and
`B_x=(N_Q,x-N_(u^c,d^c),x)/3`. The latter is a one-body occupation current
defined by the native representation. No baryon or meson spectral projector
occurs in H. Every displayed term is positive. Thus the full carrier has
`H>=min(m_u,m_d)N_colourful+N_charged_lepton`. The massless neutral factors
prevent an assertion of a unique vacuum or positive global gap.

The constants `U=2048`, `Js=Jw=1/4`, `t=1/8`, `eta=1` are declared architecture.
Three empirical central values, proton, neutron and charged-pion masses, fix
`mbar`, `delta`, and A in the retained electron energy unit. The input file
gives the primary sources and exposes every fitted value. These three mass
agreements are calibration results. Neutral-pion, vector-meson, Delta and
nuclear comparisons are excluded from this fit. Their interpretation as
physical particle masses requires a separate relativistic species dictionary.

## Exact native reduction

Twice each native colour root generator has integer real or imaginary
entries. With real root currents X, imaginary currents iY and native Cartan
currents D1 and D8,

`48 C^2=12 sum(X^2-Y^2)+12D1^2+4D8^2`.

The code constructs these matrices from the native exterior representation,
including the actual conjugate-colour frames. The spin and weak Casimirs use
`4J^2=Z^2+2Z+4Jplus^*Jplus`. Kinetic, number, and nine times B-density-square
matrices are integer occupation operators. Their exact colour-singlet
intertwiners are the epsilon and quark-conjugate contractions. The full
three-quark occupation sector has dimension 2024 and complete colour-singlet
dimension 120; the full quark-conjugate sector has dimension 576 and complete
colour-singlet dimension 64. Every complementary occupation and colour
sector remains in the source. Charge-block restrictions exhaust these colour
spaces; none selects a fitted energy by installing a projector in H.

The complete algebraic matrices are enclosed by exact dyadic Gram/residual
products and outward entry bounds. Polar-unitary comparison gives enclosures
for every ordered eigenvalue in every charge block. Floating spin-resolved
vectors and occupation probabilities are separately labelled diagnostics.
The near-degenerate even and odd baryon states are retained; the computation
does not assert a numerically isolated four-state band.

## Unique positive calibration

Let E_B(A) be the charge-one, spin-one-half lowest interaction energy with
`mbar=delta=0`. Let E_pi(A) be the charge-one, spin-zero analogue. Their exact
B-density-square operators lie respectively in `[5/9,1]` and `[0,2/9]`.
For every positive increment h the min-max principle therefore gives

`5h/9 <= E_B(A+h)-E_B(A) <= h`,
`0 <= E_pi(A+h)-E_pi(A) <= 2h/9`.

Consequently `E_B(A)-3E_pi(A)/2` is strictly increasing with slope bounded
below by 2/9. Equating this to the observed mean-nucleon energy minus three
halves of the charged-pion energy has at most one root. Explicit endpoints
at A zero and 4096 bracket that root. At zero use the baryon co-location
trial and pion positivity; at4096 use the baryon lower bound5A/9 and the
A-independent pion co-location trial. The receipt retains these rational
endpoint inequalities and a positive lower bound for mbar. Then

`mbar=(m_pi-E_pi(A))/2`, `delta=m_n-m_p`.

The producer retains nearby rational coefficients and verifies the actual
calibrated spectral residuals. The pion co-location trial bounds its
interaction energy well below the supplied pion energy; both bare masses are
strictly positive. Experimental uncertainty and physical model error are not
included in the numerical spectral enclosures.

For comparison, a free number-mass source has mean nucleon energy `3mbar`
and charged-meson energy `2mbar`. That ratio cannot fit both supplied targets.
The B-density square supplies an independent interacting baryon energy: it
vanishes for a colocated quark-conjugate pair and equals A for three colocated
quarks. This is the declared mechanism that separates those scales.

## Complete records and physical tests

On each complete admitted occupation sector, prepare its maximally mixed
state, apply a rank-one spectral filter and retain its full complement.
The success probability is the reciprocal of the full occupation dimension.
Every original occupation probability on both flags is stored, including
zeros. The orthogonal full-Fock complement has an identity failure outcome.
An idle opposite-charge reference is retained. The full nonselective filter
and occupation instrument leave the maximally mixed input unchanged, so their
total mean work is zero. These are t=0 records; rounded eigenvectors are not
asserted to be exact stationary preparations at arbitrary times.

The fitted neutron-proton difference exceeds the electron energy by
`0.00078233355931 GeV`. This opens the ground beta channel kinematically.
A charged-current matrix element, parity violation and a lifetime are not
proved by a mass difference. The complete Delta-like sectors are retained,
including their possible two-electron channels in any weak extension.

Neutral-pion, vector and Delta assignments fail their withheld mass
comparisons on this source. These are candidate-specific negative results.
The calibrated nucleon and charged-pion energies do not establish a complete
physical spectrum, continuum confinement, physical elapsed time or a common
closed observer implementation.
