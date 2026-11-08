# A complete coordinate detector on one cyclic source response

Use the unchanged open eight-site spin source

    H_B=J sum_<xy>(Sx.Sy-1/4)-B sum_x epsilon_x Sx^1,
    J=1, B=1/4,

its true unique ground Omega and the local encoded-Z probe A at vertex zero.
The source dimension is 256 and ||A||=1. Its combined graph-reflection group
G=(Z/2Z)^3 is the exact group in `native_joint_spin_spectroscopy`: a coordinate
reflection combines with the global spin flip when it exchanges sublattices.
In the Marshall basis these actions are literal configuration permutations.
They commute with H_B and fix its positive ground. Their characters are
chi_k(x)=(-1)^(k.x), with exact projectors

    Q_k=1/8 sum_(t in G) chi_k(t) U_t.

These labels are graph reflections on this open cube. They are not small
translation momenta. No sphere-coordinate readout or other field carrier is
identified with them by a shared dimension count.

## The source structure factor fixes the coordinate frame

Put a=A Omega and a_x=U_x a. The character weights and normalized modes are

    w_k=||Q_k a||^2,
    psi_k=Q_k a/sqrt(w_k),          when w_k>0.

The supported psi_k are orthonormal, directly by Q_k Q_l=delta_kl Q_k.
Their phases are fixed by the same probe vector a, without an independent
phase choice for each spectral eigenvector. The orbit Gram matrix is

    C_xy=<a_x,a_y>=sum_k chi_k(x)chi_k(y) w_k.

Its character eigenvalues are 8w_k. Thus the original orbit vectors need
not be orthogonal. Replacing their Gram matrix by the identity would change
the source responses. On the supported characters define

    eta_x=1/sqrt(8) sum_(k:w_k>0) chi_k(x) psi_k,
    E_x=|eta_x><eta_x|,
    P=sum_(k:w_k>0)|psi_k><psi_k|,     E_dark=I-P.           (1)

Character orthogonality proves sum_x E_x=P. These positive effects and the
dark effect form a complete measurement on the entire source Hilbert space,
with every dark and zero-probability outcome retained. The transformation
law is U_t eta_x=eta_(x+t), so U_t E_x U_t*=E_(x+t), while E_dark is invariant.

If all eight w_k are positive, the eta_x are orthonormal and the effects
in (1) are a projective measurement, with rank-eight code and rank-248 dark
complement. If some w_k vanish, the vectors are a tight frame: their Gram
is 1/8 sum_(k:w_k>0) chi_k(x)chi_k(y), rather than delta_xy. In that case the
same effects are a positive-operator measurement, not an eight-vector
projective basis. Dividing by zero or removing the dark complement is not
an alternative construction.

The map from a_x to eta_x is the inverse square root of the orbit Gram
on its support, equivalently the character multiplier 1/sqrt(8w_k).
It is a spectral whitening operation. Nothing here makes it a local
source operation. Equal weights would give an orthogonal untranslated
probe orbit; the actual source has unequal weights, with largest more
than ten times smallest. The resulting coordinate detector is consequently
not identified with merely measuring the literal probe's original site.

A destructive complete instrument has branches <eta_x| and I-P on orthogonal
output flags; it retains the untouched reference on every branch. This
specifies a finite source detector and its back-action. An arbitrary-reference
native circuit implementing all these branches is not reported executed.
The source's finite coherent-control framework supplies a possible compilation
route, whose complete pulse, precision and locality costs require their own
admission. The result here is the exact effect construction and its certified
finite source probabilities.

## Certification from the complete native source witness

The existing 256-vector integer witness is admitted by its exact source
matrix residual, Gram defect and isolated-ground certificate. Let e bound
the norm error of its normalized ground vector. Applying the exact local A
and exact Q_k preserves that error. Integer projection vectors give exact
proposed weights wtilde_k, and

    |w_k-wtilde_k|<=2e.

For every k the resulting strict lower bound exceeds 1/100, so the full
eight-character projective construction applies to the actual ground.
Write wmin=min_k(wtilde_k-2e). Normalization of nonzero vectors gives

    ||psi_k-psitilde_k||<=2e/sqrt(wmin)<=2e/wmin=:d.        (2)

The two eight-column isometries differ in operator norm by at most
sqrt(8)d<3d. Fourier transformation does not change this bound. Their
code projections differ by at most 6d. For the destructive instrument,
the detected-row block differs by at most 3d and the dark block by at most
6d. The complete dilation difference is at most sqrt(45)d<9d. Hence the
half-diamond distance of these complete instruments is at most 9d,
including arbitrary entangled source/reference inputs.

On the same ideal probe vector, the exact coordinate probabilities are

    Pr(x)=1/8 |sum_k chi_k(x)sqrt(w_k)|^2,
    Pr(dark)=0.                                         (3)

Exact square-root intervals evaluate (3) for the proposed weights, and the
complete instrument bound plus the input trace-distance e bounds each
actual probability by an additional 9d+e. No outcome is renormalized.
The proposal detector on the actual probe has dark probability at most e^2,
because the proposed probe belongs exactly to its code. These two statements
distinguish the ideal ground-defined detector from its certified finite
representation.

## The dynamics test is separate and fails for this code

Covariance and a complete coordinate measurement do not make a particle
sector. Since H_B commutes with Q_k, the component of H_B psi_k in the code
is only its scalar projection onto psi_k. Therefore

    ||(I-P)H_B psi_k||^2
      =<psi_k,H_B^2 psi_k>-<psi_k,H_B psi_k>^2.            (4)

The right side is the energy variance within one spatial character. The
full source matrix and integer projected vectors evaluate it exactly for
the proposed modes. The source bound ||H_B||<=14 and (2) give variance
error at most 6*14^2*d. The resulting certified lower bound in character111
exceeds 1/100. Thus the cyclic rank-eight probe code is not an invariant
subspace of this source Hamiltonian. Time evolution leaves it, with the
complement included in the detector. An energy origin shift leaves (4) unchanged.

This is a direct source test of a proposed finite position sector. It
prevents interpreting the Fourier construction alone as an autonomous
one-particle space. An invariant branch projection with nonzero readout
weight, or an explicitly controlled approximate invariant limit, is an
additional requirement. The spectral-tube readout supplies a separate test
for that requirement; no invariant branch follows from rank alone.

## Finite uncertainty and the continuous identification boundary

Within the full code define the coordinate-bit effect

    X_j=sum_x x_j |eta_x><eta_x|

and the discrete character observable

    R_j=sum_k k_j |psi_k><psi_k|=(I-U_(e_j))/2.

In coordinate representation these are (I-Z_j)/2 and (I-XPauli_j)/2.
Consequently [X_j,R_j]=iY_j/2, its trace is zero, and its square is -I/4.
For every normalized state with finite variances the Robertson inequality
gives

    Delta X_j Delta R_j >= |<Y_j>|/4.                    (5)

There is no constant positive lower bound in (5): a coordinate eigenstate
has zero coordinate variance and zero expected commutator. A Y_j eigenstate
saturates the positive value 1/4. Extending X_j by zero on the dark complement
and retaining the full group character observable leaves this conclusion
intact. The finite trace identity Tr[X,P]=0 also excludes [X,P]=i hbar I
for nonzero hbar on any nonzero finite-dimensional full space.

The labels R_j are discrete characters of reflections. Calling pi R_j a
momentum merely chooses a logarithm of a finite involution. It supplies no
continuous translation generator, physical length unit or particle position.
On an even periodic source the corresponding character-frame construction
is covariant under its genuine discrete translations. The canonical spin
continuum separately has strongly continuous envelope translations and
their spectral generators. Connecting a local, operational position
measurement to those same translations, retaining a stable particle
subspace and controlling the limit is the remaining identification.
The [continuum construction](CONTINUUM_POSITION.md) proves that the same
field-generated cyclic subspace has full absolutely continuous spatial
spectrum, and gives a covariant position observable by its positive Gram
whitening. Operational equality to the native carrier position record is
a separate requirement. The field quadratures Q and P are not substituted
for a particle coordinate/momentum pair.
