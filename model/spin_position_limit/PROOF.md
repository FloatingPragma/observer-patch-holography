# A finite-mode route from source instruments to field position

Use the ordered-spin source, selected ground, scaling clock and continuum
representation of `spin_canonical_field`. The field-generated position in
`spin_covariant_localization/CONTINUUM_POSITION.md` is defined on
`H_Q=closure{Q(F)Omega}`. Its source-vector map is

    T(Fhat/sqrt(|k|))=Q(F)Omega,
    m^2 I <= C=T*T <=3I,
    V=T C^(-1/2).

Here `V` maps scalar Fourier space unitarily onto this actual source
subspace and intertwines its existing spatial translations. The following
finite-mode construction connects this operator to complete finite
instruments. It does not identify the load centroid with excitation
position or produce a sharp particle shell.

## 1. Finite compressions converge to the source polar map

Choose a countable orthonormal basis `f_j` of scalar Fourier space from
the complex span of smooth Fourier-cut tests supported away from zero.
Such a basis follows by Gram--Schmidt from a countable dense set. The
mathematical statement allows measurable bins and tests. An executable
instance uses a computable dense test family and bins whose Fourier
matrix integrals admit certified approximations; arbitrary measurable
sets need not have computable effect matrices. The
corresponding source tests have `Fhat_j=sqrt(|k|) f_j`, so they are in the
existing test class. Let `R_n` project onto its first n vectors, and put

    T_n=T R_n,  C_n=R_n C R_n restricted to R_n,
    V_n=T_n C_n^(-1/2),  Pi_n=V_n V_n*.

Each `C_n` has spectrum in `[m^2,3]`. Extend it to the full Fourier space
by `C_n + m^2(I-R_n)`. These extensions converge strongly to C and stay
in that same positive compact spectral interval. Polynomial
approximation of `x^(-1/2)` therefore proves strong convergence of the
inverse square roots. Consequently

    V_n R_n -> V strongly,  Pi_n -> Pi_Q strongly.                 (1)

The increasing ranges are precisely `span{Q(F_j)Omega:j<=n}`. This is
convergence inside the original field representation, not the
introduction of an independent wave-function state.

Let `E_0(B)` be Fourier position's spectral measure. For a finite
measurable partition `{B_alpha}` of space, define finite matrices

    A_alpha,n=R_n E_0(B_alpha)R_n,
    E_alpha,n=V_n A_alpha,n V_n*,  E_dark,n=I-Pi_n.                (2)

The matrices are positive and sum to `R_n`. Hence (2) is a complete
positive-operator measurement on the original field space, with no
renormalization by its code probability. Equation (1) proves
`E_alpha,n -> E_Q(B_alpha)` strongly and convergence of the dark effect.
No boundary-zero assumption is needed here: the fixed set enters through
one bounded projection, not through a weak convergence of spatial labels.

Choose the complete Lüders instrument with branch operators

    L_alpha,n=V_n sqrt(A_alpha,n) V_n*,  L_dark,n=I-Pi_n.          (3)

Continuous functional calculus similarly gives strong convergence of
(3) to the position projections and `I-Pi_Q`. These operators are
selfadjoint contractions. Their complete channels converge in trace norm
on every fixed density matrix, including any fixed untouched reference:
first verify this for finite-rank trace-class operators, then use
contraction and trace-norm density. This is pointwise convergence, not
uniform diamond-norm convergence on the infinite-dimensional space.
At every finite n the dark branch remains necessary.

## 2. Bounded source probes give the finite Gram matrices

Unbounded microscopic moments cannot be passed through a bounded-history
limit without further estimates. Instead, for a real test F use

    A_p,s(F)=sin(s Q_p(F))/s,       s!=0.

For complex F take the complex linear combination of its real and
imaginary tests. Each operator is a finite linear combination of bounded
Weyl operators with known coefficients. At fixed nonzero s, the common
bounded-history limit therefore determines every matrix element

    <A_p,s(F_i)Omega_p, A_p,s(F_j)Omega_p>,
    <A_p,s(F_i)Omega_p, B_p Omega_p>,
    <A_p,s(F_i)Omega_p, U_p(t) A_p,s(F_j)Omega_p>.                (4)

Here `B_p` is any fixed bounded source history and `U_p(t)` uses the
same scaling clock. Finite-volume approximation and removal of the
selection field follow the inherited selected-state limit at each fixed
p and s. A box-size or selection-field stopping modulus is not assumed.

After this source limit, `s->0` gives `Q(F)Omega` in norm by the field's
Stone domain and the bound `|sin(sx)/s|<=|x|`. For each fixed n the
limiting Gram matrix is C_n and has lower bound m^2. The finite source
Gram matrices are therefore invertible for sufficiently accurate
approximations, and their inverse square roots converge. All matrices
in (2)--(3) are finite; their overlaps on source history vectors are
continuous functions of these Gram entries and the cross entries (4).
Thus their complete finite outcome laws, terminal overlaps on the named
history-vector dictionary, and any fixed
finite sequence of these instruments interspersed with the same time
evolution converge to those of (3). Products of finite-rank branch maps
expand into finitely many kernels (4); dark branches expand as identity
minus the code projector. All branches remain in this expansion. No trace distance between different
microscopic Hilbert spaces is asserted.

Finally let n increase. Equations (1)--(3) give the continuum position
instrument on `H_Q`, with its complement on the full history space.
The resulting order is: finite source approximation at fixed p,s,n;
the common scaling subsequence at fixed s,n; s to zero; n to infinity.
Compiler error must tend to zero before each comparison. A diagonal
choice exists by these finite lists of convergent quantities. This
existence does not provide a computable rate for the infinite ground or
scaling subsequence.

## 3. A finite probability error certificate

The construction can reject poorly conditioned finite data. Let C and
C' be n-dimensional positive Gram matrices with common lower bound
alpha>0; suppose `||C'-C||<=d_C`. The resolvent integral for the inverse
square root gives

    ||C'^(-1/2)-C^(-1/2)|| <= d_C/(2 alpha^(3/2)).                 (5)

For any declared coordinate effect `0<=A<=I`, set
`H=C^(-1/2) A C^(-1/2)` and similarly H'. Equation (5) yields

    ||H||,||H'||<=1/alpha,  ||H'-H||<=d_C/alpha^2.                (6)

For the overlaps z and z' of a subnormalized source input with the probe
vectors, assume `||z||<=M` and `||z'-z||<=d_z`. Its original outcome
weight is `z*H z`. Expanding this quadratic form gives

    |z'*H' z' - z*H z|
      <= (2M+d_z)d_z/alpha + M^2 d_C/alpha^2.                    (7)

If the coordinate matrix A is approximated by A' with
`||A'-A||<=d_A`, the additional event bound is `M^2 d_A/alpha`.
The admitted approximate matrices must themselves be positive and sum
to the identity on the finite code; entrywise rounding alone does not
certify a complete instrument. A complete instrument implementation error
adds its own total-variation bound, without division by a successful
detection probability. A common
lower bound can be certified from `lambda_min(C)>=gamma` and
`d_C<gamma`, by taking `alpha=gamma-d_C`. Scalar probabilities alone do
not certify a complete quantum instrument; the source compiler separately
bounds all input/reference outputs.

The raw measurement precision must track the sine amplitude. For a complex
source test the known linear-combination normalization can be chosen as
`Lambda=2/|s|`. If every normalized complex Gram entry has error at most
`epsilon_G`, then `d_C<=n Lambda^2 epsilon_G`. If every normalized cross
entry has error at most `epsilon_z`, then
`d_z<=sqrt(n) Lambda epsilon_z<=n Lambda epsilon_z`. These errors must
include preparation, control and sampling. Division by known coefficients
is different from renormalizing by an observed success probability. Fixed
raw accuracy is insufficient as s decreases.

## 4. Operational scope

The finite realization consumes the existing observer's declared
coherent memory, addressed source controls and pointer readout. A
compiler for arbitrary finite source-register instruments is required
for each instance of (3), with every pointer, export and failure branch
retained. Classical native repair alone is spin-blind on the protected
sector and cannot supply this operation. Knowledge of a matrix does not
by itself constitute its measurement.

The generic compiler in `spin_operational_position` supplies such
finite-dimensional implementations on its admitted coherent-service domain.
An executable instance consumes certified finite vector and effect data,
including their preparation and approximation errors. Their algebraic
existence can be combined with (1)--(7) to obtain an ordered-limit
readout. The host interpreter's resource cap, exponential circuit cost,
unknown effective ground/scaling modulus, and growing-time accuracy
remain distinct limits on running a prescribed continuum experiment.
The source state, translation action and clock are not replaced in this
comparison.

The argument establishes convergence for bounded statistics and complete
fixed histories. It supplies no uniform second-moment tail bound for
measured position or momentum. Consequently finite-detector variances
cannot be substituted into a continuum uncertainty inequality without
that additional control. The field-vector subspace H_Q is not asserted
to be invariant under time evolution: its dark complement is retained
through every experiment. A larger time-invariant response localization
is a separate construction. Neither construction proves Lorentz
covariance, a particle pole, nonzero scattering, a mass, or the
identification of a measured particle species.
