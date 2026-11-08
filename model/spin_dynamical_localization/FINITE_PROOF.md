# A complete time-invariant source-cube coordinate detector

Use the unchanged finite source of `native_ordered_ground`: the open eight
site spin cube with `J=1`, selection field `B=1/4`, and Hamiltonian `H_B`.
Its unique normalized ground is `Omega_B`, and `A=Z_0` is the unitary local
probe. Every time below uses this same finite `H_B`. This example does not
replace the zero-field continuum clock by a finite-field gap.

## The actual source has an invariant response enclosure

In the Marshall basis let `C` complement all eight configuration bits.
It is a selfadjoint involution, and exact matrix entries show
`[C,H_B]=0`. The irreducible source has a unique strictly positive ground;
`C` permutes its entries, so `C Omega_B=Omega_B`. Since `C A=-A C`, the
whole orbit `exp(-itH_B)A Omega_B` lies in the odd subspace

    P_odd=(I-C)/2,       rank P_odd=128.                     (1)

This subspace reduces the actual Hamiltonian. Its even complement is an
exact time-invariant dark outcome.

The source group `G=(Z/2Z)^3` is the combined graph-reflection group of
`native_joint_spin_spectroscopy`. In the Marshall frame its actions `U_g`
are configuration permutations. They commute with `C` and `H_B`. Let
`Q_k=(1/8)sum_g chi_k(g)U_g`. For `g!=0`, both `U_g` and `C U_g` have
sixteen fixed configurations. The identity has trace256 and `C` has trace0.
Consequently

    rank Q_k P_odd = (1/16) sum_g chi_k(g)
                         [Tr U_g-Tr C U_g] =16             (2)

for every character `k`. The certificate also constructs sixteen exact
orthogonal columns in every block. Their total dimension is128, so the
blocks exhaust (1). These are finite reflection characters, not small
continuous translation momenta.

## Preserve the source-defined first position mode

Set `a=A Omega_B`, `w_k=||Q_k a||^2`, and
`psi_k=Q_k a/sqrt(w_k)`. The full source witness certifies `w_k>1/100`.
Choose a real orthonormal frame `B_k` of the rank16 block in (2), with first
axis having positive overlap with `psi_k`. The finite algorithm picks an
axis using the certified ground approximation and checks that this sign
also holds for the actual ground. In this frame write `alpha_k=B_k*psi_k`.
It is a real unit vector with `alpha_k,0>0`. The orthogonal matrix

    F_k=[I-2 v_k v_k*/||v_k||^2] diag(-1,1,...,1),
    v_k=e_0+alpha_k,       ||v_k||^2=2(1+alpha_k,0)            (3)

satisfies `F_k e_0=alpha_k`. Thus `b_k,j=B_k F_k e_j` is a complete
character basis whose first vector is exactly the source mode `psi_k`.
Fourier transform the character index, retaining all multiplicities:

    eta_x,j=1/sqrt(8) sum_k chi_k(x) b_k,j,
    E_x=sum_(j=0)^15 |eta_x,j><eta_x,j|,
    E_dark=I-P_odd.                                        (4)

Orthogonality proves that (4) is a complete projective measurement, with
eight rank16 position effects and a rank128 dark effect. The actual group
law gives `U_g E_x U_g*=E_(x+g)`. Its restriction to the rank8 first-mode
code is exactly the source coordinate measurement of
`spin_covariant_localization`: the `j=0` vectors and probabilities agree.
The extra axes retain time-evolved response that leaves that first-mode
code. They are not additional particle species.

The source certifies strictly positive energy variance in character111.
Since `H_B` preserves characters, this proves that the minimal invariant
span of the eight first modes has dimension at least9. Equation (1) gives
an upper bound128. The certificate does not claim the exact minimal rank
or that (1) contains no other source response. It gives a complete
invariant enclosure with no lost probability. The choice of extra axes
is fixed by the stated source algorithm; uniqueness or covariance under
additional cubic rotations is not inferred.

A complete destructive instrument maps each `eta_x,j` to the classical
position flag `x` and retained multiplicity output `j`, with the even
source vector retained on the dark branch. Its sum of branch squares is
the identity. This is an explicit finite instrument definition. No native
pulse tape implementing it is reported executed.

## Certified source and frame error

Let `e` be the existing certificate's norm error for its normalized proposed
ground. Project that vector onto even parity and normalize. Projection and
the normalization inequality give ground error `e_s<=2e`. Apply the exact
unitary `A` and exact `Q_k` to obtain the proposed odd response modes. For
their rational weights `wtilde_k`,

    |w_k-wtilde_k|<=2e_s,
    ||psi_k-psitilde_k||<=2e_s/sqrt(w_min)<=2e_s/w_min=:d,
    w_min=min_k(wtilde_k-2e_s)>1/100.                        (5)

The proposed first-axis squared overlap is at least1/16, because it is the
largest of sixteen squared coefficients summing to one. The certificate
also checks the lower bound after subtracting `d`, so the actual and
proposed first overlaps are positive. Both vectors `e_0+alpha_k` in (3)
therefore have norm at least `sqrt(2)`. Normalizing those vectors costs at
most `sqrt(2)d`; the difference of their rank-one projectors is at most
twice that. Thus the corresponding Householder matrices differ in norm by
at most `4sqrt(2)d<6d`.

Taking the orthogonal direct sum over characters and the fixed finite
Fourier transform preserves this norm bound. The actual and proposed
complete decoder dilations differ in norm by at most `6d`; their half
diamond distance is at most `6d`, including arbitrary untouched references.
Both use the same even dark branch. This is a bound on complete
instruments, not on a renormalized selected position outcome.

The exact source norm bound is `||H_B||<=14`. On the normalized proposed
input the order28 Taylor polynomial for `exp(-itH_B)` has vector error

    epsilon_T <= 3^ceil(14|t|) (14|t|)^29/29!.              (6)

Here `exp(b)<=3^ceil(b)` for `b>=0`. Every polynomial coefficient and
source-matrix multiplication is rational. Directed rational square-root
intervals evaluate the normalized source frames and decoder coefficients.
All interval operations round outward to100 binary digits. A floating
point search proposes a matrix entry for the locality control; its
admission uses only these rational intervals.

For each complete position effect, replacing the normalized evolved input
by its Taylor polynomial costs at most `epsilon_T(2+epsilon_T)` in the
probability. Add `e_s+6d` for the state and complete instrument errors.
This gives certified source probabilities at times0 and1/8, with exact
zero dark probability at both times and all intermediate times. For
example the position0 intervals are

    t=0:   [0.86746856452, 0.86746856490],
    t=1/8: [0.86273312785, 0.86273312823].                   (7)

The complete eight intervals enclose total probability one. These are
coordinate readout probabilities of the actual finite source, not mass
ratios or a relativistic scattering prediction.

## The explicit finite detector is not a one-site readout

The exact position0 effect has source-basis matrix entry

    <1|E_0|254> in [-0.39057283702,-0.39057283627].           (8)

The two configurations differ at every site. An operator supported on a
proper subset of the eight sites has zero matrix element between them,
because it acts as the identity on at least one differing site. Equation
(8) therefore proves that this finite detector cannot be supported on a
proper subset of the cube. It is also incompatible with a diagonal
literal occupation readout. A growing-volume locality theorem does not
follow from this one-cube calculation. Operational approximation and
locality in the continuum retain their separate requirements.
