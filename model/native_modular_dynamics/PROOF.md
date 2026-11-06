# Native count reconstruction of modular histories and sector references

The source is the binary total-two native carrier with independent uniform
selection of its thirty seams and fair placement coins. Every no-movement
attempt is counted. Joint configuration readout is the current-only census
and signed-current dictionary in `native_thermal_return`. The constructed
objects are statistics reconstructed from its complete public records, not
additional native quantum channels or a physical clock identification.

## The continuing reference and its modular normalization

Let P be the symmetric 66-state complete-attempt transition matrix obtained
by counting all sixty seam/coin branches per input. The existing source
certificate proves `2I/3 <= P <= I`; the token graph is connected. For a
supplied positive loop length L, native return counts reconstruct

    rho_L=P^L/Tr(P^L),       H=-log P.

This reference is faithful. Its modular automorphisms are

    alpha_s^L(A)=rho_L^(is) A rho_L^(-is)
               =P^(iLs) A P^(-iLs).

Thus all these references use the same generator H, and their modular
parameters compare exactly by `t=Ls`: `alpha_s^L=alpha_(Ls/M)^M`.
Here t parametrizes `U(t)=P^(it)=exp(-itH)`. With the conventional Hamiltonian
Heisenberg action `A -> exp(i tau H) A exp(-i tau H)`, the same modular
action has `tau=-Ls`. The scalar trace normalizations cancel. The states rho_L generally differ;
clock compatibility is not equality of their expectation values. Each
rho_L is invariant under its modular flow, so ordered correlation functions
are stationary under a common shift of every t. This is a continuing regime
of the reconstructed algebra and its reference. The native transition P,
conditional return selection, modular automorphism, recorded attempt count
and physical duration are distinct objects.

## Finite native words for ordered modular correlations

Put X=I-P, so `||X||<=1/3`. For real t define c_0=1 and

    c_(k+1)=c_k (k-it)/(k+1),
    V_m(t)=sum_(k=0)^m c_k(t) X^k.

The binomial series converges in operator norm to `U(t)=P^(it)`. If
`N=ceil(|t|)>=1`, then `|c_k(t)|<=binom(N+k-1,k)`: induction uses
`|k-it|<=k+N`. Consequently the explicit rational tail

    epsilon_m(t)=(1-1/3)^(-N)
                 -sum_(k=0)^m binom(N+k-1,k) (1/3)^k

bounds `||V_m(t)-U(t)||`. At t=0 the error is zero. Rational t gives exact
rational real and imaginary polynomial coefficients. Expanding `(I-P)^k`
expresses V_m as a finite polynomial in P. No matrix logarithm or floating
eigensolver supplies a premise of this approximation certificate.

For coordinate projections Q_a and ordered times t_1,...,t_r, trace
stationarity gives

    phi_L(alpha_(t_1)(Q_a1) ... alpha_(t_r)(Q_ar))
      = (rho_L U(t_1-t_r))_(ar,a1)
        product_(j=1)^(r-1) U(t_(j+1)-t_j)_(a_j,a_(j+1)),

where `alpha_t(A)=U(t) A U(-t)`. Replacing the r relative-time unitaries by
V_m gives a finite sum of exact rational native closed-word moments. A term
with powers j_0,...,j_(r-1) has unnormalized probability

    (1/66) (P^(L+j_0))_(ar,a1)
            product_(k=1)^(r-1) (P^j_k)_(a_k,a_(k+1)).

It is the event of the stated intermediate configurations and the final
return, with all other outcomes retained. The normalizer is the independent
L-step return probability `r_L=Tr(P^L)/66`. Zero powers implement
simultaneous constraints: inconsistent repeated configurations have zero
probability. The longest term has `L+r m` body attempts. The first census,
full signed-current word and final return flag make these events readable;
on return the final configuration label also identifies the initial label.
Private census and warm-up tapes are retained but excluded from this common
public marginal. A weighted combination of these probabilities is a
reconstructed correlation, not a single acceptance probability or a native
implementation of U(t).

If the r approximation errors are epsilon_j, the correlation error is at
most `product_j(1+epsilon_j)-1`. This follows by telescoping the operator
product, using unit norms for exact unitaries and coordinate projections,
and the norm-one positive reference functional. The relative-time formula
needs r approximants rather than 2r; it also makes common time-translation
invariance exact at every truncation degree.

## Stability under record and matrix errors

Let A_j be the sum of the absolute values of the coefficients of V_m in
the P-power basis, bounded by the sum of the absolute real and imaginary
parts. Set `A=product_j A_j` and `M=product_j(1+epsilon_j)`. If every native
numerator-event probability and the normalizing return probability has
absolute error at most delta, then, for `delta<1/66`,

    |estimated polynomial correlation - polynomial correlation|
      <= delta (A+M)/(1/66-delta).

Each numerator is an unconditional event probability; zero events,
nonreturns and census failures are included. The bound accounts for signed
and complex coefficient amplification and for normalization. A separate
finite-sample guarantee is needed when estimating those probabilities from
samples. The package's systematic delta comes from the existing bounded
census/preparation theorem. No claim of practical sample complexity follows.

A separate matrix bound applies to a Hermitian contraction estimate
`0<=P_hat<=I` with entry error at most e. Then `eta=66e` bounds its operator
error. If `B_j=sum_k k |u_jk|`, polynomial telescoping gives
`||V_j(P_hat)-V_j(P)||<=eta B_j`. Also

    ||rho_hat_L-rho_L||_1 <= 2*66*L*eta/(1-66*L*eta)

when `66*L*eta<1`, because `Tr(P^L)>=1`. Combining these two explicit
bounds with the product norms controls the reconstructed correlation.
A raw empirical matrix is not asserted Hermitian, positive or contractive;
those are hypotheses of this second estimate, not consequences of entry
accuracy. Exact native-count evaluation uses rational arithmetic instead.

## Same-reference sectors and the compatibility boundary

For a nonzero coordinate projection Q, compression of rho_L gives the
faithful sector reference

    p_Q=Tr(Q rho_L),       sigma_Q=Q rho_L Q/p_Q.

Nested coordinate compressions condition associatively, with their
acceptance probabilities multiplying. These probabilities and all sector
matrix entries are native count moments. They provide compatible reference
restrictions, even when their modular flows are incompatible.

Suppose the inclusion of the Q-corner intertwined its normalized modular
flow with the full rho_L flow. Applying the intertwining to the corner
identity Q would force `rho_L^(is) Q rho_L^(-is)=Q` for every s. Equivalently
`[Q,rho_L]=0`, hence `[Q,P]=0`, since positive functional calculus recovers
P from rho_L. Every nonempty proper coordinate subset cuts an edge of the
connected token graph, so `[Q,P]!=0`. Therefore no such proper coordinate
corner intertwines the full modular flow. This covers only coordinate
compressions, not arbitrary tensor regional algebras or all observer maps.
The diagonal coordinate algebra likewise cannot be invariant under this
continuous flow: a continuous automorphism group on its finite set of
minimal projections must fix every projection, contrary to the same test.

There is an explicit record-sector completion. Choose a finite coordinate
partition `(Q_a)` and retain its classical label. Its completed reference is

    rho_bar=sum_a Q_a rho_L Q_a=sum_a p_a sigma_a.

The matrix entries are obtained from the same source moments by setting
different-label entries to zero. This corresponds to sector-filtered
reconstruction, with the label and its probability retained. Each sigma_a
is faithful on its corner. On the declared algebra `direct_sum_a Q_a M Q_a`,
the rho_bar modular action restricts exactly to the sigma_a modular action:
the p_a phases cancel. This supplies a common compatible reference for that
sector algebra. The reference has changed from rho_L, and the sector
algebra includes all reconstructed within-sector matrix statistics; it is
not merely a local occupation-bit algebra. No physical pinching apparatus
is inferred from the matrix operation alone.

The completed modular generator is generally not the restriction of H.
Compression does not commute with powers: leaving and re-entering a sector
contributes to `Q P^2 Q-(Q P Q)^2`. The exact positive two-step difference
is a retained control. Consequently the uncompleted formula `t=Ls` cannot
be transferred to independently compressed references without another
proof. Tracial references and one-dimensional corners have trivial modular
action; the certificate retains these controls as well as nontrivial
within-sector commutators. None of these conclusions identifies modular
flow with native stochastic evolution or physical time.

## Finite interaction diagnosis beyond a fixed-hopping comparison

The native effective-field package supplies the exact local generator
`J=I-P=G/60` and positive logarithm polynomials
`H_m=sum_(k=1)^m J^k/k`, with
`0<=H-H_m<=r_m I`, where `r_m=1/[2(m+1)3^m]`.
The finite reflection-positivity and source-Poissonization theorems are
inherited; they are not additional results of this package. Poissonizing
native attempts with a unit-rate total-attempt Poisson clock gives the distinct local generator J, while discrete
attempt reconstruction gives H. Changing that clock changes the operator.

Every quadratic, number-conserving one-body Hamiltonian in the fixed native
occupation dictionary has diagonal entries on the two-particle sector

    K_({i,j},{i,j})=c+epsilon_i+epsilon_j.

This includes arbitrary one-body hoppings and on-site potentials. Pairing
terms vanish when restricted between fixed two-particle basis states, so
allowing them before this compression does not remove the diagonal test.
For four distinct ports, the rectangular contrast therefore vanishes:

    K_{01,01}+K_{23,23}-K_{02,02}-K_{13,13}=0.

The native generator instead has `G_({i,j},{i,j})=10-2 A_ij`, where A is
the carrier adjacency. Its displayed contrast is -4, or -1/15 for J.
It follows that J cannot be any such quadratic operator, even after
adjusting all hoppings, diagonal potentials and a scalar energy origin.
If K is in this class, the operator norm distance obeys
`||J-K||>=1/60`, because each of the four diagonal errors is bounded by
the operator norm. All 66 native states are retained.

The same obstruction applies to the complete reconstructed H, rather than
only its first logarithm term. Compute the four diagonal entries of H_16
exactly. Each positive remainder diagonal lies in `[0,r_16]`, so the true
contrast lies within `2 r_16` of the computed contrast. The resulting
rational interval excludes zero. Its distance from zero, divided by four,
is an explicit lower bound on the distance from H to every quadratic in
the stated class. The executable imports the existing logarithm producer,
checks its generator against all native transition entries, and reports
that enclosure. This is stronger than deleting the density interaction at
fixed hopping. It makes no claim about arbitrary re-encodings of states,
auxiliary fields, chiral particles or an interacting relativistic continuum.
