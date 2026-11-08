# Transfer-ground responses from native records

The fair source determines a positive transfer generator and its unique
ground vector. Complete native record experiments determine finite
approximations to that ground's positive-energy response, with explicit
preparation and sampling costs. This is a reconstruction of correlations.
It does not execute a unitary on the native working carrier or prepare a
physical vacuum.

## The reference and the field dictionary

Use the connected fair binary source with V ports, M seams, independent
uniform seam choices and fresh fair endpoint coins. The C/E mixture is
fixed and the repair probability is one. The conserved global occupation q
is part of the readable source specification. Its population transfer is

    P_q=I-G_q/(2M),       G_q=sum_e(I-S_e).

The owning [exchange theorem](../native_exchange_clock/PROOF.md) supplies
positivity of P_q and the exact source count normalization. On the
configuration transfer Hilbert space the normalized constant vector Omega
is the unique ground of

    H_M=-2M log(P_q),      H_M Omega=0.                       (1)

Connectedness gives uniqueness, including the trivial endpoint sectors.
The native working density tau_q=I/binomial(V,q), the pure transfer state
|Omega><Omega| and a supplied closed-loop reference P_q^L/Tr(P_q^L) are
different objects. Equation (1) selects a ground of the reconstructed
operator. It asserts no preparation of that pure density on the working
carrier and chooses no Higgs ground.

Take 2<=q<=V-2. For two disjoint source seams define the diagonal field
F=(n_a-n_b)(n_c-n_d), with norm one, mean zero and variance

    v_q=4 binomial(V-4,q-2)/binomial(V,q)>0.

Its vector F Omega/sqrt(v_q) lies in primitive degree two. The owning slice
isometry identifies it with u=h/2 on the two-occupation sector, independently
of q. For several such fields the matrix

    W_ij(t)=<F_i Omega,exp(-it H_M)F_j Omega>/v_q              (2)

has a positive matrix spectral measure supported on nonnegative energies.
For arbitrary field labels i_a, times t_a and complex coefficients z_a,
the sum of conjugate(z_a)z_b W_(i_a,i_b)(t_b-t_a) is the squared norm of
sum_a z_a exp(-it_a H_M)F_(i_a)Omega/sqrt(v_q). This proves positive type
on the full finite transfer Hilbert space. It is also the normalized
ground correlation <Omega,F_i(t)F_j Omega>/v_q under the conventional
Heisenberg action. The matrix need not have orthogonal field vectors;
operator-norm estimates for a collection carry their Gram factor.

For every integer lag n, the complete four-E experiment gives

    C_ij(n)/v_q=<h_i,P_2^n h_j>/4.                           (3)

The first two E records supply outgoing contrasts and the last two supply
incoming contrasts. This reflection includes the measurement backaction.
All coin/raw flags, zero maps and exports are retained. The
[correlation proof](../native_correlation_sectors/PROOF.md) proves (3),
and the [preparation proof](../native_correlation_preparation/PROOF.md)
bounds its signed score after actual native warmup from arbitrary input,
including an untouched reference. The preparation archive is retained but
traced from this comparison; it is not used for adaptive continuation.

## Blocking controls the continuation cost

On primitive degrees at most r, the source has 0<=G<=2r Delta, where Delta
is the maximum port degree. Suppose M>=2r Delta, and choose

    b=floor[M/(2r Delta)],       Q=P^b,       u=2Mt/b.

The scalar Bernoulli inequality gives

    Q >= (1-r Delta/M)^b I >= (1-b r Delta/M)I >= I/2.        (4)

Therefore exp(-it H_M)=Q^(iu) on this subspace. Set c_0=1 and

    c_(k+1)=c_k(k-iu)/(k+1),
    V_m(Q)=sum_(k=0)^m c_k(I-Q)^k=sum_(j=0)^m a_j Q^j,
    a_j=(-1)^j sum_(k=j)^m binomial(k,j)c_k.                  (5)

For N=max(1,ceil|u|), the coefficient bound
|c_k|<=binomial(N+k-1,k) follows by induction. Thus

    ||V_m(Q)-exp(-it H_M)|| <= epsilon_m
      =2^N-sum_(k=0)^m binomial(N+k-1,k)/2^k.                (6)

At t=0 the polynomial is exactly I and the error is zero. Every coefficient
and the error bound are rational complex or rational when t is rational.
No floating logarithm or analytic continuation of fitted poles is used.
The certificate uses A_m=sum_j(|Re a_j|+|Im a_j|), which bounds the sum
of complex coefficient moduli. Since b>=M/(4r Delta), |u|<=8r Delta|t|.
Consequently a fixed-degree, fixed-time continuation has a truncation bound
and a coefficient-amplification bound independent of capacity. The required
native lag b m grows with M. Small error need not have practical sampling
cost, and finite native attempt ceilings are not removed.

Substitution of (3) into (5) gives an explicit open-record statistic

    W_hat_ij(t)=sum_(j=0)^m a_j C_ij(bj)/v_q.                (7)

It differs from (2) by at most epsilon_m. It is a signed complex combination
of separate experiments, not one acceptance probability or an implemented
quantum effect. The finite closed-loop reconstruction in
[native_modular_dynamics](../native_modular_dynamics/PROOF.md) supplies the
binomial method. Here the reference is the transfer ground, the experiments
are the source-attached open four-E records at the same q, and blocking
keeps the coefficient bounds uniform in M on bounded correlation degree.

The local exchange continuation has the additional controlled error

    |W_hat_ij(t)-<h_i,exp(-itG_2)h_j>/4|
       <=epsilon_m+|t| r^2 Delta^2/(M-r Delta),               (8)

with r=2 for these probes. This uses the owning generator comparison on
the same reducing degree subspace. It supplies a local exchange response
at the counted time t. It does not prove a spatial continuum, a Lorentz
cone, or an identification of real physical time with t.

## Preparation, failures and finite samples

Suppose every normalized lag statistic in (7) has systematic error at most
eta. Then its contribution to the reconstruction error is at most A_m eta.
The owning warmup budget proves such a bound for each four-E experiment,
uniformly in its lag and in arbitrary initial input. The reference-valued
version bounds the trace norm of the signed reference operator relative to
the same untouched reference marginal. It is not preparation of the entire
ground or working density. Combining different experiments does not create
a single joint quantum instrument.

Each complete score lies in [-1,1]. Give every one of the m+1 lag settings
L independent repetitions. Chebyshev's inequality and a union bound imply
that, with failure probability at most p, the complex sampling error in
(7) is at most s whenever

    L >= (m+1) A_m^2/(p v_q^2 s^2).                          (9)

Repeated runs must have the stated independent fresh randomness; no
independence of successive scores from one trajectory is inferred. Equation
(9) remains valid around the actual prepared expectation, so the total
error to (2) is epsilon_m+A_m eta+s. A bounded numerical example is not a
claim that the displayed number of experiments has been executed.

If a host stop is assigned score zero and its probability is at most
p_stop per setting, add A_m p_stop/v_q. This requires a proved stop bound.
The expanded interpreter's RESOURCE_STOP is retained; merely querying a
lazy kernel word does not make that bound small. The packet reports complete
native program digests only when warmup+b m+4<=2^32. Nonadmitted programs
have no digest and remain explicit failures of that sufficient budget.
There is no lower bound on the resources of every possible protocol in this
nonadmission result. Complete raw-law and array errors require their own
bounds; none is hidden in exact count arithmetic.

## An intervention condition that separates the two references

For the same contrast field and the transfer ground, define the normalized
linear phase-kick susceptibility

    chi_H(t)=i<Omega,[F(t),F]Omega>/v_q=-2 Im W(t).            (10)

It is the derivative at epsilon=0 of the future F expectation after the
state kick exp(i epsilon F), divided by v_q. Put
mu_k=<h,G_2^k h>/4. The scalar inequality sin x>=x-x^3/6 for x>=0 gives

    chi_G(t)>=2t mu_1-t^3 mu_3/3,
    chi_H(t)>=2t mu_1-t^3 mu_3/3
                    -2t r^2 Delta^2/(M-r Delta),            (11)

for t>0. The exact source moments give positive lower bounds at the
retained times. Since ||F||<=1, the second derivative of the kicked
expectation has norm at most four. For epsilon>0 its finite change is at
least

    epsilon v_q chi_H_lower(t)-2 epsilon^2.                 (12)

The retained witness makes (12) strictly positive. All quantities in its
lower bound are exact fractions; independent dense spectral calculations
check the orientation and the finite-kick result.

On the native working carrier, any unitary kick leaves tau_q invariant.
There is also a stronger statement for the bare record interface: every
fine C/E history is a weighted partial injection of configuration basis
vectors. Its effect is diagonal. For an arbitrary initial state and untouched
reference, tracing the terminal working carrier gives, in each retained
record block, a weighted sum of diagonal reference blocks rho_(x,x).
A diagonal ideal phase kick changes none of these blocks. Complete record and
reference laws are therefore unchanged by that kick, including zero
outcomes, record-adaptive selectors and retained host/FIFO stops. The
terminal quantum working operator is not asserted unchanged; a subsequent
coherent analyzer is outside this bare-record conclusion.

No exact implementation of this phase kick is inferred from that algebraic
statement. If a declared implementation has complete half-diamond error
delta from the ideal phase channel, its subsequent bare-record score can
change by at most 2 delta. Processor failures must be included in delta.
Such an allowance can be compared with the finite margin in (12); the
reconstruction does not silently set implementation error to zero.

Equations (10)-(12) exclude an identity dictionary that identifies the
working state with the transfer ground, maps that same phase intervention
to exp(i epsilon F), and maps the native subsequent record readout to the
transfer Heisenberg field. Such a dictionary would assign both zero and a
strictly positive change to one intervention. This does not exclude a
different encoding, a different reference preparation, an observer with
additional coherent controls, or a different physical field dictionary.
It states a testable obligation for that dictionary: its state and
intervention map must reproduce the recovered positive-energy response,
with errors smaller than the explicit margin. Reconstructing (7) alone
does not meet that obligation.

The centered supplied lattice gauge-matter branch, its ground and its
continuum field net obey their own hypotheses. They are not identified with
Omega by this construction. Neither relativistic particle identification
nor a canonical physical-premise closure is asserted.
