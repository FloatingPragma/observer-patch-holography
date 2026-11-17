# Thermal return records from native current readout

The construction concerns the binary total-two carrier and the declared
independent uniform seam scheduler with fair placement coins. Joint access
to signed currents and an adaptive seam controller are admitted. Every
census uses the existing native controller; no initial configuration is
supplied to its policy. The result supplies a finite approximation to a
closed-history experiment, including its rejected and unresolved outcomes.
It does not identify a laboratory temperature or a physical unit of time.

## Two current censuses implement full-configuration closure

The current controller in `model/native_clock/configuration.py` has two
relevant modes. On a binary input configuration x, a successful restoration
returns x and reports x. Its total success probability is s_R(x). A
successful terminal measurement reports its input y with probability s_M(y);
its final work configuration may differ from y. Fine successful words are
partial injections with singleton input support. A conclusive signed-current
word therefore retrodicts the input, independently of an initial label.

The exact controller enumeration gives all 66 values of both success laws.
Its uniform bounds at four seed scans and 32 retries are

    1-s_R(x) <= 35/2^31,       1-s_M(y) <= 33/2^31.

The maximum native costs are 2648 and 2520 attempts respectively. These
bounds include unsuccessful probes, routing attempts and all timeouts.
The original labels and private tapes are retained as records. On classical
input laws, restoration is an approximate nondisturbing census; no such
statement is made about arbitrary coherent input states.

After the first successful restoration, the reported x gives both initial
loads in the incident-clock region. Every incident signed current updates
these two bits by the exact continuity equation; a remote repair changes
neither bit. Thus every admitted intermediate regional field reading is
computed from this initial census and the body's current records. No extra
measurement or initial-bit disclosure is required for the local field word.
The second census reports the body's endpoint y. Equality of the two
reported configurations implements full-state return. Equality of only
local loads is a different experiment, as the retained control verifies.

Write M_w(y,x) for the subprobability kernel of a completed body word w,
and N_v(y,x) for its waiting-timeout prefixes. For any initial classical law
p, the exact coarse output probabilities are

    P(initial census failure) = sum_x p(x)(1-s_R(x)),
    P(timeout,v) = sum_xy p(x)s_R(x) N_v(y,x),
    P(final census failure,w) = sum_xy p(x)s_R(x) M_w(y,x)(1-s_M(y)),
    P(return,w) = sum_x p(x)s_R(x) M_w(x,x)s_M(x),
    P(different,w) = sum_(x != y) p(x)s_R(x) M_w(y,x)s_M(y).

Their sum is one. The implementation returns every fine census and body
record; the certificate above is its explicitly declared public marginal.
No census failure is interpreted as a configuration, and no probability is
renormalized before a requested return conditioning. Independent fresh coin
and seam choices are required in the respective stages. The terminal
census disturbance occurs after the compared endpoint, which its currents
retrodict; it therefore does not change the comparison criterion.

## Error and resource accounting

Couple each ideal census with its native implementation until a failure.
For every initial law, the complete public marginal total variation from the same
body with ideal censuses is at most

    epsilon_census = max_x(1-s_R(x)) + max_y(1-s_M(y)).

The simpler upper bound is 17/2^29 at the stated settings. The comparison
includes both new failure flags. Successful record masses are dominated by
those of the ideal census experiment, so forgetting their failures and
normalizing generally changes the successful distribution.

There are nine incident seams out of thirty. With h attempts allowed per
incident interval and L intervals, the lost completion mass is

    delta_cap = 1-(1-(7/10)^h)^L.

Coupling with the uncapped experiment gives complete-law error at most
`epsilon_census + delta_cap`. For exactly uniform preparation, successful
return-word masses are submeasures of the ideal uncapped ones. If the ideal
return probability is p_return, their conditional total variation is at most
`min(1,(epsilon_census+delta_cap)/p_return)`. This follows by expressing the
ideal conditional measure as the mixture of the normalized retained
submeasure and the lost submeasure. It does not assume a state-independent
census success probability.

A deterministic thirty-seam sweep removes collisions from every total-two
configuration. The existing native mixing theorem then bounds the binary
preparation error after n independent uniform-seam attempts by

    epsilon_warm <= min(1,33(19079/20000)^n).

All warm-up records remain stored but are ignored by the subsequent
controller and by this output marginal. No product independence from that
tape is asserted. The complete public marginal error from arbitrary full-sector
classical inputs is bounded by

    delta = epsilon_warm + epsilon_census + delta_cap.

For `p_return > delta`, the general normalization estimate
`2 delta/(p_return-delta)` bounds the conditional error; cap it at one.
If that denominator is nonpositive, only the trivial bound one is asserted.
The exact ideal return probability is also computed, and `p_return >= 1/66`
is an analytic uniform lower bound. These bounds charge conditioning rather
than treating a rare return as a free preparation.

The attempt cap for the complete protocol is

    30 + n + 2648 + L h + 2520

with the covering sweep and warm-up requested, and excludes them otherwise.
There is no physical-duration interpretation of this cap. The controller
makes no extra field probes; a different observable requiring a destructive
probe needs a separate compiled restoration and its own error and cost.

## Energy and inverse temperature follow in a specified order

Let A and B be the actual remote and incident subtransition matrices in
column convention, D=I-A, and K=B D^{-1}. Native endpoint enumeration gives
A and B symmetric and nonnegative, `D>0`, and the incident triangle plus
positive exchange kernels gives

    0 < (3/38) I <= T=D^{-1/2} B D^{-1/2} <= I.

The exact source certificate verifies these premises. The upper bound uses
`D-B=I-(A+B)>=0`. The Hamiltonian is consequently defined from the native
transfer, before any return condition or temperature choice:

    H=-log T,                 exp(-n H)=T^n.

A finite polynomial approximation to H has the source-derived operator tail
`(1-alpha)^(m+1)/((m+1)alpha)` when `T>=alpha I`. No target energy spectrum is
used to obtain H. The admitted local field projectors commute with D, so
similarity and trace cyclicity give, for uniform initial preparation,

    P(return)=Tr(K^L)/66=Tr(T^L)/66,
    P(word | return)=Tr(F_L T ... F_1 T)/Tr(T^L),
    rho_L=T^L/Tr(T^L)=exp(-L H)/Tr(exp(-L H)).

Thus the supplied closed-history length selects the internal inverse
temperature `beta=L` for the independently reconstructed H. This faithful
finite Gibbs functional satisfies the finite matrix KMS identity. The
certificate rechecks native short-history counts, the exact inverse and
KMS boundary identities, and the full-return versus local-return control.
The finite compiled census protocol approximates these closed cylinders
with the complete-outcome error bounds above.

This is a reconstruction of a declared ensemble of classical source
histories. It does not prepare rho_L as a coherent native quantum working
state, uniquely identify a full matrix density from local records, implement
every complex effect, or identify stochastic native waiting with the unitary
continuation of H. In particular the unconditioned stationary open observer
has the ground functional of this transfer, whose finite-beta KMS boundary
fails by at least `min(beta,1)*25/594` on the retained local field. The return
condition changes the ensemble and its failures remain part of the record.

## Separate complete-attempt variant

The same two-census adapter also permits L global native attempts, with each
seam selected uniformly. This variant has no incident waiting or timeout.
The first census resolves the entire initial configuration. Each signed
current then updates its two endpoint loads, so the observer records the
full configuration after every attempt without additional probes. This
requires sharing those global signed currents, unlike the regional body.

Let P be the symmetric complete-attempt population kernel, in input-row
convention. For a full configuration word `(x_1,...,x_L)`, the ideal native
path weight from x_0 is

    (1/66) product_(j=1)^L P(x_(j-1),x_j).

The actual successful weight multiplies this by `s_R(x_0)s_M(x_L)`;
initial and terminal census failures complete its probability law. The
same marginal comparison has error at most epsilon_census. Its return
probability before census error is `Tr(P^L)/66`. The finite certificate
checks every actual seam/coin branch and all words, including zero-weight
words, for L=1 or 2. For arbitrary L, the displayed product and the same
coupling proof apply. Native counting and the reconstruction theorem in
`model/native_quantum_reconstruction` interpret these words using the
separate Hamiltonian `-log P`, rather than the incident-clock `-log T`.

With a covering sweep and n warm-up attempts the cost is
`30+n+2648+L+2520`; the complete public marginal error is bounded by
`epsilon_warm+epsilon_census`. Its conditional bound charges the same return
normalization. A finite runtime record contains every private tape, but the
comparison theorem concerns the explicitly defined common public marginal;
it asserts no product independence from the stored preparation or census
records. Neither clock variant makes an operational all-complex-effect or
physical-temperature identification.
