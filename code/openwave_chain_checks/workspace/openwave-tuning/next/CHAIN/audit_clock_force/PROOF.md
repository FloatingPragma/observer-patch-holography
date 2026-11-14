# Independent clock and force audit

The frozen finite quantum clock and the frozen finite classical force
measurement satisfy their stated mathematical claims. Their input admission
interfaces have failures on invalid densities, rates and durations. The
machine report separates these findings. No source identity, preparation or
clock equivalence between the two constructions is proved.

The clock source file has SHA256
`eee29493c0cc8a2bfdaf2868fdce82e10d055535f8ee0e6a0b4b818ef9fb38db`.
The force manifest has SHA256
`f4f2f6b490aef1a87eb0fd6b5c07c429d5f394baf4414e0f39b7dc309babba1e`.
`AUDIT.json` binds all transitive files read by the finite-member audit.
The producer's `verify.py` and `reference_review.py` are not executed as
independent evidence.

## Clock instrument

The independent reconstruction uses the declared 384 ordered mask/record
states, direct bit-count CAR signs, the six half-record updates, periodic
electric energy and cyclic face translation. It reconstructs every entry of
the actual Hamiltonian. The global interacting ground projector and the
four fair-reservoir branches are reconstructed from this Hamiltonian and a
direct seven-mode signed swap. Their preparation agrees with the pinned
array and preserves every initial/final reservoir label.

For each basis input, seam, coin and selected C/E law, the integer repair
fixes a current, target mask, record update and CAR sign. Grouping precisely
by the published record yields 12 C, 18 E and 30 MIX Kraus operators. Each
reconstructed operator agrees with the adapter. The Gram sum is the identity
on all 384 states. This checks amplitudes and coherences as well as the
classical column sums. The three seams have weight 1/3 and each coin has
weight 1/2. MIX also has law weight 1/2.

Let J be the complete native map and C(rho)=-i[H,rho]. For uniformization
rate Lambda at least lambda, the marked instrument consists of
sqrt(1-lambda/Lambda) I and sqrt(lambda/Lambda) K_r. At ordered times
0<t1<...<tm<T, the complete path density is

    exp(-Lambda T) Lambda^m A_b rho A_b* dt1 ... dtm,

where A_b contains every unitary interval and every selected marked
operator. Summing marks at fixed times preserves trace. Integrating the
ordered simplex supplies T^m/m!, and summing m gives total trace one.
Every native zero-current record therefore belongs to the instrument.

The infinitesimal generator is C+lambda(J-Id). If z counts native attempts,
x counts moving attempts and y counts auxiliary idles, the generator is

    C + lambda z (J_zero+x J_move) - lambda Id
      + (Lambda-lambda)(y-1) Id.

At x=1, trace preservation of J closes the generating-function equation.
Its trace is exp(T[lambda(z-1)+(Lambda-lambda)(y-1)]). Native attempts and
idles are independent Poisson counts. Moving attempts have a state-dependent
law and cannot replace native attempts in this formula.

The zero- and one-attempt states are

    rho0 = exp(-lambda T) U(T) rho U(T)*,
    rho1 = lambda exp(-lambda T) integral_0^T
           U(T-s) J(U(s)rho U(s)*) U(T-s)* ds.

The remaining state is exp(TL)rho-rho0-rho1, the positive sum of all paths
with at least two attempts. This is an unlimited-process coarse outcome.
The audit reconstructs rho0 spectrally and rho1 by independent ordered-time
quadrature, with 18/24-node agreement. Direct density ODEs using reconstructed
H and K independently check the complete actual prepared output and its
tail. All twelve reservoir/count states are retained. The tail probability
is approximately 0.0360514933 at T=0.37 and lambda=0.8.

For each conditional history the coherent durations sum to T, the commit
count is N, and the declared interval-service count is N+I+1. A concrete
zero-current mark plus idle has one native commit and three coherent
services. Its conditional state and ordered-time density agree with direct
operator multiplication. The idle is an identity operation; the zero-current
native mark is a measured branch and can alter coherence.

Multiplying H and lambda by a gives aL and therefore a common time scaling.
Multiplying lambda alone leaves C unchanged. The audit rejects equality of
these two interventions on an actual 36-state number sector. Independently,
the zero-attempt conditional state changes an actual native terminal
probability by about 0.001995954. Hence deterministic T=cN is incompatible
with the retained positive-probability zero-attempt branch. Deleting native
no-ops removes about 0.756915923 of the probability of one native attempt
on the declared charged preparation.

The number-block semigroup solver is restricted to block-diagonal densities.
Full-history and coherent-service maps retain cross-number coherences. No
native compiler or physical duration for the coherent H service is supplied.

## Force instrument and finite caps

The second construction has two classical one-load copies on 96 ports,
288 seams per copy and a fair literal native coin. The audit checks actual
native outputs for all 11,520 first-event marks across the two endpoint
preparations and five coupling values. It evaluates the 512-round field
operator by the independent spectral polynomial

    G_K = [1-(I-L/12)^K] L^+ on the neutral subspace.

The spectral pair energies agree with the frozen local-stencil table. No
9,216-state generator or semigroup is constructed.

For a state s, let r_e(s) be each scheduled seam rate and Lambda_s their
sum over both copies. Each fine copy/seam/coin mark has intensity r_e/2.
The literal rule transfers the load on exactly one of its two coin branches
when the seam is incident. Its other outcomes are retained holds. Therefore
at zero coupling the total attempt rate is 576 nu and the transfer rate is
6 nu.

For cap one at horizon T, the no-event probability is exp(-Lambda_s T).
Every fine mark has integrated probability g_s r_e/2, where
g_s=(1-exp(-Lambda_s T))/Lambda_s. All 1,152 marks carry the cap flag and
terminal position, including holds. Along with the no-event output these
are 1,153 outcomes with total probability one. Their fine waiting density
is exp(-Lambda_s t) r_e/2 dt for 0<t<T.

The cap-two proof uses the native rate kernel A, including its holding
diagonal, with column sum Lambda. From an initial x, only 13 first-step
positions y occur. Put

    h_xy = integral_0^T exp(-Lambda_x t) exp(-Lambda_y(T-t)) dt.

The zero-event probability is exp(-Lambda_x T), the one-event non-cap
probability at y is A_yx h_xy, and the capped probability at z is

    sum_y A_zy A_yx (g_x-h_xy)/Lambda_y.

Summing z uses the exact column sum Lambda_y. The total cap probability is
sum_y A_yx(g_x-h_xy). Combining all three classes gives one. The audit
evaluates these expressions on the actual 13-state first-step support and
checks h_xy against independent waiting-time quadrature. Actual frozen
executions at caps one and two are replayed, including their waits, marks,
native currents, holds, terminal flags and complete path densities.

## Force observable and controlled arithmetic transfer

With rho=delta_x-delta_y, the supplied controller energy is
U=(kappa/2)rho^T G_K rho. An oriented charge-sigma hop changes rho by
-sigma B_e. Expanding the same quadratic form gives

    Delta U = -kappa sigma B_e^T G_K rho
              + (kappa/2) B_e^T G_K B_e.

The finite self-response term is required. On the measured seam its value
is about 0.156256764, while the force is about 0.034922890. Omitting it
therefore gives a different force even on this one finite member.

The scheduler supplies r=nu exp(-beta Delta U/2), so the moving rate is
w=r/2. Reverse rates have the detailed-balance ratio exp(-beta Delta U).
For the equally weighted endpoint preparations x=0 and x=1, with y=9,
the unconditional score j_e/g_s has expectation

    J = [w(0,y to 1,y)-w(1,y to 0,y)]/2
      = (nu/2) sinh(beta F/2),
    F = U(0,y)-U(1,y).

The score uses all capped outputs and is zero on no event. The measured
drift is approximately 0.0349512914 per supplied source time. The zero-field
second-current rate is nu/2, so D0=nu/4 and

    |J-beta D0 F| <= nu |beta F/2|^3 cosh(|beta F/2|)/12.

The audit checks this relation at five couplings. Removing cap outcomes
removes about 0.438 probability; removing holding marks removes about
0.433. Conditioning on the no-cap result gives drift zero. Relabelling the
same law with doubled nu changes its no-event probability from about
0.5621 to 0.3160 and is rejected.

The arithmetic variant uses an independently executed integer recurrence
with 40 fractional bits, zero start and 512 rounds for each candidate pair.
All 22 exact rational energies and residual certificates agree. Let
theta=(23/24)^512, eta=2 theta+(408/2^40)(1-theta), delta=eta+2 theta,
and ell=beta kappa delta. The two energy functions differ by at most
kappa delta. Their seam log-rate ratio has magnitude at most ell.

At common positions, only twelve incident seam rates can differ. The
frozen bound is eight nu per incident seam. Couple shared mark intensities
and identical fair coins until the first disagreement. The total mismatch
intensity is at most 96 nu (exp(ell)-1). Stopping both histories at the same
finite attempt cap preserves the coupling, including all cap and hold
outcomes. The complete marked path total variation therefore obeys

    TV <= min(1,96 nu T ell/(1-ell)),  0<=ell<1.

The outward rational bound is 670801/10^15. With frozen exposure score and
g_s>=67/100000, its expectation error is at most 670801/335000000000.
These comparisons use exact-real exponentials and the ideal frozen
stencil. They contain no general bound for floating library evaluation.
The static field theorem supplies no dynamic trajectory refinement or
continuum particle-force theorem.

## Admission failures and consumer boundary

`MALFORMED_INPUTS.json` gives sparse reproducible fixtures. A trace-one
Hermitian matrix with a negative eigenvalue passes the frozen clock's
semigroup and trajectory admission. A matrix containing NaN passes its
semigroup admission and returns NaN. `require_block_diagonal` checks only
off-number entries; trajectory checks trace without positivity or
finiteness. The semigroup APIs must preserve valid subnormalized branch
inputs when these guards are hardened.

The force `first_event` callable accepts a negative horizon or one negative
seam rate and returns negative probabilities whose sum is approximately
one. Normalization alone is insufficient admission. These malformed inputs
are outside the frozen finite member's hypotheses. Claims of general
defensive input validation fail; the valid-member mathematical results
pass. A consumer must validate domains or use the exact pinned member.

The source domains are distinct: interacting 384-state C10 with charged
reservoir reference versus two mobile classical 96-port copies. Names such
as native, clock, field or force do not supply a source/state/reference/
clock transfer. Neither adapter supplies physical seconds, a C10-to-force
transfer, mechanical inertia, Maxwell/Lorentz dynamics, or a derived
electromagnetic coupling. Architecture choices are retained as declared
inputs. No public source, acceptance configuration or official row is changed.
