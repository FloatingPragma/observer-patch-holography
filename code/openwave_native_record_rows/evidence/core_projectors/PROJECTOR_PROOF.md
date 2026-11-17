# The compact core is an actual source spectral band

This derivation follows the sealed native bound-core response run. It is
retrospective mathematical work. The new archive analysis was separately
preregistered before execution; the earlier source run was not blind to this
later theorem. No new native endpoint execution is claimed here.

Keep the exact source, graph, sectors, normalizations and response instrument
of `../20261114_codex_bound_core_response/`. Let T_a translate both occupied
sites by a on V=(Z/LZ)^d, L even and at least4. For chi(a)=(-1)^(sum a),

    Pi_K = (1/|V|) sum_a chi(a) T_a

is an orthogonal projector: translation composition gives Pi_K^2=Pi_K and
inversion gives selfadjointness. It selects the whole total-momentum
(pi,...,pi) sector, including all relative-displacement orbits. Shorter
unordered-pair orbits have their actual stabilizer signs; none are silently
counted as generic full-length orbits. The native source transfer commutes
with translations because its edge/coin apparatus is uniform.

Let D_adj be the diagonal projector onto adjacent unordered pairs. It
commutes with every T_a. Extend a pair function by zero on coincident
arguments. On a function in the Pi_K sector,

    f(x+e_i,y) = -f(x,y-e_i),
    f(x-e_i,y) = -f(x,y+e_i).

All off-diagonal hopping terms cancel, even adjacent to the forbidden
coincidence diagonal: both paired terms are then zero. A nonadjacent pair
has4d mixed boundary edges, while an adjacent pair has4d-2. Therefore on the
ENTIRE Pi_K sector the native generator obeys the exact operator identity

    W J_2 Pi_K = (2d I-D_adj) Pi_K.

Consequently

    Pi_bound = D_adj Pi_K
             = (2d I-WJ_2) Pi_K

is the actual lower spectral projector in that sector. It has energy
(2d-1)/W, while its orthogonal K-sector complement has energy2d/W. The
adjacent orbit for each axis has length|V|, giving rankd. In configuration
coordinates its diagonal is1/|V| on each adjacent pair and zero elsewhere.
The normalized source state functional is therefore

    rho_bound=Pi_bound/d,
    Tr(rho_bound A)=Tr(Pi_bound A)/d.

This is positive, normalized and source selected by exact spectral data;
it is not an arbitrary classical adjacent-pair state relabeled "bound".
The complete archived native one-event matrices determine J_2. The analyzer
checks the integer multiples of these operator identities on every matrix
entry, as well as symmetry, idempotence, rank and the diagonal measure.

Translation is an observable relabeling of decoded configurations. A twisted
return word with final load T_a times its public initial load reads
Tr(T_a P^m) from complete retained histories. The signed sum with chi(a)
implements the momentum projector as a source count functional, keeping
negative scores and all nonreturns. Positivity belongs to the reconstructed
projector, not to each score separately. The one-event matrix is obtained
from the frozen full native branch archive; no target field kernel enters.

For the specified core-current experiment, the first recorded outcome is the
position bond b=(x,i). All subsequent comparison preparation, frozen center
and native current outcomes depend on b and the recorded edge/coin. Its
effects on the original pair space are diagonal. Thus every complete event
outcome effect E_r satisfies

    Tr(rho_bound E_r)
       = (1/W) sum_adjacent_b E_r(b,b).

The right side is exactly the executed uniformly prepared bond experiment,
including each zero current and comparison outcome. This is a complete
instrument-statistics identity on the actual spectral functional. It neither
claims a coherent working-state preparation nor requires that additional
service for this finite diagonal test. An arbitrary later off-diagonal probe
would be a different task and is not implemented by the classical mixture.

The intrinsic center is measured before the response, so the source state
is not falsely described as center-of-mass localized. Its measured offsets
are +/-e_i/2, giving R^2=1/4 and R_parallel^2=1/(4d). The actual compared
current response and independently executed point denominator, with fixed
total charge2, give g_eff(k)/g_eff(0)=cos(kappa/2) along diagonal k. The
half-suppression onset is |k|R_parallel=pi/3. The same band, source law,
record instrument, charge normalization and attempt clock supply all parts
of this statement. Universal beta functions, SI charge, coherent preparation
and a fixed physical-radius continuum are outside this finite claim.
