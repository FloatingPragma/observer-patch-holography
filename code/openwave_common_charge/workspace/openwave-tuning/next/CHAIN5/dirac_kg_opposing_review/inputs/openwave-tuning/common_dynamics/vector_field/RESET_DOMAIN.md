# Correlated native reset inputs

The reset domain is the complete 78-dimensional, twelve-port, total-two
native sector, with arbitrary quantum references and correlations between
carriers. This is stronger than a classical mixture assumption. It follows
from the actual fine source instrument, rather than from the minimum of
the diagonal success probabilities alone.

`model/native_clock/configuration.py:source_audit` checks both C and E laws
on every seam and all 78 input columns. Each fine map is a scaled partial
injection. Signed currents reconstruct its predecessor. The configuration
controller declares success only when its retained current tape identifies
one initial configuration x. Its deterministic binary restoration then
routes the final configuration to the requested fixed target t. Therefore
every conclusive fine reset Kraus operator has the form

    K_w = a_w |t><x(w)|.

The success effect is consequently

    E_success = sum_x s_x |x><x|,
    s_x = sum_{w:x(w)=x} |a_w|^2.

Its least eigenvalue is exactly min_x s_x. The exact dynamic program
`success_probabilities(4,32,full=True,mode='reset',target=t)` supplies those
78 probabilities for each of the two targets used in COMPOSITION.json.
This diagonal conclusion uses the rank-one tape property; it would not
follow from evaluating diagonal matrix entries of an arbitrary effect.

For a correlated input rho_CR, the successful fine branch is

    |a_w|^2 |t><t| tensor <x(w)|rho_CR|x(w)>.

Keep w and every failure tape in the full instrument. For comparison,
define the ideal fine success branch using weight |a_w|^2/s_x; all s_x are
positive. This ideal instrument also retains a resolved initial label and
its reference correlations. The successful difference is a positive map
with trace at most epsilon=1-min_x s_x, while the retained failure output
has trace at most epsilon. Their orthogonal flags give half-diamond
distance at most epsilon for arbitrary untouched references. After the
private resolved label and tape are marginalized at this comparison
boundary, the ideal output is the fixed blank tensor the original reference
marginal. Conditioning on success does not preserve that marginal exactly.

Each later reset acts on a conditional state that may remain correlated
with every earlier tape. The same operator bound still applies. Summing
the conditional first-failure bounds gives

    P(any timeout) <= 3203(1-p_blank) + 3200(1-p_ready).

Thus the COMPOSITION union bound requires no independence of initially
unknown carriers. An exact product of unconditional success factors is
valid only under the additional independence/configuration assumptions of
that particular preparation. In a general correlated preparation, the
complete law instead uses conditional success factors and retains every
first-failure record. Arbitrary-reference bounds for the subsequent
initialized coherent evolution are a separate operator-composition claim.

The pointer's fresh ready carrier is prepared with this same reset
instrument before its native readout word. The added clock reset and
pointer in CLOCK_JOIN use the identical argument. The proof does not
silently prepare a coherent uniform affine field state; the field source
here is the declared definite F=0 preparation.
