# Correlation degrees of the native connected transfer

The unchanged fair native repair supplies an exhaustive orthogonal
decomposition of each fixed-number population-transfer space. Primitive
correlation degree r has dimension binomial(V,r)-binomial(V,r-1). Its block
embeds unitarily into every occupation q with r<=q<=V-r. This identifies
which transfer modes enter each correlation order; equal eigenvalues can
occur in different orders.

Two disjoint native seams give the degree-two contrast
phi=(n_a-n_b)(n_c-n_d). Four complete E reads measure its lag covariance:
the first two use outgoing contrasts d-2j and the final two use incoming d.
All backaction, 625 raw slots per two-read pair and zero outcomes are retained.
The normalized covariance is h*P2^t h/4, independent of the global occupation
on the same graph for 2<=q<=V-2. It depends on graph size and attempt clock.
The variance is 4*binomial(V-4,q-2)/binomial(V,q); small variance amplifies
preparation and statistical errors.

[The proof](PROOF.md) gives the finite algebra and complete instrument join.
The classical Boolean-slice decomposition is standard; the source attachment
is the literal native transfer, the exact reflected record probe and its
finite controls. Existing positivity, logarithm, exclusion duality and
preparation results are cited with their own scopes. No physical particle,
gauge representation, chirality or relativistic action is selected. The
algebraic intertwiners do not identify individual raw quantum branches across
occupations, and the transfer reference is not a prepared coherent vacuum.

The evaluator checks actual federation levels 0,1,2 without enumerating their
many-body carriers. The full-history interpreter retains its separate
resource and FIFO stops. Exact sparse kernels establish the expectation,
not execution of an exponentially expanded history tree.

```sh
python -m model.native_correlation_sectors.packet --write
python -m model.native_correlation_sectors.packet --check
python -m pytest -q -W error model/native_correlation_sectors
```
