# Arithmetic and small gain terms

All inequalities below concern the actual finite dimensions, 128 repetitions,
80 matching applications per macrostep and stored angles. `check_math.py`
checks their numerical premises with rational enclosures.

Let u=2^-53 and δ bound the error of the computed reciprocal square root
`a=fl(1/sqrt(12))`. The integer second frame column has squared norm 12 and
coefficients at most 2. Comparing direct floating division with multiplication
by a bounds its error by

    sqrt(12) δ + gamma_2 sqrt(12)|a|,
    gamma_2=2u/(1−2u).

With δ≤u this is below `3.5δ+2.001u`. Tensoring the two nearly orthonormal
frames adds their two errors and their product; the full pair-frame error
is below `10δ+6u`. The producer allowance `f=100δ+16u` exceeds twice this
error. The first columns are dyadic and exact.

An ideal single matching is a direct sum of two-point unitary mixes. Each
pair matching is the tensor mix except for disjoint 2×2 contact squares,
which are replaced by one diagonal native phase. Using the absolute-value
coefficient norm ≤sqrt(2) per factor, the independent four-term pair
evaluation needs less than `256u+16τ` relative Frobenius allowance per
matching, including complex products, additions and phase error τ. This
already exceeds the direct real-operation count. The producer uses the
larger `512u+8τ`; all actual τ are independently enclosed and every tighter
independent finite residual bound lies below the producer certificate.
The target reference phase and residual norm rounding are also charged.

For accumulated relative gain G<1 and actual frame error d≤f/2,

    (norm+d)G+d ≤ norm·G+f.

Thus the apparent omitted frame×gain cross term is inside the explicit frame
margin; exact unitarity of rounded centers is unnecessary. The stored whole-
product comparison separately covers its initial frame and all matching
rounding. Both source/reference branches are included when forming the final
2×2 Gram. The `2·stored` allowance bounds their state discrepancy.

A normalized complex dot on N entries needs at most gamma_(8N) times its
Cauchy-Schwarz magnitude, including its real and imaginary accumulation.
The producer's `4 gamma_(16N)` allowance exceeds the errors of the overlap,
norm and complete reduced-reference record together at the checked norms.
The explicit geometric gain guard is retained afterwards. This is distinct
from the much larger native-word and prepared-input errors.

The native matrices are compared to exact dyadic-angle targets with a
48-term alternating rational Taylor enclosure, plus a separate rational
square-root interval. No assumption about a system trigonometric library's
correct rounding supplies the certificate. The source-phase step-duration
rounding is visible in the one-step residual and in the reported duration
check; it is not silently identified with the nominal clock.
