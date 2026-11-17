# Finite-word probability and conditioning certificates

The native endpoint table gives the position mean
P=I−L/10 with optical eigenvalues3/5 and2/5. For refinement epsilon=2^−k,
the retained active/blank choice register gives P_epsilon=I−epsilon²L/10.
In each optical mode the native current-reversal word has eigenvalues
z=p±i sqrt(1−p²), where p=1−epsilon² lambda/10 and lambda=4,6. Its invariant
four-dimensional sector is the direct sum of these two signed pairs. The
choice decoder restricts there to D=(I2,0). The finite source preparation,
all native incidence means and complete re-encoding preserve this sector.
The source's separate raw/folded-register proof justifies this restriction;
the unused padding sector does not acquire population.

For a counter of size Q, the actual outcome r has word multiplier

F_r(z)=Q^−1 sum(t=0,…,Q−1) z^t exp(−2 pi i r t/Q)
      =(1−z^Q)/(Q(1−z exp(−2 pi i r/Q))).

The accepted counter outcomes1,…,Q/2−1 give the Schur matrix
G_ij=sum_r F_r(z_i) conjugate(F_r(z_j)). Consequently the exact completely
positive selected map is B[G entrywise (B† rho B)]B†. All four signed
populations and their coherences are retained. Completeness of the full
counter follows from finite Fourier orthogonality. The omitted outcomes
are its complementary instrument, not a normalization adjustment.

The code propagates the actual unnormalized cooling vector, accepted counter
and decoder through every optical path block. Four native target eraser
outcomes have operator P_slot^4/9. The accepted plane filter is scalar on
this invariant source sector. The same re-encoding map acts on every path
block. The selected returned clock is
tau_s=(rho00+rho11+s rho01+s rho10)/2. Its trace, all marker/selector flags and
the plane-filter factors give the complete selected preparation probability.
An independently prepared pair has the product cost. This calculation
includes the terminal selected path factor.

The exact conditional state tau_s/Tr(tau_s), rather than an ideal ray, gives
signed compact populations w_plus,w_minus. The pair's unwanted population is
1−w_plus²−w_minus². Opposite-sign zero-frequency products are included in that
complement. Its Ramsey moment at stride n is
[sum_j w_j z_j^(n 2^k)]². This is the full product-state observable justified
by the separately executed independent-source admission and factorization;
agreement of marginals alone would not justify the product.

Every arithmetic operation encloses its exact value. The Fourier root uses
pi=16 atan(1/5)−4 atan(1/239). The tangent subtraction identity gives
4 atan(1/5)−atan(1/239)=pi/4 in the specified positive branch. After N terms,
the arctangent series remainder is bounded by
|x|^(2N+1)/((2N+1)(1−x²)). The sine/cosine series use the first omitted term
divided by1−x² as an absolute remainder bound. No floating normalization of a
root is used. Matrix powers, geometric sums, conditioning denominators and
probabilities use outward interval operations. The installed mpmath1.3.0
directed interval arithmetic and square root are the arithmetic substrate;
the interval-complex conjugation is explicitly (real,−imaginary).

Exact dyadic endpoints are the machine-readable certificate. Rounded decimal
strings are display values. A narrow interval around a known source value
certifies arithmetic for the declared finite model, not measurement precision
for an apparatus with additional perturbations.

The record-only frequency reader applies the invariant-band cosine inequality
to these binary laws with its separately certified population bound. Its
inverse cosine is enclosed using sqrt, six exact half-angle reductions and
the bounded arctangent series. Every compatible alias at the second stride
is intersected with the first-stride interval. The source eigenvalue enters
only the subsequent prediction comparison. The single calibration multiplies
the resulting interval by511 sqrt(5)/2 after decoding.

## Propagation of complete implementation error

Let a complete instrument, including failure flags, produce the selected
subnormalized state rho with p=Tr(rho)>0. Let a perturbed complete instrument
produce sigma, and suppose ||rho−sigma||_1≤e<p. Then its selected probability
q is positive and |p−q|≤e. The triangle inequality gives

||rho/p−sigma/q||_1
 ≤ ||rho−sigma||_1/p + ||sigma||_1 |1/p−1/q|
 ≤ 2e/p.

Thus conditional trace distance is at most e/p. A certified lower bound p_min
may replace p. If complete source channels have half-diamond errors delta_j,
telescoping composition and contraction give e≤2 sum_j delta_j for any input,
including retained references and adaptive classical flags. Conditioning
therefore contributes at most2 sum_j delta_j/p_min to any selected binary
probability or selected invariant-band population.

For two independently perturbed selected preparations with trace-distance
bounds t1,t2, the product state differs by trace distance at most t1+t2.
Its binary measurement probabilities and unwanted-band population acquire
that same allowance. A measurement-channel half-diamond allowance adds to
the binary-probability allowance. These terms must be propagated into the
reader's eta and delta before applying the interval theorem. Source branches
before selection and the selected denominator cannot be omitted.

The exact declared mathematical apparatus has zero additional channel
perturbation. The outward computation certifies its finite probabilities and
impurity. A claimed perturbed implementation must supply the complete error
budget and conditioning bound above; a small floating residual cannot replace
that obligation. Increasing the actual accepted source probability directly
improves this conditional error budget without changing the theorem.
