# Native charge-sector phases and a common unit convention

On the complete two-load configuration graph, the magnetic compression is
`T_k=I-L_k/D`, where `D=2M` and the active character is `k` modulo six.
The Dirichlet form is a sum of squared covariant edge differences, so
`L_k` is positive semidefinite. The unchanged zero and parity operators have
integer entries. The extra order-three voltage affects one configuration
edge, between `(0,1)` and `(1,5)`.

Write its off-diagonal change as `alpha E_ij+conj(alpha) E_ji`. For
characters 1 and 5 the parent is `L_3`, with `2 Re(alpha)=-3` and
`|alpha|²=3`; for characters 2 and 4 the parent is `L_0`, with
`2 Re(alpha)=3` and `|alpha|²=3`. Multilinearity of the determinant gives

`p_changed(x)=p_parent(x)-2 Re(alpha) C_ij(x)-|alpha|² p_delete_ij(x)`.

The two-entry term has negative sign because it exchanges the selected
columns. The off-diagonal cofactor is calculated from a single adjugate
column: `B_0=I`, `B_r=L B_(r-1)+c_r I`, where `c_r` are the coefficients
of the parent's characteristic polynomial. Cayley-Hamilton proves this
recurrence. Thus two integer parent polynomials, their two deleted-pair
minors and sparse adjugate columns determine all six sectors. Independent complete-matrix
Faddeev calculations over three prime fields check the implementation.

Exact polynomial root isolation encloses each lowest eigenvalue and its
multiplicity; a Sturm count excludes an omitted lower positive root.
Hermiticity and the Dirichlet form exclude negative roots. For the prepared
cyclic span, `cos(omega)=1-lambda/D`. Directed rational cosine bounds
enclose the principal phase. The full event carrier includes its fixed
phase complements, so these prepared-sector lower phases are not a gap of
the whole carrier.

The displayed unit convention uses one walk-event duration `1/sqrt(D)` in
source phase-time units. It reuses the single electron anchor attached to
the source-clock transition `1/3072`, with no extra mass target. The
displayed energy is this anchor scale times `sqrt(D) omega`. The clock and
time convention are explicit inputs. They do not establish a stationary
joint Hamiltonian, a spatial momentum, or a rest-particle pole. These are
the undressed source/active-memory spectra; a spectrum for the joint
electric-link program requires its actual coupled operator.

For characters 2 and 4, every parity edge phase is trivial. Only the
distinguished edge contributes to the Dirichlet energy of the normalized
constant vector. Its squared phase difference is exactly three. For
`N=choose(P,2)` configurations, the Rayleigh principle therefore gives
`lambda_min <= 3/N`. Since `arcsin(z) <= z/sqrt(1-z²)` for `0<=z<1`,

`D omega_min² <= (6/N)/(1-3/(2ND))`.

This exact bound is evaluated on 240, 960 and 3840 ports under the same
clock rule. It tends to zero in this carrier enlargement, ruling out a
positive uniform prepared-sector gap under that rule. The proof does not
enumerate the large two-load vector. It uses the literal voltage formula
and the unique distinguished configuration edge.

Under the uniform fixed-two-load reference, the probability that both loads
occupy the first twelve-port cell is `66/N`. The corresponding bounded
indicator has variance `(66/N)(1-66/N)`, which also vanishes. This is a
statement about that state and fixed graph-length convention, not an
impossibility theorem for other densities, observables or architectures.
No four-dimensional quantum field or Yang-Mills limit follows here.
