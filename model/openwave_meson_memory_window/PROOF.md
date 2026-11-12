# A controlled interval for native meson emission

## Complete source and family reduction

The source, filled neutral reference, current, parent preparation and fixed
box coefficient are those of `openwave_meson_ring`. No mass or current is
changed. Its exact spin-zero even-parity reduction has six parent
coordinates and `36 L^2` emitted coordinates. Charge conservation admits
zero or one charged electron and its neutral hole; this is an exact source
sector, not a pair truncation. All other charge, spin, parity, preparation
failure and source-reference records remain in the complete instrument.

Write the full admitted Hamiltonian as

    H = [[A, G C†], [G C, D]],   G=1/64,   C†C=24 I6.

Inside the emitted domain, charge permits only the diagonal native weak
current. Consequently quark family f and lepton family ell are conserved
there. Each of the nine blocks is one universal `4 L^2` matrix plus the
scalar `m_f+m_ell-2`. The injected vector is `C_f1 psi`, where the family
matrix is the actual native finite-word current and `||psi||^2=8`.
The producer checks the complete native matrices and injected columns on
rings of two and four sites. The algebraic identity is independent of L.
No heavy family is discarded by evaluating a light-channel lower bound.

The universal matrix is specified directly by four positions
`(x_q,x_a,x_e,x_h)` and two total-spin-zero vectors. Hadron hopping is -1/8;
its same-site energies are 9/4 and 11/4, and its different-site energy is
127/24. The emitted native contact is

    3 I + 1_(x_q=x_e) Ve + 1_(x_q=x_h) Vh,
    Ve=[[-1/2,-sqrt(3)/2],[-sqrt(3)/2,-3/2]],
    Vh=[[-1/2,-sqrt(3)/2],[-sqrt(3)/2,1/2]].

The electron and canonical hole each have nearest-neighbour hopping -q/4.
The hole spin is the conjugate representation; its physical momentum is
conjugate to its own coordinate. Reflection uses integer orbit columns
with Gram matrix 2I. The parent even orbital matrix is
`[[9/4,-1/4],[-1/4,127/24]]`. Put
`c=sqrt((1+73/sqrt(5473))/2)` and
`d=sqrt((1-73/sqrt(5473))/2)`. Its normalized ground vector is `(c,d)`.
The local injection is `(c,sqrt(3)c)` on equal hadron sites and
`(d,sqrt(3)d)` on unequal sites. Reflection supplies the factor sqrt(2).
These formulas prove the entire reduced construction, including its
normalization. The parent energy is
`m_mu-1+(181-sqrt(5473))/48-3G`.

The frozen box uses the exact rational midpoints of the charged analyzer's
certified effect intervals, as the source compiler specifies. Thus `m_mu`
and `q=(m_mu-2)/2-11G/4` are exact rationals in this member. The family word
is enclosed independently through its complete three-dimensional native
spectral certificate. Polar-frame and isolated-eigenvector bounds certify
`|C_01|>19/100`; no measured mixing angle enters this bound.

## Cubic feedback error

Let `(a(t),b(t))=exp(-iHt)(a0,0)` and `lambda=G sqrt(24)`. Define

    a_free(t)=exp(-iAt)a0,
    b1(t)=-iG integral_0^t exp(-iD(t-s)) C a_free(s) ds.

Unitarity and the exact off-diagonal norm give `||b(t)||<=lambda t`.
Variation of constants in the parent equation then gives
`||a(t)-a_free(t)||<=lambda^2 t^2/2`. Substitution into the emitted equation
proves

    ||b(t)-b1(t)|| <= lambda^3 t^3/6.

This bound includes all returned amplitude and all subsequent emissions.
It does not use a norm bound on D, an energy cutoff or an installed decay
rate. It holds in the infinite-line Hilbert space because D is bounded
for this finite-species nearest-neighbour source. The six-dimensional
memory equation and its exact reconstruction therefore admit a useful
finite-time error estimate without eliminating any emitted family.

## Numerical and volume certification

The producer evaluates the universal first-emission vector by a segmented
Taylor polynomial of the inhomogeneous equation. Every algebraic constant
is enclosed at60 decimal places and rounded to a dyadic denominator2^40.
The integer reflection-orbit assembly introduces no square-root rounding.
The operator error is bounded by
`2 error(q)+error(parent_energy)+2^-41+G error(sqrt(3))`;
the injection error is four times the largest injection-coefficient error.
All other coefficients are exact dyadic numbers. Rows have at most17
nonzero entries including the inhomogeneous column, and the absolute
row and column sums of the unscaled augmented matrix are below207.

Steps are at most1/200 and degree24. The receipt gives the geometric
majorant for every Taylor tail. For binary64 arithmetic, gamma=2^-44
bounds each complex sparse row, division, scalar multiplication and
addition, including construction of the scaled step matrix. The sum of
intermediate Taylor norms is below16; error propagation within a step is
below8. Hence `24*256*gamma` bounds one step's roundoff, and a factor2
bounds propagation across the whole interval. The exact inhomogeneous
semigroup has norm at most `exp(G||psi||t)<2`. The receipt adds coefficient,
roundoff and truncation errors. The final norm lower bound uses a2^60
integer vector and integer square root, avoiding an unbounded floating
norm comparison. No floating propagation diagnostic is used as proof.

The hopping part alone has norm below103. The on-site interaction picture
preserves outgoing support. Two contact-local initial sources therefore
have identical Dyson words on the ring and line before
`R=(L-2)/2` hopping steps. The inherited tail
`2 sum_(k>=R) (103t)^k/k!` bounds their vector difference. It applies to
common spatial instruments with an outside-line complement. It does not
transfer a finite-ring Fourier instrument to an infinite-line momentum
measurement. Combining this tail, the cubic feedback bound, the numerical
bound and `|C_01|>19/100` gives a strictly positive full-source infinite-line
light-channel probability at t=1/100 and t=1/10.

## Complete observations and work

The preparation success is `c^2/4>6/25`. Its failed quantum branch remains
an outcome. The terminal instrument distinguishes no pair, the selected
light pair and every other emitted family; all source complements are
retained. Multiplication by6/25 gives an original-input lower bound.
The full source energy is conserved. Pair/no-pair dephasing removes only
the off-diagonal contact, so its mean work is
`-2G Re<a,C†b>` and its absolute value is at most `2 lambda^2 t`.
The fine family instrument has its exact full-source Lüders work and is
not asserted to have this same two-outcome bound.

Finite rings recur. The large two-site conversion peak at source time720
fails volume stability and supplies no decay width. The controlled
finite interval proves native outgoing emission with explicit error; it
supplies neither a physical strange-particle identity nor an irreversible
lifetime.
