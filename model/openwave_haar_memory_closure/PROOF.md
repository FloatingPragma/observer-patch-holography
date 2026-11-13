# Computed interior self-energy and controlled retained memory

All statements use the actual refined full-Haar source at(alpha,beta)=(16,1/4),
with no new physical-time conversion. P is its complete201-state electric frame
and Q=1-P. The source theorem proves D=QHQ>=80Q. The exact all-channel action
constructs B=QHP, including every output representation and fusion multiplicity;
the row/column Schur bound supplies a rational b²>=||B||². A=PHP is the actual
finite compression. This is an orthogonal Hilbert-space decomposition, not a
replacement of the full source by its finite frame.

## An actual computed correction with an infinite-channel remainder

For real z<80 the full Feshbach operator is
F(z)=A-z-B*(D-z)^(-1)B. Its zeros reproduce the full eigenvalue problem below80.
The positive self-energy obeys0<=Sigma(z)<=b²/(80-z). For0<=z<80,
||Sigma(z)-Sigma(0)||<=b² z/[80(80-z)].

Split D=D0-V, where D0=16 E_Q+9 and V=(1/4)Q S Q. Then D0>=89 and ||V||<=9.
Consequently the normalized Neumann ratio is at most9/89<1, on the *entire*
omitted Hilbert space. The executable computes the first term
G=B* D0^(-1) B exactly: every off-block spin-network row contributes its full
weighted outer product divided by16 E+9. This is a finite sum because B has a
finite-dimensional source and Wilson multiplication has finite fusion support.
The infinite Q space is retained in the remainder bound

||Sigma(0)-G|| <= b²·9/[89(89-9)].

Thus ||Sigma(z)-G|| is bounded by this remainder plus
b² z/[80(80-z)]. G is positive semidefinite by its exact weighted-Gram
construction; it is nonzero and contains off-diagonal generated interactions.
The output gives its exact sparse-matrix digest, trace, selected diagonal
values and the full norm remainder. This computes a genuine correction to the
retained generator. It does not assert that the correction is diagonal, purely
a Wilson coupling, independent of energy, or a complete RG transformation.

## Exact memory and explicit tail

Let U(t)=P exp(-tH)P on P. Eliminating Q from the block evolution gives

U'(t)=-A U(t)+integral_0^t K(t-s)U(s)ds,
K(t)=B* exp(-tD) B,
||K(t)||<=b² exp(-80t).

The transfer theorem supplies the positive self-adjoint full semigroup; all
block operations are valid because B is bounded and D has the stated positive
self-adjoint realization. The integrated kernel tail beyond T is at most
b² exp(-80T)/80. Duhamel and the contractivity of U and exp(-tA) give

||U(t)-exp(-tA)|| <= b²[t/80-(1-exp(-80t))/80²]
                  <= min(b²t²/2,b²t/80).

The executable encloses these quantities by rational exponential bounds, without
sampling a time grid or substituting a fitted memory rate.

## A same-clock projection control

For the existing boundary-only Haar isometry J let P_b=JJ*, Q_b=1-P_b and
B_b(t)=J*exp(-tH)J. Exactly,
B_b(s+t)-B_b(s)B_b(t)=J*exp(-sH)Q_b exp(-tH)J.
At equal times this is a positive Gram operator. Take the coarse Haar-constant
probe Omega (its fine image is Haar constant), not the interacting ground state.
The earlier exact spatial Haar integration gives J*S J=0. Together with the
executed full-source moments <S²>=9,<S³>=0,<S⁴>=963/4 and E S Omega=3 S Omega,

||Q_b H Omega||=3 beta,
||H² Omega||²=81 alpha² beta²+3888 alpha beta³+(6999363/4)beta⁴
             =9321795/1024 <96².

For nonnegative H, |exp(-tH)-I+tH|<=t² H²/2 as a spectral-vector estimate.
Hence the boundary-only vacuum defect lies between
max(0,3 beta t-48t²)² and(3 beta t+48t²)². At t=1/128 its positive lower bound is
9/1048576, so the original boundary-only compression is not an autonomous
semigroup in this same clock. Scalar energy shifts cannot remove a nonzero
defect.

The201-mode P retains Omega and H Omega exactly. Therefore its vacuum leakage
is at most t²||H²Omega||/2 and its equal-time vacuum semigroup defect is at most
t⁴||H²Omega||²/4. At t=1/1024 this upper bound is less than1/225 of the certified
boundary-only lower bound, an actual controlled improvement from retaining
interior modes. The exact discarded second-step norm is
||Q H²Omega||²=beta⁴·303/2=303/512>0. No exact autonomy follows. The O(t⁴)
improvement is for this specified vacuum probe; arbitrary retained states can
still have O(t²) defects, and no interacting-vacuum or operator-norm improvement
is inferred from this one vector.

These results close a finite interior-memory premise with an actual generated
operator and all-representation remainder. They do not establish a uniform
iteration over volume and spacing, a derived physical coupling flow, Lorentz
invariance, four-dimensional Yang–Mills identification or the Clay mass gap.
