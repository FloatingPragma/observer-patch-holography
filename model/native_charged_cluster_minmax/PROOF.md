# Min–max bounds for the charged Wilson cube

The source is the finite orthogonal Wilson cube admitted by
`field_record_addenda/native_orthogonal_wilson_charged_band/model.py`.
Every frozen input hash is checked by `band.admit()`. The construction uses
\(g=1/64\), total photon capacity three, the fixed compensated charge frame,
reference cost \(4|R|\), and native-complement completion \(c=21/2\).
The helper's separate \(23/2\) completion is evaluated as a separate row.
No physical unit or continuum interpretation enters the calculation.

Put \(\lambda_1=\sqrt5/2\),
\(\lambda_2=\sqrt{17-4\sqrt3}/2\),
\(v=1/16+1/(12\sqrt6)\), and \(q_*^2=(13+5\sqrt6)/48\).
For charge \(a=0,1\), let \(C_0=8,C_1=67/8\),
\(\nu_a=24gq_*+C_ag^2\), and \(d_a=(12v+C_a)g^2\).
The frozen bridge supplies a signed code-vacuum shift upper bound
\(\Delta_U\). The free lowest positive two-dimensional projector is
constructed exactly as a polynomial in the literal 32-mode matrix.
Its rank and idempotence are checked using rational matrix identities and
the field \(\mathbb Q(\sqrt5)\).

On this free eigenspace with photon vacuum, the odd phase has zero
compression. The finite-coordinate inequality
\(|1-\langle e^{igQ_e}\rangle|\le g^2v/2\), together with twelve bond
norms bounded by two, gives \(\|PVP\|\le d_1\).
The min–max principle bounds the two trial energies above by
\(\lambda_1+d_1\); the third code energy is at least
\(\lambda_2-\nu_1\). The native complement is reducing and has lower bound
\(c-\Delta_U\). Consequently the full gap above the rank-two cluster is

\[
\min(\lambda_2-\nu_1,c-\Delta_U)-\lambda_1-d_1
>0.17175185.
\]

All endpoints in `results.json` are outward rational bounds. The code
gap is greater than \(0.19242959\). The multiplicity is a
fixed-charge source statement, without a physical spin interpretation.
The exact corner rotation lift has cube \(-I\); its traces on the free
two-space are \(1,-1\). The corner mirror composed with time reversal has
square \(+I\), so this symmetry calculation supplies no Kramers assertion.

For the charged current-square bound, diagonalize the positive real photon
covariance into Hermitian mode currents \(J_\mu\). The frozen bridge proves
\(C h^* C^\dagger=-h\), \(C J_e^* C^\dagger=J_e\) and
\(\operatorname{Tr}J_e=0\), so the filled negative sea projector \(P_-\)
has \(\operatorname{Tr}(P_-J_\mu)=0\). For a normalized added positive
orbital \(f\), the fermion Wick identity gives

\[
\langle d\Gamma(J_\mu)^2\rangle_{P_-+|f\rangle\langle f|}
=\operatorname{Tr}[P_-J_\mu(1-P_-)J_\mu]
 +\langle f|J_\mu^2-2J_\mu P_-J_\mu|f\rangle.
\]

The subtracted Pauli-blocking term is nonnegative. Summing the mode
covariances therefore bounds the charged current square by the neutral
one plus the operator norm of \(T_4/4+T_6/(2\sqrt6)\). Its exact Hermitian
row-norm bound is \((15/16)/4+(5/8)/(2\sqrt6)\). The estimate is uniform
over the charged free two-space and justifies the executed
`charged_A_squared_upper = neutral_A_squared_upper + cpbound` bound.

For energy intervals, the checker separately bounds the odd and even
photon-parity residuals. The current-square row bounds are \(15/16,5/8\).
Let \(r_a^2\) be the sum of their squared residual bounds, as evaluated in
`audit.compute()`, and \(q_c=\lambda_{2,L}-\nu_1-\lambda_{1,U}>0\).
Schur-complement lower bounds give

\[
L_0=-d_0-r_0^2/(2-\nu_0),\qquad
L_1=\lambda_{1,L}-d_1-r_1^2/q_c.
\]

The upper bounds are \(d_0,\lambda_{1,U}+d_1\). In the charged lower-bound
argument, an eigenvalue above \(\lambda_1-d_1\) satisfies the claimed
bound directly; below it the complement denominator is at least \(q_c\).
Subtracting the actual neutral energy and the exact reference cost four
therefore puts the charged transition, divided by \(m\), in
\([L_1-d_0,\lambda_{1,U}+d_1-L_0]\), enclosed by
\([1.11134603,1.12278560]\).
The raw reference cost is four, rather than \(4m\).

The [perturbative coefficients](self_energy/PROOF.md) are a separate
calculation. Their unproved compressed-matrix transfer allowance does not
enter these min–max bounds.
