# Exact momentum flatness on the periodic two-cell source

The source is the complete periodic eight-cell Hamiltonian
\(H=H_{\rm loc}+2N_{\rm nonblank}+2^{-34}T\) of
`openwave_neutral_stationary`. Every fermion, scalar and 119-state link is
retained. Its rank-1024 neutral band and actual ground energy are inherited
on this exact member. No continuum, larger-volume or antiperiodic transfer
is made.

## Three exact complete-carrier symmetries

For coordinate axis \(i\), let \(Q_i\) be fermion parity on the four cells
with \(x_i=1\). Let \(R_i\) exchange the two oppositely anchored links on
each parallel axis-\(i\) pair and apply the full representation inversion
\(J|r,a,b\rangle=|\bar r,b,a\rangle\) to both. All matter and scalar cells
are fixed. The exact coordinate and current identities are

\[
 J U_{ab}J=U_{ba}^{\dagger},\qquad
 J L_0^A J=L_1^A,\qquad J^2=1.
\]

The executor checks all 59 coordinate identities and all 24 endpoint-current
identities on the full 119-state link. Swapping the reversed edge restores
each endpoint current to the same cell, so \(R_i\) preserves every Gauss
constraint. It commutes with \(H_{\rm loc}\) and the link penalty. Every
native local source word has even fermion degree, so \(Q_i\) commutes with
the complete local source too.

The periodic side-two hopping terms all have coefficient \(i/4\).
Inversion/exchange maps the axis-\(i\) forward word to the other edge's
adjoint word, with opposite coefficient. The Pauli matrix is Hermitian;
the spin and native coordinate transpositions retain every entry. Therefore

\[
 R_iT_jR_i=(-1)^{\delta_{ij}}T_j,
 \qquad Q_iT_jQ_i=(-1)^{\delta_{ij}}T_j.
\]

The second identity follows because precisely an axis-\(i\) hop crosses the
selected plane. Hence

\[
 C_i=Q_iR_i,\qquad C_i^2=1,\qquad[C_i,C_j]=0,
 \qquad[C_i,H]=0.                                      \tag{1}
\]

All 72 axis/edge word cases are checked. Occupied full-H column controls
include nonblank links, multiple fermion occupations and quantum scalar
states. Omitting the matter parity or representation inversion produces a
nonzero commutator. Substituting the antiperiodic boundary source also fails.

Let \(T_j^{\rm lat}\) be the complete periodic translation, including its
fermionic occupation signs. Translation along \(i\) exchanges the two
matter planes, so

\[
 T_j^{\rm lat} C_i (T_j^{\rm lat})^*
   =(-1)^{N\delta_{ij}}C_i.                             \tag{2}
\]

The link permutation/inversion commutes with translations. Equation (2)
is checked on arbitrary full occupation/link/scalar keys. The neutral
raw-lepton-seven, baryon-zero sector has odd parity by the native identity
\((-1)^N=(-1)^{L+B_3}\).

## Exact dispersion obstruction

In the odd sector, \(C_i\) flips the \(i\)-th translation character and
preserves all other characters and the energy. The three commuting \(C_i\)
thus unitarily identify the full Hamiltonian at all eight momenta. Every
odd-sector eigenvalue has multiplicity divisible by eight, and every
corresponding momentum splitting is exactly zero. This is an operator
identity to all orders in the hopping coefficient; its remainder is zero.

Inside the known atomic band, the neutral local factor has even parity and
the seven ground factors have odd parity. On its location projector \(P_y\),

\[
 C_iP_y=(-1)^{y_i}P_y.                                  \tag{3}
\]

The executor verifies the complete 1024-row graded translation and symmetry
algebra. Each of the eight simultaneous \(C_i\) sectors has rank 128. Let
\(E\) be the isolated interacting-band projector. The canonical graph
isometry from \(P\) to \(E\) is a function of their projectors and therefore
intertwines both (1) and the translations. Consequently its exact effective
operator is a direct sum of eight unitarily equivalent 128-row blocks.
Their internal spin energies need not be equal; their multiplicity across
momentum is exact.

The ordinary second-order formula is
\(K_2=-PTQ[Q(H_0-a)Q]^{-1}QTP\). Returning all links to blank after two
hops uses one edge twice. Each endpoint then has even CAR degree, so its
parity cannot exchange the neutral location with an odd ground factor.
This gives zero off-diagonal position blocks at second order. The generic
canonical graph error \(6(|\lambda|\,360)^3/d^2\) cannot by itself prove
exact flatness; (1)-(3) close the higher-order question without computing
all charged resolvents or importing a dispersion law.

## Complete spatial records

Define the dressed location projectors

\[
 D_y=E\prod_{i=1}^3\frac{1+(-1)^{y_i}C_i}{2}.
\]

They sum to \(E\), have rank 128 and commute exactly with \(H\). A local
atomic input in \(P_y\), followed by any spectral filter of \(H\), stays in
the same simultaneous symmetry sector. Every different atomic-position
effect \(P_x\), \(x\ne y\), has probability exactly zero at every real
source time. For a state in \(D_y\), the inherited projector angle
\(\|E-P\|\leq\eta=v/d\) gives

\[
 \Pr(P_y)\geq1-\eta^2,\qquad
 \Pr(I-P)\leq\eta^2,
 \qquad\Pr(P_x)=0\quad(x\ne y).                        \tag{4}
\]

The complete virtual-link/profile complement \(I-P\) is retained and can
vary in time. These exact forbidden position records are stronger than a
finite-time small-spreading estimate, but exclude using this member as a
moving neutral massive particle. A larger periodic geometry does not fix
every matter coordinate under the link reflection used above; it requires
its own source isolation and propagation calculation. No nonzero dispersion
on that different source is inferred here.
