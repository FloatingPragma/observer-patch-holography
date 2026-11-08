# Kinetic winding and the free observer field

The original Weyl hopping has a complete eight-node momentum-space winding
catalogue. On the centered full-ground family
`c=v/a`, `0<=U=o(c^(2/7))`, fixed `v>0`, the existing
[free-limit proof](../ballistic_free_boundary/README.md) maps these nodes to
eight massless Weyl fields on the same state and clock. This is an internal
source-to-field map for that family. It assigns no Standard Model particle,
and the default interacting ratio `U/c=8` is outside its hypotheses.

## Exact winding

The literal source link `i c sigma_j/2` and its adjoint give

\[
 h(k)=-c\sum_{j=1}^3\sigma_j\sin k_j.
\]

All zeros on the momentum torus are `k_b=pi b`, `b in {0,1}^3`.
Writing `k=k_b+q`, their normalized Jacobians and winding degrees are

\[
 A_b=\operatorname{diag}(-(-1)^{b_1},-(-1)^{b_2},-(-1)^{b_3}),
 \qquad \chi_b=\det A_b=-(-1)^{|b|}.
\]

For an outward-oriented sphere `|q|=r`, `0<r<=1`, the homotopy
`A_b((1-t) sin(q)+t q)`, `0<=t<=1`, never vanishes: each nonzero
coordinate keeps the sign of the corresponding coordinate of `A_b q`.
The normalized source symbol therefore has the same degree as the linear
map `A_b`. The bound `|sin q| >= (1-r^2/6)|q|` gives the explicit boundary
gap `c r(1-r^2/6)`. A continuous traceless Hermitian perturbation of norm
strictly below this bound preserves the degree by the straight homotopy.
The statement requires a nonzero symbol on that sphere. A scalar identity
term can move the energy of a degeneracy and does not preserve a zero-energy
particle interpretation merely by preserving the eigenprojectors.

For `P_s=(I+s n.sigma)/2`, `s=+/-1`, fix the convention

\[
 C(P_s)=\frac{1}{2\pi i}\int\operatorname{Tr}(P_s\,dP_s\wedge dP_s),
 \qquad A=-i u^\dagger du.
\]

The Pauli identity gives
`Tr(P_s [d_theta P_s,d_phi P_s])=s i sin(theta)/2` for the unit-degree
sphere. Hence `C(P_s)=s chi_b`. Both signs are checked directly by a
gauge-invariant triangulated eigenvector flux for the nonlinear source
symbol, separately from the exact projector identity.

| Source node, in units of pi | Kinetic degree | Positive-band Chern number | Free Weyl chirality |
| --- | ---: | ---: | ---: |
| (0,0,0) | -1 | -1 | -1 |
| (0,0,1) | +1 | +1 | +1 |
| (0,1,0) | +1 | +1 | +1 |
| (0,1,1) | -1 | -1 | -1 |
| (1,0,0) | +1 | +1 | +1 |
| (1,0,1) | -1 | -1 | -1 |
| (1,1,0) | -1 | -1 | -1 |
| (1,1,1) | +1 | +1 | +1 |

The total degree is zero. These are momentum-space band defects, not an
exhaustive catalogue of real-space defects or branch-commuting repair charges.

## Observer-field meaning on the proved free family

The exact special-unitary source frame in
[the geometric Weyl construction](../weyl_geometric_modular/README.md)
conjugates `sum_j (A_b)_jj sigma_j p_j` to `chi_b sigma.p`.
With `k=k_b+a p`, `c=v/a`, the kinetic error is at most
`v a^2 |p|^3/6`. The free-limit theorem controls the interaction and the full
vacuum hierarchy in the stated family; the geometric construction supplies
the compatible Lorentz spin representation and common stress tensor.
Positive-energy particles have `E=v|p|` and helicity `chi_b/2` in units
`hbar=1`; holes of the filled negative-energy band have opposite charge and
helicity `-chi_b/2`. Their physical momentum is the negative of the removed
mode's momentum. A hole is the antiparticle of that Weyl field, not a ninth
independent field. Four fields of each chirality can be paired into four
massless Dirac fields. No color, weak representation or measured particle
name follows from this pairing.

## Source mass channels and the limit of the protection

The exhaustive onsite anticommuting family in
[the native mass construction](../native_mass_channels/README.md) is

\[
 B(x)=m_0(-1)^{x_1+x_2+x_3}I+
       \sum_j m_j(-1)^{x_j}\sigma_j.
\]

In the complete sixteen-row spin/valley space, the first term exchanges
`b` with `b xor (1,1,1)`, and term `j` exchanges it with `b xor e_j`.
Each exchange joins opposite winding. If `Gamma_chi` has blocks `chi_b I`,
then every kinetic matrix commutes with `Gamma_chi` and every mass matrix
anticommutes with it. All four mass matrices also anticommute with all three
kinetic matrices. These exact identities are independently joined to the
literal onsite Fourier transform. The terms break the indicated one-step
translations while preserving two-step translations, so their gapping does
not contradict the isolated-node winding protection.

The [source mass-shell proof](../source_mass_shell/README.md) gives masses
`|m0 +/- sqrt(m1^2+m2^2+m3^2)|` when these coefficients are supplied.
It does not select the coefficients or prove spontaneous order in the full
interacting source. In particular, a kinetic winding alone does not prove
a positive-residue interacting pole, a Green-function winding, a physical
mass or protection under every native repair branch. Those require the
actual state's spectral response and specified repair dynamics.
