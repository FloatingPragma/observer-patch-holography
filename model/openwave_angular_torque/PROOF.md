# Angular momentum of a finite localized matter and Maxwell construction

Let P_N project the three-dimensional harmonic-oscillator functions onto all
states with total excitation at most N, with N at least two. The projection
commutes with ordinary spatial rotations. Its coordinate wavefunctions are
Gaussian times polynomials, so every retained matter mode is localized.
Define p_i=P_N(-i partial_i)P_N and
L_i=P_N(-i epsilon_ijk x_j partial_k)P_N. The latter equals
-i epsilon_ijk a_j^* a_k within complete shells. It is not the product of
separately truncated coordinates and momenta. The exact relations are
[L_i,L_j]=i epsilon_ijk L_k and [L_i,p_j]=i epsilon_ijk p_k.

Use the exact native sixteen-row Clifford representation. The source unitary
in the stationary magnetic construction decomposes it into four complete
Dirac blocks. All four blocks and both spectral signs are retained. In each
block J=L+Sigma/2 generates actual rotations of the coordinate wavefunction
and its native spinor. This construction supplies an exact continuous
rotation action because the regulator consists of complete angular shells.
A finite cubic site set has no such action by position permutations: a
continuous map from the connected rotation group to a finite permutation
group is constant. Exact cubic covariance alone therefore does not determine
a continuous orbital generator.

The finite transverse field is

    A_b(x) = exp(-|x|^2/(2R^2)) (b cross x)/2,

with all three real amplitudes b and all three canonical momenta pi.
It is divergence free and square integrable together with its curl. Direct
Gaussian integration gives

    integral |A_b|^2 = c_0 |b|^2,  c_0=pi^(3/2)R^5/4,
    integral |curl A_b|^2 = c_1 |b|^2, c_1=5 pi^(3/2)R^3/8.

Thus the retained Maxwell Hamiltonian is
H_f=kappa |pi|^2/(2c_0)+c_1 |b|^2/(2kappa). Its rotational Noether charge is
b cross pi. This follows from the actual Maxwell mode metric and the vector
transformation of b; no spin label is assigned to the field.

The matter source is

    D(b)=beta m+alpha_i P_N(p_i-q A_b,i)P_N,
    K(b)=|D(b)|+Omega N,
    E(psi,b,pi)=<psi,K(b)psi>+H_f.

The oscillator energy trap is a declared primitive. It is applied to the
positive kernel, and is distinct from the scalar mass well of the cubic
source. The field regulator is also distinct from its periodic cochain
field. The same native Clifford algebra is used, without a claim that these
changed source Hamiltonians share their spectra, ground states or clock.
The anticommutation of beta with the kinetic matrix implies |D|>=m.

For a signed spectral resolution D u_a=lambda_a u_a, the actual derivative
of |D| is

    (delta |D|)_ab = (lambda_a+lambda_b)
                    delta D_ab/(|lambda_a|+|lambda_b|).

This includes degenerate eigenvalues and opposite spectral signs. Define
g_j=<psi,partial_bj K psi> and the retained current j=-g. A derivative of
the energy therefore supplies the current, including every retained mode.
Rotational covariance U(R)K(b)U(R)^*=K(Rb) implies

    d< J >/dt = b cross g,
    bdot=kappa pi/c_0, pidot=-c_1 b/kappa-g,
    d(b cross pi)/dt = -b cross g.

Consequently the complete angular momentum <J>+b cross pi and total energy
are conserved exactly by this finite semiclassical matter/field law.
The separate canonical spin can exchange angular momentum with the spatial
wavefunction and the field; it is not the conserved generator.

The executable checks the full matrix Lie identities, a noncubic rotation,
current derivatives on states with both Dirac spectral signs, and the
instantaneous torque identity. It evolves every matter coordinate and all
six field phase-space coordinates at two tolerances. The reported drifts
are numerical integration diagnostics, separate from the exact conservation
theorem. Complete spectral J_z probabilities are supplied on independent
preparations at each observation time. Removing field reaction, or replacing
L by the cross product of truncated coordinates, fails the stated controls.

For prescribed fields the matter law has its full number-conserving Fock
extension with the empty vacuum and spectator blocks retained. Its coupled
state-dependent Maxwell law is semiclassical mean field, not a linear
completely positive map on arbitrary entangled inputs. The construction
does not include longitudinal electrostatic backreaction, a radiative
continuum, or a physical electron and unit calibration.
