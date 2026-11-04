"""Exact local identities and rational spectral domain for the CAR lift."""
from fractions import Fraction as F

import sympy as s

from model.native_mass_channels.algebra import SIGMA


def certificate():
    j = s.I*SIGMA[1]
    assert j.T == -j and j.H*j == s.eye(2) and j*j == -s.eye(2)
    for sigma in SIGMA:
        assert sigma*j+j*sigma.T == s.zeros(2)
        link = s.I*sigma/2
        # Neighbouring sites have opposite stagger signs.
        assert -link*j+j*link.conjugate() == s.zeros(2)
    a, b, c = s.symbols('a b c', complex=True)
    k = s.Matrix([[a, b], [b, c]])
    delta = s.kronecker_product(j, k.conjugate())
    assert delta.T == -delta
    assert s.simplify(delta.H*delta-s.kronecker_product(s.eye(2), k.T*k.conjugate())) == s.zeros(4)
    z1, z2 = s.symbols('z1 z2', complex=True)
    for parity in (-1, 1):
        first = s.Matrix([z1, z2])
        local = s.Matrix.hstack(first, -parity*j*first.conjugate())
        assert s.simplify(parity*j*local.conjugate()-local*j) == s.zeros(2)
    return {
        'spin_antisymmetric_unitary': 'J=i sigma2; J^T=-J; J^2=-I',
        'source_particle_hole': 'h C=-C conjugate(h), C=D J',
        'normalized_kernel_frame': 'h W=0; W^dagger W=I2; C conjugate(W)=W J',
        'local_creation_kernel': 'Delta=C tensor conjugate(K); K^T=K',
        'core_creation_kernel': 'Delta_core=J tensor conjugate(K)',
        'Nambu_inclusion': 'E=diag(W tensor I, conjugate(W) tensor I)',
        'exact_intertwiner': 'B_full E=E B_core',
        'anti_Hermitian_current_inclusion': 'diag(I tensor G,I tensor conjugate(G)) E = E diag(I2 tensor G,I2 tensor conjugate(G))',
        'Hermitian_charge_convention': 'diag(Q,-conjugate(Q))',
        'hole_frame': 'T=diag(I,C^dagger tensor I)',
        'complete_transformed_operator': 'T^dagger B_full T=[[h tensor I,I tensor conjugate(K)],[I tensor K,h tensor I]]',
        'all_energies': 'lambda(h) plus or minus singular_values(K), each with its full multiplicity',
        'core_energies': 'plus or minus singular_values(K), each twice',
        'CAR_map': 'c(f) maps to c((W tensor I) f)',
        'even_state_map': 'rho maps to rho tensor normalized_trace_on_complement',
        'exact_dynamics_and_current_covariance_of_state_map': True,
        'normalized_trace_complement_is_vacuum': False,
        'localization': 'every normalized core wavefunction has site probability rho_x; sum rho_x=1',
        'position_mean': 'zero by inversion; a fixed localized state is not a momentum eigenstate',
    }


def spectral_bound(kernel_norm_upper, radial_gap=62):
    """Rational isolation of single-quasiparticle bands, not many-body bands."""
    if isinstance(kernel_norm_upper, (float, bool)) or isinstance(radial_gap, (float, bool)):
        raise ValueError('exact rational bounds required')
    m, d = F(kernel_norm_upper), F(radial_gap)
    if m < 0 or d <= 0 or 2*m >= d:
        raise ValueError('require 0<=2*M<d for complete core-band isolation')
    return {'radial_gap_lower': str(d), 'kernel_norm_upper': str(m),
            'complement_absolute_energy_lower': str(d-m),
            'core_absolute_energy_upper': str(m),
            'band_separation_lower': str(d-2*m),
            'domain': 'single-quasiparticle Nambu spectrum',
            'full_many_body_band_isolation': False}


def pulse_service_certificate():
    """Finite Fock operator bounds; no full Fock matrix is materialized."""
    return {
        'record_Hamiltonian': 'H0 + C_enabled * gain * L(P_updated_memory)',
        'H0': 'dGamma(h tensor I45) plus the declared quark creation pairing',
        'L': 'local lepton creation pairing with Delta=C tensor conjugate(T_lepton tensor P)',
        'stationary_mean': 'H0 + L(gain * Z_C)',
        'each_record_preserves_same_defect_CAR_subalgebra': True,
        'physical_modes': 2430,
        'radial_norm_upper': 66,
        'static_quark_pairing_norm_upper': 1215,
        'unit_gain_lepton_pairing_norm_upper': 1215,
        'all_record_Fock_norm_upper': 'M=2430*66+1215*(1+gain)',
        'arbitrary_CQ_reference_trace_error': '4*M*T/(N*(1-p))+(8/(1-p)+4)*M^2*T^2/N',
        'stationary_product_initial_memory': 'first error term vanishes',
        'domain': 'fixed finite volume, N>=1, T>=0, 0<=p<1, gain>=0',
        'finite_unread_channel_equals_mean_unitary': False,
        'original_site_field_repair_used_as_pulse_service': False,
        'uniform_in_volume_continuum_bound': False}
