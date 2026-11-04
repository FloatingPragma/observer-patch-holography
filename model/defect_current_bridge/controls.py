"""Complete finite source checks and an independent complex-kernel control."""
import numpy as np
from scipy import sparse as sp

from model.mcpr_native_memory_bridge.construction import design
from model.mcpr_native_memory_bridge.certificate import mass_norm_bound
from . import algebra, construction as c


def _bounded(errors, bound=1e-9):
    if any(not np.isfinite(value) or value >= bound for value in errors.values()):
        raise ArithmeticError('full finite operator residual exceeded its diagnostic bound')
    return {'absolute_Frobenius_residual_upper': str(bound),
            'checked_identities': sorted(errors), 'every_check_passed': True}


def source_case(branch='forward', inverse_branch=0):
    source = design(branch, inverse_branch)
    mass_bound = mass_norm_bound(branch)
    kernel = c.native_kernel(source['mass_matrix_over_v'])
    value = c.assemble(kernel)
    errors = c.residuals(value)
    action_errors, native_errors = c.current_errors(value)
    errors['all12_actual_action_currents'] = max(action_errors)
    errors['all12_native_port_currents'] = max(native_errors)
    sv = np.linalg.svd(kernel.toarray(), compute_uv=False)
    energies = np.linalg.eigvalsh(value['core'].toarray())
    errors['complete_core_spectrum'] = float(np.max(np.abs(np.sort(energies)-np.sort(np.r_[np.repeat(sv, 2), -np.repeat(sv, 2)]))))
    lepton = np.sort(np.linalg.eigvalsh(source['mass_matrix_over_v']))
    errors['memory_MCPR_same_eigenvalues'] = float(np.max(np.abs(lepton-np.sort(np.array(source['masses_over_v'], float)))))
    if sv.max() >= 1 or len(sv) != 45 or np.sum(sv < 1e-12) != 3:
        raise ArithmeticError('complete source kernel and rational norm admission changed')
    if np.sum(energies < -1e-12) != 84 or np.sum(energies > 1e-12) != 84 or np.sum(abs(energies) < 1e-12) != 12:
        raise ArithmeticError('full core Nambu multiplicities changed')
    density = value['source']['density']
    points = np.array(value['source']['points'])
    errors['position_mean'] = float(np.linalg.norm(density@points))
    return {'P_branch': branch, 'inverse_memory_branch': inverse_branch,
            'operator_diagnostics': _bounded(errors),
            'all_current_count': len(action_errors), 'all_native_port_count': len(native_errors),
            'radial_rows': 54, 'internal_rows': 45, 'physical_annihilator_modes': 2430,
            'all_Nambu_rows': 4860, 'defect_annihilators': 90, 'defect_Nambu_rows': 180,
            'complement_annihilators': 2340, 'kernel_rank': 42,
            'core_negative_Nambu_rows': 84, 'core_zero_Nambu_rows': 12, 'core_positive_Nambu_rows': 84,
            'single_quasiparticle_bounds': algebra.spectral_bound(1),
            'exact_source_mass_bound': mass_bound,
            'energy_unit': 'declared v(P), no measured clock calibration',
            'lepton_eigenvalues': 'exactly the three retained MCPR m_i(P)/v(P) in the analytic construction',
            'up_eigenvalues': ['6^-6/sqrt(2)', '6^-3/sqrt(2)', '1/sqrt(2)'],
            'down_eigenvalues': ['6^-6/sqrt(2)', '6^-4/sqrt(2)', '6^-2/sqrt(2)'],
            'neutrino_eigenvalues': ['0', '0', '0'],
            'all_massless_core_directions_retained': True}


def independent_control():
    """Direct full216-row diagonalization with a complex symmetric kernel."""
    k = np.array([[1/7, (2+3j)/11], [(2+3j)/11, (-1+2j)/13]])
    value = c.assemble(k)
    full = np.linalg.eigvalsh(value['bdg'].toarray())
    radial = np.linalg.eigvalsh(value['h'].toarray())
    singular = np.linalg.svd(k, compute_uv=False)
    expected = np.sort(np.array([energy+sign*mass for energy in radial for mass in singular for sign in (-1, 1)]))
    errors = c.residuals(value)
    errors['full216_eigenvalue_enumeration'] = float(np.max(abs(full-expected)))
    bad = c.assemble(k, staggered=False)
    e, h = bad['inclusion'], bad['bdg']
    leakage = c.norm(h@e-e@(e.conj().T@h@e))
    if leakage <= 1e-3:
        raise ArithmeticError('unstaggered pairing must fail core preservation')
    # Local site dephasing fails preservation of the source zero-mode subspace.
    _, w, _, source = c.radial()
    origin = source['P_rows']
    projector = sp.csr_matrix((np.ones(2), (origin, origin)), shape=(54, 54))
    repair_leakage = c.norm(projector@w-w@(w.conj().T@projector@w))
    if repair_leakage <= 1e-3:
        raise ArithmeticError('local site resolution must leak out of the core')
    return {'complex_symmetric_control': _bounded(errors),
            'independent_full_Nambu_rows': 216,
            'unstaggered_pairing_leakage_lower': '0.001',
            'origin_site_resolution_leakage_lower': '0.001',
            'finite_diagonalization_is_an_exact_interval_certificate': False}
