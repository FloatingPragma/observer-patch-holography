"""Reconstruct the concrete fixed-repair coefficient and its action boundary."""
import hashlib
import json
from pathlib import Path

import sympy as s

from tools import strict_json
from ..common_source import native_repair_mass as repair
from ..common_source import native_repair_mass_attachment as attachment
from ..common_source.chiral_candidate import action
from ..conditional_standard_model.spectrum import ActionInputs, _zero
from . import native, spectrum

ROOT = Path(__file__).resolve().parents[2]
RECEIPT = 'model/common_source/native_repair_mass_results.json'


def certificate():
    retained = strict_json.load(ROOT/RECEIPT)
    artifact = ROOT/'model/common_source'/retained['artifact']['name']
    if hashlib.sha256(artifact.read_bytes()).hexdigest() != retained['artifact']['sha256']:
        raise ValueError('the retained native detector artifact is stale')
    exact = retained['exact_kernel']
    u, p0, y = (s.Matrix([[s.sympify(x) for x in row] for row in exact[key]])
                 for key in ('native_U', 'native_seed_projector', 'Y'))
    records = repair.native_records()
    p = sum(row[-1] for row in records if row[-2])
    if len(records) != 128 or p != s.Rational(exact['p']) or p != s.Rational(1, 12):
        raise ArithmeticError('the concrete native event probability changed')
    if (not _zero(u.T*u-s.eye(3)) or not _zero(u**5-s.eye(3))
            or not _zero(p0*p0-p0) or s.simplify(s.trace(p0)) != 1):
        raise ArithmeticError('the retained native rotation/projector architecture changed')
    weights = s.Matrix([(1-p)*p**k/(1-p**5) for k in range(5)])
    transition = s.Matrix([[s.Rational(x) for x in row] for row in exact['transition']])
    if transition*weights != weights or weights != s.Matrix(list(map(s.Rational, exact['stationary_weights']))):
        raise ArithmeticError('the complete native-memory stationary weights changed')
    reconstructed = sum((weights[k]*u**k*p0*(u.T)**k for k in range(5)), s.zeros(3))
    if not _zero(reconstructed-y):
        raise ArithmeticError('the native stationary readout coefficient changed')
    spectral = attachment.spectral_certificate(y)
    data = ActionInputs(s.sqrt(2), 1, 1, 1, 1, y, y, y,
                        ('unchanged fixed-repair coefficient on all three Higgs channels',))
    zero = spectrum.MajoranaInputs((0, 0, 0), (s.eye(3),), 1, ('no Weinberg term in this source law',))
    old_kernel = spectrum.full_kernel(matter=data, neutrinos=zero)
    p45 = s.kronecker_product(native.dictionary()['permutation'], s.eye(3))
    native_kernel = s.ImmutableSparseMatrix(p45.T*old_kernel*p45)
    serial = [[int(i), int(j), str(s.expand(x))] for (i, j), x in sorted(native_kernel.todok().items())]
    digest = hashlib.sha256(json.dumps(serial, separators=(',', ':')).encode()).hexdigest()
    prior = retained['same_current_Higgs_and_heldout_test']
    target = prior['same_current_Higgs_attachment']
    if digest != target['neutral_K45_sparse_sha256']:
        raise ArithmeticError('the actual fixed-repair K45 differs from its retained action')
    neutrino_indices = action.ordering()['massless_neutrino_internal_family_indices']
    if any(any(old_kernel[i, j] != 0 for j in range(45)) for i in neutrino_indices):
        raise ArithmeticError('the fixed-repair source acquired a Majorana coefficient')
    intervals = spectral['ascending_ratios_after_one_exact']
    windows = prior['separate_heldout_comparison']['retained_broad_charged_lepton_ratio_windows']
    disjoint = []
    for interval, window in zip(intervals, windows, strict=True):
        lo, hi = map(s.Rational, interval)
        left, right = map(s.Rational, window)
        disjoint.append(bool(hi < left or right < lo))
    if disjoint != [True, True]:
        raise ArithmeticError('the retained fixed-law ratio exclusion changed')
    return {'readout_source': RECEIPT,
            'readout_producer': 'model/common_source/native_repair_mass_native_readout.py',
            'action_attachment': 'model/common_source/native_repair_mass_attachment.py',
            'raw_native_records_recomputed': len(records), 'native_advance_probability': str(p),
            'retained_rotation_and_projector_bound_by_parent_manifest': True,
            'stationary_Y_reconstructed_exactly': action._serialize(y),
            'stationary_weights': list(map(str, weights)),
            'complete_native_frame_K45_sha256': digest,
            'all_three_charged_Yukawa_matrices_equal_this_Y': True,
            'CKM_in_common_mass_basis': 'I3', 'Majorana_coefficient': 'zero3x3',
            'minimal_neutrino_masses': ['0', '0', '0'],
            'source_Y_positive_eigenvalue_intervals_exact': spectral['positive_eigenvalue_intervals_exact'],
            'charged_mass_ratios_after_lightest_exact': intervals,
            'retained_comparison_windows': windows, 'each_window_disjoint_exact': disjoint,
            'actual_detector_retained_fixed_probe_settings': retained['native_fixed_basis_readout']['complete_domains']['fixed_probe_settings'],
            'actual_detector_all_native_syndrome_and_pointer_outcomes_retained': retained['native_fixed_basis_readout']['execution']['all128_pointer_labels_and_native_syndrome_outputs_retained'],
            'expensive_detector_reexecuted_by_this_package': False,
            'retained_detector_artifact_sha256': retained['artifact']['sha256'],
            'reading_fixed_specs_determines_this_declared_model': True,
            'fixed_model_matches_the_retained_charged_ratio_windows': False,
            'source_repair_law_selects_the_separate_empirical_kappa': False,
            'all_other_interacting_source_laws_excluded': False}
