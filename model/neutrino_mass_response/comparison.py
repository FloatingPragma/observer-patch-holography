"""An explicit empirical unit example, separate from the native calculation."""
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
INPUT = 'model/common_source/calibrated_neutrino_parameters.json'


def calculate():
    packet = json.loads((ROOT/INPUT).read_text())
    split = s.Rational(packet['oscillation_inputs']['delta_m31_squared_eV2'])
    lightest = s.Rational(packet['unmeasured_architecture_choices']['lightest_mass_eV'])
    scale = s.sqrt(lightest**2+split)/10**9
    assert scale.is_positive
    return {'general_unit_law': 'every source energy gap Delta maps to s_GeV Delta for one positive common s_GeV',
        'source_fixes_s_GeV': False,
        'physical_rest_pole_attachment_proved': False,
        'two_positive_unit_Dirac_control_gaps': ['1', '5'],
        'conditional_rest_masses_if_the_attachment_is_supplied': ['s_GeV', '5 s_GeV'],
        'comparison_calibration': {'role': 'empirical illustrative choice of one common source energy unit',
            'input_file': INPUT, 'input_sha256': hashlib.sha256((ROOT/INPUT).read_bytes()).hexdigest(),
            'declared_external_source': packet['source'],
            'consumed_mass_squared_splitting_eV2': str(split),
            'consumed_lightest_mass_architecture_eV': str(lightest),
            'unit_choice': 's_GeV=sqrt(m_lightest^2+Delta_m31^2)/10^9',
            's_GeV_exact': str(scale),
            'conditional_GeV_values_exact': [str(scale), str(5*scale)],
            'conditional_GeV_values_numeric': [float(scale), float(5*scale)],
            'active_or_sterile_species_assigned_to_either_value': False,
            'independent_mass_prediction': False,
            'normal_ordering_and_zero_lightest_mass_are_source_predictions': False,
            'oscillation_fit_reproduced': False}}
