"""An explicit independent energy/length attachment, including its failures."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
from model.openwave_universal_energy_geometry.source import profile
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def comparison():
    inputs=json.loads((HERE/'inputs.json').read_text())
    old=json.loads((ROOT/'model/openwave_string_positive_cone/results.json').read_text())['common_length_attachment']
    a=F(old['a_fm']);u=F(old['electron_unit_GeV']);hc=F(inputs['hbar_c_GeV_fm'])
    g=1/(4*np.pi*float(profile()['kappa_exact']))
    value=g*float(a/(u*hc));target=float(F(inputs['Newton_G_natural_GeV_minus2']))
    speed=float(u*a/hc)
    assert value/target>1e39 and speed<.01
    return {'declared_common_unit_map':'one geometry energy unit equals the conditional electron unit; one geometry edge equals the independently strong-calibrated length',
            'energy_unit_GeV_exact':str(u),'length_unit_fm_exact':str(a),
            'length_input_ancestry':'The existing strong length is attached using the deuteron quadrupole and an illustrative string-tension target. The common hyperfine source has seven empirical strong inputs; this dimensional test transfers no Hamiltonian equivalence.',
            'exact_strong_mass_parent':'model/openwave_common_hyperfine_spectrum/results.json',
            'separate_dimensional_targets':['deuteron quadrupole','illustrative string tension1GeV/fm'],
            'source_stiffness_or_q_tuned_to_G':False,
            'dimensionless_Newton_response':g,
            'conditional_G_natural_GeV_minus2':value,
            'external_G_natural_GeV_minus2':target,
            'conditional_to_observed_strength_ratio':value/target,
            'conditional_photon_speed_over_observed_c':speed,
            'formula':'G_nat = g_dimensionless * a_fm/(u_GeV*hbar_c_GeV_fm); c_model/c = u_GeV*a_fm/hbar_c_GeV_fm',
            'physical_G_strength_match':False,'physical_common_light_speed_match':False,
            'common_strong_geometry_Hamiltonian_transfer_proved':False,
            'external_reference_inputs':inputs}
