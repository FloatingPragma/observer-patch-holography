"""One source-bound packet for neutral transport and coupled sterile controls."""
import hashlib
import json
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from . import comparison, neutral, self_energy, sterile, sterile_records, transport

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def bindings():
    seeds = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.update({self_energy.PARENT, comparison.INPUT,
        'model/defect_mass_join/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/neutral_composite/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/ground_composite_response/PROOF.md'})
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    pins = bindings()
    result = {'schema': 'oph.neutrino_mass_response.v1', 'source_bindings': pins,
        'native_neutral_and_declared_composite_Majorana_operators': neutral.certificate(),
        'complete_few_site_neutral_transport': transport.calculate(),
        'complete_native_neutral_normal_selfenergy': self_energy.calculate(),
        'coupled_full32_sterile_source': sterile.operator_certificate(),
        'complete_unit_Dirac_lepton_sterile_control': sterile_records.calculate(False),
        'complete_unit_sterile_Majorana_lepton_control': sterile_records.calculate(True),
        'conditional_one_scale_GeV_comparison': comparison.calculate(),
        'scope': {'native_source_definitions_use_empirical_neutrino_inputs': False,
            'same_ground_active_and_sterile_addition_and_removal_measures_executed': True,
            'all_physical_Gauss_sectors_of_each1280control_compared': True,
            'full32_all_unit_ground_or_spatial_neutral_poles_computed': False,
            'tau_antineutrino_is_a_sterile_species': False,
            'three_observed_generations_or_mixing_matrix_derived': False,
            'extra_cyclic_or_sterile_modes_discarded': False,
            'neutrino_rest_mass_or_Lorentz_limit_proved': False,
            'source_fixes_a_GeV_scale': False}}
    if bindings() != pins:
        raise ValueError('the neutral source changed during calculation')
    return result


def write():
    value = certificate()
    (HERE/'results.json').write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    return value
