"""Use the existing edge theorem with a certified source tail premise.

The producer's complete finite receipt is checked before it is consumed.
No original weight is replaced by a conditional success probability.
"""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from model.spin_shell_readout.threshold import bracket,non_atomic_control

ROOT=Path(__file__).resolve().parents[4]
SOURCE=ROOT/'model/coherent_defect_wave/chiral_joint/two_cell/results.json'


def source_receipt():
    result=json.loads(SOURCE.read_text())
    if result['schema']!='oph.coherent_defect_wave.chiral_joint.two_cell.v1':
        raise ValueError('wrong physical source receipt')
    for name,digest in result['input_sha256'].items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError('stale physical source input: '+name)
    return result


def attach_branch(source,kind):
    if kind not in ('add','remove'):raise ValueError('charged success branch required')
    preparation=source['original_attempt_readout']['preparation']
    p0=F(preparation['original_ground_success_probability'])
    fail=F(preparation['original_preparation_failure_probability'])
    if p0+fail!=1:raise ValueError('incomplete source preparation')
    branch=next(row for row in source['conditional_readout']['records'] if row['kind']==kind)
    if not branch['unreachable_conserved_blocks_have_exact_zero_effect_on_prepared_ground']:
        raise ValueError('uncontrolled lower conserved-sector response')
    windows=branch['complete_energy_windows']
    first=windows[0]
    if first['spectral_rank']!=1:
        raise ValueError('aggregate first-window weight is not a proved first-atom weight')
    lo,hi=map(F,first['energy_transfer_interval'])
    lower,upper=(p0*F(q) for q in first['probability_interval'])
    if not 0<lower<=upper<=1:raise ValueError('first atom has no certified original weight')
    if any(F(row['energy_transfer_interval'][0])<=hi for row in windows[1:]):
        raise ValueError('first atom is not separated from all other reachable energies')
    # The actual first atom gives nu([theta,theta+e])>=lower*e for all0<e<=1.
    # All lower source energies are excluded by the complete spectrum and
    # conserved-sector selection, not inferred from missing detector counts.
    epsilon=F(1,1000000)
    conclusion=bracket(lo-epsilon,(F(0),F(0)),hi,(lower,F(1)),epsilon,lower,1)
    return {'source_branch':kind,'reference':'the same unique two-cell physical ground',
            'source_energy_clock':source['conditional_readout']['source_energy_clock'],
            'first_atom_source_energy_interval':list(map(str,(lo,hi))),
            'first_atom_original_probability_interval':list(map(str,(lower,upper))),
            'source_tail_premise':{'amplitude':str(lower),'power':1,'epsilon_max':'1',
                'reason':'complete source spectrum plus certified positive first-atom weight'},
            'existing_edge_theorem_conclusion':{k:str(v) if isinstance(v,F) else v for k,v in conclusion.items()},
            'all_other_attempts_retained_as_complement':True,
            'complement_probability_interval':[str(1-upper),str(1-lower)],
            'finite_threshold_is_relativistic_mass':False,
            'continuum_uniform_tail_modulus_proved':False}


def certificate_data():
    source=source_receipt()
    rows=[attach_branch(source,kind) for kind in ('add','remove')]
    controls=[non_atomic_control(k) for k in (3,6,9)]
    assert all(row['atom_weight']==0 for row in controls)
    assert all(row['bracket']['threshold_lower_strict']<row['threshold']<=row['bracket']['threshold_upper']
               for row in controls)
    return {'actual_finite_source_edges':rows,
            'existing_nonatomic_control_levels':[row['kind'] for row in controls],
            'nonatomic_example_is_OPH_source_output':False,
            'one_isolated_pole_required_for_every_physical_species':False,
            'continuum_em_charge_or_GeV_mass_identified':False,
            'source_receipt_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest()}
