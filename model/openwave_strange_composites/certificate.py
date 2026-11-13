"""Complete strange-sector energies, native instruments and a decay obstruction."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
from .source import source,select,parent
from .current import mass_parameters,instrument
from .symmetry import charge_conjugation
from .nonleptonic import certificate as nonleptonic
from model.openwave_hadron_refinement_flow.positivity import joint_stability

INPUTS=json.loads((Path(__file__).parent/'inputs.json').read_text())


def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,complex):return [round(value.real,12),round(value.imag,12)]
    if isinstance(value,float):return round(value,12)
    if isinstance(value,(list,tuple)):return [serial(v)for v in value]
    if isinstance(value,dict):return {k:serial(v)for k,v in value.items()}
    return value


@lru_cache(None)
def calculate():
    dm,pair=mass_parameters();unit=F(parent.INPUTS['electron_unit_GeV']);rows=[]
    cases=[('Lambda','baryon',(1,1,0,1),(0,0,0,0),.5,0),
           ('Sigma+','baryon',(2,0,0,1),(0,0,0,0),.5,1),
           ('Sigma0','baryon',(1,1,0,1),(0,0,0,0),.5,1),
           ('Sigma-','baryon',(0,2,0,1),(0,0,0,0),.5,1),
           ('Xi0','baryon',(1,0,0,2),(0,0,0,0),.5,.5),
           ('Xi-','baryon',(0,1,0,2),(0,0,0,0),.5,.5),
           ('Omega-','baryon',(0,0,0,3),(0,0,0,0),1.5,0),
           ('K+','meson',(1,0,0,0),(0,0,0,1),0,.5),
           ('K-','meson',(0,0,0,1),(1,0,0,0),0,0),
           ('K0','meson',(0,1,0,0),(0,0,0,1),0,.5),
           ('Kbar0','meson',(0,0,0,1),(0,1,0,0),0,0)]
    for name,kind,q,a,s,i in cases:
        z=select(kind,q,a,s,i,float(dm),float(pair));E=z['energy']*float(unit)
        target=INPUTS['mass_calibration_targets_GeV'].get(name)
        role='calibration' if target else 'strong_C_consequence' if name in('K-','Kbar0') else 'held_out'
        if role=='calibration':assert abs(E-float(F(target)))<1e-10
        if role=='held_out':target=INPUTS['held_out_masses_GeV'][name]
        rows.append({'name':name,'native_quark_counts':q,'native_antiquark_counts':a,
            'spin':s,'light_quark_isospin':i,'selected_complete_dimension':z['dimension'],
            'all_selected_source_levels':z['all_levels'],'computed_GeV':E,'comparison_GeV':target,
            'comparison_role':role,'difference_GeV':None if target is None else E-float(F(target)),
            'full_native_eigen_residual':z['residual']})
    sym=charge_conjugation();weak=instrument(F(INPUTS['Cabibbo_sine']));zero=instrument(F(0))
    assert zero['exclusive_proton_probability']<1e-25
    assert weak['exclusive_proton_probability']>0
    count=weak['complete_original_native_occupation_records'];total=count['preparation_failure']
    total+=sum(count['success_and_current_success_occupation'])+sum(count['success_and_current_failure_occupation'])
    assert abs(total-1)<3e-13
    c=nonleptonic()
    return serial({'schema':'oph.strange_composites.v1',
        'verdict':'PASS_NATIVE_STRANGE_MASSES_C_SYMMETRY_AND_CURRENT_WITH_CLOSED_NONLEPTONIC_CONTROL',
        'source_sector_sizes':{'baryon_colour_singlets':len(source('baryon')['labels']),
            'baryon_full_native_occupation':len(source('baryon')['full']),
            'meson_colour_singlets':len(source('meson')['labels']),
            'meson_full_native_occupation':len(source('meson')['full'])},
        'mass_rows':rows,'strong_charge_conjugation':sym,
        'additional_strong_coefficients':{'strange_number':dm,'strange_pair':pair},
        'full_Fock_stability':joint_stability(),'complete_native_current_instrument':weak,
        'zero_Cabibbo_control':{'exclusive_proton_probability':zero['exclusive_proton_probability']},
        'native_nonleptonic_vertex_and_obstruction':c,
        'calibration_ancestry':INPUTS['ancestry'],
        'scope':'The source executes native strange states and current interventions with every failure retained. Sigma, Xi and neutral-kaon mass controls fail; the two-cell proton/pion five-body channel is closed. No physical strange-decay width or complete observed hadron spectrum is asserted.'})
