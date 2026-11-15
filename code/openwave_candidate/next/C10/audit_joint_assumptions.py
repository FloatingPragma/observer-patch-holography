"""Adversarial scope controls distinct from the joint numerical reconstruction."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(HERE))
from charged_reservoir import BoundarySource


def run():
    paths=[Path(__file__),HERE/'joint_source.py',HERE/'JOINT_SOURCE.json',HERE/'joint_source_states.npz',HERE/'JOINT_SOURCE_INDEPENDENT_REVIEW.json',HERE/'charged_reservoir.py',HERE/'modular_field_adapter.py']
    pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    data=json.loads((HERE/'JOINT_SOURCE.json').read_text());archive=np.load(HERE/'joint_source_states.npz');h=archive['hamiltonian'];m=BoundarySource(6);numbers=np.array([mask.bit_count() for mask,_ in m.states]);bare=m.hamiltonian
    bareground=min(float(np.linalg.eigvalsh(bare[np.ix_(numbers==q,numbers==q)])[0]) for q in range(7))
    actual=data['references']['global_allN']['reference_energy'];reference_shift=actual-bareground
    ids=np.flatnonzero(numbers==3);e,v=np.linalg.eigh(h[np.ix_(ids,ids)]);pair=next((i,i+1) for i in range(len(e)-1) if abs(e[i+1]-e[i])<1e-9)
    columns=v[:,list(pair)];psi=columns@np.ones(2)/np.sqrt(2);rho=np.outer(psi,psi)
    p=columns@columns.T;luders=p@rho@p;fine=sum((np.outer(x,x)@rho@np.outer(x,x) for x in columns.T),start=np.zeros_like(rho))
    disturbance=float(np.linalg.norm(luders-fine,ord='nuc')/2);interference=float(psi@luders@psi-psi@fine@psi)
    rows=data['references']['global_allN']['complete_reservoir_response'];boundary_sensitivity=0.;initial_number_prob={}
    for row in rows:
        boundary_sensitivity+=sum((3-x['source_number'])**2*x['probability'] for x in row['complete_source_number_energy_records'])
    field=h-bare;field_min=float(np.linalg.eigvalsh(field)[0])
    checks={'zero_field_reference_not_transferable':reference_shift>1,
       'declared_field_positive_semidefinite':field_min>-1e-12,
       'terminal_POVM_does_not_fix_continuing_instrument':abs(disturbance-.5)<1e-12 and abs(interference-.5)<1e-12,
       'boundary_zero_energy_choice_observable_in_work':abs(boundary_sensitivity-.5)<1e-12,
       'supplied_pulse_not_closed_energy_conserving':data['references']['global_allN']['mean_reservoir_apparatus_work']>1,
       'fixed_background_conjugation_failure_retained':data['conjugation_image_states_in_same_background']==0}
    if not all(checks.values()):raise ValueError(checks)
    if any(hashlib.sha256((BASE/p).read_bytes()).hexdigest()!=value for p,value in pins.items()):raise ValueError('sourcechanged')
    return {'schema':'oph.c10.joint_assumptions_adversarial_review.v1','checks':checks,
      'zero_field_ground':bareground,'interacting_ground':actual,'wrong_reference_energy_shift':reference_shift,
      'actual_degenerate_energy_group':{'source_number':3,'indices':list(pair),'energies':[float(e[i]) for i in pair],'Luders_vs_fine_readout_trace_distance':disturbance,'subsequent_interference_probability_difference':interference},
      'boundary_energy_sensitivity':'Adding supplied H_B=eta*B^2 changes mean pulsework by eta/2 on this reference; declaredmodelsets eta=0.',
      'mean_boundary_squared_change':boundary_sensitivity,
      'claims':[
        'All probability and energy readouts are complete terminal effects; archived states precede the terminal energy PVM. Continuing degenerate measurement needs explicit instrument choice.',
        'Source+reservoir charge changes cancel perbranch; fairreservoir settings mixinitial totalcharge sectors, not a neutralvacuum prep.',
        'ZeroextraBenergy, reservoirlevel.23 and idealpulse are apparatuschoices. Pulse duration, rate andphysicalclock are not derived.',
        'Periodic half-record energy differs from integerhardwallE2b; no certified regulator/continuum transfer follows.',
        'Fullfieldmixescharacter/uniformcodes, so only full384state interactingreferences apply.',
        'FixedG=(111000) background is not chargeconjugationinvariant; conjugated background is required.',
        'Physicalspecies/spin/electricunit/radiation/Maxwell3D/calibratedmass/vacuum are not identified.'],
      'reviewer_identity':'Same Codex worker as numericreconstruction, distinctadversarialmethod; not a second independent reviewer or opposite-team pass.',
      'source_sha256':pins,'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('JOINT_ASSUMPTIONS_ADVERSARIAL_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
