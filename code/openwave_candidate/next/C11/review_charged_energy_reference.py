"""Independent exact-linear pushforward of archived interacting charged records."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,math
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]


def execute():
    source_path=BASE/'next/C10/JOINT_SOURCE.json';display_path=HERE/'CHARGED_ENERGY_UNITS.json'
    source=json.loads(source_path.read_text());display=json.loads(display_path.read_text())
    for doc in (source,display):
        for relative,pin in doc['source_sha256'].items():
            if hashlib.sha256((BASE/relative).read_bytes()).hexdigest()!=pin:raise ValueError('custody drift: '+relative)
    reference=source['references']['global_allN'];ground=reference['reference_energy'];branches=reference['complete_reservoir_response']
    addition=[b for b in branches if b['source_charge_change']==1][0]
    order=sorted(addition['complete_source_number_energy_records'],key=lambda r:(-r['probability'],r['energy']))
    selected=order[0];anchor_source=selected['energy']-ground;anchor=display['calibration']['supplied_anchor_GeV']
    exact=lambda x:Fraction.from_float(float(x))
    scale=exact(anchor)/exact(anchor_source);count=0;zeros=0;maximum=0.;work_source=Fraction(0);work_display=Fraction(0)
    expected_peaks=[]
    for original,converted in zip(branches,display['complete_response']):
        assert original['reservoir_initial']==converted['reservoir_initial'] and original['reservoir_final']==converted['reservoir_final']
        assert original['probability']==converted['original_unconditional_probability']
        raw=original['complete_source_number_energy_records'];bins=converted['all_number_energy_bins'];assert len(raw)==len(bins)
        for row,unit in zip(raw,bins):
            for key in ('source_number','energy','degeneracy','probability'):assert row[key]==unit[key]
            excitation=row['energy']-ground
            maximum=max(maximum,abs(float(scale*exact(excitation))-unit['calibrated_excitation_GeV']),abs(excitation-unit['source_excitation']))
            count+=1;zeros+=int(row['probability']==0)
        work_source+=exact(original['unnormalized_apparatus_work']);work_display+=exact(converted['calibrated_unnormalized_apparatus_work_GeV'])
        peak=sorted(raw,key=lambda r:(-r['probability'],r['energy']))[0]
        expected_peaks.append((original['source_charge_change'],peak['energy']-ground,peak['probability']))
    assert count==1392 and maximum<1e-15
    assert abs(float(scale*work_source-work_display))<1e-15
    add=next(x for x in expected_peaks if x[0]==1);remove=next(x for x in expected_peaks if x[0]==-1)
    ratio=exact(remove[1])/exact(add[1]);assert abs(float(ratio)-display['scale_free_removal_to_addition_peak_ratio'])<1e-15
    assert display['calibration']['empirical_anchors_used']==1 and not display['calibration']['anchor_agreement_is_prediction']
    assert display['physical_particle_masses_predicted']==[] and anchor==.000511
    # Exact invariance of the reference-subtracted energy under a common zero shift.
    shift=Fraction(17,3);shift_invariant=all((exact(row['energy'])+shift)-(exact(ground)+shift)==exact(row['energy'])-exact(ground) for b in branches for row in b['complete_source_number_energy_records'])
    paths=[Path(__file__),HERE/'charged_energy_units.py',display_path,source_path,BASE/'next/C10/CANONICAL_SOURCE_CONTRACT.json']
    return {'schema':'oph.c11.charged_energy_reference_review.v1','reviewer':'Second distinct Codex worker; no producer import/replay','outcome':'PASS_WITH_SINGLE_SUPPLIED_CALIBRATION_SCOPE','checks':{'all1392_unconditional_bins_preserved':count==1392,'no_probability_or_zero_filtering':True,'one_common_linear_scale':maximum<1e-15,'complete_work_mean_same_scale':abs(float(scale*work_source-work_display))<1e-15,'exact_stored_number_peak_ratio':True,'common_energy_zero_shift_cancels':shift_invariant,'anchor_not_empirical_prediction':True},'bins_checked':count,'exact_zero_probability_bins_retained':zeros,'maximum_energy_pushforward_error':maximum,'work_pushforward_error':abs(float(scale*work_source-work_display)),'exact_ratio_of_stored_binary_source_excitation_numbers':str(ratio),'ratio_float':float(ratio),'anchor':{'supplied_GeV':anchor,'selected_source_excitation':anchor_source,'unconditional_peak_probability':selected['probability'],'selection_independent_of_target_mass':'sort original addition bins by descending unconditional weight, then lower energy','exact_linear_scale_of_stored_numbers':str(scale)},'scope':'Exact rational arithmetic checks the linear map of stored finite binary64 records, not an exact analytic mass theorem. One511keV displayanchor is supplied; all other numbers are conditional source excitations under that same scale. No physicalspecies, restmomentum, antiparticle symmetry or empiricalparticle mass follows.','producer_runs':0,'large_diagonalizations':0,'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
if __name__=='__main__':
 r=execute();HERE.joinpath('CHARGED_ENERGY_UNITS_REFERENCE_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
