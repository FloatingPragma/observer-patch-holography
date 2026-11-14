"""Independent complete-bin calibration audit; no source or producer rerun."""
from pathlib import Path
import json,hashlib,math
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
def run():
    paths=[Path(__file__),HERE/'charged_energy_units.py',HERE/'CHARGED_ENERGY_UNITS.json',BASE/'next/C10/JOINT_SOURCE.json',BASE/'next/C10/CANONICAL_SOURCE_CONTRACT.json']
    pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    saved=json.loads((HERE/'CHARGED_ENERGY_UNITS.json').read_text());joint=json.loads((BASE/'next/C10/JOINT_SOURCE.json').read_text());reference=joint['references']['global_allN'];ground=reference['reference_energy'];responses=reference['complete_reservoir_response']
    addition=next(x for x in responses if x['source_charge_change']==1);removal=next(x for x in responses if x['source_charge_change']==-1)
    peak=lambda rows:sorted(rows,key=lambda x:(-x['probability'],x['energy']))[0]
    add=peak(addition['complete_source_number_energy_records']);remove=peak(removal['complete_source_number_energy_records']);anchor=add['energy']-ground;scale=.000511/anchor
    count=0;error=0.;work=0.;prob=0.
    for raw,display in zip(responses,saved['complete_response']):
        if len(raw['complete_source_number_energy_records'])!=len(display['all_number_energy_bins']):raise ValueError('bin omitted')
        for original,shown in zip(raw['complete_source_number_energy_records'],display['all_number_energy_bins']):
            for k,v in original.items():
                if shown[k]!=v:raise ValueError('original bin metadata/probability changed')
            excitation=original['energy']-ground
            error=max(error,abs(excitation-shown['source_excitation']),abs(scale*excitation-shown['calibrated_excitation_GeV']))
            # An arbitrary additive source energy convention cancels exactly.
            error=max(error,abs(((original['energy']+7)-(ground+7))-excitation))
            count+=1
        error=max(error,abs(scale*raw['unnormalized_apparatus_work']-display['calibrated_unnormalized_apparatus_work_GeV']))
        work+=display['calibrated_unnormalized_apparatus_work_GeV'];prob+=display['original_unconditional_probability']
    ratio=(remove['energy']-ground)/anchor
    if count!=1392 or error>1e-12 or abs(prob-1)>1e-12 or abs(ratio-saved['scale_free_removal_to_addition_peak_ratio'])>1e-12:raise ValueError('complete calibration audit failed')
    if saved['physical_particle_masses_predicted']!=[] or saved['calibration']['anchor_agreement_is_prediction'] is not False or saved['full_criterion_met'] is not False:raise ValueError('calibration promoted')
    if any(hashlib.sha256((BASE/p).read_bytes()).hexdigest()!=h for p,h in pins.items()):raise ValueError('source changed')
    return {'schema':'oph.c11.charged_units.independent_review.v1','reviewer':'Independent Codex worker, not opposite-team Fable audit','outcome':'PASS_AS_EXPLICIT_DISPLAY_CALIBRATION','complete_bins_checked':count,'original_response_probability':prob,'maximum_bin_or_work_error':error,'source_addition_anchor':anchor,'GeV_per_source_unit':scale,'scale_free_removal_to_addition_peak_ratio':ratio,'calibrated_mean_work_GeV':work,'scope':'One supplied 0.000511 GeV anchor; all 1392 bins/zeros and original unconditional probabilities retained. Source energy shifts cancel in excitation. No particle/electron identity, physical mass prediction or charge-conjugation symmetry inferred.','source_sha256':pins,'full_criterion_met':False,'official_rows_changed':[]}
if __name__=='__main__':
    r=run();Path(__file__).with_name('CHARGED_ENERGY_UNITS_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
