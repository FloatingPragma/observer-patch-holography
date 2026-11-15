"""Executed all-character Z6 energetic selection, retaining branch distinctions."""
from pathlib import Path
import json,hashlib
import numpy as np
from character_reference import blocks,FermionicSharedRecords,Construction,fable,BASE,FABLE


def run():
    m=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,6)
    masks,sectors=blocks(m);report=[];herr=kerr=0.;spectra=[]
    for chi,(theta,v,h,ops) in enumerate(sectors):
        levels=np.linalg.eigvalsh(h);spectra.extend(levels)
        single=np.linalg.eigvalsh(fable.Gauged(3,m.construction.seams,.23,phase=theta).hamiltonian())
        herr=max(herr,float(np.max(abs(m.hamiltonian@v-v@h))))
        for law in ops:
            for fullops,smallops in zip(m.operators[law],ops[law]):
                for (_,k),(_,l) in zip(fullops,smallops):kerr=max(kerr,float(np.max(abs(k@v-v@l))))
        gap=float(min(levels[levels>levels[0]+1e-10])-levels[0])
        report.append({'character':chi,'cycle_angle':float(2*np.pi*chi/6),'CAR_ground':float(levels[0]),'CAR_fixed_number_gap':gap,'CAR_spectrum':levels.tolist(),'one_load_lowest_signed_energy':float(single[0]),'one_load_positive_band_edge':float(min(single[single>0])),'negative_filled_sum_check':float(single[single<0].sum())})
    full,fullvectors=np.linalg.eigh(m.hamiltonian);ground=float(full[0]);multiplicity=int(sum(abs(full-ground)<1e-10))
    groundvector=fullvectors[:,0];energygroups=[]
    for q,value in enumerate(full):
        if energygroups and abs(value-energygroups[-1][0])<1e-9:energygroups[-1][1].append(q)
        else:energygroups.append([float(value),[q]])
    support={}
    for law in m.operators:
        joint=[];marginal=np.zeros(len(energygroups))
        for ops in m.operators[law]:
            for meta,k in ops:
                branch=k@groundvector/np.sqrt(3);probabilities=abs(fullvectors.conj().T@branch)**2
                distribution=[float(probabilities[ids].sum()) for _,ids in energygroups];marginal+=distribution
                joint.append({**meta,'event_probability':float(np.vdot(branch,branch).real),'all_energy_group_probabilities':distribution})
        supported=[value-ground for (value,_),prob in zip(energygroups,marginal) if value>ground+1e-9 and prob>1e-12]
        support[law]={'all_energy_groups':[{'energy':value,'excitation':value-ground,'degeneracy':len(ids),'probability':float(prob)} for (value,ids),prob in zip(energygroups,marginal)],'complete_joint_native_energy_records':joint,'first_numerically_supported_excitation':min(supported),'support_probability_threshold':1e-12,'mean_excitation_work':float(sum((value-ground)*prob for (value,_),prob in zip(energygroups,marginal)))}
    minima=[r['character'] for r in report if abs(r['CAR_ground']-ground)<1e-10]
    negmin=min(r['one_load_lowest_signed_energy'] for r in report);posmin=min(r['one_load_positive_band_edge'] for r in report)
    paths=[Path(__file__),Path(__file__).with_name('character_reference.py'),Path(__file__).with_name('fermionic_shared_records.py'),FABLE]
    return {'schema':'oph.c9.character_minima.v1','parameters':{'graph':'triangle','mass':.23,'modulus':6,'fixed_CAR_number':3,'full_dimension':m.n,'native_C_E_MIX':'same fulloccupation literalbalancing/directCARphase; no potential'},'checks':{'all_H_intertwiners':bool(herr<1e-12),'all_C_E_MIX_intertwiners':bool(kerr<1e-12),'full_spectrum_is_character_union':bool(np.max(abs(np.sort(spectra)-full))<1e-12),'all_CAR_sector_grounds_equal_negative_filled_sums':bool(all(abs(r['CAR_ground']-r['negative_filled_sum_check'])<1e-12 for r in report))},'characters':report,'global_ground':ground,'global_ground_multiplicity':multiplicity,'global_fixed_number_gap':float(min(full[full>ground+1e-10])-ground),'CAR_minimizing_characters':minima,'one_load_lowest_signed_energy_minimizers':[r['character'] for r in report if abs(r['one_load_lowest_signed_energy']-negmin)<1e-10],'one_load_positive_band_edge_minimizers':[r['character'] for r in report if abs(r['one_load_positive_band_edge']-posmin)<1e-10],'native_ground_energy_readout':support,'native_support_scope':'Complete fullH eigenspace PVM after nativeevent, supplied apparatus; allzero/degenerate groups retained. Characters algebraically invariant; numerical support withinχ uses declared1e-12 threshold. Spectralgap is not automatically nativeaccessible or physicalmass.','errors':{'H_intertwiner':herr,'native_intertwiner':kerr},'scope':'Energetic comparison of allcharacters in declared fixednumberCAR source, no assumedresult; minimumpositiveone-loadbranch is distinct from global signed one-loadground. Native sectorinvariance does not prepare/cool/select a character. No binding, physicalspin/species, physicalvacuum or CARderivation.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
if __name__=='__main__':
 r=run();Path(__file__).with_name('CHARACTER_MINIMA.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:r[k] for k in ['checks','CAR_minimizing_characters','one_load_lowest_signed_energy_minimizers','one_load_positive_band_edge_minimizers','global_ground','global_ground_multiplicity','global_fixed_number_gap']}))
