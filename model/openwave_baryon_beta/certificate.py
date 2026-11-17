"""Complete probability laws and independent finite propagation checks."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from scipy.linalg import expm
from .decay import construction,eigensystem,evolve,records,G
from model.openwave_calibrated_hadrons.source import parameters,INPUTS
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def pins():
    seeds=set(SOURCE_FILES)
    for name in ('openwave_baryon_beta','openwave_baryon_weak_formfactors','openwave_calibrated_hadrons'):
        seeds.update(str(p.relative_to(ROOT))for p in (ROOT/'model'/name).iterdir()if p.suffix in ('.py','.md','.tex')or p.name=='inputs.json')
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(seeds))}

def calculate():
    d=construction();e=eigensystem(d);h=d['H'].toarray();rows=[records(d,t,e)for t in(0,16,32,64,128,256,512,1024)]
    probe=expm(-16j*(h-e[2]*np.eye(len(h))))@d['initial']
    propagation_error=float(np.linalg.norm(probe-evolve(d,16,e)));assert propagation_error<2e-8
    zero=construction(0.);ez=eigensystem(zero);z=evolve(zero,512,ez)
    zero_emission=float(np.sum(abs(z[zero['ne']!=0])**2));assert zero_emission<1e-20
    energy=lambda v:float(np.vdot(v,d['H']@v).real)
    prep_work=d['success']*energy(d['initial'])+energy(d['failure'])-energy(d['seed'])
    hf=d['native']['h'];charge=np.rint(np.diag(d['native']['q']).real).astype(int)
    edelta=float(np.linalg.eigvalsh(hf[np.ix_(charge==2,charge==2)])[0])
    en=float(np.vdot(d['neutron'][:,0],hf@d['neutron'][:,0]).real)
    unit=float(INPUTS['electron_unit_GeV']);endpoint=(float(parameters()['delta'])-1)*unit
    reverse_margin=en-edelta-2;assert endpoint>0 and reverse_margin<0
    ref_current=float(np.vdot(d['initial'],d['V']@d['initial']).real);assert abs(ref_current+1)<1e-11
    # Each colour-singlet label expands into disjoint native occupation words.
    b=d['native']['B'].tocsc();colour=[]
    for hi in range(120):
        start,end=b.indptr[hi:hi+2];colour.append([[str(d['native']['basis'][int(i)]),float(abs(a)**2)]for i,a in zip(b.indices[start:end],b.data[start:end])])
    assert max(abs(sum(x[1]for x in row)-1)for row in colour)<1e-14
    return {'schema':'native-baryon-beta-v1','source_sha256':pins(),
        'construction':{'native_mode_count':60,'three_quark_raw_dimension':2024,'lepton_number_four_dimension':70,
            'complete_colour_singlet_product_dimension':d['source_dimension'],'complete_charge_zero_dimension':d['charge_zero_dimension'],
            'exact_twice_spin_projection_one_dimension':d['source_spin_dimension'],'contact_coefficient':G,
            'lepton_kinetic_coefficient':d['kinetic'],'kinetic_selection':'(calibrated neutron-proton gap minus electron energy)/2; declared outgoing-box scale',
            'electron_unit_GeV':unit,'outgoing_neutral_interpretation':'Hole in filled native neutral doublet reference; massless for this member.'},
        'calibrated_threshold':{'electron_endpoint_GeV':endpoint,'double_electron_Delta_endpoint_GeV':reverse_margin*unit,
            'endpoint_domain':'Strong hadron source plus free outgoing leptons; not an interacting pole endpoint.',
            'target_independent_mass_prediction':False,'on_shell_two_electron_channel_allowed':False,
            'virtual_two_electron_channel_retained':True},
        'filled_neutral_reference_contact_shifts':{
            'exact_neutron_reference_current_expectation':'-1','exact_proton_reference_current_expectation':'1',
            'neutron_source_expectation_shift_GeV':-G*unit,'proton_source_expectation_shift_GeV':G*unit,
            'reference_expectation_gap_GeV':(float(parameters()['delta'])-2*G)*unit,
            'scope':'Reference expectation energies include GV. The fitted masses belong to H_h. Neither these expectation values nor the strong/free endpoint are interacting physical pole masses or a continuum decay endpoint.'},
        'complete_instrument':{'spin_charge_product_indices':d['keep'].tolist(),'lepton_masks':[str(x)for x in d['masks']],
            'native_colour_mask_weights_by_hadron_label':colour,'failure_probability':1-d['success'],
            'complement':'All other Fock occupations retained as orthogonal zero-probability outcomes on the prepared source; preparation failure retains its full quantum state.',
            'probability_rule':'For each position or momentum run, original outcome probability = preparation success times stored product probability times native colour weight. Both experiments separately include failure.'},
        'records':rows,'preparation':{'seed_product_index':d['seed_index'],'success_probability':d['success'],
            'mean_work_electron_units':prep_work,'seed_energy_electron_units':energy(d['seed'])},
        'checks':{'independent_pade_vs_spectral_state_error_at_16':propagation_error,
            'zero_contact_emission_at_512':zero_emission,
            'max_energy_drift':max(abs(r['source_energy']-rows[0]['source_energy'])for r in rows)},
        'scope':'Complete finite coherent conversion, including virtual two-electron states and recurrence. Calibrated masses and declared box/contact. No observed neutron lifetime, continuum outgoing law, or successful axial-ratio prediction.'}

def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group();g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=parser.parse_args()
    before=pins();out=calculate();assert before==pins();path=HERE/'results.json'
    if args.check:assert equivalent(json.loads(path.read_text()),out)
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'dimensions':out['construction'],'threshold':out['calibrated_threshold'],'checks':out['checks'],
        'at_512':{k:v for k,v in out['records'][6].items()if k not in ('position_probability','momentum_probability')}},indent=2))

if __name__=='__main__':main()
