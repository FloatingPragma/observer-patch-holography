"""Reconstruct every retained paired chain before assigning conditional GeV."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from model.correlated_hadron_spectrum import fit
from model.quenched_hadron_binding import analysis as old
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
HBARC=.1973269804  # GeV fm: same decimal conversion in retained HAD/had_site.py.raw.
R0_BENCHMARK_FM=.5  # Explicit supplied convention, not an empirical measurement.


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update(('model/correlated_hadron_spectrum/results.json','model/correlated_hadron_spectrum/PROOF.md',
                  'model/quenched_hadron_binding/recorded_tables.json','model/quenched_hadron_binding/results.json',
                  'model/quenched_hadron_binding/inputs/HAD/had_site.py.raw'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def authenticate(value):
    for path,digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:raise ValueError('changed parent: '+path)


def coordinates(chain,kappas,anchors):
    inv={'rho':anchors['rho']/chain['mrho_phys'],
         'string':anchors['string']/chain['sqrt_sigma'],
         'r0_benchmark':chain['r0']*HBARC/R0_BENCHMARK_FM}
    values={'kappa_c':chain['kappa_c'],'kappa_phys':chain['kappa_phys'],
        'r0_sqrt_sigma':chain['r0']*chain['sqrt_sigma'],
        'mrho_over_sqrt_sigma':chain['mrho_phys']/chain['sqrt_sigma'],
        'mrho_r0':chain['mrho_phys']*chain['r0']}
    for name,scale in inv.items():
        values['inverse_a_GeV/'+name]=scale
        values['a_fm/'+name]=HBARC/scale
        values['r0_fm/'+name]=chain['r0']*HBARC/scale
        for channel,key in [('rho','mrho_phys'),('nucleon','mN_phys')]:values['hadron_GeV/'+name+'/'+channel]=scale*chain[key]
        values['hadron_GeV/'+name+'/pion']=scale*np.sqrt(chain['mpi2_phys'])
    for left,right in [('rho','string'),('rho','r0_benchmark'),('r0_benchmark','string')]:
        values['scale_ratio/'+left+'/'+right]=inv[left]/inv[right]
        values['scale_difference_GeV/'+left+'/'+right]=inv[left]-inv[right]
    for label,kappa in [('input_'+str(k),k) for k in kappas]+[('ratio_tuned',chain['kappa_phys']),('critical',chain['kappa_c'])]:
        bare=.5*(1/kappa-1/chain['kappa_c'])
        values['a_mq/'+label]=bare
        for scale,value in inv.items():values['bare_Wilson_GeV/'+scale+'/'+label]=bare*value
    return values


def summary(central,replicas):
    names=list(central);matrix=np.array([[r[k] for k in names] for r in replicas]);cov=old.jk_cov(matrix)
    if len(replicas)!=32 or not np.isfinite(matrix).all():raise ValueError('all 32 paired finite replicas required')
    return {'names':names,'central':central,'replicas':matrix.tolist(),
            'statistical':{k:float(old.jk_err(matrix[:,j])) for j,k in enumerate(names)},'covariance':cov.tolist()}


def calculate():
    pins=bindings();parent=json.loads((ROOT/'model/correlated_hadron_spectrum/results.json').read_text());authenticate(parent)
    dump=json.loads(fit.TABLE.read_text())['lattices'][fit.RULE['lattice']]
    archived=json.loads(fit.OLD_RESULT.read_text());anchors={'rho':archived['measured']['m_rho_MeV']/1000,
        'string':archived['convention']['sqrt_sigma_MeV']/1000}
    potential,(sigma,sigma_j,r0,r0_j)=old.potential(dump)
    selected=[i for i,k in enumerate(dump['kappas']) if k in parent['kappas_used']]
    choices=[('primary','SP','primary',selected),('alt_window','SP','alt',selected),('smeared_sink','SS','primary',selected)]
    if len(selected)>=4:choices.append(('lightest_three','SP','primary',selected[1:]))
    variants={};chains={}
    for label,sink,window,indices in choices:
        central,replicas=fit.mass_chain(dump,parent['fits'],sigma,sigma_j,r0,r0_j,sink,window,indices)
        for key,value in central.items():
            if abs(value-parent['chiral_variants'][label][key][0])>2e-12:raise ArithmeticError('frozen chiral chain changed')
        variants[label]=summary(coordinates(central,dump['kappas'],anchors),[coordinates(r,dump['kappas'],anchors) for r in replicas])
        chains[label]={'central':central,'paired_replicas':replicas}
    primary=variants['primary'];cross=[]
    for h,chain in chains.items():
        for p,pdata in potential['fits'].items():
            changed=dict(chain['central'],sqrt_sigma=float(np.sqrt(pdata['sigma'])),r0=float(pdata['r0']))
            cross.append({'hadron_fit':h,'potential_fit':p,'values':coordinates(changed,dump['kappas'],anchors)})
    shifts={k:max(abs(r['values'][k]-v) for r in cross) for k,v in primary['central'].items()}
    # A paired contrast uses the full covariance. Shuffling breaks the source join.
    controls=[];names=primary['names'];matrix=np.asarray(primary['replicas']);cov=np.asarray(primary['covariance'])
    for left,right in [('rho','string'),('rho','r0_benchmark'),('r0_benchmark','string')]:
        a,b=[names.index('inverse_a_GeV/'+x) for x in (left,right)]
        measured=float(old.jk_err(matrix[:,a]-matrix[:,b]));reconstructed=float(np.sqrt(max(0,cov[a,a]+cov[b,b]-2*cov[a,b])))
        if abs(measured-reconstructed)>1e-13:raise ArithmeticError('paired contrast covariance differs')
        controls.append({'scales':[left,right],'paired_difference_statistical_GeV':measured,
            'independent_error_countercontrol_GeV':float(np.sqrt(cov[a,a]+cov[b,b])),
            'cyclically_mispaired_countercontrol_GeV':float(old.jk_err(matrix[:,a]-np.roll(matrix[:,b],1))),
            'covariance_identity_error':abs(measured-reconstructed),
            'central_difference_GeV':primary['central']['scale_difference_GeV/'+left+'/'+right],
            'crossed_fit_allowance_GeV':shifts['scale_difference_GeV/'+left+'/'+right]})
    # Critically, tuning kappa_phys and kappa_c must be repeated together.
    cs=chains['primary'];tuned=np.array([.5*(1/r['kappa_phys']-1/r['kappa_c']) for r in cs['paired_replicas']])
    unpaired=np.array([.5*(1/r['kappa_phys']-1/cs['central']['kappa_c']) for r in cs['paired_replicas']])
    small=primary['central']['a_mq/ratio_tuned'];identity=cs['central']['mpi2_phys']/(2*cs['central']['slope_mpi2'])
    if abs(small-identity)>1e-14:raise ArithmeticError('tuned bare mass identity failed')
    result={'schema':'oph.paired_hadron_quark_scale.v1','source_bindings':pins,'all_checks_pass':True,
        'source_lattice':parent['source_lattice'],'all_input_kappas':dump['kappas'],'chiral_fit_kappas':parent['kappas_used'],
        'frozen_parent_fits_reused':48,'paired_stream_count':32,'paired_chains':chains,
        'conditional_coordinates':variants,'crossed_fit_choices':cross,'crossed_fit_sensitivity':shifts,
        'calibration_ancestry':{'rho_anchor_GeV':anchors['rho'],'rho_status':'imported empirical anchor; held fixed in conditional statistics',
            'string_anchor_GeV':anchors['string'],'string_status':'supplied conventional benchmark, not an exact measured physical truth',
            'r0_benchmark_fm':R0_BENCHMARK_FM,'r0_status':'explicit benchmark parameter; not measured or inferred by this calculation',
            'hbar_c_GeV_fm_decimal':HBARC,'hbar_c_role':'same retained unit conversion decimal; rounding uncertainty not propagated',
            'light_ratio':dump['ratio_phys_input'],'light_ratio_role':'empirical pion/rho ratio selects kappa_phys after fixed fitting',
            'external_anchor_uncertainties_propagated':False,'anchor_covariances_supplied':False},
        'parameterized_r0_scale':'a^-1 = (r0/a) hbar_c / r0_physical; benchmark results scale as 0.5 fm / r0_physical',
        'paired_covariance_controls':controls,
        'tuned_mass_control':{'a_mq_from_inverse_kappas':small,'a_mq_from_pion_squared_over_twice_slope':identity,
            'paired_statistical':float(old.jk_err(tuned)),'fixed_kappa_c_countercontrol_statistical':float(old.jk_err(unpaired))},
        'scope':{'bare_Wilson_action_coordinate_definition':'a m_q = (1/kappa - 1/kappa_c)/2',
            'kappa_critical_mass_zero_is_definition':True,'renormalized_MSbar_mass':False,'identified_up_down_or_strange_mass':False,
            'rest_particle_pole_mass':False,'empirical_targets_select_fit_windows':False,'replicas_removed':False,
            'fit_sensitivity_is_statistical_confidence_interval':False,'total_uncertainty_or_significance_claimed':False,
            'sea_quark_continuum_finite_volume_and_autocorrelation_systematics_bounded':False,
            'physical_source_clock_attachment':False}}
    if bindings()!=pins:raise ArithmeticError('source changed during paired evaluation')
    return result


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    value=calculate();data=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if a.write:(HERE/'results.json').write_text(data)
    elif (HERE/'results.json').read_text()!=data:raise ValueError('paired scale receipt differs')
    print('PASS: all 32 paired scales, bare Wilson coordinates, fit allowances and covariance controls.')
