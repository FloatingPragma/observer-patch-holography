"""Recompute scheme outputs from every retained raw loop record."""
from pathlib import Path
import numpy as np
from .coupling import jackknife,gauge_control,fundamental_casimir
from model.colour_deconfinement.centre import colour_factor_identity


def run():
    rows=[];by_side={}
    for path in sorted((Path(__file__).parent/'data').glob('*.npz')):
        with np.load(path,allow_pickle=False) as z:
            side=int(z['side']);x=z['loops'];beta=float(z['beta'])
            assert x.shape==(int(z['samples']),5) and np.all(np.isfinite(x))
            assert np.max(abs(x[:,0]-z['plaquette']))<2e-13
            gauge=gauge_control(z['final_links'],side)
            assert gauge<2e-12
            blocks=[jackknife(x,side,b) for b in (5,10,20,40,80)]
            row={'raw_file':path.name,'side':side,'beta':beta,'bare_g_squared':6/beta,
                'seed':int(z['seed']),'start':str(z['start']),'warmup_sweeps':int(z['warmup']),
                'measurement_stride':int(z['stride']),'measured_records':len(x),
                'block_analyses':blocks,'first_half':jackknife(x[:200],side,20),
                'second_half':jackknife(x[200:],side,20),'final_configuration_gauge_error':gauge}
            rows.append(row);by_side.setdefault(side,[]).append(x.copy())
    if len(rows)!=4 or set(by_side)!={6,8}:raise ValueError('four complete cold/hot production ensembles required')
    combined=[]
    for side,series in sorted(by_side.items()):
        x=np.concatenate(series);report=jackknife(x,side,20)
        difference=np.array(report['successive_scale_differences']);error=np.array(report['difference_standard_errors'])
        assert np.all(difference>5*error)
        combined.append({'side':side,'independent_chains':len(series),**report,
            'scale_difference_to_block_error':(difference/error).tolist()})
    return {'schema':'oph.finite_own_ensemble_scale_coupling.v1',
        'source':'Declared four-dimensional periodic pure-colour Wilson ensemble from five-level seam trace cost',
        'five_level_trace_action_factor':'3/5','fundamental_colour_Casimir':str(fundamental_casimir()),
        'random_seam_action_identity_error':colour_factor_identity(),
        'scheme':'g_C^2(R)=-log[W(R,R)W(R-1,R-1)/W(R,R-1)^2]/k_R(N)',
        'scale_distances_in_lattice_units':[.5,1.5,2.5],
        'tree_normalization':'k_R=f(R,R)+f(R-1,R-1)-2f(R,R-1), from the finite free Green covariance',
        'rows':rows,'combined_independent_starts':combined,
        'imported_one_loop_running_curve':False,'physical_GeV_scale_selected':False,
        'continuum_beta_function_derived':False,'particle_core_onset_identified':False,
        'statistical_scope':'Blocked Monte Carlo diagnostics, not rigorous mixing or confidence bounds. Finite volume and cutoff effects remain.',
        'verdict':'PASS'}
