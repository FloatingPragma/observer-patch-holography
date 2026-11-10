"""Complete scalar/link-loop energy and translation instruments on one ground."""
from functools import lru_cache
import sympy as s
from .model import SETTINGS,operator_column
from ..periodic_spectrum.spectrum import resolved_block,neutral_translation,_kernel,sector_certificate
from ..periodic_spectrum.basis import coordinate_rows
from ..periodic_spectrum.readout import GROUND,ground_approximation,approximation,clusters,SCALE
from ..neutral_composite.spectral_bounds import projector_error

TARGET=(0,2,0,0)


@lru_cache(maxsize=None)
def unresolved_transfer(setting,outcome='success'):
    sourcekeys,source,_,_,_=resolved_block(*GROUND)
    keys,b,_,_,_=resolved_block(*TARGET);positions={key:i for i,key in enumerate(keys)};entries={}
    for col,key in enumerate(sourcekeys):
        for dest,value in operator_column(key,setting,outcome).items():
            assert dest in positions,('omitted scalar/loop Gauss output',dest)
            entries[positions[dest],col]=value
    raw=s.ImmutableSparseMatrix(len(keys),len(sourcekeys),entries);image=raw*source
    rows=coordinate_rows(*TARGET[:3]);restricted=image.extract(rows,range(image.cols))
    result=s.ImmutableSparseMatrix(restricted.rows,restricted.cols,
        {(i,j):s.expand(s.conjugate(s.I**((keys[rows[i]][0]&63).bit_count()))*value)
         for (i,j),value in restricted.todok().items()})
    assert (b*result-image).applyfunc(s.expand)==s.zeros(*image.shape)
    return result


@lru_cache(maxsize=None)
def translation_transfer(setting,sign=1):
    if sign not in (-1,1):raise ValueError('translation involution sign')
    t=unresolved_transfer(setting);c,_=neutral_translation(2,0)
    f,rows=_kernel(c-sign*s.eye(c.rows));projected=(t+sign*c*t)/2
    result=s.ImmutableSparseMatrix(projected.extract(rows,range(t.cols)))
    assert f*result==projected
    return (*TARGET[:3],sign),result


def transfer_checks():
    g0=resolved_block(*GROUND)[2];g=resolved_block(*TARGET)[2];rows=[]
    for setting in SETTINGS:
        t=unresolved_transfer(setting);f=unresolved_transfer(setting,'failure')
        assert (t.H*g*t+f.H*g*f-g0).applyfunc(s.expand)==s.zeros(g0.rows)
        metric=s.zeros(g0.rows);dims=[]
        for sign in (-1,1):
            target,ts=translation_transfer(setting,sign);gs=resolved_block(*target)[2]
            metric+=ts.H*gs*ts;dims.append(gs.rows)
            if setting!='scalar_number' and sign==-1:assert ts==s.zeros(*ts.shape)
        assert metric==t.H*g*t
        rows.append({'setting':setting,'source_dimension':g0.rows,'full_target_dimension':g.rows,
                     'translation_dimensions':dims,'complete_source_Kraus_identity':True,
                     'exact_Gauss_transfer_and_translation_completeness':True,
                     'negative_translation_success_exactly_zero':setting!='scalar_number'})
    return rows


def outward(value,error,maximum=s.Rational(1,3)):
    denominator=1<<27
    return [str(s.Rational(s.floor(max(0,value-error)*denominator),denominator)),
            str(s.Rational(s.ceiling(min(maximum,value+error)*denominator),denominator))]


def complete_readout():
    w,norm,ground_error=ground_approximation();g=resolved_block(*TARGET)[2]
    ground=list(map(s.Rational,sector_certificate(0)['global_ground_interval']))
    total_low=s.Integer(0);total_high=s.Integer(0);settings=[]
    for setting in SETTINGS:
        t=unresolved_transfer(setting);x=t*w
        success=s.factor((x.H*g*x)[0]/(3*norm));success_interval=outward(success,ground_error/3)
        failure_interval=outward(s.Rational(1,3)-success,ground_error/3)
        total_low+=s.Rational(failure_interval[0]);total_high+=s.Rational(failure_interval[1])
        phases=[];positive=[]
        for sign in (-1,1):
            target,ts=translation_transfer(setting,sign);gs=resolved_block(*target)[2]
            groups,intervals=clusters(target);zero=ts==s.zeros(*ts.shape)
            if not zero:
                z,_,_,_=approximation(target);q=s.ImmutableMatrix(z.tolist())/SCALE
                amplitudes=q.T*(gs*(ts*w))
            records=[]
            for indices,lo,hi in groups:
                if zero:window=['0','0']
                else:
                    probability=s.factor(sum(s.conjugate(amplitudes[j])*amplitudes[j] for j in indices)/(3*norm))
                    error=(ground_error+projector_error(target,indices))/3
                    window=outward(probability,error)
                delta=[lo-ground[1],hi-ground[0]]
                row={'energy_indices':indices,'rank':len(indices),
                     'source_energy_interval':[str(lo),str(hi)],
                     'energy_transfer_interval':list(map(str,delta)),
                     'setting_success_energy_phase_probability_interval':window,
                     'original_including_ground_preparation_probability_interval':
                         [str(s.Rational(v)/3424) for v in window],
                     'elastic_ground_projector':sign==1 and indices==[0]}
                records.append(row);total_low+=s.Rational(window[0]);total_high+=s.Rational(window[1])
                if delta[0]>0 and s.Rational(window[0])>0:positive.append({'translation_sign':sign,**row})
            assert sum(row['rank'] for row in records)==gs.rows
            phases.append({'translation_sign':sign,'phase':'pi' if sign==-1 else '0',
                'target_dimension':gs.rows,'success_map_exactly_zero':zero,'complete_energy_records':records})
        settings.append({'setting':setting,'choice_probability':'1/3',
            'total_success_probability_interval':success_interval,'retained_failure_probability_interval':failure_interval,
            'failure_unnormalized_state':'F rho_ground F/3; full unresolved target retained',
            'failure_energy_and_translation_readout_performed':False,
            'original_including_ground_preparation_failure_probability_interval':
                [str(s.Rational(v)/3424) for v in failure_interval],
            'zero_translation_phase_probabilities':{'pi/2':'0','-pi/2':'0'},
            'complete_translation_branches':phases,
            'lowest_certified_positive_inelastic_record':min(positive,key=lambda r:s.Rational(r['energy_transfer_interval'][0])) if positive else None})
    assert total_low<=1<=total_high
    return {'source_reference':'same complete neutral periodic Gauss ground',
            'ground_preparation_success_probability':'1/3424',
            'retained_ground_preparation_failure_probability':'3423/3424',
            'retained_ground_preparation_failure_unnormalized_state':'(I_neutral-P_ground)/3424',
            'conditional_on_recorded_ground_preparation':True,'success_probabilities_renormalized':False,
            'elastic_projectors_retained':True,'settings':settings,
            'complete_conditional_probability_enclosure':[str(total_low),str(total_high)],
            'complete_original_probability_enclosure':[str(s.Rational(3423,3424)+total_low/3424),str(s.Rational(3423,3424)+total_high/3424)],
            'Higgs_photon_W_Z_gluon_particle_identification':False,
            'physical_Lorentz_spin_or_mass':False,'GeV_units':False}
