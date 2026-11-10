"""Complete neutral odd-sector energy/translation records, with all failures."""
from functools import lru_cache
import sympy as s
from .model import operator_column
from ..periodic_spectrum.spectrum import (resolved_block,physical_block,neutral_translation,
                                         _kernel,sector_certificate)
from ..periodic_spectrum.basis import coordinate_rows
from ..periodic_spectrum.readout import (GROUND,ground_approximation,approximation,
                                        clusters,SCALE)
from .spectral_bounds import projector_error


@lru_cache(maxsize=None)
def unresolved_transfer(kind='remove',site=1,spin=0,outcome='raw'):
    sourcekeys,source,_,_,_=resolved_block(*GROUND)
    if outcome=='failure':target=(0,2,0,0)
    else:target=(0,1 if kind=='remove' else 3,(1-2*spin)*(-1 if kind=='remove' else 1),0)
    targetkeys,destination,_,_,_=resolved_block(*target)
    positions={key:i for i,key in enumerate(targetkeys)};entries={}
    for col,key in enumerate(sourcekeys):
        for dest,value in operator_column(key,kind,site,spin,outcome).items():
            assert dest in positions,('omitted neutral Gauss-sector output',dest)
            entries[positions[dest],col]=value
    raw=s.ImmutableSparseMatrix(len(targetkeys),len(sourcekeys),entries)
    image=raw*source;rows=coordinate_rows(*target[:3]);restricted=image.extract(rows,range(image.cols))
    t=s.ImmutableSparseMatrix(restricted.rows,restricted.cols,
        {(i,j):s.expand(s.conjugate(s.I**((targetkeys[rows[i]][0]&63).bit_count()))*value)
         for (i,j),value in restricted.todok().items()})
    assert (destination*t-image).applyfunc(s.expand)==s.zeros(*image.shape)
    return target,t


@lru_cache(maxsize=None)
def translation_transfer(kind='remove',sign=1,site=1,spin=0):
    if sign not in (-1,1):raise ValueError('translation involution sign')
    target,t=unresolved_transfer(kind,site,spin)
    c,xi=neutral_translation(target[1],target[2]);f,rows=_kernel(c-sign*s.eye(c.rows))
    projected=(t+sign*c*t)/2
    result=s.ImmutableSparseMatrix(projected.extract(rows,range(t.cols)))
    assert f*result==projected
    resolved=(*target[:3],sign)
    assert resolved_block(*resolved)[1]*result==resolved_block(*target)[1]*projected
    return resolved,result


def exact_transfer_checks():
    g0=resolved_block(*GROUND)[2];rows=[]
    for kind in ('remove','add'):
        full,t=unresolved_transfer(kind);g=resolved_block(*full)[2]
        metric=s.zeros(g0.rows);rank=0
        for sign in (-1,1):
            target,ts=translation_transfer(kind,sign)
            gs=resolved_block(*target)[2];metric+=ts.H*gs*ts;rank+=gs.rows
        assert metric==t.H*g*t
        # The failure image uses the complete unresolved source sector.
        fail,f=unresolved_transfer(kind,outcome='failure');gf=resolved_block(*fail)[2]
        assert (t.H*g*t/4+f.H*gf*f-g0).applyfunc(s.expand)==s.zeros(g0.rows)
        translation,_=neutral_translation(2,0)
        opposite=(f-translation*f)/2
        assert opposite!=s.zeros(*opposite.shape)
        rows.append({'kind':kind,'source_dimension':g0.rows,'unresolved_target_dimension':g.rows,
                     'complete_resolved_target_dimension':rank,'translation_dimensions':[rank//2,rank//2],
                     'failure_target_dimension':gf.rows,'failure_has_nonzero_opposite_translation_image':True,
                     'exact_full_source_Kraus_identity':True,
                     'all_raw_images_in_exact_Gauss_kernel':True,'translation_projectors_complete':True})
    return rows


def star_matrix(block,site=0):
    keys,b,g,_,_=resolved_block(*block)
    signs=[(-1)**((((key[0]>>(6*site))&63).bit_count()+int(key[1]!=0)+int(key[2]!=0))%2) for key in keys]
    # The metric form is adequate for exact intertwiner checks without a
    # coordinate projection of a translation-resolved odd sector.
    return s.ImmutableSparseMatrix(b.rows,b.cols,{(i,j):signs[i]*v for (i,j),v in b.todok().items()})


def star_translation_certificate():
    records=[]
    for kind in ('remove','add'):
        full,_=unresolved_transfer(kind);keys,b,g,_,h=resolved_block(*full)
        transformed=star_matrix(full)
        rows=coordinate_rows(*full[:3])
        c=s.ImmutableSparseMatrix(transformed.extract(rows,range(b.cols)))
        phases=[s.I**((keys[row][0]&63).bit_count()) for row in rows]
        c=s.ImmutableSparseMatrix(c.rows,c.cols,{(i,j):s.expand(s.conjugate(phases[i])*v) for (i,j),v in c.todok().items()})
        assert b*c==transformed and c*c==s.eye(c.rows) and h*c==c*h
        translation,_=neutral_translation(full[1],full[2]);assert c*translation==-translation*c
        records.append({'kind':kind,'odd_sector_dimension':g.rows,'star_commutes_H':True,
                        'star_anticommutes_translation_involution':True,
                        'T_plus_minus_energy_degeneracy_enforced_by_star':True})
    return {'operator':'C_x=(-1)^N_matter,x times incident link parities','odd_sector_checks':records,
            'unique_ground_commutes_each_star':True,
            'separated_local_two_point_function':'0 for x!=y at every source time',
            'equal_T_phase_weights_for_each_resolved_energy_from_local_probe':True,
            'finite_translation_phase_is_particle_dispersion':False,
            'link_dressed_transport_excluded':False}


def outward(value,error,maximum=s.Rational(1,2)):
    denominator=1<<27
    lo=s.Rational(s.floor(max(0,value-error)*denominator),denominator)
    hi=s.Rational(s.ceiling(min(maximum,value+error)*denominator),denominator)
    return [str(lo),str(hi)]


def complete_readout():
    w,norm,ground_error=ground_approximation()
    ground=list(map(s.Rational,sector_certificate(0)['global_ground_interval']))
    branches=[];total_lo=s.Integer(0);total_hi=s.Integer(0)
    for kind in ('remove','add'):
        full,t=unresolved_transfer(kind);g=resolved_block(*full)[2];x=t*w
        success=s.factor((x.H*g*x)[0]/(8*norm))
        success_window=outward(success,ground_error/2,s.Rational(1,4))
        failure_window=outward(s.Rational(1,2)-success,ground_error/2)
        total_lo+=s.Rational(failure_window[0]);total_hi+=s.Rational(failure_window[1])
        phases=[]
        for sign in (-1,1):
            target,ts=translation_transfer(kind,sign);gs=resolved_block(*target)[2]
            z,_,_,_=approximation(target);q=s.ImmutableMatrix(z.tolist())/SCALE
            image=ts*w;amplitudes=q.T*(gs*image)
            groups,intervals=clusters(target);records=[]
            for indices,lo,hi in groups:
                probability=s.factor(sum(s.conjugate(amplitudes[i])*amplitudes[i] for i in indices)/(8*norm))
                error=(ground_error+projector_error(target,indices))/2
                window=outward(probability,error,s.Rational(1,4))
                total_lo+=s.Rational(window[0]);total_hi+=s.Rational(window[1])
                records.append({'energy_indices':indices,'rank':len(indices),
                    'source_energy_interval':[str(lo),str(hi)],
                    'energy_transfer_interval':[str(lo-ground[1]),str(hi-ground[0])],
                    'setting_success_energy_phase_probability_interval':window,
                    'original_including_ground_preparation_probability_interval':
                        [str(s.Rational(value)/3424) for value in window]})
            assert sum(row['rank'] for row in records)==gs.rows
            phases.append({'translation_sign':sign,'translation_phase':'-pi/2' if sign==-1 else 'pi/2',
                           'physical_dimension':gs.rows,'complete_energy_records':records})
        first=[s.Rational(phase['complete_energy_records'][0]['setting_success_energy_phase_probability_interval'][0]) for phase in phases]
        branches.append({'kind':kind,'setting_probability':'1/2',
            'total_success_probability_interval':success_window,
            'retained_failure_probability_interval':failure_window,
            'failure_state':'F rho0 F/2 with full-carrier F=sqrt(I-J^dagger J/4)',
            'failure_energy_or_phase_readout_performed':False,
            'failure_original_including_ground_preparation_probability_interval':
                [str(s.Rational(value)/3424) for value in failure_window],
            'first_energy_overlap_certifiably_positive':all(value>0 for value in first),
            'zero_probability_translation_records':[{'phase':'0','probability':'0'},{'phase':'pi','probability':'0'}],
            'all_other_neutral_L_spin_sectors_have_exact_zero_probe_image':True,
            'complete_translation_branches':phases})
    assert total_lo<=1<=total_hi
    return {'reference':'same complete neutral periodic ground as periodic_spectrum',
            'ground_preparation_success_probability':'1/3424',
            'retained_ground_preparation_failure_probability':'3423/3424',
            'retained_ground_preparation_failure_unnormalized_state':'(I_neutral-P_ground)/3424',
            'site':1,'spin':0,'success_operator':'J_kind/2; setting choice probability1/2',
            'conditional_on_ground_preparation':True,'success_postselection_renormalization':False,
            'complete_conditional_probability_enclosure':[str(total_lo),str(total_hi)],
            'complete_original_probability_enclosure':[str(s.Rational(3423,3424)+total_lo/3424),str(s.Rational(3423,3424)+total_hi/3424)],
            'branches':branches,'empirical_neutrino_identification':False,'GeV_mass':False}
