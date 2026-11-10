"""Certified complete charged energy windows on the actual two-cell ground.

All probabilities refer to the prepared neutral ground.  Its preparation
success/failure is a separate retained instrument, never a vacuum selection.
"""
from functools import lru_cache
from math import isqrt
from collections import defaultdict
import sympy as s
from ..model import fermion
from .basis import coordinate_rows
from .spectrum import (physical_block,rational_block,certified_block,
                       sector_certificate,_integer_matrix,_sparse_integer_product)

BITS=44
SCALE=1<<BITS
GROUND=(0,2,0)


def sqrt_upper(value,bits=80):
    value=s.Rational(value)
    if value<0:raise ValueError('nonnegative square required')
    scale=1<<bits
    numerator=int(value.p)*scale*scale;denominator=int(value.q)
    n=isqrt(numerator//denominator)
    if n*n*denominator<numerator:n+=1
    return s.Rational(n,scale)


@lru_cache(maxsize=None)
def approximation(block):
    _,embedding,g,k,generator=physical_block(*block)
    d,z=rational_block(*block)
    hd,hi=_integer_matrix(generator);gd,gi=_integer_matrix(g)
    hz=_sparse_integer_product(hi,z);gz=_sparse_integer_product(gi,z)
    residual=SCALE*hz-hd*z*d[None,:]
    gr=_sparse_integer_product(gi,residual)
    residual_squares=[s.Rational(sum(residual[i,j]*gr[i,j] for i in range(g.rows)),
                                gd*hd**2*SCALE**4) for j in range(g.rows)]
    assert all(r>=0 for r in residual_squares)
    gram=z.T@gz
    gram_den=gd*SCALE**2
    for i in range(g.rows):gram[i,i]-=gram_den
    gram_error=s.Rational(max(sum(abs(v) for v in row) for row in gram),gram_den)
    assert gram_error<s.Rational(1,1<<23)
    return z,tuple(s.Rational(int(e),SCALE) for e in d),residual_squares,gram_error


@lru_cache(maxsize=1)
def ground_approximation():
    sector=sector_certificate(0)
    assert sector['ground_is_certifiably_simple']
    assert sector['ground_window_candidate_labels']==[{'L':2,'twoSz':0,'eigen_index':0}]
    _,_,g,_,_=physical_block(*GROUND)
    z,d,residuals,_=approximation(GROUND)
    w=s.ImmutableMatrix([s.Rational(int(v),SCALE) for v in z[:,0]])
    norm=s.factor((w.T*g*w)[0]);assert norm>0
    other_lows=[]
    for b in sector['complete_conserved_blocks']:
        for j,interval in enumerate(b['energy_intervals']):
            if (b['L'],b['twoSz'],j)!=(2,0,0):other_lows.append(s.Rational(interval[0]))
    separation=min(other_lows)-d[0];assert separation>0
    trace_error=sqrt_upper(residuals[0]/norm)/separation
    assert trace_error<s.Rational(1,10**8)
    return w,norm,trace_error


def probe_column(key,mode,kind):
    mask,link,hl,hr=key
    if kind in ('add','remove'):
        out=fermion(mask,mode,kind=='add')
        if out is None:return {}
        # The even-source fermionic swap gives a fixed minus removal phase.
        return {(out[0],link,hl,hr):out[1]*(-1 if kind=='remove' else 1)}
    present=bool(mask&(1<<mode))
    if kind=='add_fail':return {key:s.Integer(-1)} if present else {}
    if kind=='remove_fail':return {} if present else {key:s.Integer(1)}
    raise ValueError('unknown complete reservoir branch')


@lru_cache(maxsize=None)
def probe_matrix(site,spin,kind):
    if site not in (0,1) or spin not in (0,1):raise ValueError('site/spin must be0 or1')
    mode=6*site+3*spin
    if kind=='add':
        boundary=[0,0];boundary[site]=-1;target=(tuple(boundary),1,1-2*spin)
    elif kind=='remove':
        boundary=[0,0];boundary[site]=1;target=(tuple(boundary),3,2*spin-1)
    elif kind in ('add_fail','remove_fail'):target=GROUND
    else:raise ValueError('unknown branch')
    source_keys,source_embedding,_,_,_=physical_block(*GROUND)
    target_keys,target_embedding,_,_,_=physical_block(*target)
    positions={key:i for i,key in enumerate(target_keys)}
    entries={}
    for j,key in enumerate(source_keys):
        for dest,value in probe_column(key,mode,kind).items():
            assert dest in positions,('probe leaves declared diagonal Gauss sector',dest)
            entries[positions[dest],j]=value
    raw=s.ImmutableSparseMatrix(len(target_keys),len(source_keys),entries)
    image=raw*source_embedding
    rows=coordinate_rows(*target)
    result=s.ImmutableSparseMatrix(image.extract(rows,range(image.cols)))
    # The spectrum's embedding includes D=i^(N_left); remove that phase
    # when reading the exact coordinate rows of its rational Gauss basis.
    result=s.ImmutableSparseMatrix(result.rows,result.cols,
        {(i,j):s.expand(s.conjugate(s.I**((target_keys[rows[i]][0]&63).bit_count()))*v)
         for (i,j),v in result.todok().items()})
    assert target_embedding*result==image
    return target,result


def energy_clusters(block):
    proof=certified_block(*block)
    intervals=[tuple(map(s.Rational,pair)) for pair in proof['energy_intervals']]
    groups=[]
    for j,(lo,hi) in enumerate(intervals):
        if groups and lo<=groups[-1][2]:
            groups[-1][0].append(j);groups[-1][2]=max(groups[-1][2],hi)
        else:groups.append([[j],lo,hi])
    return groups,intervals


def projector_error(block,indices,intervals):
    z,energies,residuals,delta=approximation(block)
    outside=[j for j in range(len(intervals)) if j not in indices]
    if not outside:return delta
    sep=min(min(abs(energies[i]-intervals[j][0]),abs(energies[i]-intervals[j][1]))
            for i in indices for j in outside)
    for i in indices:
        assert all(not lo<=energies[i]<=hi for j,(lo,hi) in enumerate(intervals) if j in outside)
    assert sep>0
    # Equal-rank spectral projection perturbation, certified by individual
    # column residuals outside the exact separated energy cluster.
    sine=sqrt_upper(sum(residuals[i] for i in indices)/(1-delta))/sep
    return min(s.Integer(1),sine)+delta


@lru_cache(maxsize=4)
def complete_experiment(site=1,spin=0):
    w,norm,ground_error=ground_approximation()
    ground_interval=list(map(s.Rational,sector_certificate(0)['global_ground_interval']))
    branches=[];all_lows=[];all_highs=[];branch_weights=[];approx_total=s.Integer(0)
    for kind in ('add','add_fail','remove','remove_fail'):
        target,t=probe_matrix(site,spin,kind)
        _,_,g,_,_=physical_block(*target)
        x=t*w
        branch_weight=s.factor((x.H*g*x)[0]/(2*norm))
        z,_,_,_=approximation(target)
        q=s.ImmutableMatrix(z.tolist())/SCALE
        amplitudes=q.T*(g*x)
        clusters,intervals=energy_clusters(target)
        records=[]
        for indices,lo,hi in clusters:
            probability=s.factor(sum(s.conjugate(amplitudes[j])*amplitudes[j] for j in indices)/(2*norm))
            assert probability.is_Rational
            error=(ground_error+projector_error(target,indices,intervals))/2
            lower=max(s.Integer(0),probability-error);upper=min(s.Rational(1,2),probability+error)
            # Export stable coarse rational outward bounds, never normalized.
            precision=1<<27
            low=s.Rational(s.floor(lower*precision),precision)
            high=s.Rational(s.ceiling(upper*precision),precision)
            records.append({'energy_indices':indices,'spectral_rank':len(indices),
                'source_energy_interval':[str(lo),str(hi)],
                'energy_transfer_interval':[str(lo-ground_interval[1]),str(hi-ground_interval[0])],
                'probability_interval':[str(low),str(high)]})
            all_lows.append(low);all_highs.append(high);approx_total+=probability
        exact_complete_weight=s.factor(branch_weight)
        branch_weights.append(exact_complete_weight)
        precision=1<<27
        branch_low=s.Rational(s.floor(max(0,exact_complete_weight-ground_error/2)*precision),precision)
        branch_high=s.Rational(s.ceiling(min(s.Rational(1,2),exact_complete_weight+ground_error/2)*precision),precision)
        branches.append({'kind':kind,'boundary':target[0],'L':target[1],'twoSz':target[2],
            'complete_energy_windows':records,
            'unreachable_conserved_blocks_have_exact_zero_effect_on_prepared_ground':True,
            'original_branch_weight_interval':[str(branch_low),str(branch_high)]})
    assert sum(branch_weights)==1
    assert sum(all_lows)<=1<=sum(all_highs)
    assert abs(approx_total-1)<s.Rational(1,10**8)
    # Complete Kraus identity holds on every physical vector in the source
    # block, independently of the reference approximation.
    _,_,gsource,_,_=physical_block(*GROUND)
    complete=s.zeros(gsource.rows)
    for kind in ('add','add_fail','remove','remove_fail'):
        target,t=probe_matrix(site,spin,kind);g=physical_block(*target)[2]
        complete+=t.H*g*t/2
    assert complete==gsource
    return {'site':site,'spin_component':spin,'source_ground_block':{'boundary':0,'L':2,'twoSz':0},
            'source_ground_is_exact_unique_physical_ground':True,
            'conditional_on_retained_ground_preparation':True,
            'actual_even_parity_fermionic_reservoir_and_compensator':True,
            'exact_complete_Kraus_metric_identity':True,
            'source_energy_clock':'same full two-cell Hamiltonian',
            'records':branches,'probabilities_renormalized':False,
            'total_probability_enclosure':[str(sum(all_lows)),str(sum(all_highs))],
            'physical_rest_mass_or_GeV_identification':False}


def unconditional_experiment(site=1,spin=0):
    """Keep the failed binary preparation and all conditional energy records."""
    from .instrument import ground_preparation
    prepared=ground_preparation()
    success=s.Rational(prepared['original_ground_success_probability'])
    failure=s.Rational(prepared['original_preparation_failure_probability'])
    conditional=complete_experiment(site,spin)
    records=[{'outcome':'preparation_failure','probability_interval':[str(failure)]*2,
              'energy_measurement_performed':False}]
    for branch in conditional['records']:
        for energy in branch['complete_energy_windows']:
            records.append({'outcome':branch['kind'],'boundary':branch['boundary'],
                'L':branch['L'],'twoSz':branch['twoSz'],**energy,
                'probability_interval':[str(success*s.Rational(v)) for v in energy['probability_interval']],
                'energy_measurement_performed':True})
    low=sum(s.Rational(r['probability_interval'][0]) for r in records)
    high=sum(s.Rational(r['probability_interval'][1]) for r in records)
    assert low<=1<=high
    return {'site':site,'spin_component':spin,'preparation':prepared,
            'record_count_including_preparation_failure':len(records),
            'records':records,'probabilities_renormalized':False,
            'unreachable_conserved_blocks_have_exact_zero_effect':True,
            'total_probability_enclosure':[str(low),str(high)]}


def reflected_experiments():
    """Numerical enclosures of the experiments related by the exact symmetry."""
    right=complete_experiment(1,0);left=complete_experiment(0,1)
    pairs=[]
    # A preserves the spin label of a mode: addition spin up becomes removal
    # spin up. Theta supplies the independent spin reversal used here.
    for a,b in (('add','remove'),('remove','add'),('add_fail','remove_fail'),('remove_fail','add_fail')):
        ra=next(r for r in right['records'] if r['kind']==a)
        rb=next(r for r in left['records'] if r['kind']==b)
        assert len(ra['complete_energy_windows'])==len(rb['complete_energy_windows'])
        for ea,eb in zip(ra['complete_energy_windows'],rb['complete_energy_windows']):
            assert ea['spectral_rank']==eb['spectral_rank']
            for key in ('energy_transfer_interval','probability_interval'):
                ia=list(map(s.Rational,ea[key]));ib=list(map(s.Rational,eb[key]))
                assert max(ia[0],ib[0])<=min(ia[1],ib[1])
        pairs.append({'right':a,'left':b,'compared_energy_windows':len(ra['complete_energy_windows'])})
    return {'exact_pairing':'antiunitary cell-reflected charge conjugation composed with spin reversal',
            'right_site_spin':[1,0],'left_site_spin':[0,1],
            'pairs':pairs,'same_endpoint_addition_removal_equality_asserted':False}


if __name__=='__main__':
    import json
    print(json.dumps(complete_experiment(),indent=2,sort_keys=True))
