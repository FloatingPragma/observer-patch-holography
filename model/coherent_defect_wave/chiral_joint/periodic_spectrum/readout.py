"""Complete energy-on-success spectroscopy of the actual periodic ground."""
from functools import lru_cache
import sympy as s
from .spectrum import resolved_block,rational_block,certified_block,sector_certificate
from .basis import coordinate_rows
from ..model import fermion
from ..two_cell.spectrum import _integer_matrix,_sparse_integer_product
from ..two_cell.readout import sqrt_upper

SCALE=1<<44
GROUND=(0,2,0,1)


@lru_cache(maxsize=None)
def approximation(block):
    _,embedding,g,k,generator=resolved_block(*block)
    d,z=rational_block(*block)
    hd,hi=_integer_matrix(generator);gd,gi=_integer_matrix(g)
    hz=_sparse_integer_product(hi,z);gz=_sparse_integer_product(gi,z)
    residual=SCALE*hz-hd*z*d[None,:];gr=_sparse_integer_product(gi,residual)
    squares=[s.Rational(sum(residual[i,j]*gr[i,j] for i in range(g.rows)),
                      gd*hd**2*SCALE**4) for j in range(g.rows)]
    assert all(r>=0 for r in squares)
    gram=z.T@gz;den=gd*SCALE**2
    for j in range(g.rows):gram[j,j]-=den
    delta=s.Rational(max(sum(abs(v) for v in row) for row in gram),den)
    assert delta<s.Rational(1,1<<23)
    return z,tuple(s.Rational(int(x),SCALE) for x in d),squares,delta


@lru_cache(maxsize=1)
def ground_approximation():
    proof=sector_certificate(0)
    assert proof['ground_certifiably_simple'] and proof['ground_window']==[
        {'L':2,'twoSz':0,'translation_sign':1,'eigen_index':0}]
    z,energies,residuals,_=approximation(GROUND);g=resolved_block(*GROUND)[2]
    w=s.ImmutableMatrix([s.Rational(int(x),SCALE) for x in z[:,0]])
    norm=s.factor((w.T*g*w)[0]);assert norm>0
    other=[]
    for row in proof['blocks']:
        for j,pair in enumerate(row['energy_intervals']):
            if (row['L'],row['twoSz'],row['translation_involution_sign'],j)!=(2,0,1,0):
                other.append(s.Rational(pair[0]))
    separation=min(other)-energies[0];assert separation>0
    error=sqrt_upper(residuals[0]/norm)/separation
    assert error<s.Rational(1,10**8)
    return w,norm,error


@lru_cache(maxsize=None)
def probe_matrix(kind,site=1,spin=0):
    if kind not in ('add','remove') or site not in (0,1) or spin not in (0,1):
        raise ValueError('specified charged success probe required')
    boundary=[0,0];boundary[site]=-1 if kind=='add' else 1
    L=1 if kind=='add' else 3
    Sz=(1-2*spin)*(1 if kind=='add' else -1)
    target=(tuple(boundary),L,Sz,0)
    sourcekeys,source,_,_,_=resolved_block(*GROUND)
    targetkeys,destination,_,_,_=resolved_block(*target)
    positions={key:i for i,key in enumerate(targetkeys)};entries={};mode=6*site+3*spin
    for j,key in enumerate(sourcekeys):
        mask,*rest=key;step=fermion(mask,mode,kind=='add')
        if step is not None:
            out=(step[0],*rest);assert out in positions
            entries[positions[out],j]=step[1]*(-1 if kind=='remove' else 1)
    raw=s.ImmutableSparseMatrix(len(targetkeys),len(sourcekeys),entries)
    image=raw*source;rows=coordinate_rows(*target[:3]);t=image.extract(rows,range(image.cols))
    t=s.ImmutableSparseMatrix(t.rows,t.cols,{(i,j):s.expand(s.conjugate(s.I**((targetkeys[rows[i]][0]&63).bit_count()))*v)
        for (i,j),v in t.todok().items()})
    assert destination*t==image
    return target,t


def clusters(block):
    proof=certified_block(*block);intervals=[tuple(map(s.Rational,pair)) for pair in proof['energy_intervals']]
    groups=[]
    for j,(lo,hi) in enumerate(intervals):
        if groups and lo<=groups[-1][2]:groups[-1][0].append(j);groups[-1][2]=max(groups[-1][2],hi)
        else:groups.append([[j],lo,hi])
    return groups,intervals


def projector_error(block,indices,intervals):
    z,energies,residuals,delta=approximation(block)
    outside=[j for j in range(len(intervals)) if j not in indices]
    if not outside:return delta
    assert all(not intervals[j][0]<=energies[i]<=intervals[j][1] for i in indices for j in outside)
    separation=min(abs(energies[i]-edge) for i in indices for j in outside for edge in intervals[j])
    assert separation>0
    sine=sqrt_upper(sum(residuals[i] for i in indices)/(1-delta))/separation
    return min(s.Integer(1),sine)+delta


@lru_cache(maxsize=None)
def charged_spectroscopy(site=1,spin=0):
    w,norm,ground_error=ground_approximation();gsource=resolved_block(*GROUND)[2]
    ground=list(map(s.Rational,sector_certificate(0)['global_ground_interval']))
    branches=[];successmetric=s.zeros(gsource.rows);sumlow=s.Rational(1,2);sumhigh=s.Rational(1,2)
    for kind in ('add','remove'):
        target,t=probe_matrix(kind,site,spin);g=resolved_block(*target)[2]
        successmetric+=t.H*g*t
        x=t*w;weight=s.factor((x.H*g*x)[0]/(2*norm))
        # Exact same-boundary conjugation and unique ground prove 1/4; the
        # full rational approximation independently agrees within its error.
        assert abs(weight-s.Rational(1,4))<=ground_error/2
        z,_,_,_=approximation(target);q=s.ImmutableMatrix(z.tolist())/SCALE
        amplitude=q.T*(g*x);groups,intervals=clusters(target);records=[]
        for indices,lo,hi in groups:
            probability=s.factor(sum(s.conjugate(amplitude[j])*amplitude[j] for j in indices)/(2*norm))
            assert probability.is_Rational
            error=(ground_error+projector_error(target,indices,intervals))/2
            precision=1<<27
            lower=s.Rational(s.floor(max(0,probability-error)*precision),precision)
            upper=s.Rational(s.ceiling(min(s.Rational(1,2),probability+error)*precision),precision)
            records.append({'energy_indices':indices,'rank':len(indices),
                'source_energy_interval':[str(lo),str(hi)],
                'energy_transfer_interval':[str(lo-ground[1]),str(hi-ground[0])],
                'original_swap_probability_interval':[str(lower),str(upper)],
                'joint_source_reference_phase_splitting_factors':
                    ['0','1/2','0','1/2'],
                'phase_labels':['0','pi/2','pi','3pi/2']})
            sumlow+=lower;sumhigh+=upper
        branches.append({'kind':kind,'target_boundary':target[0],'L':target[1],'twoSz':target[2],
            'exact_total_branch_probability':'1/4','complete_energy_windows':records})
    assert successmetric==gsource
    assert sumlow<=1<=sumhigh
    # The two no-energy-readout alternatives have the complementary n and
    # I-n effects and each original weight 1/4. Their quantum state is kept.
    failures=[{'kind':kind,'probability':'1/4','energy_readout_performed':False,
               'post_swap_quantum_state_and_all_apparatus_records_retained':True}
              for kind in ('add_fail','remove_fail')]
    return {'source':'actual two-oriented-link periodic source; all scalar/link registers quantum',
            'source_reference':'unique complete neutral Gauss ground, L2 Sz0 T+1',
            'original_ground_preparation_success':'1/3424','original_ground_preparation_failure':'3423/3424',
            'conditional_on_retained_ground_preparation':True,
            'source_ground_trace_distance_error_upper_bound':'1/100000000',
            'site':site,'spin':spin,'charged_success_branches':branches,'retained_failures':failures,
            'complete_probability_enclosure':[str(sumlow),str(sumhigh)],
            'conditional_outcome_count_including_zero_phase_records':
                4*sum(len(row['complete_energy_windows']) for row in branches)+len(failures),
            'original_attempt_outcome_count_including_preparation_failure':
                4*sum(len(row['complete_energy_windows']) for row in branches)+len(failures)+1,
            'original_attempt_probability_rule':
                'multiply each post-ground record by 1/3424; retain ground-filter failure 3423/3424',
            'zero_weight_energy_phase_records':
                'zero events retained; phase factors specify joint weights without conditioning on a zero event',
            'protocol':'energy and joint-reference phase are read only after successful charge exchange',
            'probabilities_renormalized':False,
            'joint_phase_energy_degeneracy_is_reference_position_structure':True,
            'source_only_particle_momentum_or_rest_mass_identified':False,
            'GeV_prediction':False}
