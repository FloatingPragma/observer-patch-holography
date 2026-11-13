"""Two-load coherent transport with the declared active memory voltage.

Exact amplitudes use Z[zeta], zeta=exp(i*pi/3), zeta^2=zeta-1.
The event register is retained. Original native records have their explicit
inverse-recoding commit policy; an unread record channel is a different map.
"""
from collections import deque
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from types import MappingProxyType
import json

from model.openwave_active_native_charge import source as charge
from model.openwave_native_process import source as native
from model.openwave_coherent_native_architecture import source as architecture

ROOTS=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))
ZERO=(0,0)

def add(x,y):return(x[0]+y[0],x[1]+y[1])
def sub(x,y):return(x[0]-y[0],x[1]-y[1])
def scale(c,x):return(c*x[0],c*x[1])
def mul(x,y):return(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def conj(x):return(x[0]+x[1],-x[1])
def norm2(x):return x[0]*x[0]+x[0]*x[1]+x[1]*x[1]
def real2(x):return 2*x[0]+x[1]
def phase(k,chi):
    if type(k)is not int or not 0<=k<6:raise ValueError('one of six memory characters required')
    return ROOTS[(k*chi)%6]
def rotate(x,e):
    a,b=x;e%=6
    if e==0:return x
    if e==1:return(-b,a+b)
    if e==2:return(-a-b,a)
    if e==3:return(-a,-b)
    if e==4:return(b,-a-b)
    return(a+b,-a)

def voltage(before,after,seam,extra=True):return charge.voltage(before,after,seam,extra)

@dataclass(frozen=True)
class Carrier:
    level:int|None
    ports:int
    seams:tuple
    pairs:tuple
    index:object
    edges:tuple
    seam_edges:tuple
    degree:tuple
    configuration_neighbours:tuple
    @property
    def states(self):return len(self.pairs)
    @property
    def events(self):return 2*len(self.seams)

@lru_cache(None)
def carrier(level=None):
    if level is None:ports=12;seams=tuple(charge.SEAMS)
    elif type(level)is int and level==0:
        g=native.graph(level);ports=g.ports;seams=tuple(g.seams)
    else:raise ValueError('local12 ports or complete level0 federation required')
    pairs=tuple(combinations(range(ports),2));index={p:i for i,p in enumerate(pairs)}
    edges=[];by_seam=[];neighbours=[[]for _ in pairs]
    for a,b in seams:
        row=[]
        for spectator in range(ports):
            if spectator in(a,b):continue
            x=tuple(sorted((a,spectator)));y=tuple(sorted((b,spectator)))
            i,j=index[x],index[y]
            if i>j:i,j=j,i;x,y=y,x
            chi=voltage(x,y,(a,b));assert voltage(y,x,(a,b))==(-chi)%6
            item=(i,j,chi);row.append(item);edges.append(item);neighbours[i].append(j);neighbours[j].append(i)
        by_seam.append(tuple(row))
    return Carrier(level,ports,seams,pairs,MappingProxyType(index),tuple(edges),tuple(by_seam),tuple(map(len,neighbours)),tuple(map(tuple,neighbours)))

@lru_cache(None)
def _phase_edges(level,k,extra):
    g=carrier(level);phase(k,0)
    if type(extra)is not bool:raise ValueError('explicit voltage control required')
    return tuple((i,j,k*(c if extra else 3*(c%2))%6)for i,j,c in g.edges)

def phase_edges(g,k,extra=True):return _phase_edges(g.level,k,extra)

def apply_T_num(g,values,k,extra=True):
    """D*T_k action on every cyclotomic integer source coordinate."""
    if len(values)!=g.states:raise ValueError('complete source vector required')
    out=[scale(g.events-d,x)for d,x in zip(g.degree,values,strict=True)]
    for i,j,e in phase_edges(g,k,extra):
        out[j]=add(out[j],rotate(values[i],e));out[i]=add(out[i],rotate(values[j],-e))
    return out

def apply_T0_integer(g,values):
    if len(values)!=g.states:raise ValueError('complete source vector required')
    out=[(g.events-d)*x for d,x in zip(g.degree,values,strict=True)]
    for i,j,_ in g.edges:out[i]+=values[j];out[j]+=values[i]
    return out

def cyclic_step(g,a,b,k,extra=True):
    t=apply_T_num(g,b,k,extra);return[scale(-g.events,x)for x in b],[add(scale(g.events,x),scale(2,y))for x,y in zip(a,t,strict=True)]

def terminal_numerators(g,a,b,k,extra=True):
    """D times complete terminal weights; the squared term uses T_0."""
    tb=apply_T_num(g,b,k,extra);squares=apply_T0_integer(g,[norm2(x)for x in b])
    return[g.events*norm2(x)+q+real2(mul(conj(x),y))for x,q,y in zip(a,squares,tb,strict=True)]

def event_swap_column(g,i,event,k,extra=True):
    if not 0<=i<g.states or not 0<=event<g.events:raise ValueError('valid source/event address required')
    phase(k,0)
    if event%2==0:return i,event,(1,0)
    a,b=g.seams[event//2];pair=g.pairs[i]
    if(a in pair)==(b in pair):return i,event,(1,0)
    out=tuple(sorted(b if p==a else a if p==b else p for p in pair))
    return g.index[out],event,phase(k,voltage(pair,out,(a,b),extra))

def apply_T(g,values,k,extra=True):
    """Floating full-source diagnostic/API; exact certificates use apply_T_num."""
    import numpy as np
    values=np.asarray(values);out=values*(np.asarray([g.events-d for d in g.degree]).reshape((-1,)+(1,)*(values.ndim-1)))
    out=out.astype(complex)
    for i,j,e in phase_edges(g,k,extra):
        p=complex(ROOTS[e][0]+ROOTS[e][1]/2,ROOTS[e][1]*(3**.5)/2)
        out[j]+=p*values[i];out[i]+=p.conjugate()*values[j]
    return out/g.events

def apply_S(g,amplitudes,k,extra=True):
    """Full source/event/reference action; no source/event coordinates dropped."""
    import numpy as np
    a=np.asarray(amplitudes,dtype=complex)
    if a.shape[:2]!=(g.states,g.events):raise ValueError('full source and event axes required')
    out=a.copy()
    for r,row in enumerate(g.seam_edges):
        for i,j,c in row:
            exponent=k*c%6 if extra else(3*k*(c%2))%6
            p=complex(ROOTS[exponent][0]+ROOTS[exponent][1]/2,ROOTS[exponent][1]*(3**.5)/2)
            out[j,2*r+1]=p*a[i,2*r+1];out[i,2*r+1]=p.conjugate()*a[j,2*r+1]
    return out

def cyclic_terminal(g,a,b,k,extra=True):
    import numpy as np
    a=np.asarray(a);b=np.asarray(b)
    return abs(a)**2+apply_T(g,abs(b)**2,0).real+2*np.real(a.conjugate()*apply_T(g,b,k,extra))

def distances(g,origin=0):
    if type(origin)is not int or not 0<=origin<g.states:raise ValueError('valid source origin required')
    d=[-1]*g.states;d[origin]=0;q=deque([origin])
    while q:
        i=q.popleft()
        for j in g.configuration_neighbours[i]:
            if d[j]<0:d[j]=d[i]+1;q.append(j)
    assert min(d)>=0
    return tuple(d)

def matching_code_projection(g,arc_state,cell=0):
    """Actual projection amplitudes for the fixed six matching words × events."""
    import numpy as np
    from model.common_source.record_code import MATCHINGS
    state=np.asarray(arc_state)
    if state.shape[0]!=g.states:raise ValueError('full two-load source required')
    groups=[[g.index[tuple(v+12*cell for v in p)]for p in row]for row in MATCHINGS]
    return np.stack([sum(state[i]for i in row)/(6**.5)for row in groups])

def interval(value,places=15):
    scale10=10**places;lo=value.numerator*scale10//value.denominator;hi=-((-value.numerator*scale10)//value.denominator)
    return[str(F(lo,scale10)),str(F(hi,scale10))]

def packet(g,k,radius,steps,extra=True):
    """Exact complete terminal law for a normalized configuration-ball packet."""
    if type(radius)is not int or radius<0 or not steps or any(type(n)is not int or n<0 for n in steps):raise ValueError('nonnegative integer preparation and times required')
    dist=distances(g);initial=[int(d<=radius)for d in dist];B=sum(initial);a=[(v,0)for v in initial];b=[ZERO]*g.states;power=1;records=[]
    adjacent={tuple(e)for e in g.seams};adj=[p in adjacent for p in g.pairs]
    first_cell=[max(p)<12 for p in g.pairs]
    for n in range(max(steps)+1):
        if n in steps:
            nums=terminal_numerators(g,a,b,k,extra);den=g.events*power*power*B
            assert min(nums)>=0 and sum(nums)==den
            tb=apply_T_num(g,b,k,extra);overlap=ZERO
            for i,v in enumerate(initial):
                if v:overlap=add(overlap,add(scale(g.events,a[i]),tb[i]))
            fidelity=F(norm2(overlap),g.events*g.events*power*power*B*B)
            digest=sha256()
            for v in nums:digest.update(str(v).encode()+b'\n')
            records.append({'steps':n,'complete_source_outcomes':g.states,'complete_vector_numerators_sha256':digest.hexdigest(),'common_denominator_sha256':sha256(str(den).encode()).hexdigest(),'total_probability_exact':'1','minimum_probability_nonnegative':True,'origin_probability_interval_exact':interval(F(nums[0],den)),'outside_preparation_ball_interval_exact':interval(F(sum(v for v,d in zip(nums,dist,strict=True)if d>radius),den)),'adjacent_pair_probability_interval_exact':interval(F(sum(v for v,s in zip(nums,adj,strict=True)if s),den)),'both_loads_in_first_cell_interval_exact':interval(F(sum(v for v,s in zip(nums,first_cell,strict=True)if s),den)),'mean_configuration_distance_interval_exact':interval(F(sum(v*d for v,d in zip(nums,dist,strict=True)),den)),'complete_arc_return_fidelity_interval_exact':interval(fidelity)})
        if n<max(steps):a,b=cyclic_step(g,a,b,k,extra);power*=g.events
    return {'character':k,'ports':g.ports,'source_states':g.states,'event_states':g.events,'configuration_edges':len(g.edges),'preparation_radius':radius,'preparation_support':B,'extra_order_three_voltage':extra,'records':records,'continuous_monitoring':False,'physical_lifetime_or_rest_mass':None}

def explicit_integer_step(g,psi,k,extra=True):
    """D*W on the full retained arc vector over the cyclotomic integers."""
    if len(psi)!=g.states or any(len(r)!=g.events for r in psi):raise ValueError('all source/event arcs required')
    reflected=[]
    for row in psi:
        total=ZERO
        for z in row:total=add(total,z)
        reflected.append([sub(scale(2,total),scale(g.events,z))for z in row])
    out=[row[:]for row in reflected]
    for r,row in enumerate(g.seam_edges):
        for i,j,c in row:
            e=k*(c if extra else 3*(c%2))%6
            out[j][2*r+1]=rotate(reflected[i][2*r+1],e);out[i][2*r+1]=rotate(reflected[j][2*r+1],-e)
    return out

def first_exit(g,k,radius=1,horizon=32,extra=True):
    """Actual absorbing first-exit instrument; no terminal-tail substitution."""
    if g.states>66:raise ValueError('explicit monitored control uses full local3960 arcs')
    d=distances(g);inside=[x<=radius for x in d];psi=[[(int(i==0),0)]*g.events for i in range(g.states)];power=1;exited=F();survival=F(1);mean_truncated=F();records=[];first_below_half=None;digest=sha256()
    for n in range(1,horizon+1):
        mean_truncated+=survival;psi=explicit_integer_step(g,psi,k,extra);power*=g.events;den=g.events*power*power
        weights=[sum(norm2(z)for z in row)for row in psi]
        for i,w in enumerate(weights):
            if not inside[i]:
                p=F(w,den);exited+=p;digest.update(f'{n}:{i}:{p}\n'.encode());psi[i]=[ZERO]*g.events
        survival=F(sum(w for w,s in zip(weights,inside,strict=True)if s),den)
        assert survival+exited==1
        if first_below_half is None and survival<F(1,2):first_below_half=n
        if n in(1,2,4,8,16,horizon):records.append({'steps':n,'survival_interval_exact':interval(survival),'cumulative_exit_interval_exact':interval(exited),'complete_total_probability_exact':'1'})
    return {'character':k,'extra_order_three_voltage':extra,'region_configuration_radius':radius,'retained_arc_dimension':g.states*g.events,'horizon':horizon,'exit_labels':horizon*sum(not x for x in inside),'terminal_survival_labels':sum(inside),'zero_probability_slots_retained':True,'all_exit_node_and_time_probabilities_sha256':digest.hexdigest(),'first_monitored_survival_below_half':first_below_half,'expected_min_exit_time_horizon_interval_exact':interval(mean_truncated),'records':records,'infinite_time_tail_assumed_zero':False,'physical_time_unit':None}

def spectrum_census(g,k):
    """Complete arc census, including all modes outside the prepared cyclic span."""
    phase(k,0);N=g.states;D=g.events;q=len(g.edges);zero=int(k==0)
    cyclic=2*N-zero;plus=q-N+2*zero;minus=N*D-q-N;dispersive=2*(N-zero)
    assert min(plus,minus)>=0 and plus+minus+dispersive==N*D
    return {'character':k,'source_states':N,'event_states':D,'full_arc_dimension':N*D,'swap_negative_dimension':q,'cyclic_dimension':cyclic,'constant_cyclic_mode_dimension':zero,'fixed_plus_one_dimension':plus,'fixed_minus_one_dimension':minus,'dispersive_dimension':dispersive,'unprepared_fixed_plus_one_dimension':plus-zero,'unprepared_fixed_minus_one_dimension':minus,'unprepared_modes_assigned_as_particles':False}

def positive_sector_bounds(g):
    """Loop/path lower bound on every nonzero-character magnetic Laplacian."""
    loops=((g.index[(0,11)],3,3),(g.index[(0,1)],2,3));ecc={b:max(distances(g,b))for b,_,_ in loops};out=[]
    for k in range(1,6):
        b,v,ell=loops[0 if k%2 else 1];c2=4 if k%2 else 3
        lower=1/(2*g.states*(F(ell,c2)+ecc[b]))
        out.append({'character':k,'loop_voltage':v,'loop_length':ell,'loop_base_eccentricity':ecc[b],'laplacian_minimum_lower_exact':str(lower),'minus_log_transfer_minimum_lower_exact':str(lower/g.events),'prepared_cyclic_phase_squared_lower_exact':str(2*lower/g.events)})
    return {'transfer_floor_exact':str(1-F(2*max(g.degree),g.events)),'nonzero_characters':out,'uniform_positive_continuum_gap_derived':False}

def full_arc_control():
    """Full retained reference operator controls, with exact structural proof."""
    import numpy as np
    g=carrier();N,D=g.states,g.events;eye=np.eye(N);checks=[]
    for k in range(6):
        arcs=np.repeat(eye[:,None,:],D,axis=1)/np.sqrt(D);a=eye.astype(complex);b=np.zeros((N,N),complex);total=0
        for step in range(3):
            effects=np.einsum('xej,xek->xjk',arcs.conjugate(),arcs,optimize=True)
            tb=apply_T(g,b,k);cross_b=b.conjugate()[:,:,None]*b[:,None,:]
            predicted=a.conjugate()[:,:,None]*a[:,None,:]+apply_T(g,cross_b,0)+a.conjugate()[:,:,None]*tb[:,None,:]+tb.conjugate()[:,:,None]*a[:,None,:]
            assert np.max(abs(effects-predicted))<2e-12 and np.max(abs(effects.sum(axis=0)-eye))<2e-12
            total+=N**3
            if step<2:
                arcs=apply_S(g,2*arcs.mean(axis=1,keepdims=True)-arcs,k);a,b=-b,a+2*tb
        checks.append({'character':k,'source_basis_columns':N,'retained_arcs':N*D,'terminal_matrix_unit_values_checked':total,'actual_numeric_residual_below_guard':True,'guard':2e-12,'guard_is_physical_implementation_error_certificate':False})
    # One exact direct-state control also includes complex cyclotomic input.
    exact=[]
    for k in range(6):
        a=[(int(i==0),int(i==1))for i in range(N)];b=[ZERO]*N;B=sum(map(norm2,a));psi=[[v]*D for v in a];power=1
        for step in range(4):
            expected=terminal_numerators(g,a,b,k);actual=[sum(norm2(z)for z in row)for row in psi]
            assert actual==expected and sum(actual)==D*power*power*B
            if step<3:psi=explicit_integer_step(g,psi,k);a,b=cyclic_step(g,a,b,k);power*=D
        exact.append({'character':k,'four_complete_exact3960_arc_comparisons':True})
    return {'floating_complete_operator_controls':checks,'exact_complex_state_controls':exact,'arbitrary_reference_identity':'S is a phase permutation; R=2AA*-I is a reflection; recurrence and terminal identities hold as operator identities.'}

def matching_code_numerators(g,a,b,k,cell=0,extra=True):
    """Six exact fixed-code weights, denominator6*D times input denominator."""
    from model.common_source.record_code import MATCHINGS
    groups=[[g.index[tuple(v+12*cell for v in p)]for p in row]for row in MATCHINGS];tb=apply_T_num(g,b,k,extra);out=[]
    for row in groups:
        sa=ZERO;sb=ZERO;st=ZERO
        for i in row:sa=add(sa,a[i]);sb=add(sb,b[i]);st=add(st,tb[i])
        square=len(g.seams)*norm2(sb)
        for e in range(len(g.seams)):
            total=ZERO
            for i in row:
                j,_,p=event_swap_column(g,i,2*e+1,k,extra)
                total=add(total,mul(conj(p),b[j]))
            square+=norm2(total)
        out.append(g.events*norm2(sa)+square+real2(mul(conj(sa),st)))
    return out

def matching_retention(k,steps=(0,1,4,16),extra=True):
    from model.common_source.record_code import MATCHINGS
    g=carrier();initial={g.index[p]for p in MATCHINGS[0]};a=[(int(i in initial),0)for i in range(g.states)];b=[ZERO]*g.states;power=1;rows=[]
    for n in range(max(steps)+1):
        if n in steps:
            nums=matching_code_numerators(g,a,b,k,extra=extra);den=36*g.events*power*power
            assert min(nums)>=0 and sum(nums)<=den
            rows.append({'steps':n,'all_six_code_outcomes_interval_exact':[interval(F(v,den))for v in nums],'code_population_exact':str(F(sum(nums),den)),'code_population_interval_exact':interval(F(sum(nums),den)),'source_code_failure_exact':str(1-F(sum(nums),den)),'source_code_failure_interval_exact':interval(1-F(sum(nums),den)),'total_with_failure_exact':'1'})
        if n<max(steps):a,b=cyclic_step(g,a,b,k,extra);power*=g.events
    return {'character':k,'initial_matching_word':0,'source_code_recovery_applied':False,'all_source_event_complements_retained':True,'records':rows}

def contact_identity():
    """All4356 entries of the native parity/contact decomposition."""
    g=carrier();one=[[0]*12 for _ in range(12)];K=[[0]*g.states for _ in range(g.states)]
    for a,b in g.seams:one[a][a]+=1;one[b][b]+=1;one[a][b]-=1;one[b][a]-=1
    for i,j,c in g.edges:
        K[i][i]+=1;K[j][j]+=1;s=(-1)**c;K[i][j]-=s;K[j][i]-=s
    adjacent=set(g.seams);entries=0
    for i,(a,b)in enumerate(g.pairs):
        for j,(c,d)in enumerate(g.pairs):
            free=one[a][c]*int(b==d)-one[a][d]*int(b==c)-one[b][c]*int(a==d)+one[b][d]*int(a==c)
            assert K[i][j]==free-2*int(i==j and(a,b)in adjacent);entries+=1
    changes={}
    for k,base in((1,3),(2,0)):
        row=[]
        for i,j,c in g.edges:
            delta=sub(scale(-1,phase(k,c)),scale(-1,phase(base,c)))
            if delta!=ZERO:row.append({'from_pair':list(g.pairs[i]),'to_pair':list(g.pairs[j]),'forward_matrix_difference_Zzeta':list(delta),'reverse_matrix_difference_Zzeta':list(conj(delta))})
        assert len(row)==1 and row[0]['from_pair']==[0,1]and row[0]['to_pair']==[1,5]
        changes[str(k)]=row
    return {'identity':'L3=dGamma(L)-2P_adj','all_matrix_entries_exact':entries,'contact_rank':len(adjacent),'native_double_occupancy_event':'identity, not the exterior-swap minus sign','free_exterior_identification':False,'added_interaction_parameter':None,'rank_two_order_three_voltage_updates':changes,'physical_strong_interaction_identified':False}

def validate_loops(g):
    result=[]
    for base,target in((g.index[(0,11)],3),(g.index[(0,1)],2)):
        witness=None
        for b in g.configuration_neighbours[base]:
            for c in g.configuration_neighbours[b]:
                if c!=base and base in g.configuration_neighbours[c]:
                    path=(base,b,c,base);total=0
                    for i,j in zip(path,path[1:]):
                        seam=tuple(sorted(set(g.pairs[i])^set(g.pairs[j])));assert seam in g.seams;total+=voltage(g.pairs[i],g.pairs[j],seam)
                    if total%6==target:witness=path;break
            if witness:break
        assert witness is not None
        result.append({'voltage':target,'configuration_path':[list(g.pairs[i])for i in witness],'length':3})
    return result

def _conjugate_row(row,k,source):
    # For real initial coefficients, all6th-root signs conjugate exactly.
    out=json.loads(json.dumps(row));out['character']=k;out['execution']='exact conjugacy transfer of the complete computed outcome vector';out['conjugate_of_character']=source;return out

def run():
    local=carrier();large=carrier(0);controls=full_arc_control();packets=[]
    for g,steps in((local,(0,1,4,16,64)),(large,(0,1,18,36,72))):
        rows={k:packet(g,k,0,steps)for k in range(4)}
        for row in rows.values():row['execution']='direct exact cyclotomic-integer computation'
        rows[4]=_conjugate_row(rows[2],4,2);rows[5]=_conjugate_row(rows[1],5,1);packets.extend(rows[k]for k in range(6))
    held_width=packet(local,1,1,(0,1,4,16,64));phase_control=packet(local,1,0,(0,1,4,16,64),False)
    monitored=[first_exit(local,k,horizon=32)for k in range(6)]
    monitored_control=first_exit(local,1,horizon=32,extra=False)
    retention=[matching_retention(k)for k in range(6)]
    assert any(a['origin_probability_interval_exact']!=b['origin_probability_interval_exact']for a,b in zip(packets[1]['records'],phase_control['records'],strict=True))
    return {'verdict':'PASS_COMPLETE_NATIVE_CHARGED_PROPAGATION_AND_LOCALIZATION_TESTS','construction':architecture.construction_descriptor(),'program':'retained coherent native event, fixed active memory voltage and uniform-event reflection','cyclotomic_basis':'Z[zeta], zeta=exp(i*pi/3), zeta^2=zeta-1','source_voltage_is_unchanged':True,'original_record_commit_policy':'inverse coherent recoding before commitment of all25 native labels; this is not unread C/E repetition','exact_kernel':'T_k=I-L_k/(2M), S_k^2=I, W_k=S_k(2AA*-I)','exact_terminal_probability':'|a_x|^2+(T_0|b|^2)_x+2Re(conj(a_x)*(T_k b)_x)','complete_source_packet_laws':packets,'held_out_preparation_width':held_width,'fixed_added_voltage_off_control':phase_control,'full_arc_controls':controls,'full_arc_complement_census':[spectrum_census(g,k)for g in(local,large)for k in range(6)],'exact_loop_controls':[{'ports':g.ports,'loops':validate_loops(g)}for g in(local,large)],'positive_nonzero_sector_bounds':[{'ports':g.ports,**positive_sector_bounds(g)}for g in(local,large)],'native_contact_identity':contact_identity(),'compact_localized_apparatus_modes':compact_modes(),'localization_protocol':{'preparation':'one adjacent pair at(0,1), no energy-band postselection','continuous_first_exit_instruments':monitored,'fixed_voltage_off_first_exit_control':monitored_control,'fixed_matching_code_retention':retention,'terminal_probability_is_not_first_exit_survival':True,'stable_localized_particle_identified':False,'physical_lifetime':None,'physical_rest_mass_GeV':None},'scope':{'complete_two_load_federation_retained':True,'all_six_charge_sectors_certified':True,'four_and_five_use_exact_conjugacy_for_real_preparations':True,'coherent_phase_character_preparation_and_reflection_declared':True,'enlarged_full_reflection_compiled':False,'physical_electric_unit_derived':False,'uniform_positive_continuum_gap':False,'full_arc_fixed_zero_phase_modes_retained':True,'physical_particle_spectrum_identified':False,'independent_GeV_mass_prediction':False,'OpenWave_physical_criterion_awarded':False}}

def _sparse_S(g,state,k,extra=True):
    out={}
    for(i,e),z in state.items():
        j,ee,p=event_swap_column(g,i,e,k,extra);assert(j,ee)not in out;out[j,ee]=mul(p,z)
    return out

def _compact_eigencheck(g,state,k,eigenvalue,extra=True):
    sums={}
    for(i,_),z in state.items():sums[i]=add(sums.get(i,ZERO),z)
    assert all(z==ZERO for z in sums.values())
    # A*psi=0, so the reflection is exactly -I on this sparse vector.
    actual=_sparse_S(g,{key:scale(-1,z)for key,z in state.items()},k,extra)
    assert actual=={key:scale(eigenvalue,z)for key,z in state.items()}
    return sum(norm2(z)for z in state.values())

def _flat_triangle(g):
    for a in range(g.states):
        for b in g.configuration_neighbours[a]:
            for c in g.configuration_neighbours[b]:
                if c!=a and a in g.configuration_neighbours[c]:
                    path=(a,b,c,a);edges=[]
                    for i,j in zip(path,path[1:]):
                        seam=tuple(sorted(set(g.pairs[i])^set(g.pairs[j])));r=g.seams.index(seam);chi=voltage(g.pairs[i],g.pairs[j],seam);chi0=voltage(g.pairs[i],g.pairs[j],seam,False);edges.append((i,j,r,chi,chi0))
                    if sum(e[3]for e in edges)%6==0 and sum(e[4]for e in edges)%6==0:return edges
    raise ArithmeticError('a flat triangle is required')

def compact_modes():
    """Exact compact apparatus modes, distinct from uniform-event packets."""
    rows=[]
    for g in(carrier(),carrier(0)):
        elementary=0
        # Explicitly check every source pair and every character for both
        # actual and voltage-off policies, with all zero coordinates implicit.
        for i in range(g.states):
            state={(i,0):(1,0),(i,2):(-1,0)}
            for k in range(6):
                for extra in(True,False):assert _compact_eigencheck(g,state,k,-1,extra)==2;elementary+=1
        # Each moving edge also supplies an independent four-arc eigenvector.
        edge_checks=0
        for r,edges in enumerate(g.seam_edges):
            for i,j,c in edges:
                for k in range(6):
                    p=phase(k,c);state={(i,2*r+1):(1,0),(j,2*r+1):p,(i,0):(-1,0),(j,0):scale(-1,p)}
                    assert _compact_eigencheck(g,state,k,-1)==4;edge_checks+=1
        holds=sum(g.events-d-1 for d in g.degree);full_minus=g.states*g.events-len(g.edges)-g.states
        assert holds+len(g.edges)==full_minus
        # A flat actual triangle produces a six-arc fixed +1 circulation.
        triangle=_flat_triangle(g);circulations=[]
        for k in range(6):
            state={};coefficient=(1,0)
            for i,j,r,c,_ in triangle:
                after=mul(coefficient,phase(k,c));state[i,2*r+1]=coefficient;state[j,2*r+1]=scale(-1,after);coefficient=after
            assert coefficient==(1,0)and _compact_eigencheck(g,state,k,1)==6
            # This selected flat triangle avoids the switched voltage edge.
            assert _compact_eigencheck(g,state,k,1,False)==6
            circulations.append({'character':k,'eigenvalue_exact':'1','norm_squared_before_normalizing_exact':'6','source_support':[list(g.pairs[i])for i,_,_,_,_ in triangle],'nonzero_event_arcs':6})
        rows.append({'ports':g.ports,'source_states':g.states,'full_event_carrier':g.states*g.events,'source_local_idle_difference_checks':elementary,'moving_edge_four_arc_checks':edge_checks,'two_idle_even_event_mode':{'event_labels':[0,2],'amplitudes_relative_exact':['1','-1'],'normalization_squared_exact':'2','W_eigenvalue_exact':'-1','uniform_event_preparation_overlap_exact':'0','single_source_position_survival_for_every_integer_time_exact':'1','holds_under_arbitrary_moving_branch_voltage':True},'complete_minus_one_band_compact_basis':{'source_local_hold_difference_dimension':holds,'independent_edge_mode_dimension':len(g.edges),'sum_equals_entire_minus_one_eigenspace':full_minus,'largest_basis_support_arcs':4},'flat_triangle_plus_one_modes':circulations})
    return {'members':rows,'additional_preparation':'Coherent event-bit Hadamard followed by relative minus phase on event2; source and memory unchanged. This is distinct from uniform A.','uniform_A_packet_population_exact':'0','source_only_position_measurements_populate_idle_difference_modes':False,'exact_infinite_step_persistence_in_stated_program':True,'global_minus_phase_is_an_observed_clock_by_itself':False,'physical_particle_identity':False,'rest_mass_GeV':None,'additional_Nature_prediction':False}
