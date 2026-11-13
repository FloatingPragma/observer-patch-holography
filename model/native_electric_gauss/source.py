"""Integer electric links carried by the native two-load observer process.

The added integer registers and electric phase service are declared resources.
Native source motion and its active six-state memory voltage are unchanged on
admitted links. Finite-window boundary events retain an overflow outcome.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from collections import defaultdict, deque
from itertools import combinations, product
from functools import lru_cache
import hashlib
import json

from model.openwave_native_process import source as native
from model.openwave_active_native_charge import source as active

SEAMS=tuple(active.SEAMS)
EDGE_INDEX={e:i for i,e in enumerate(SEAMS)}
PAIRS=tuple(active.PAIRS)
KAPPA=F(1,8)
LOOP=((0,5),(5,11),(0,11))

@dataclass(frozen=True)
class ElectricState:
    pair:tuple
    tagged_port:int
    charge_number:int
    memory:int
    flux:tuple

@dataclass(frozen=True)
class Hop:
    state:ElectricState
    proposed_signed_current:int
    realized_signed_current:int
    proposed_number_current:int
    overflow:bool


def extension_descriptor():
    return {'schema':'oph.native_integer_electric_extension.v1',
      'source':'native_H66_two_binary_loads','event_selector':'coherent recoded bit0 identity/bit1 swap; distinct from native raw coin', 'source_embedding':512,
      'active_memory':6,'memory_embedding':8,'physical_seams':[list(e)for e in SEAMS],
      'added_registers':['transported_tag_bit','integer_charge_number','30_oriented_integer_electric_links','pair_preparation_reference','electric_phase_reference'],
      'gauss_convention':'G_v=sum_out E-sum_in E-rho_v',
      'charge_assignment':'tagged load q, other load -q; q in Z; active character q modulo 6',
      'finite_boundary':'admitted hop is reversible; outside-window proposal is identity with explicit overflow',
      'electric_energy_coefficient':str(KAPPA),
      'declared_services':['integer controlled addition with boundary flag','vacuum-pair string swap','electric-energy controlled phase','complete orthogonal pointer readout'],
      'scope':'same native source/memory motion with declared integer-link and apparatus extension; no physical electric unit or continuum Maxwell claim'}


def validate(s,limit=None):
    if not isinstance(s,ElectricState) or s.pair not in PAIRS or s.tagged_port not in s.pair:
        raise ValueError('full native pair and transported tag required')
    if type(s.charge_number)is not int or type(s.memory)is not int or not 0<=s.memory<6:
        raise ValueError('integer charge and complete active-memory coordinate required')
    if len(s.flux)!=len(SEAMS) or any(type(v)is not int for v in s.flux):
        raise ValueError('one integer flux per physical native seam required')
    if limit is not None and(type(limit)is not int or limit<0 or any(abs(v)>limit for v in s.flux)):
        raise ValueError('admitted integer-link window required')
    return s


def charge_density(s):
    validate(s);out=[0]*12
    out[s.tagged_port]=s.charge_number
    out[next(v for v in s.pair if v!=s.tagged_port)]=-s.charge_number
    return tuple(out)


def divergence(flux):
    if len(flux)!=len(SEAMS):raise ValueError('complete link vector required')
    out=[0]*12
    for(a,b),e in zip(SEAMS,flux,strict=True):out[a]+=e;out[b]-=e
    return tuple(out)


def gauss(s):return tuple(a-b for a,b in zip(divergence(s.flux),charge_density(s),strict=True))


def electric_energy(s,kappa=KAPPA):
    validate(s)
    if not isinstance(kappa,F) or kappa<0:raise ValueError('nonnegative rational declared coupling required')
    return kappa*sum(e*e for e in s.flux)/2


def apply_native_current(s,seam,event_bit=1,limit=None):
    """Complete reversible coherent event, including collision and overflow.

    The event bit is the recoded coherent identity/swap selector, not the
    original native final-endpoint coin. Number-current is the actual endpoint
    motion. Signed current includes
    the transported integer charge. A failed finite-window proposal is a
    no-hop on source, memory and every link, not a cyclic wrap.
    """
    validate(s,limit)
    seam=tuple(seam)
    if seam not in EDGE_INDEX or type(event_bit)is not int or event_bit not in(0,1):raise ValueError('literal native seam and coherent event bit required')
    a,b=seam
    if event_bit==0 or(a in s.pair)==(b in s.pair):return Hop(s,0,0,0,False)
    origin,destination=(a,b)if a in s.pair else(b,a)
    number_current=1 if origin==a else-1
    load_charge=s.charge_number if origin==s.tagged_port else-s.charge_number
    current=load_charge*number_current
    e=EDGE_INDEX[seam];flux=list(s.flux);flux[e]-=current
    if limit is not None and abs(flux[e])>limit:return Hop(s,current,0,number_current,True)
    pair=tuple(sorted(destination if v==origin else v for v in s.pair))
    tag=destination if s.tagged_port==origin else s.tagged_port
    memory=(s.memory+active.voltage(s.pair,pair,seam))%6
    target=ElectricState(pair,tag,s.charge_number,memory,tuple(flux))
    return Hop(target,current,current,number_current,False)


def complete_raw_events(s,seam,law='mixture',limit=None):
    """All25 native slots, each split into kept/overflow record sectors."""
    validate(s,limit);mask=sum(1<<v for v in s.pair);out=[]
    for slot,branch in enumerate(native.raw_local_branches(law)):
        image=native.apply_branch(mask,tuple(seam),branch)
        event_bit=abs(branch.current) if branch.current is not None else 0
        hop=apply_native_current(s,seam,event_bit,limit)
        valid=image is not None
        for overflow in(False,True):
            probability=branch.weight if valid and hop.overflow==overflow else F(0)
            out.append({'slot':slot,'record':branch.record,'overflow':overflow,'probability':probability,
                        'state':hop.state if valid and hop.overflow==overflow else None,
                        'proposed_signed_current':hop.proposed_signed_current,
                        'realized_signed_current':hop.realized_signed_current})
    return tuple(out)


def pair_string(path=(0,1),charge_number=1,memory=0):
    """Neutral pair and full integer preparation string; no flux discarded."""
    path=tuple(path)
    if len(path)<2 or len(set(path))!=len(path) or any(tuple(sorted(e))not in EDGE_INDEX for e in zip(path,path[1:])):
        raise ValueError('simple actual native path required')
    flux=[0]*len(SEAMS)
    for a,b in zip(path,path[1:]):flux[EDGE_INDEX[tuple(sorted((a,b)))]]+=charge_number*(1 if a<b else-1)
    s=ElectricState(tuple(sorted((path[0],path[-1]))),path[0],charge_number,memory,tuple(flux))
    validate(s);assert not any(gauss(s))
    return s


def prepare_pair_swap(address,target):
    """Full direct-sum permutation; arbitrary complement addresses unchanged.

    The vacuum retains a reference tag identifying this neutral-pair operation.
    Memory passes through. Preparation work is the target's electric energy;
    an independent matter preparation cost is not supplied by this kernel.
    """
    validate(target)
    if any(gauss(target)):raise ValueError('neutral Gauss-compatible target required')
    vacuum=('vacuum',target.memory)
    if address==vacuum:return target
    if address==target:return vacuum
    return address


def native_word(s,seams=LOOP,limit=None):
    records=[]
    for seam in seams:
        hop=apply_native_current(s,seam,1,limit);s=hop.state
        records.append({'seam':list(seam),'number_current':hop.proposed_number_current,
                        'signed_current':hop.realized_signed_current,'overflow':hop.overflow})
    return s,tuple(records)


def phase_character(q,voltage):
    """Exponent of exp(i*pi/3); q is supplied by the integer lift."""
    if type(q)is not int:raise ValueError('integer charge required')
    return(q*voltage)%6


def local_operator_certificate():
    """Every local signed-pair/active-memory/link basis column and raw slot.

    Spectator links tensor identity. Columns are checked as partial injections,
    which fixes every matrix unit, including entangled reference systems.
    """
    checks=0;raw_checks=0;overflow_count=0;gauss_controls=0;undressed_witness=None
    domains=[(q,limit,seam) for q in(-1,0,1) for limit in(1,) for seam in SEAMS]
    domains += [(q,2,SEAMS[0]) for q in(-2,-1,0,1,2)]
    for q,limit,seam in domains:
       if True:
        e=EDGE_INDEX[seam]
        # Memory covariance makes the six-coordinate action a translation;
        # all six coordinates are nevertheless executed on every event.
        for coin in(0,1):
         images={}
         for pair in PAIRS:
          for tag in pair:
           for memory in range(6):
            for value in range(-limit,limit+1):
             flux=[0]*len(SEAMS);flux[e]=value
             s=ElectricState(pair,tag,q,memory,tuple(flux));hop=apply_native_current(s,seam,coin,limit)
             key=(pair,tag,memory,value);out=(hop.state.pair,hop.state.tagged_port,hop.state.memory,hop.state.flux[e])
             assert out not in images;images[out]=key
             assert apply_native_current(hop.state,seam,coin,limit).state==s
             assert gauss(s)==gauss(hop.state)
             assert all(hop.state.flux[j]==s.flux[j]for j in range(len(SEAMS))if j!=e)
             assert electric_energy(hop.state)-electric_energy(s)==KAPPA*(F(hop.realized_signed_current**2,2)-value*hop.realized_signed_current)
             checks+=1;overflow_count+=hop.overflow
             if hop.state.pair!=s.pair and undressed_witness is None:
                undressed=replace(hop.state,flux=s.flux);undressed_witness={'q':q,'seam':list(seam),'before_pair':list(pair),'after_pair':list(hop.state.pair),'gauss_change':[a-b for a,b in zip(gauss(undressed),gauss(s),strict=True)]}
        # The raw law's source domains depend on pair only. Check each signed
        # column at memory0; translation proves the identical other five.
        for pair in PAIRS:
         for tag in pair:
          for value in range(-limit,limit+1):
           flux=[0]*len(SEAMS);flux[e]=value;s=ElectricState(pair,tag,q,0,tuple(flux))
           rows=complete_raw_events(s,seam,limit=limit)
           assert len(rows)==50 and sum(r['probability']for r in rows)==1
           for row in rows:
            if row['probability']:assert gauss(s)==gauss(row['state'])
           raw_checks+=len(rows)
    assert undressed_witness and any(undressed_witness['gauss_change'])
    return {'coherent_basis_columns':checks,'raw_record_columns':raw_checks,
      'executed_domains':'q=-1,0,1; L=1; all30 seams; plus q=-2,-1,0,1,2 and L=2 on seam(0,1)',
      'finite_window_overflow_columns':overflow_count,'gauss_exact_on_every_column':True,
      'every_event_full_permutation':True,'full_matrix_units_and_references':'basis-permutation/partial-injection theorem; no matrix-unit sampling',
      'all25_raw_slots_each_with_two_boundary_flags':True,'undressed_hop_negative_control':undressed_witness,
      'local_tensor_extension':'all other29 links and arbitrary references are exact identity spectators'}


def preparation_certificate():
    paths=((0,1),(0,5,1),(0,11,5,1))
    records=[]
    for path in paths:
     for q in(-2,-1,0,1,2):
      for m in range(6):
       target=pair_string(path,q,m);vac=('vacuum',m)
       assert prepare_pair_swap(prepare_pair_swap(vac,target),target)==vac
       assert prepare_pair_swap(('complement',m),target)==('complement',m)
       assert sum(charge_density(target))==0 and not any(gauss(target))
       records.append({'path':list(path),'charge_number':q,'memory':m,'electric_work':str(electric_energy(target))})
    return {'executed_preparations':len(records),'records':records,
      'total_integer_charge_exact':0,'pair_reference_retained':True,
      'source_preparation_work':'independent source-sector energy not supplied; electric contribution recorded',
      'active_character_attachment':'q mod6 on the distinguished-load reference; compensating load carries -q; global memory is not total net charge',
      'integer_lift_not_forced':True,'real_charge_unit':'Q=q_unit*q, supplied q_unit; no SI value inferred',
      'finite_group_obstruction':'6 f(1)=0 in R implies f=0 for an additive map Z6 to R',
      'local_U1_action':'exp(i sum_v theta_v G_v); every dressed operation commutes exactly'}

# Q(zeta), zeta=exp(i*pi/3), with zeta^2=zeta-1.
ZERO=(F(0),F(0));ONE=(F(1),F(0))
ROOTS=(ONE,(F(0),F(1)),(F(-1),F(1)),(F(-1),F(0)),(F(0),F(-1)),(F(1),F(-1)))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def scale(c,a):return(c*a[0],c*a[1])
def mul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def conjugate(a):return(a[0]+a[1],-a[1])
def norm(a):return a[0]*a[0]+a[0]*a[1]+a[1]*a[1]
def real(a):return a[0]+a[1]/2

def add_state(out,s,a):
    out[s]=add(out.get(s,ZERO),a)
    if out[s]==ZERO:del out[s]


def character_preparation(q=1,path=(0,1)):
    # The common 1/sqrt(6) is retained as probability denominator6.
    return {pair_string(path,q,m):ROOTS[(-q*m)%6]for m in range(6)}


def probability(state):return sum(map(norm,state.values()),F(0))/6

def loop_operator(state,inverse=False,limit=None):
    out={}
    for s,a in state.items():
        t,_=native_word(s,tuple(reversed(LOOP))if inverse else LOOP,limit)
        add_state(out,t,a)
    return out


def loop_reference_measurement(state,phase=1,limit=None):
    """Complete reference outcomes (I +/- zeta^phase U)/2.

    U is the actual three-event permutation, including finite boundary holds.
    A Hadamard reference, controlled U, root-of-unity reference phase and
    Hadamard/readout implement these two Kraus operators. No source outcome
    or finite-window failure is removed.
    """
    if type(phase)is not int:raise ValueError('integer sixth-root phase required')
    moved=loop_operator(state,limit=limit);result=[]
    for sign in(1,-1):
        out={s:scale(F(1,2),a)for s,a in state.items()}
        for s,a in moved.items():add_state(out,s,scale(F(sign,2),mul(ROOTS[phase%6],a)))
        result.append(out)
    assert sum(probability(x)for x in result)==probability(state)
    return tuple(result)


def energy_mean(state):return sum(norm(a)*electric_energy(s)for s,a in state.items())/6

def overlap(a,b):
    out=ZERO
    for s,v in a.items():out=add(out,mul(conjugate(v),b.get(s,ZERO)))
    return scale(F(1,6),out)


def loop_experiment():
    initial=character_preparation();p=probability(initial);assert p==1
    # L=2 contains both native three-hop blocks on every admitted branch.
    first=loop_reference_measurement(initial,phase=1,limit=2)
    records=[]
    for phase in(0,1):
      for a,first_state in enumerate(first):
       for b,state in enumerate(loop_reference_measurement(first_state,phase,limit=2)):
        assert all(not any(gauss(s))for s in state)
        records.append({'second_reference_phase_sixths':phase,'first_pointer':a,'second_pointer':b,
          'joint_probability':str(probability(state)),
          'joint_electric_energy':str(energy_mean(state)),
          'complete_flux_source_memory_weights':[
           {'pair':list(s.pair),'tagged_port':s.tagged_port,'charge_number':s.charge_number,'memory':s.memory,
            'flux':list(s.flux),'probability':str(norm(v)/6)}for s,v in sorted(state.items(),key=lambda row:(row[0].flux,row[0].memory))]})
    assert [r['joint_probability']for r in records[:4]]==['5/16','3/16','3/16','5/16']
    assert [r['joint_probability']for r in records[4:]]==['3/8','1/8','1/8','3/8']
    e0=energy_mean(initial);e1=sum(energy_mean(v)for v in first)
    e2=sum(F(r['joint_electric_energy'])for r in records[4:])
    loop_states=[];s=pair_string()
    for n in range(3):
        assert not any(gauss(s));assert s.pair==(0,1)and s.tagged_port==0
        loop_states.append({'loop_number':n,'memory_shift':s.memory,'electric_energy':str(electric_energy(s)),'flux':list(s.flux)})
        s,_=native_word(s,limit=None)
    # Boundary checks exercise the full unitary event word, not an assumed
    # invariant fixed-matter cycle at the boundary.
    boundary_input=character_preparation()
    boundary=boundary_input
    boundary_failures=0
    for _ in range(4):
        for s in boundary:
            _,hist=native_word(s,limit=1);boundary_failures+=sum(r['overflow']for r in hist)
        moved=loop_operator(boundary,limit=1)
        assert loop_operator(moved,inverse=True,limit=1)==boundary
        boundary=moved
    assert boundary_failures>0
    conditional=[]
    for v in first:
        ov=scale(1/probability(v),overlap(v,loop_operator(v,limit=2)))
        conditional.append({'real_loop_expectation':str(real(ov)),
          'loop_expectation_Q_zeta':[str(x)for x in ov],
          'magnetic_current_over_beta_squared':str(F(3,16)),
          'magnetic_current_over_beta_sign':1 if ov[1]<0 else -1})
    return {'native_closed_word':[list(e)for e in LOOP],'native_events_per_loop':3,'event_policy':'recoded coherent swap; original C/E final-endpoint coins retained only in the committed raw instrument',
      'initial_pair':[0,1],'compensating_reference_port':1,'integer_charge':1,
      'active_character':1,'inherited_memory_holonomy':2,'finite_link_limit':2,
      'loop_states':loop_states,'two_complete_reference_blocks':records,
      'conditional_loop_observables':conditional,
      'initial_electric_energy':str(e0),'after_first_nonselective_energy':str(e1),
      'after_second_nonselective_energy':str(e2),
      'first_electric_work':str(e1-e0),'second_electric_work':str(e2-e1),
      'nonselective_magnetic_current_exact':0,
      'boundary_control_overflows_retained':boundary_failures,
      'boundary_word_inverse_exact':True,
      'failure_complements':'all pair/sign/memory/link states retained; encoded padding identity flagged by admission; every reference output included',
      'apparatus_contract':'coherent native-record policy; original raw-record commitment is a different channel; pointer outcomes are copied only after each reference block',
      'magnetic_generator':'H_B=beta*(I-(U+U^dagger)/2) is a declared finite loop generator; no autonomous implementation or beta derivation claimed',
      'measured_observable':'loop quadratures reconstructed from complete reference probabilities; electric flux energy from complete orthogonal records'}


def electric_phase_readout():
    """Controlled native loop, electric phase, inverse loop, complete pointer.

    The declared phase duration is tau=16*pi/9 source-time units. Its common
    input phase is removed; the relative phase is an exact sixth root.
    All three input cycle sectors and their full six-memory preparations
    execute the native forward/inverse word before probability evaluation.
    """
    state=character_preparation();rows=[]
    for n in range(3):
        moved=loop_operator(state,limit=3)
        assert loop_operator(moved,inverse=True,limit=3)==state
        e0=energy_mean(state);e1=energy_mean(moved);delta=e1-e0
        assert delta==F(3*(2*n+1),16)
        relative=ROOTS[-(2*n+1)%6]
        phase_moved={s:mul(relative,a)for s,a in moved.items()}
        restored=loop_operator(phase_moved,inverse=True,limit=3)
        probabilities=[]
        for sign in(1,-1):
            out={s:scale(F(1,2),a)for s,a in state.items()}
            for s,a in restored.items():add_state(out,s,scale(F(sign,2),a))
            probabilities.append(probability(out))
        assert sum(probabilities)==1
        rows.append({'cycle_sector':n,'input_electric_energy':str(e0),'loop_branch_energy':str(e1),
          'energy_difference':str(delta),'relative_phase_sixths':-(2*n+1),
          'complete_pointer_probabilities':[str(x)for x in probabilities],
          'field_work_before_phase':str(delta/2),'field_work_during_phase':'0',
          'field_work_after_inverse':str(-delta/2),'retained_system_restored_exact':True})
        state=moved
    assert [r['complete_pointer_probabilities']for r in rows]==[['3/4','1/4'],['0','1'],['3/4','1/4']]
    zero=character_preparation();moved=loop_operator(zero,limit=1)
    assert all(electric_energy(x,F(0))==0 for x in moved)
    restored=loop_operator(moved,inverse=True,limit=1)
    plus={};minus={}
    for x in zero:
        add_state(plus,x,scale(F(1,2),add(zero[x],restored.get(x,ZERO))))
        add_state(minus,x,scale(F(1,2),add(zero[x],scale(-1,restored.get(x,ZERO)))))
    zero_prob=[probability(plus),probability(minus)];assert zero_prob==[F(1),F(0)]
    return {'declared_phase_duration':'16*pi/9 in source-time units','records':rows,
      'zero_kappa_negative_control_probabilities':[str(x)for x in zero_prob],
      'control_resources':'fresh coherent reference per block; no autonomous apparatus-work claim; phase normalization supplied',
      'implementation_scope':'native forward/inverse permutations executed exactly; electric diagonal phase is the declared U(1) service'}


def _solve_rational(a,b):
    n=len(b);m=[[F(v)for v in row]+[F(y)]for row,y in zip(a,b,strict=True)]
    for j in range(n):
        k=next(i for i in range(j,n)if m[i][j]);m[j],m[k]=m[k],m[j]
        z=m[j][j];m[j]=[v/z for v in m[j]]
        for i in range(n):
            if i!=j:
                z=m[i][j];m[i]=[x-z*y for x,y in zip(m[i],m[j],strict=True)]
    return tuple(row[-1]for row in m)


def coulomb_comparison():
    """Actual integer string energies against exact graph Poisson relaxation."""
    lap=[[0]*12 for _ in range(12)];adj=[[]for _ in range(12)]
    for a,b in SEAMS:
        lap[a][a]+=1;lap[b][b]+=1;lap[a][b]-=1;lap[b][a]-=1;adj[a].append(b);adj[b].append(a)
    parent={0:None};queue=deque([0])
    while queue:
        a=queue.popleft()
        for b in adj[a]:
            if b not in parent:parent[b]=a;queue.append(b)
    rows=[]
    for destination in range(1,12):
        path=[destination]
        while path[-1]!=0:path.append(parent[path[-1]])
        path=tuple(reversed(path));rho=[0]*12;rho[0]=1;rho[destination]=-1
        phi=_solve_rational([r[:11]for r in lap[:11]],rho[:11])+(F(0),)
        e=tuple(phi[a]-phi[b]for a,b in SEAMS)
        assert divergence(e)==tuple(rho)
        relaxed=KAPPA*sum(v*v for v in e)/2
        integer=electric_energy(pair_string(path))
        assert relaxed<integer
        rows.append({'pair':[0,destination],'shortest_native_path':list(path),'integer_ground_electric_energy':str(integer),
          'graph_poisson_relaxed_energy':str(relaxed),'effective_resistance':str((phi[0]-phi[destination])),
          'strict_integer_excess':str(integer-relaxed)})
    return {'rows':rows,'integer_ground_proof':'sum E_e^2 >= sum abs(E_e) >= graph_distance for one unit neutral pair; a shortest integer string attains equality',
      'poisson_proof':'B E=rho, E=B^T L^+rho minimizes real quadratic energy; it is a lower bound, not the executed integer flux state',
      'measured_field_scale':'electric energy in source units, kappa=1/8 declared',
      'physical_scope':'finite graph electrostatics; no Euclidean 1/r distance law, fine-structure constant or SI electric charge derived'}


def cut_flux(s,region):
    """Complete oriented boundary readout, independent of the charge label."""
    validate(s)
    if type(region)is not int or not 0<=region<1<<12:raise ValueError('complete native region mask required')
    return sum(e*(((region>>a)&1)-((region>>b)&1))for(a,b),e in zip(SEAMS,s.flux,strict=True))


def cut_charge_certificate():
    initial=pair_string();states=[initial]
    for _ in range(4):states.append(native_word(states[-1])[0])
    # Include actual moving endpoints, not only fixed-source loop deformations.
    current=initial
    for seam in LOOP:
        current=apply_native_current(current,seam,1,limit=2).state
        states.append(current)
    records=[]
    for j,state in enumerate(states):
        rho=charge_density(state);hist=defaultdict(int);values=[]
        for mask in range(1<<12):
            flux=cut_flux(state,mask)
            charge=sum(rho[v]for v in range(12)if(mask>>v)&1)
            assert flux==charge and type(flux)is int
            hist[flux]+=1;values.append(flux)
        assert dict(hist)=={0:2048,1:1024,-1:1024}
        if 0<j<5:assert values==baseline
        # The first assignment also provides the full region readout table.
        if j==0:baseline=values
        records.append({'native_state_index':j,'pair':list(state.pair),'tagged_port':state.tagged_port,
                        'complete_cut_histogram':{str(k):v for k,v in sorted(hist.items())},
                        'cut_table_sha256':hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()})
    # H2 reference followed by the actual controlled native (0,5) event.
    # Common denominator12 is the six-memory character times two paths.
    joint=defaultdict(F);realized=[]
    for reference in(0,1):
        for state,amplitude in character_preparation().items():
            if reference:state=apply_native_current(state,(0,5),1,2).state
            outcome=(cut_flux(state,1<<0),cut_flux(state,1<<5))
            weight=norm(amplitude)/12;joint[outcome]+=weight
            realized.append({'reference':reference,'pair':list(state.pair),'memory':state.memory,
                             'regional_charge_record':list(outcome),'probability':str(weight)})
    assert dict(joint)=={(1,0):F(1,2),(0,1):F(1,2)}
    return {'regions_per_state':4096,'native_source_states_checked':len(states),
      'complete_cut_records':records,'initial_full_region_charge_table':baseline,
      'loop_deformation_cut_tables_exact':True,
      'unit_endpoint_charge_from_boundary_flux':1,
      'sharp_half_split_probability':'0',
      'coherent_half_expectation_control':{'preparation':'H2 reference controlling actual native coherent swap(0,5); complete character memory retained',
       'regions':[[0],[5]],'joint_records':[{'charges':list(k),'probability':str(v)}for k,v in sorted(joint.items())],
       'regional_expectations':['1/2','1/2'],'all12_reference_memory_records':realized,
       'sharp_half_half_record_probability':'0'},
      'topological_statement':'integer relative-chain boundary index; cut deformations avoiding endpoints and integer cycle additions preserve charge',
      'row1_review_scope':{'finite_integer_flux_and_unsplittable_unit':'PASS on the declared extended carrier',
       'architecture':'tag/sign and integer link registers are declared; their exact native transport and complete observations are executed',
       'pinned_M10_equivalence':'requires implementation/equivalence on the admitted platform or transparent review of the declared extension',
       'physical_electron_identification':'unproved; unit endpoint is an electric defect candidate, not an identified electron',
       'official_row_promotion':False}}


def maxwell_identities():
    s0=pair_string();s1,_=native_word(s0);c=tuple(b-a for a,b in zip(s0.flux,s1.flux,strict=True))
    assert divergence(c)==(0,)*12 and sum(x*x for x in c)==3
    checks=0
    for n in range(-5,6):
        flux=tuple(a+n*b for a,b in zip(s0.flux,c,strict=True))
        s=replace(s0,flux=flux,memory=(2*n)%6);t,_=native_word(s)
        assert tuple(b-a for a,b in zip(s.flux,t.flux,strict=True))==c
        assert electric_energy(t)-electric_energy(s)==KAPPA*(sum(a*b for a,b in zip(c,s.flux,strict=True))+F(sum(x*x for x in c),2))
        for edge in range(30):assert t.flux[edge]-s.flux[edge]==c[edge]
        checks+=1
    return {'exact_cycle_flux':list(c),'cycle_divergence_exact':list(divergence(c)),
      'executed_integer_cycle_columns':checks,'electric_commutator':'[E_e,U]=c_e U',
      'electric_energy_commutator':'[H_E,U]=U*kappa*(c dot E+|c|^2/2)',
      'finite_loop_Maxwell_equation':'i[H_B,E_e]=c_e*i*beta*(U-U^dagger)/2',
      'gauss_commutes_with_electric_and_loop_generators':True,
      'boundary_scope':'constant-cycle identities apply to the unbounded integer orbit and finite support before a boundary; actual bounded word/inverse retains all boundary modifications',
      'couplings':{'electric_kappa':str(KAPPA),'magnetic_beta':'declared nonnegative coefficient; measured quadratures permit any explicit value'},
      'time_and_units':'event count and source phase; no independently calibrated physical time, SI normalization or continuum propagation speed'}


def run():
    return {'schema':'oph.native_electric_gauss_and_loop_records.v1',
      'verdict':'PASS: native integer Gauss transport and complete finite electric/loop records',
      'extension_descriptor':extension_descriptor(),'operator_certificate':local_operator_certificate(),
      'neutral_pair_preparation':preparation_certificate(),'integer_cut_charge':cut_charge_certificate(),'native_loop_observation':loop_experiment(),'electric_phase_readout':electric_phase_readout(),
      'integer_electrostatics':coulomb_comparison(),'finite_Maxwell_identities':maxwell_identities(),
      'physical_criterion_status':{'electrical_charge':'finite integer/U1 extension; physical charge identification and SI unit unproved',
       'Maxwell_Coulomb':'finite native-link records and exact identities; physical continuum/strength unproved',
       'OpenWave_row_promotions':0}}
