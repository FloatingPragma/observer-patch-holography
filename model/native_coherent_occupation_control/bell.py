"""Exact native boundary parity resources with complete local-number accounting."""
from fractions import Fraction as F
import numpy as np
from model.native_program import interpreter as native
from model.federated_observer import source
from model.native_fair_local_preparation import coupling


def layout(level=0):
    if type(level)is not int or not 0<=level<=2:
        raise ValueError('bounded source geometry level required')
    graph=source.federation(level)
    edge=next(e for e in graph.seams if e[0]//12!=e[1]//12)
    spare=tuple(next(12*(a//12)+i for i in range(12) if 12*(a//12)+i!=a) for a in edge)
    return graph,edge,spare


def code_masks(edge,spare):
    if len(set((*edge,*spare)))!=4:
        raise ValueError('two disjoint dual rails required')
    return tuple((1<<(edge[0] if a else spare[0]))|(1<<(edge[1] if b else spare[1]))
                 for a,b in ((0,0),(0,1),(1,0),(1,1)))


def branch_run(level=0,choi=False):
    graph,edge,spare=layout(level);masks=code_masks(edge,spare)
    command={'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':'E','repeat':1}
    program=native.program(federation={'level':level,'total':2},biased=False,
                           commands=[command],capacity=1,overflow='export_oldest')
    density={(x,i if choi else 0,y,j if choi else 0):F(1,4)
             for i,x in enumerate(masks) for j,y in enumerate(masks)}
    native._validate_density(density,program['carrier'])
    rows=native._outcomes(program,command,density)
    if len(rows)!=25:raise ArithmeticError('complete original source address count changed')
    trace=lambda rho:sum(v for (x,a,y,b),v in rho.items() if x==y and a==b)
    accepted=[];nonzero=[]
    for row in rows:
        if row.density:nonzero.append(row)
        if row.record.law=='E' and row.record.raw==0:
            expected={(masks[i],i if choi else 0,masks[j],j if choi else 0):F(1,8)
                      for i in (0,3) for j in (0,3)}
            if row.record.label!=0 or row.record.current!=0 or row.density!=expected:
                raise ArithmeticError('native E zero record is not the coherent even projector')
            accepted.append(row)
    if len(accepted)!=2 or sum(trace(row.density) for row in rows)!=1:
        raise ArithmeticError('the complete Bell instrument lost probability or a coin')
    if sum(trace(row.density) for row in accepted)!=F(1,2):
        raise ArithmeticError('exact Bell acceptance changed')
    if not choi:
        rejected=[r for r in nonzero if r not in accepted]
        if len(rejected)!=4 or any(len(r.density)!=1 or trace(r.density)!=F(1,8) for r in rejected):
            raise ArithmeticError('rejected resource states are not fully known basis branches')
        moved=sum(trace(r.density) for r in rejected if r.record.current!=0)
        if moved!=F(1,4):raise ArithmeticError('restoration exposure changed')
    return {'level':level,'ports':graph.ports,'edge':edge,'spares':spare,
            'complete_raw_slots':len(rows),'nonzero_raw_slots':len(nonzero),
            'Bell_success_mass':F(1,2),'per_coin_success_mass':F(1,4),
            'arbitrary_encoded_reference_checked':choi,
            'source_resources':program['declared_resources']},rows


def restore_basis(mask,edge,target_first_cell_number,attempts):
    """Exact finite endpoint restoration, retaining failed terminal branches."""
    if type(attempts)is not int or not 0<=attempts<=32:
        raise ValueError('bounded native restoration count required')
    if (type(mask)is not int or mask<0 or type(edge)is not tuple or len(edge)!=2
            or any(type(x)is not int or x<0 for x in edge) or edge[0]>=edge[1]
            or edge[0]//12==edge[1]//12 or type(target_first_cell_number)is not int
            or not 0<=target_first_cell_number<=12):
        raise ValueError('an exact recorded boundary state and cell-total target required')
    firstcell=edge[0]//12
    count=lambda x:((x>>(12*firstcell))&4095).bit_count()
    if count(mask)==target_first_cell_number:
        return (((),mask,F(1),'RESTORED'),)
    if count(mask)-target_first_cell_number!=((mask>>edge[0])&1)-((mask>>edge[1])&1):
        raise ValueError('the boundary endpoint move cannot restore this cell total')
    paths=[((),mask,F(1),'ACTIVE')]
    for _ in range(attempts):
        out=[]
        for history,state,weight,status in paths:
            if status!='ACTIVE':out.append((history,state,weight,status));continue
            for coin in (0,1):
                a,b=edge;na=(state>>a)&1;nb=(state>>b)&1
                j=(na-nb-(2*coin-1)*abs(na-nb))//2
                image=state if not j else state^((1<<a)|(1<<b))
                out.append(((*history,(coin,j,na-nb)),image,weight/2,
                            'RESTORED' if count(image)==target_first_cell_number else 'ACTIVE'))
        paths=out
    return tuple((h,x,w,'RESTORATION_TIMEOUT' if s=='ACTIVE' else s) for h,x,w,s in paths)


def retry_budget(rounds,restoration,initial_tree=True):
    if (type(rounds)is not int or not 1<=rounds<=16 or type(restoration)is not int
            or not 0<=restoration<=32 or type(initial_tree)is not bool):
        raise ValueError('bounded positive rounds, restoration count and tree flag required')
    abort=F(1,2**(restoration+2));continuation=F(1,2)-abort
    success=F(1,2)*sum((continuation**j for j in range(rounds)),F())
    abort_mass=abort*sum((continuation**j for j in range(rounds-1)),F())
    exhausted=F(1,2)*continuation**(rounds-1)
    if success+abort_mass+exhausted!=1:raise ArithmeticError('complete bounded retry law lost mass')
    pulses=44*rounds;attempts=rounds+(rounds-1)*restoration+22*initial_tree
    return {'rounds':rounds,'restoration_attempt_ceiling':restoration,
            'Bell_success_mass_given_two_eligible_cells':success,
            'restoration_timeout_mass':abort_mass,'round_exhaustion_mass':exhausted,
            'maximum_native_attempts':attempts,'maximum_local_source_pulses':pulses,
            'maximum_native_record_slots':attempts+1,
            'reused_local_work_bits_per_cell':12,'simultaneous_source_cells':2,
            'round_counter_bits':rounds.bit_length(),'restoration_counter_bits':restoration.bit_length(),
            'complete_local_service_error':coupling.complete_bound(pulses),
            'ideal_source_return_interface_required':True,
            'monolithic_native_grammar_execution_claimed':False}


def service_words():
    _,edge,spare=layout();beam=coupling.words()['beam'][1];out=[]
    for endpoint,other in zip(edge,spare,strict=True):
        positions=(endpoint%12,other%12)
        for inverse in (False,True):
            word=coupling.occupation.compiler.inverse(beam) if inverse else beam
            word=coupling.occupation.compiler.remap(word,positions)
            declaration=native.source_word(word,12)
            out.append({'cell':endpoint//12,'local_ports':positions,'inverse':inverse,
                        'source_word_digest':native.digest(declaration),
                        'resources':declaration['declared_resources']})
    return tuple(out)


def beam_execution_certificate():
    """Fresh full-column beam/inverse executions with all first failures."""
    coupling.synthesis.certify_primitive.cache_clear()
    beam=coupling.words()['beam'][1];rows=[]
    for inverse in (False,True):
        word=coupling.occupation.compiler.inverse(beam) if inverse else beam
        unique=coupling.synthesis.certify_word(word)
        run=coupling.recovery.execute_complete(word,2,np.eye(4,dtype=complex))
        ideal=coupling.ideal_matrix('beam')
        if inverse:ideal=ideal.conj().T
        bound=coupling.complete_bound(len(word))
        if np.linalg.norm(run.success-ideal,2)>float(F(bound['success_operator_error_upper'])):
            raise ArithmeticError('complete local beam execution exceeds its certified envelope')
        rows.append({'inverse':inverse,'pulses':len(word),'complete_input_columns':4,
                     'certified_primitive_types':unique,'error':bound,
                     'all_first_failure_outputs_retained':True})
    return tuple(rows)


def charge_block_check(level=0):
    """Actual boundary branches have a definite, recorded cell-number shift."""
    _,rows=branch_run(level,choi=True);_,edge,_=layout(level)
    cells=(edge[0]//12,edge[1]//12)
    for row in rows:
        for x,_,y,_ in row.density:
            left=tuple(((x>>(12*c))&4095).bit_count() for c in cells)
            right=tuple(((y>>(12*c))&4095).bit_count() for c in cells)
            if left!=right:raise ArithmeticError('a complete native branch created cell-charge coherence')
            if row.record.kind=='repair' and row.record.current is not None:
                if left!=(1-row.record.current,1+row.record.current):
                    raise ArithmeticError('the retained current does not determine the cell charge shift')
    return True
