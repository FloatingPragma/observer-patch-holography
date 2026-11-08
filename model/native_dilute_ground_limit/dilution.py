"""Exact finite dilution bounds for the source Dicke ground."""
from fractions import Fraction as F
from math import isqrt
from model.native_program import interpreter as native


def integer(value,lower=0):
    if type(value)is not int or value<lower:raise ValueError('exact nonnegative integer required')
    return value


def population(ports,total,support):
    integer(ports,1);integer(total);integer(support)
    if total>ports or support>ports or support>64:
        raise ValueError('valid source sector and at most 64 local ports required')
    if ports-total<support:return F(1)
    empty=F(1)
    for j in range(support):empty*=F(ports-total-j,ports-j)
    occupied=1-empty
    if not 0<=occupied<=min(F(1),F(support*total,ports)):
        raise ArithmeticError('finite occupation probability violates its exact union bound')
    return occupied


def sqrt_upper(value,bits=128):
    if type(value)not in (int,F)or value<0 or type(bits)is not int or not 1<=bits<=192:
        raise ValueError('nonnegative exact value and bounded precision required')
    value=F(value);scale=1<<bits;n=isqrt(value.numerator*scale*scale//value.denominator)
    if F(n*n,scale*scale)<value:n+=1
    return F(n,scale)


def history_bound(ports,total,supports):
    if type(supports)is not tuple or not 1<=len(supports)<=64:
        raise ValueError('one to 64 declared complete intervention supports required')
    probabilities=tuple(population(ports,total,s)for s in supports)
    bound=min(F(1),2*sum((sqrt_upper(p)for p in probabilities),F()))
    return {'occupied_probabilities':probabilities,'complete_terminal_halftrace_upper':bound,
            'all_classical_history_TV_upper':bound,'arbitrary_clock_rescaling_allowed':True,
            'branchwise_local_number_preservation_required':True,
            'source_ground_input_required':True,'postselection_used':False}


def observable_bounds(ports,total,left,right):
    p=population(ports,total,left);r=population(ports,total,right)
    return {'norm_one_mean_error_upper':2*p,'norm_one_variance_upper':4*p,
            'norm_one_connected_correlation_upper':4*sqrt_upper(p*r),
            'connected_spectral_weight_upper':4*p,'uniform_over_controlled_times':True}


def complete_number_preserving_instrument(labels,branches,ports):
    """Bounded exact real-rational checker for (4), including zero slots.

    General complex instruments satisfy the proof; this finite checker
    certifies its real-rational subset with explicit nonnegative weights.
    """
    integer(ports)
    if ports>4 or type(labels)is not tuple or not 1<=len(labels)<=64:
        raise ValueError('bounded distinct complete labels and at most four ports required')
    if (any(type(x)is not str or not x for x in labels)or len(set(labels))!=len(labels)
            or type(branches)is not tuple or len(branches)!=len(labels)):
        raise ValueError('typed complete branch menu required')
    d=1<<ports;gram=[[F()for _ in range(d)]for _ in range(d)];vacuum=[]
    for label,branch in zip(labels,branches,strict=True):
        if type(branch)is not tuple or len(branch)!=3 or branch[0]!=label:
            raise ValueError('every declared fine label including zeros required')
        _,weight,matrix=branch
        if type(weight)not in (int,F)or not 0<=weight<=1 or type(matrix)is not tuple or len(matrix)!=d:
            raise ValueError('bounded exact weighted local Kraus matrix required')
        if any(type(row)is not tuple or len(row)!=d or any(type(x)not in (int,F)for x in row)for row in matrix):
            raise ValueError('complete exact square local matrix required')
        for i in range(d):
            for j in range(d):
                if i.bit_count()!=j.bit_count()and matrix[i][j]:
                    raise ValueError('branch changes local occupation or imports a charged ancilla')
                gram[i][j]+=weight*sum(matrix[k][i]*matrix[k][j]for k in range(d))
        vacuum.append(F(weight)*matrix[0][0]**2)
    if gram!=[[F(i==j)for j in range(d)]for i in range(d)]or sum(vacuum)!=1:
        raise ValueError('complete trace preservation and scalar vacuum law required')
    return tuple(vacuum)


def native_endpoint_certificate():
    states=native.local.basis();dimension=len(states)
    if dimension!=66 or any(sum(x)!=2 for x in states):raise ArithmeticError('literal H66 source changed')
    initial={(i,0,j,0):F(1,66)for i in range(66)for j in range(66)}
    command={'op':'repair','selector':{'kind':'fixed','seam':[0,1]},'law':'E','repeat':1}
    program=native.program(shell='H66',biased=False,commands=[command],capacity=1)
    history=native.execute(program,initial)
    if history.status!='COMPLETE'or history.attempted_steps!=1 or len(history.branches)!=25:
        raise ArithmeticError('public native complete endpoint execution changed')
    accepted=[row for row in history.branches if row.record.law=='E'and row.record.raw==0]
    expected={(i,0,j,0):F(1,132)for i,x in enumerate(states)for j,y in enumerate(states)
              if x[0]==x[1]and y[0]==y[1]}
    if len(accepted)!=2 or any(row.density!=expected for row in accepted):
        raise ArithmeticError('full coherent even-corner endpoint map differs')
    unequal=F(2*2*10,12*11)
    if sum(native.trace(row.density)for row in history.branches)!=1:
        raise ArithmeticError('complete native output lost mass')
    if sum(native.trace(row.density)for row in accepted)!=1-unequal:
        raise ArithmeticError('exact native ground parity probability differs')
    return {'ports':12,'total':2,'raw_slots':25,'unequal_probability':unequal,
            'complete_classical_vacuum_comparison_TV':unequal,
            'terminal_halftrace_linear_term':unequal/2,
            'terminal_halftrace_sqrt_radicand':unequal-F(3,4)*unequal**2,
            'public_native_resources':history.resources,'full_coherent_accepted_operators_checked':True}


def density_closure_certificate():
    """All twelve exact density columns on the literal H66 source."""
    from model.primitive_observer import transfer
    source=transfer.native_counts();transfer.exact_audit(source)
    states=source.states;counts=source.counts;V=12;q=2;D=66
    if len(states)!=D or tuple(states)!=tuple(native.local.basis()):
        raise ArithmeticError('native ground density chart changed')
    L=[[0 for _ in range(V)]for _ in range(V)]
    for a,b in native.local.SEAMS:
        L[a][a]+=1;L[b][b]+=1;L[a][b]-=1;L[b][a]-=1
    for i,state in enumerate(states):
        for port in range(V):
            left=60*state[port]-sum(counts[i][j]*states[j][port]for j in range(D))
            right=sum(state[k]*L[k][port]for k in range(V))
            if left!=right:raise ArithmeticError('exact ground density intertwining failed')
    residue=F(q*(V-q),V*(V-1))
    for a in range(V):
        for b in range(V):
            covariance=F(sum(x[a]*x[b]for x in states),D)-F(q*q,V*V)
            if covariance!=residue*(F(a==b)-F(1,V)):
                raise ArithmeticError('exact ground density spectral residue failed')
    return {'source_ports':V,'source_total':q,'complete_source_dimension':D,
            'exact_density_intertwining_entries':D*V,'exact_covariance_entries':V*V,
            'mean_zero_unit_l2_density_residue':residue,
            'transfer_frequency':'-2M log(1-lambda/(2M))',
            'same_ground_and_native_Laplacian_dictionary':True,
            'square_root_frequency_substituted':False}


def visibility_requirement(ports,total,norm,weight):
    """Necessary support for a requested connected spectral weight; no sufficiency."""
    integer(ports,1);integer(total)
    if (total>ports or type(norm)not in (int,F)or norm<=0
            or type(weight)not in (int,F)or not 0<=weight<=norm*norm):
        raise ValueError('valid source sector, positive norm and possible bounded weight required')
    threshold=F(weight)/(4*F(norm)**2)
    support=F(ports)*threshold/total if total else F()if weight==0 else None
    integer_support=None if support is None else -(-support.numerator//support.denominator)
    return {'minimum_scaled_occupation':threshold,'minimum_integer_support':integer_support,
            'positive_weight_impossible_in_empty_sector':total==0 and weight>0,
            'necessary_condition_only':True,'sufficient_particle_identification_claimed':False}
