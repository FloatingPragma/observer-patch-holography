"""Exact cyclic reduction of the declared coherent native event hardware."""
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from types import SimpleNamespace

from model.openwave_native_process import source as native
from model.openwave_coherent_native_architecture import source as architecture


def distance(g,origin=0):
    values=[-1]*g.ports;values[origin]=0;q=deque([origin])
    while q:
        a=q.popleft()
        for b in g.neighbours[a]:
            if values[b]<0:values[b]=values[a]+1;q.append(b)
    if min(values)<0:raise ValueError('connected graph required')
    return values


def integer_transfer(g,values):
    """D T = D I-L, D=2M; exact integer sparse action."""
    if len(values)!=g.ports:raise ValueError('complete port vector required')
    D=2*len(g.seams);out=[D*x for x in values]
    for a,b in g.seams:
        z=values[a]-values[b];out[a]-=z;out[b]+=z
    return out


def cyclic_step(g,a,b):
    """Numerators in ψ=A a+S A b, each over D^n sqrt(B)."""
    D=2*len(g.seams);tb=integer_transfer(g,b)
    return [-D*x for x in b],[D*x+2*y for x,y in zip(a,tb,strict=True)]


def port_numerators(g,a,b):
    """Full source probabilities, including both event-subspace components."""
    D=2*len(g.seams);tb=integer_transfer(g,b);tb2=integer_transfer(g,[x*x for x in b])
    return [D*x*x+y+2*x*z for x,y,z in zip(a,tb2,tb,strict=True)]


def interval(value,places=15):
    scale=10**places;lo=value.numerator*scale//value.denominator
    hi=-((-value.numerator*scale)//value.denominator)
    return [str(F(lo,scale)),str(F(hi,scale))]


def packet(g,radius,steps):
    """Uniform real graph-ball preparation; exact terminal port instrument."""
    if type(radius)is not int or radius<0 or any(type(n)is not int or n<0 for n in steps):
        raise ValueError('nonnegative integer width and coherent steps required')
    dist=distance(g);initial=[int(d<=radius) for d in dist];B=sum(initial)
    a=initial[:];b=[0]*g.ports;D=2*len(g.seams);power=1;out=[]
    for n in range(max(steps)+1):
        if n in steps:
            num=port_numerators(g,a,b);den=D*power*power*B
            if min(num)<0 or sum(num)!=den:raise ArithmeticError('complete walk norm/positivity failed')
            digest=sha256()
            for x in num:digest.update(str(x).encode()+b'\n')
            at_origin=F(num[0],den);outside=F(sum(x for x,d in zip(num,dist,strict=True) if d>radius),den)
            mean=F(sum(x*d for x,d in zip(num,dist,strict=True)),den)
            # Complete encoded-state return to Aψ, not compressed norm.
            ta=integer_transfer(g,b)
            overlap=sum((D*a[i]+ta[i])*initial[i] for i in range(g.ports))
            fidelity=F(overlap*overlap,D*D*power*power*B*B)
            out.append({'coherent_steps':n,'complete_port_probability_vector_sha256':digest.hexdigest(),
                        'denominator_sha256':sha256(str(den).encode()).hexdigest(),
                        'source_origin_probability_interval':interval(at_origin),
                        'outside_preparation_ball_probability_interval':interval(outside),
                        'mean_graph_distance_interval':interval(mean),
                        'full_encoded_return_fidelity_interval':interval(fidelity),
                        'total_probability_exact':'1','negative_outcomes':0})
        if n<max(steps):a,b=cyclic_step(g,a,b);power*=D
    return {'ports':g.ports,'seams':len(g.seams),'radius':radius,'preparation_ports':B,
            'held_out_preparation_width':radius==2,'records':out}


def small_graph():
    from model.observer import quantum_mixing
    seams=tuple(quantum_mixing.seams()); neighbours=[[]for _ in range(12)]
    for a,b in seams:neighbours[a].append(b);neighbours[b].append(a)
    return SimpleNamespace(ports=12,seams=seams,neighbours=neighbours)


def full_arc_spectrum(g):
    """Exact complement dimensions, using connectedness and positive T."""
    distance(g);N=g.ports;M=len(g.seams);D=2*M
    if D<=2*max(map(len,g.neighbours)):raise ValueError('strictly positive population transfer required')
    fixed_plus=M-N+2;fixed_minus=N*D-M-N;dispersive=2*(N-1)
    assert fixed_plus>=1 and fixed_minus>=0 and fixed_plus+fixed_minus+dispersive==N*D
    return {'ports':N,'seams':M,'full_source_event_dimension':N*D,
            'cyclic_dimension':2*N-1,'cyclic_constant_plus_one_dimension':1,
            'dispersive_dimension':dispersive,'full_W_plus_one_multiplicity':fixed_plus,
            'full_W_minus_one_multiplicity':fixed_minus,
            'unprepared_complement_plus_one_multiplicity':fixed_plus-1,
            'unprepared_complement_minus_one_multiplicity':fixed_minus,
            'source_swap_minus_one_multiplicity':M,
            'classification':'unprepared architectural event modes; no particle assignment'}


def explicit_arc_step(g,psi,phase=-1):
    """Exact full arc execution, including signed rational amplitudes."""
    E=2*len(g.seams)
    if len(psi)!=g.ports or any(len(row)!=E for row in psi):raise ValueError('all source/event arcs required')
    means=[sum(row,F())/E for row in psi]
    reflected=[[means[x]+phase*(v-means[x]) for v in row]for x,row in enumerate(psi)]
    out=[row[:]for row in reflected]
    for e,(a,b) in enumerate(g.seams):out[a][2*e+1],out[b][2*e+1]=reflected[b][2*e+1],reflected[a][2*e+1]
    return out


def arc_control():
    g=small_graph();D=2*len(g.seams)
    # All12 input columns independently execute on all720 arcs; their
    # complete Gram checks every144 source operator columns/reference.
    columns=[[[F(int(x==j))for _ in range(D)]for x in range(g.ports)]for j in range(g.ports)]
    a=[[int(x==j)for x in range(g.ports)]for j in range(g.ports)];b=[[0]*g.ports for _ in range(g.ports)]
    checked=0;terminal_checked=0
    for n in range(5):
        power=D**n
        transfers=[integer_transfer(g,z)for z in b]
        for j in range(g.ports):
            nums=port_numerators(g,a[j],b[j]);den=D*power*power
            actual=[sum(v*v for v in row)/D for row in columns[j]]
            assert actual==[F(z,den)for z in nums]
            for k in range(g.ports):
                explicit=[sum((x*y for x,y in zip(row1,row2,strict=True)),F())/D
                          for row1,row2 in zip(columns[j],columns[k],strict=True)]
                pair_transfer=integer_transfer(g,[x*y for x,y in zip(b[j],b[k],strict=True)])
                cross=[D*a[j][x]*a[k][x]+pair_transfer[x]
                       +a[j][x]*transfers[k][x]+transfers[j][x]*a[k][x]
                       for x in range(g.ports)]
                assert explicit==[F(z,D*power*power)for z in cross]
                terminal_checked+=g.ports
                assert sum(explicit)==int(j==k);checked+=1
        if n<4:
            columns=[explicit_arc_step(g,x)for x in columns]
            for j in range(g.ports):a[j],b[j]=cyclic_step(g,a[j],b[j])
    basis=[[F(int(x==0))for _ in range(D)]for x in range(12)]
    first=[explicit_arc_step(g,basis,z)for z in (-1,1)];assert first[0]==first[1]
    second=[explicit_arc_step(g,first[i],z)for i,z in enumerate((-1,1))]
    return_prob=[sum(v*v for v in wave[0])/D for wave in second]
    assert return_prob[1]==1 and return_prob[0]!=1
    # Exact full controlled-swap permutation has one exchanged pair per seam.
    swap=list(range(12*D))
    for e,(x,y)in enumerate(g.seams):swap[x*D+2*e+1],swap[y*D+2*e+1]=swap[y*D+2*e+1],swap[x*D+2*e+1]
    exchanged=sum(i!=x for i,x in enumerate(swap))//2
    spectrum=full_arc_spectrum(g)
    assert exchanged==30 and spectrum['full_W_plus_one_multiplicity']==20
    assert spectrum['full_W_minus_one_multiplicity']==678 and spectrum['dispersive_dimension']==22
    return {'complete_arc_carrier_dimension':720,'source_basis_columns':12,'matrix_unit_Gram_checks':checked,
            'matrix_unit_terminal_port_effect_checks':terminal_checked,
            'full_spectrum_census':spectrum,'actual_swap_exchanged_pairs':exchanged,
            'explicit_steps':4,'same_first_coherent_event':True,
            'reflection_two_step_origin_probability':str(return_prob[0]),
            'identity_continuation_two_step_origin_probability':str(return_prob[1]),
            'all_different_histories_kept_coherent':True}


def run():
    controls=arc_control();rows=[];dispersion=[];spectra=[]
    for level in (0,1,2):
        g=native.graph(level);scale=isqrt(2*len(g.seams));steps=(0,1,scale,2*scale)
        for radius in (0,1,2):rows.append(packet(g,radius,steps))
        D=2*len(g.seams)
        spectra.append(full_arc_spectrum(g))
        dispersion.append({'ports':g.ports,'laplacian_upper':12,
                           'rescaled_squared_phase_minus_laplacian_upper_exact':str(F(144,2*D-12)),
                           'relative_upper_factor_exact':str(F(1,1-F(6,D)))})
    return {'verdict':'PASS_COHERENT_NATIVE_PORT_PROPAGATION','construction':architecture.construction_descriptor(),
            'named_program':'coherently recoded native event followed by uniform-event reflection',
            'exact_recurrence':'a_next=-b; b_next=a+2*T*b; psi=A*a+S*A*b',
            'exact_port_probability':'abs(a_x)^2+(T abs(b)^2)_x+2 Re(conj(a_x)*(Tb)_x)',
            'dispersion':'cos(omega)=1-lambda_L/(2M)',
            'finite_dispersion_bounds':dispersion,
            'complete_arc_spectrum_census':spectra,
            'zero_mode':'constant port amplitude has omega=0',
            'coherent_clock_convention':'steps/sqrt(2M); dimensionless declared comparison, not laboratory time',
            'full_arc_control':controls,'packets':rows,
            'first_arrival_measurement':False,'terminal_instrument':'all source port projectors, with complete sum I on the one-load carrier',
            'preparation_and_event_reflection':'coherent uniform preparation and its reflection are declared hardware services',
            'fine_graph_full_reflection_compiled':False,
            'native_quantum_channel_intertwiner':'original C/E recovered only at inverse-recode then original-record commit intervention',
            'original_unread_quantum_channel_equals_swap_trace':False,
            'Klein_Gordon_physical_mass_term_derived':False,'physical_clock_or_speed':False,
            'new_GeV_rest_masses':False,'physical_OpenWave_row_award':False}
