"""Exact unread response derived from complete native branch maps.

The output is a marginal of the recorded process, not its whole observer
history. Source marginals suffice for the fixed independent seam schedule;
the complete native records remain available to the observer.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json

from model.openwave_native_process import source as native


def population_column(g, vertex):
    if type(vertex) is not int or not 0 <= vertex < g.ports:
        raise ValueError('one source port required')
    denominator = 2*len(g.seams)
    out = {j: F(1, denominator) for j in g.neighbours[vertex]}
    out[vertex] = 1-F(len(g.neighbours[vertex]), denominator)
    return out


def coherence_multiplier(g, left, right):
    if left == right or not (0 <= left < g.ports and 0 <= right < g.ports):
        raise ValueError('two distinct source ports required')
    adjacent = int(right in g.neighbours[left])
    return 1-F(3*(len(g.neighbours[left])+len(g.neighbours[right]))-2*adjacent,
               4*len(g.seams))


def apply_population(g, values):
    if len(values) != g.ports:
        raise ValueError('complete source population required')
    out = list(values); denominator = 2*len(g.seams)
    for a,b in g.seams:
        current = (values[a]-values[b])/denominator
        out[a] -= current; out[b] += current
    return out


def coarse(values):
    if len(values)%12:
        raise ValueError('complete twelve-port cells required')
    return [sum(values[i:i+12], F()) for i in range(0,len(values),12)]


def lift(values):
    return [value/12 for value in values for _ in range(12)]


def endpoint_matrix_unit_certificate():
    """Independent endpoint formulas checked on all16 operator columns.

    A third spectator port supplies the one-load off-endpoint column. The
    fourth endpoint occupation tests the full binary local complement too.
    """
    matrices = {}
    for law in ('C','E','mixture'):
        cols=[]
        for x,y in product(range(4),repeat=2):
            expected=defaultdict(F)
            for kind,cw in (('C',F(1,2) if law=='mixture' else F(law=='C')),
                            ('E',F(1,2) if law=='mixture' else F(law=='E'))):
                for c in (0,1):
                    # Independent binary averaging formula, including00/11.
                    dx=(x&1)-((x>>1)&1);dy=(y&1)-((y>>1)&1)
                    jx=(dx-(2*c-1)*abs(dx))//2
                    jy=(dy-(2*c-1)*abs(dy))//2
                    if (jx==jy if kind=='C' else dx==dy):
                        xx=x+jx; yy=y+jy  # removing bit0 and adding bit1 addsj.
                        expected[xx,yy]+=cw/2
            expected={k:v for k,v in expected.items() if v}
            actual=native.unread_column(x,y,(0,1),law)
            assert actual==expected,(law,x,y,actual,expected)
            assert len(native.branch_column(x,y,(0,1),law))==25
            cols.append([x,y,[[a,b,str(v)] for (a,b),v in sorted(actual.items())]])
        matrices[law]=cols
    # All endpoint/spectator one-particle columns certify global tensor lift.
    for i,j in product(range(3),repeat=2):
        actual=native.unread_column(1<<i,1<<j,(0,1))
        if i==j:
            expected=({(1<<i,1<<i):F(1)} if i==2 else
                      {(1,1):F(1,2),(2,2):F(1,2)})
        elif {i,j}=={0,1}:
            expected={}
        elif 2 in (i,j):
            expected={(1<<i,1<<j):F(1,4)}
        assert actual==expected
    return {'complete_endpoint_operator_columns':48,'one_load_endpoint_spectator_columns':9,
            'raw_slots_per_seam':25,'zero_hold_and_padding_retained':True,
            'all_branch_maps_by_law':matrices}


def graph_certificate(level):
    g=native.graph(level);M=len(g.seams)
    # Accumulate every actual seam-population column from its two active
    # local basis columns. Spectator columns contribute identity exactly.
    summed=[defaultdict(F,{i:F(1)}) for i in range(g.ports)]
    for a,b in g.seams:
        for i in (a,b):
            actual=native.unread_column(1<<i,1<<i,(a,b))
            summed[i][i]-=F(1,M)
            for (x,y),value in actual.items():
                assert x==y and x.bit_count()==1
                summed[i][x.bit_length()-1]+=value/M
    columns=[]
    for i in range(g.ports):
        column={j:v for j,v in summed[i].items() if v}
        assert column==population_column(g,i)
        assert sum(column.values())==1 and min(column.values())>=0
        columns.append([[j,str(v)] for j,v in sorted(column.items())])
    # Exact endpoint-incidence partition for every off-diagonal basis column.
    counts=Counter()
    for i in range(g.ports):
        for j in range(i+1,g.ports):
            a=int(j in g.neighbours[i]);di=len(g.neighbours[i]);dj=len(g.neighbours[j])
            unaffected=M-di-dj+a
            one_endpoint=di+dj-2*a
            numerator=4*unaffected+one_endpoint
            assert numerator==4*M-3*(di+dj)+2*a
            counts[(di,dj,a,numerator)]+=1
            counts[(dj,di,a,numerator)]+=1
    assert sum(counts.values())==g.ports*(g.ports-1)
    # Every cell-projection column is compared with its uniform-lift closure.
    coarse_columns=[]
    for i in range(g.ports):
        row=[F() for _ in range(g.cells)]
        for j,value in population_column(g,i).items():row[j//12]+=value
        coarse_columns.append(row)
    max_error=F(); witness=None
    for cell in range(g.cells):
        average=[sum((coarse_columns[12*cell+p][c] for p in range(12)),F())/12
                 for c in range(g.cells)]
        for p in range(12):
            error=sum(abs(a-b) for a,b in zip(coarse_columns[12*cell+p],average,strict=True))
            if error>max_error:max_error,witness=error,[cell,p]
    assert max_error==F(11,12*M)
    # A two-step uniform-input counterexample needs neither eigenfit nor
    # an assumed coarse evolution. The exact memory term is CT(I-SC)TS.
    initial=lift([F(i==0) for i in range(g.cells)])
    first=apply_population(g,initial)
    actual2=coarse(apply_population(g,first))
    closed2=coarse(apply_population(g,lift(coarse(first))))
    memory=[a-b for a,b in zip(actual2,closed2,strict=True)]
    assert sum(memory)==0 and any(memory)
    return {'level':level,'ports':g.ports,'cells':g.cells,'seams':M,
            'complete_population_columns':g.ports,'complete_offdiagonal_columns':g.ports*(g.ports-1),
            'population_operator':'T=I-L/(2M)',
            'population_columns_sha256':sha256(json.dumps(columns,separators=(',',':')).encode()).hexdigest(),
            'coherence_multiplier_classes':[{'degrees':[a,b],'adjacent':bool(c),'multiplier_exact':str(F(v,4*M)),
                                            'ordered_column_count':count}
                                           for (a,b,c,v),count in sorted(counts.items())],
            'cell_projection_one_step_l1_max_exact':str(max_error),
            'intensity_scaled_projection_residual_exact':str(2*M*max_error),
            'cell_projection_counterexample_cell_and_port':witness,
            'two_step_uniform_cell_memory_l1_exact':str(sum(map(abs,memory))),
            'two_step_cell_memory_vector_exact':list(map(str,memory)),
            'source_stationary_state':'I/ports; no target density or energy is used',
            'positive_population_spectrum_lower_bound_exact':str(1-F(g.max_degree,M))}


def run():
    return {'schema':'oph.native_response.v1','verdict':'PASS: complete native finite response and exact coarse-memory obstruction',
        'member':{'parent':'model/openwave_native_process','total_load':1,
                  'seam_schedule':'independent uniform existing seams','C_weight_exact':'1/2','E_weight_exact':'1/2',
                  'repair_probability_exact':'1','coin_probability_exact':'1/2',
                  'stationary_input':'full one-load maximally mixed density',
                  'response':'the channel is linear, so its exact derivative at the stationary input equals itself'},
        'endpoint_certificate':endpoint_matrix_unit_certificate(),
        'graphs':[graph_certificate(level) for level in (0,1,2)],
        'full_channel':{'diagonal':'I-L/(2M)',
            'offdiagonal':'Phi(|i><j|)=[1-(3(deg(i)+deg(j))-2A_ij)/(4M)]|i><j|, i!=j',
            'reference':'All identities tensor with the identity on any idle reference.',
            'record_boundary':'The formula is the exact unread-source marginal. Complete raw records are retained by the underlying instrument.',
            'memory_identity':'C T^(n+1) S = A C T^n S + sum_{j=0}^{n-1} C T Q (Q T Q)^(n-1-j) Q T S C T^j S; A=CTS,Q=I-SC'},
        'cosine_symbol':{'M10_1_object':'sum_u[1-cos(a u.q)]/(2 a^2) over twelve icosahedral directions',
            'native_transfer_proved':False,
            'reason':'The executed source is a finite port-pair graph channel. No translation representation, plane-wave intertwiner or second-order wave-time law is supplied. Even the natural cell-mass projection fails an exact basis-column intertwiner.',
            'reuse':'The exact orbit expansion remains valid for its declared cosine operator; it is not used as the native response.'},
        'scope':{'exact_native_operator_response':True,'full_record_history_replaced_by_unread_source':False,
                 'continuum_dispersion_or_physical_particle_identification':False,'new_Nature_validation':False}}
