"""Finite lazy native M6 words for local oscillators, routes and edge probes."""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
from providers import native
import geometry as g


def register_bits(cap):
    if type(cap) is not int or cap < 1:
        raise ValueError('positive local number cap required')
    return cap.bit_length()


def remap(word, wires):
    for kind, targets, angle in word:
        yield kind, tuple(wires[q] for q in targets), angle


def pauli_word(label, angle, wires):
    return remap(native.pauli(label,native.angle(0,0,angle)),wires)


def pair_word(left, right, angle, wires, controls=(), quadrature='X'):
    """Exact commuting Pauli expansion of one controlled two-level rotation.

    X is |r><l|+|l><r|; Y is i|r><l|-i|l><r|. Other columns,
    including local padding and false predicates, have identity action.
    """
    k=len(wires); controls=tuple(controls)
    if not 0<=left<right<1<<k or quadrature not in ('X','Y'):
        raise ValueError('ordered code pair and X/Y quadrature required')
    if set(wires)&{q for q,v in controls} or len(set(wires))!=k:
        raise ValueError('disjoint target and predicate registers required')
    differing=left^right
    for mask in range(1<<k):
        nY=(mask&differing).bit_count()
        if nY%2 != (quadrature=='Y'):
            continue
        sign=(-1)**((nY//2)+(left&mask).bit_count())
        label=''.join(('Y' if mask>>j&1 else 'X') if differing>>j&1
                      else ('Z' if mask>>j&1 else 'I') for j in range(k))
        for control_mask in range(1<<len(controls)):
            control_sign=(-1)**sum(v for j,(q,v) in enumerate(controls) if control_mask>>j&1)
            labels=label+''.join('Z' if control_mask>>j&1 else 'I' for j in range(len(controls)))
            coefficient=F(2*sign*control_sign,1<<(k+len(controls)))
            yield from pauli_word(labels,F(angle)*coefficient,tuple(wires)+tuple(q for q,v in controls))


def signed_prefix_sum(cap, mask):
    """Sum sign(n), sum n sign(n) for 0<=n<=cap, in O(log cap)."""
    states={(True):(1,0)}
    for bit in reversed(range(max(1,cap.bit_length()))):
        nxt={}
        for tight,(count,weighted) in states.items():
            maximum=(cap>>bit)&1 if tight else 1
            for value in range(maximum+1):
                key=tight and value==maximum
                sign=-1 if value and mask>>bit&1 else 1
                c,w=nxt.get(key,(0,0))
                nxt[key]=(c+sign*count,w+sign*(weighted+count*value*(1<<bit)))
        states=nxt
    return tuple(sum(v[j] for v in states.values()) for j in (0,1))


def diagonal_linear_word(cap, slope, intercept, wires, controls=()):
    """exp(i[slope*n+intercept]) on n<=cap; identity on padding."""
    k=len(wires); controls=tuple(controls)
    for mask in range(1<<k):
        count,weighted=signed_prefix_sum(cap,mask)
        coefficient=(F(slope)*weighted+F(intercept)*count)/(1<<k)
        for cmask in range(1<<len(controls)):
            sign=(-1)**sum(v for j,(q,v) in enumerate(controls) if cmask>>j&1)
            label=''.join('Z' if mask>>j&1 else 'I' for j in range(k))
            label+=''.join('Z' if cmask>>j&1 else 'I' for j in range(len(controls)))
            yield from pauli_word(label,coefficient*sign/(1<<len(controls)),tuple(wires)+tuple(q for q,v in controls))


def diagonal_values_word(values, wires, controls=()):
    k=len(wires)
    if len(values)!=1<<k:
        raise ValueError('all padded diagonal columns required')
    for mask in range(1<<k):
        coefficient=sum((F(x)*(-1)**((j&mask).bit_count()) for j,x in enumerate(values)),F())/(1<<k)
        for cmask in range(1<<len(controls)):
            sign=(-1)**sum(v for j,(q,v) in enumerate(controls) if cmask>>j&1)
            label=''.join('Z' if mask>>j&1 else 'I' for j in range(k))
            label+=''.join('Z' if cmask>>j&1 else 'I' for j in range(len(controls)))
            yield from pauli_word(label,coefficient*sign/(1<<len(controls)),tuple(wires)+tuple(q for q,v in controls))


def translation_word(cap, gamma, slices, bits, wires, controls=(), quadrature='X'):
    if type(slices) is not int or slices<1 or len(wires)!=register_bits(cap):
        raise ValueError('positive subdivision and complete local register required')
    for _ in range(slices):
        for n in range(1,cap+1):
            coefficient=F(gamma)*g.midpoint(g.interval_sqrt(n,bits))/slices
            yield from pair_word(n-1,n,coefficient,wires,controls,quadrature)


def harmonic_free_word(cap, time, slices, bits, wires):
    # Compression of the full p^2 includes its full diagonal at n=cap.
    for _ in range(slices):
        yield from diagonal_linear_word(cap,-F(time)*g.A/(2*slices),-F(time)*g.A/(4*slices),wires)
        for n in range(cap-1):
            angle=F(time)*g.A*g.midpoint(g.interval_sqrt((n+1)*(n+2),bits))/(4*slices)
            yield from pair_word(n,n+2,angle,wires)


def layout(coexact_cap, harmonic_cap):
    caps=[coexact_cap]*190+[harmonic_cap]*3
    modes=[]; offset=5
    for cap in caps:
        modes.append(tuple(range(offset,offset+register_bits(cap))))
        offset+=register_bits(cap)
    return {'spin':(0,1),'cell':(2,3,4),'modes':modes,'clock':offset,
            'probe':offset+1,'memory_bits':offset+2,'caps':caps}


def preparation_contract(coexact_cap, harmonic_cap, timeout=F(1,10**8)):
    info=layout(coexact_cap,harmonic_cap)
    wires=sum(map(len,info['modes']))
    return {'encoding':'190 coexact local number registers and 3 harmonic Hermite registers',
            'coexact_cap':coexact_cap,'harmonic_cap':harmonic_cap,'field_qubits':wires,
            'prepared_codeword':'all field bits zero',
            'native_reset':native.reset_repetition_plan(wires,F(timeout)),
            'coexact_state':'vacuum of frozen diagonal Lambda comparison source',
            'harmonic_state':'product Gaussian Var(q)=Var(p)=1/2; exact auxiliary Hermite vacuum',
            'charge_sector':'independent label q; physical E=B^T L^+ q+U p_co+H p_h',
            'padding':'local n>cap is retained with zero generator; initial padding weight zero',
            'vacuum_transfer':'transfer.geometry_transfer prices comparison vacuum versus actual K vacuum',
            'conditioning_applied':False}


def spin_diagonalizer(axis, wires=(0,1)):
    basis=[]
    for j,letter in enumerate(('XX','XY','XZ')[axis]):
        if letter=='X':basis.append(('Y',(j,),native.angle(F(-1,4))))
        if letter=='Y':basis.append(('X',(j,),native.angle(F(1,4))))
    yield from remap(native.negate_template(basis)+native.cnot(0,1),wires)


def route_axis_word(axis, info, slices, bits):
    diagonalizer=list(spin_diagonalizer(axis,info['spin']))
    yield from diagonalizer
    for cell in range(8):
        for sign in (-1,1):
            predicate=tuple((q,(cell>>j)&1) for j,q in enumerate(info['cell']))+((info['spin'][1],int(sign<0)),)
            for mode,(cap,wires) in enumerate(zip(info['caps'],info['modes'])):
                gamma=g.midpoint(g.coefficient(cell,axis,sign,mode,bits))
                yield from translation_word(cap,gamma,slices,bits,wires,predicate)
    target=[g.current(cell,axis,1)['destination'] for cell in range(8)]
    assert target==[g.current(cell,axis,-1)['destination'] for cell in range(8)]
    # The size-two destination map is a single cell-coordinate bit flip.
    difference={i^j for i,j in enumerate(target)}
    assert len(difference)==1
    mask=difference.pop(); assert mask and mask&(mask-1)==0
    yield from remap(native.primitive_templates()['flip'],(info['cell'][mask.bit_length()-1],))
    yield from native.negate_template(diagonalizer)


def free_field_word(info, time, slices, bits, include_coulomb=True):
    for mode,(cap,wires) in enumerate(zip(info['caps'],info['modes'])):
        if mode<190:
            omega=g.midpoint(g.frequency(mode,bits))
            yield from diagonal_linear_word(cap,-F(time)*omega,-F(time)*omega/2,wires)
        else:
            yield from harmonic_free_word(cap,time,slices,bits,wires)
    if include_coulomb:
        values=[-F(time)*row['coulomb'] for row in g.data()['sectors']]
        yield from diagonal_values_word(values,info['cell'])


def macro_word(coexact_cap,harmonic_cap,slices,bits,h=F(1,8),mass=F(1)):
    info=layout(coexact_cap,harmonic_cap)
    yield 'Z',(info['spin'][0],),native.angle(0,0,-F(h)*F(mass))
    for axis in range(3):
        yield from route_axis_word(axis,info,slices,bits)
    yield from free_field_word(info,F(1),slices,bits)
    yield 'Z',(info['clock'],),native.angle(0,0,-F(h)/16)


def plane_preparation_word(info):
    """After all-zero reset, prepare all8 root cells and the common clock+."""
    for wire in (*info['cell'],info['clock']):
        yield from remap(native.primitive_templates()['hadamard'],(wire,))


def edge_characteristic_word(edge,time,info,slices,bits):
    """Controlled exp(i time E_edge); the probe reads its full Fourier law."""
    control=((info['probe'],1),)
    values=[F(time)*row['particular'][edge] for row in g.data()['sectors']]
    yield from diagonal_values_word(values,info['cell'],control)
    for mode,(cap,wires) in enumerate(zip(info['caps'],info['modes'])):
        gamma=F(time)*g.midpoint(g.edge_coefficient(edge,mode,bits))
        yield from translation_word(cap,gamma,slices,bits,wires,control,'Y')


def edge_ramsey_word(edge,time,info,slices,bits,setting='X'):
    if setting not in ('X','Y'):
        raise ValueError('X/Y Fourier probe setting required')
    yield from remap(native.primitive_templates()['hadamard'],(info['probe'],))
    yield from edge_characteristic_word(edge,time,info,slices,bits)
    if setting=='Y':
        yield 'Z',(info['probe'],),native.angle(F(1,4))
    yield from remap(native.primitive_templates()['hadamard'],(info['probe'],))


def pointer_word(probe,work0,work1):
    """From two blank work bits, prepare READY code2 and copy to C/E pointer."""
    if len({probe,work0,work1})!=3:
        raise ValueError('fresh distinct native pointer wires required')
    yield from remap(native.primitive_templates()['flip'],(work1,))
    yield from remap(native.cnot(2,0)+native.cnot(0,1),(work0,work1,probe))


def pulses(word,bits):
    for kind,targets,angle in word:
        yield native.certified_pulse(kind,targets,angle,bits)


def translation_error(cap,gamma,slices,bits):
    # Adjacent ladder pairs are the only noncommuting pairs. Their commutator
    # norm is |gamma|^2 sqrt(n(n+1)); use an exact rational majorant.
    trotter=F(gamma)**2*cap*(cap+1)/(4*slices)
    algebraic=abs(F(gamma))*cap/F(1<<bits)
    return trotter+algebraic


def harmonic_free_error(cap,time,slices,bits):
    normsum=abs(F(time))*g.A*((2*cap+1)+cap*(cap+1))/4
    return normsum**2/(2*slices)+abs(F(time))*g.A*cap/F(1<<(bits+2))


def resource_bound(coexact_cap,harmonic_cap,slices,macros=1):
    """Finite integer bound without materializing any whole-field matrix."""
    info=layout(coexact_cap,harmonic_cap)
    pairs=sum(cap*(1<<register_bits(cap))//2 for cap in info['caps'])
    maxk=max(map(len,info['modes']))+4
    per_pauli=8*maxk+1
    kick_terms=3*16*16*slices*pairs
    harmonic_terms=3*slices*((1<<register_bits(harmonic_cap))+(harmonic_cap-1)*(1<<register_bits(harmonic_cap))//2)
    phases=sum(1<<register_bits(cap) for cap in info['caps'][:190])+8
    total=macros*(per_pauli*(kick_terms+harmonic_terms+phases)+200)
    return {'native_pulses_upper':total,'load_export_services_upper':2*total,
            'field_qubits':sum(map(len,info['modes'])),'memory_bits':info['memory_bits'],
            'enumerated':False,'growth':'O(macros*slices*193*R^2*log R), fixed 16-route predicates',
            'reset_and_pointer_services':'separate complete instruments; excluded from unitary pulse count',
            'physical_dimension_expression':f'({coexact_cap+1})^190 * ({harmonic_cap+1})^3',
            'padding_count_expression':f'2^{sum(map(len,info["modes"]))} - ({coexact_cap+1})^190*({harmonic_cap+1})^3'}


def universal_pulse_error(bits, rational_angle_upper):
    coefficients,norms=native.source_constants(bits)
    pi=native.constants(bits)['pi']; unit=F(1,1<<bits)
    theta=max(abs(pi[0]),abs(pi[1]))+F(rational_angle_upper)+8*unit
    theta_error=(pi[1]-pi[0])+8*unit
    return max(sum((norm*(theta*(hi-lo)+theta_error*max(abs(lo),abs(hi))+8*unit)
                    for (lo,hi),norm in zip(row,norms)),F()) for row in coefficients.values())


def static_holding_contract(time=F(8), tolerance=F(1,10**8)):
    """Common 190+3 field preparation/holding, with a separate charge word.

    A charge-sector Ramsey experiment can cancel this common field unitary
    in its ideal reduced law. Reset timeout and native export errors are
    retained; full-state and edge comparisons also use the stated tails.
    """
    import transfer
    time=F(time); tolerance=F(tolerance)
    if time<0 or not 0<tolerance<1:
        raise ValueError('nonnegative exact holding time and tolerance required')
    co=1; harmonic=32
    while F(transfer.static_harmonic_box(time,harmonic)['state_error_upper_exact'])>tolerance:
        harmonic*=2
    normsum=time*g.A*((2*harmonic+1)+harmonic*(harmonic+1))/4
    slices=max(1,(3*normsum**2/(2*tolerance)).__ceil__())
    k=register_bits(harmonic); per=8*k+1
    pulses=190*2*9+3*slices*((1<<k)+(harmonic-1)*(1<<k)//2)*per
    precision=max(80,(pulses*tolerance.denominator//tolerance.numerator+1).bit_length()+80)
    angle=max(F(2),time*(harmonic+1)**2)
    while pulses*universal_pulse_error(precision,angle)>tolerance:
        precision*=2
    native_error=pulses*universal_pulse_error(precision,angle)
    rounding=190*time*F(3,2)/(1<<precision)+3*time*g.A*harmonic/F(1<<(precision+2))
    return {'schema':'oph.native_field.static_holding.v1','time_exact':str(time),
            'coexact_cap':co,'harmonic_cap':harmonic,'subdivisions':slices,'precision_bits':precision,
            'preparation':preparation_contract(co,harmonic,tolerance),
            'native_pulses_upper':pulses,'native_load_export_services_upper':2*pulses,
            'native_export_operator_error_exact':str(native_error),
            'Trotter_operator_error_exact':str(3*normsum**2/(2*slices)),
            'algebraic_rounding_operator_error_exact':str(rounding),
            'harmonic_state_transfer':transfer.static_harmonic_box(time,harmonic),
            'coexact_vacuum_truncation_error_exact':'0',
            'pulse_iterator':'pulses(free_field_word(layout(coexact_cap,harmonic_cap),time,subdivisions,precision_bits,include_coulomb=False),precision_bits)',
            'common_unitary':'all190 half-frequency phases and all3 harmonic p^2 words retained',
            'charge_energy':'supplied separately by consuming sector word; no charge dependence in this common factor',
            'reduced_sector_clock_law':'common ideal field unitary cancels for every initial common field state; cutoff/splitting error does not affect this reduced identity',
            'full_state_scope':'true actual-K vacuum comparison additionally requires spectral geometry transfer; no cancellation claimed for E_edge records'}
