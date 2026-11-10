"""Interval-certified full three-color five-block collision at N=1.

Exact small-cycle check of the declared Hamiltonian and projected packet
intervention. The physical low sector is exactly reducing; this check
asserts no three-dimensional continuum.
"""
from fractions import Fraction as F
from functools import lru_cache
import itertools
import math
import mpmath as mp


def exact_endpoint(interval,upper=False):
    sign,man,exponent,_=interval._mpi_[int(upper)]
    return F((-1)**sign*man)*(F(2)**exponent)


@lru_cache(maxsize=1)
def certificate():
    iv=mp.iv;iv.dps=65
    V=5;D=iv.mpf(2);coupling=iv.mpf(1)/8
    root=iv.sqrt(5)
    cs=[iv.mpf(1),(root-1)/4,-(root+1)/4,-(root+1)/4,(root-1)/4]
    sn1=iv.sqrt(10+2*root)/4;sn2=iv.sqrt(10-2*root)/4
    sn=[iv.mpf(0),sn1,sn2,-sn2,-sn1]
    phase=[iv.mpc(cs[j],sn[j]) for j in range(V)]
    omega=[iv.sqrt(D*D-4*D*coupling*c) for c in cs]
    u=[(iv.sqrt(D/w)+iv.sqrt(w/D))/2 for w in omega]
    v=[(iv.sqrt(D/w)-iv.sqrt(w/D))/2 for w in omega]
    zeta=[v[j]/u[j] for j in range(V)]
    Z=[[sum((zeta[j]*cs[(j*(x-y))%V] for j in range(V)),iv.mpf(0))/V for y in range(V)] for x in range(V)]
    determinant=iv.mpf(1)
    for z in zeta:determinant*=1-z*z
    vacuum=determinant**(iv.mpf(3)/4)
    @lru_cache(None)
    def pairs(sites):
        if not sites:return iv.mpf(1)
        if len(sites)%2:return iv.mpf(0)
        first=sites[0]
        return sum((Z[first][sites[j]]*pairs(sites[1:j]+sites[j+1:]) for j in range(1,len(sites))),iv.mpf(0))
    words=list(itertools.product(range(4),repeat=V));ids={w:i for i,w in enumerate(words)}
    def packet(sign):
        output=[]
        for word in words:
            value=iv.mpc(0)
            for x in range(V):
                if word[x]!=1:continue
                for y in range(V):
                    if word[y]!=2:continue
                    remaining=list(word);remaining[x]=remaining[y]=0
                    gaussian=iv.mpf(1)
                    for color in (1,2,3):
                        gaussian*=pairs(tuple(i for i,label in enumerate(remaining) if label==color))
                    value+=phase[(sign*(x-y))%V]*gaussian
            output.append(vacuum*value/(V*u[1]*u[1]))
        return output
    def inner(left,right):return sum((iv.mpc(x.real,-x.imag)*y for x,y in zip(left,right)),iv.mpc(0))
    incoming=packet(1);outgoing=packet(-1)
    pin=inner(incoming,incoming).real;pout=inner(outgoing,outgoing).real
    incoming=[x/iv.sqrt(pin) for x in incoming]
    outgoing=[x/iv.sqrt(pout) for x in outgoing]
    overlap=inner(incoming,outgoing)
    outgoing=[y-overlap*x for x,y in zip(incoming,outgoing)]
    outnorm=iv.sqrt(inner(outgoing,outgoing).real)
    outgoing=[x/outnorm for x in outgoing]
    rows=[{} for _ in words]
    for origin,word in enumerate(words):
        rows[origin][origin]=F(2*sum(label>0 for label in word))
        for x in range(V):
            y=(x+1)%V
            for color in (1,2,3):
                if word[x] in (0,color) and word[y] in (0,color):
                    target=list(word);target[x]=color-word[x];target[y]=color-word[y]
                    dest=ids[tuple(target)]
                    rows[dest][origin]=rows[dest].get(origin,F(0))-F(1,8)
    for row,entries in enumerate(rows):
        for col,value in entries.items():assert rows[col][row]==value
    ivrows=[[(j,iv.mpf(c.numerator)/c.denominator) for j,c in row.items() if c] for row in rows]
    def apply(vector):return [sum((c*vector[j] for j,c in row),iv.mpc(0)) for row in ivrows]
    matched_free_generator=inner(outgoing,apply(incoming))
    assert matched_free_generator.real.a>0 or matched_free_generator.real.b<0
    order=100
    clocks=[F(0),F(1,4),F(1,2),F(3,4),F(1)]
    term=incoming[:]
    evolved={clock:incoming[:] for clock in clocks}
    powers={clock:iv.mpf(1) for clock in clocks}
    for n in range(1,order+1):
        term=[-iv.j*x/n for x in apply(term)]
        for clock in clocks[1:]:
            powers[clock]*=iv.mpf(clock.numerator)/clock.denominator
            evolved[clock]=[x+powers[clock]*y for x,y in zip(evolved[clock],term)]
    amplitude=inner(outgoing,evolved[F(1)])
    norm_bound=iv.mpf(95)/8
    tail=iv.exp(norm_bound)*norm_bound**(order+1)/math.factorial(order+1)
    imag=amplitude.imag
    if imag.b<0:component=-imag
    else:
        assert imag.a>0
        component=imag
    assert component.a>tail.b
    amplitude_lower=exact_endpoint(component-tail)
    probability=(amplitude_lower/2)**2
    series=[]
    for clock in clocks:
        value=iv.mpf(clock.numerator)/clock.denominator
        remainder=iv.exp(norm_bound*value)*(norm_bound*value)**(order+1)/math.factorial(order+1)
        amp=inner(outgoing,evolved[clock])
        magnitude=iv.sqrt(amp.real**2+amp.imag**2)
        probability_lo=max(F(0),exact_endpoint(magnitude-remainder))**2
        probability_hi=min(F(1),max(F(0),exact_endpoint(magnitude+remainder,True))**2)
        if clock==0:probability_lo=probability_hi=F(0)
        moment_error=2*remainder+remainder**2
        marginals=[]
        for site in range(V):
            colors=[]
            for color in (1,2,3):
                expectation=sum((psi.real**2+psi.imag**2 for word,psi in zip(words,evolved[clock]) if word[site]==color),iv.mpf(0))
                lower=max(F(0),exact_endpoint(expectation-moment_error))
                upper=min(F(1),exact_endpoint(expectation+moment_error,True))
                colors.append({'component':color-1,'ideal_interval_exact':[str(lower),str(upper)],
                               'display_midpoint':float((lower+upper)/2)})
            empty=sum((psi.real**2+psi.imag**2 for word,psi in zip(words,evolved[clock]) if word[site]==0),iv.mpf(0))
            ground_lower=max(F(0),exact_endpoint(empty-moment_error))
            ground_upper=min(F(1),exact_endpoint(empty+moment_error,True))
            marginals.append({'block':site,'components':colors,
                              'ground_ideal_interval_exact':[str(ground_lower),str(ground_upper)],
                              'ground_display_midpoint':float((ground_lower+ground_upper)/2),
                              'higher_native_sector_ideal_probability':'0'})
        series.append({'clock_Delta_times_t':str(clock),
                       'ideal_hit_probability_interval_exact':[str(probability_lo),str(probability_hi)],
                       'ideal_complement_probability_interval_exact':[str(1-probability_hi),str(1-probability_lo)],
                       'ideal_hit_display_midpoint':float((probability_lo+probability_hi)/2),
                       'occupation_marginals':marginals,
                       'numerical_status':'outward interval arithmetic plus explicit operator Taylor remainder',
                       'uniform_marginals_follow_exact_translation_symmetry':True,
                       'classical_trajectories_claim':False})
    from .free_comparator import compare
    comparator=compare(incoming,outgoing,words,clocks,cs,exact_endpoint)
    for actual,free in zip(series,comparator,strict=True):actual.update(free)
    return {'schema':'oph.five_block_native_collision.v1','macro_graph':'five-site cycle, degree two',
            'native_eight_cell_blocks':5,'physical_source_cells':40,'replicas_per_site':1,
            'exact_reducing_source_dimension':4**5,'native_components_retained':3,
            'D_over_Delta':'2','t_over_Delta':'1/8','clock_Delta_times_T':'1',
            'input':'projected free alpha:k=2pi/5, beta:-k pair, normalized',
            'output':'projected momentum-swapped pair, orthogonalized to input',
            'full_native_complement_retained_by_invariance':True,
            'free_unprojected_channel_exactly_zero':True,
            'same_admitted_free_channel_exactly_zero':False,
            'same_admitted_free_generator_matrix_element_interval':str(matched_free_generator),
            'N1_Hamiltonian_identity':'H_N1=P_N H_free P_N on the admitted carrier',
            'free_comparator':'full quadratic Fock evolution, same admitted input and orthogonalized detector',
            'free_comparator_evolution_oscillator_cutoff':False,
            'input_Gaussian_projection_weight_interval':str(pin),
            'output_Gaussian_projection_weight_interval':str(pout),
            'initial_overlap_before_orthogonalization_interval':str(overlap),
            'Taylor_order':order,'Hamiltonian_norm_upper':'95/8',
            'amplitude_real_interval':str(amplitude.real),'amplitude_imaginary_interval':str(imag),
            'Taylor_remainder_norm_upper_exact':str(exact_endpoint(tail,True)),
            'exact_source_amplitude_lower':str(amplitude_lower),
            'conditional_hit_probability_lower_after_preparation_margin':str(probability),
            'conditional_hit_probability_lower_decimal':str(float(probability)),
            'complete_time_series':series,
            'actual_interacting_ground_input_claim':False,'three_dimensional_continuum_claim':False}


if __name__=='__main__':
    import json
    print(json.dumps(certificate(),indent=2))
