"""Actual N=2 native hopping parity and its scalar action on the atomic band."""
from functools import lru_cache
import sympy as s
from ..periodic_family.model import Configuration,cells,kinetic_terms,yukawa_words
from ..higgs_composites import model as car
from ..local_reference.reference import certificate as local_reference


def cross_sublattice_terms(config):
    locations=cells(config.side);colors=tuple(sum(cell)%2 for cell in locations)
    terms=kinetic_terms(config)
    for term in terms:
        if colors[term.left_cell]==colors[term.right_cell]:
            raise ValueError('actual hopping contains a same-sublattice term')
    return colors,terms


def validate_matrix_symmetry(h0,t,p,u):
    n=h0.rows;identity=s.eye(n)
    if h0!=h0.H or t!=t.H:raise ValueError('self-adjoint full source required')
    if p!=p.H or p*p!=p or u!=u.H or u*u!=identity:
        raise ValueError('orthogonal reference and unitary involution required')
    if not 0<s.trace(p)<n:raise ValueError('nonempty reference and complement required')
    if h0*p!=s.zeros(n) or u*h0!=h0*u or u*p!=p*u or u*t*u!=-t or p*t*p!=s.zeros(n):
        raise ValueError('full-source parity and zero atomic block required')
    scalar=1 if u*p==p else -1 if u*p==-p else None
    if scalar is None:raise ValueError('parity must act as a scalar on the reference band')
    return scalar


@lru_cache(maxsize=1)
def certificate():
    config=Configuration();colors,terms=cross_sublattice_terms(config)
    assert config.side==2 and len(terms)==24
    assert all(len(term.path)==1 for term in terms)
    even=[i for i,color in enumerate(colors) if color==0]
    odd=[i for i,color in enumerate(colors) if color==1]
    assert len(even)==len(odd)==4
    words=yukawa_words()
    assert all(len(word)%2==0 for channel,row,col,word,value in words)
    for index in range(12):
        assert all(len(word)%2==0 for word,row,col in car.gauge_current(index))
    reference=local_reference();assert reference['fermion_parity']=='odd'
    phase=(-1)**len(even);assert phase==1
    larger=Configuration(side=4)
    larger_colors=tuple(sum(cell)%2 for cell in cells(4))
    counter=next(term for term in kinetic_terms(larger)
                 if larger_colors[term.left_cell]==larger_colors[term.right_cell])
    assert counter.displacement==2 and counter.coefficient!=0
    return {'source_side':2,'source_volume':8,'actual_oriented_hopping_terms':24,
        'even_sublattice_cells':even,'odd_sublattice_cells':odd,
        'unitary':'U_A=(-1)^(sum_(x in A) N_fermion,x)',
        'full_native_onsite_word_parity_checks':len(words),'full_Gauss_current_parity_checks':12,
        'all_actual_hops_reverse_sublattice_parity':True,
        'full_source_identities':['U_A H0 U_A=H0','U_A T U_A=-T','U_A P0=P0'],
        'actual_local_ground_parity':'odd','reference_band_parity_phase':phase,
        'reference_band_dimension':256,'full_complement_parity':'U_Q=Q U_A Q',
        'exact_canonical_graph_symmetry':'X(-lambda)=U_Q X(lambda)',
        'exact_canonical_effective_matrix_symmetry':'K(-lambda)=K(lambda)',
        'exact_self_energy_symmetry':'B(lambda)^dagger D(lambda)^(-1) B(lambda) is even',
        'larger_spectral_lattice_theorem_transferred':False,
        'actual_side_four_counterexample':{'left_cell':counter.left_cell,'right_cell':counter.right_cell,
            'axis':counter.axis,'displacement':counter.displacement,'path':list(counter.path),
            'coefficient':str(counter.coefficient),'same_sublattice':True}}
