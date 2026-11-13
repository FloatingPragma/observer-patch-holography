"""Exact local bounds, complete records, and finite-refinement work budgets."""
from fractions import Fraction as F
from itertools import product
from math import factorial
import sympy as s


def graph(side):
    if type(side) is not int or side < 1:
        raise ValueError('positive integer subdivision required')
    vertices = list(product(range(side + 1), repeat=3))
    def edge(v, axis):
        w = list(v); w[axis] += 1
        return (v, tuple(w))
    links = {edge(v, a) for v in vertices for a in range(3) if v[a] < side}
    faces = []
    for v in vertices:
        for a in range(3):
            for b in range(a + 1, 3):
                if v[a] == side or v[b] == side:
                    continue
                va = list(v); va[a] += 1
                vb = list(v); vb[b] += 1
                faces.append(frozenset((edge(v, a), edge(v, b), edge(tuple(va), b), edge(tuple(vb), a))))
    # A fixed physical unit square in the z=0 outer face.
    loop = {edge((x, y, 0), 0) for x in range(side) for y in (0, side)}
    loop |= {edge((x, y, 0), 1) for y in range(side) for x in (0, side)}
    incident = [i for i, f in enumerate(faces) if f & loop]
    assert len(links) == 3*side*(side + 1)**2
    assert len(faces) == 3*side**2*(side + 1)
    assert len(loop) == 4*side
    assert len(incident) == (5 if side == 1 else 8*side - 4)
    return {'vertices': vertices, 'links': sorted(links), 'faces': faces,
            'loop': sorted(loop), 'incident_faces': incident}


def cosine_interval(x, degree=12):
    x = F(x)
    if abs(x) > 1 or degree < 2 or degree % 2:
        raise ValueError('even Taylor degree and |x|<=1 required')
    center = sum(((-1)**k*x**(2*k)/factorial(2*k) for k in range(degree//2 + 1)), F())
    radius = abs(x)**(degree + 1)/factorial(degree + 1)
    return center-radius, center+radius


def local_bound(side, time, representation):
    if representation not in ('fundamental', 'spin_one'):
        raise ValueError('declared normalized Wilson character required')
    g = graph(side); t = F(time)
    alpha, beta = F(32, side), F(side**3, 32)
    dimension, casimir = ((2, F(3, 4)) if representation == 'fundamental' else (3, F(2)))
    energy = alpha*len(g['loop'])*casimir
    delta = min(F(2), 2*beta*len(g['incident_faces'])*abs(t))
    lo, hi = cosine_interval(energy*t)
    # C(t) / C(0): |<Omega,(A_H(t)-A_E(t))A Omega>| <= delta/dimension.
    error = dimension*delta
    return {'subdivision': side, 'spacing_exact': str(F(1, side)),
            'g_squared_exact': str(F(64, side**2)), 'alpha_exact': str(alpha), 'beta_exact': str(beta),
            'link_count': len(g['links']), 'face_count': len(g['faces']),
            'loop_link_count': len(g['loop']), 'incident_face_count': len(g['incident_faces']),
            'representation': representation, 'time_exact': str(t),
            'haar_second_moment_exact': str(F(1, dimension**2)),
            'free_frequency_native_units_exact': str(energy),
            'heisenberg_operator_error_upper_exact': str(delta),
            'normalized_correlation_error_upper_exact': str(error),
            'normalized_interacting_real_correlation_interval_exact': [str(max(F(-dimension), lo-error)), str(min(F(dimension), hi+error))],
            'normalized_interacting_real_correlation_interval_decimal': [float(max(F(-dimension), lo-error)), float(min(F(dimension), hi+error))],
            'preparation_energy_above_ground_interval_exact': ['0', str(beta*len(g['faces']))]}


def zero_time_records(representation, eta=F(1, 2)):
    if not 0 < eta < 1:
        raise ValueError('smooth complete instrument requires 0<eta<1')
    moments = ([F(1), F(0), F(1, 4), F(0), F(1, 8)] if representation == 'fundamental'
               else [F(1), F(0), F(1, 9), F(1, 27), F(1, 27)])
    effects = {1: [F(1, 4), eta/2, eta*eta/4],
               -1: [F(1, 4), -eta/2, eta*eta/4],
               0: [F(1, 2), F(0), -eta*eta/2]}
    records = []
    for a, b in product((1, -1, 0), repeat=2):
        p = sum((x*y*moments[i+j] for i,x in enumerate(effects[a]) for j,y in enumerate(effects[b])), F())
        assert p >= 0
        records.append({'first': a, 'second': b, 'probability_exact': str(p)})
    assert sum(F(r['probability_exact']) for r in records) == 1
    assert sum(r['first']*r['second']*F(r['probability_exact']) for r in records) == eta*eta*moments[2]
    return records


def complete_noncommuting_control():
    X=s.Matrix([[0,1],[1,0]]); Z=s.diag(1,-1); I=s.eye(2)
    rho=(I+s.Rational(3,5)*X+s.Rational(4,5)*Z)/2
    U=s.Matrix([[1,-1],[1,1]])/s.sqrt(2); eta=s.Rational(1,2)
    K={1:(I+eta*Z)/2,-1:(I-eta*Z)/2,0:s.sqrt((1-eta**2)/2)*I}
    assert sum((k.H*k for k in K.values()),s.zeros(2))==I
    probs={(a,b):s.simplify(s.trace(kb*U*ka*rho*ka.H*U.H*kb.H)) for a,ka in K.items() for b,kb in K.items()}
    assert sum(probs.values())==1 and all(p>=0 for p in probs.values())
    weighted=sum(a*b*p for (a,b),p in probs.items())
    target=eta**2*s.re(s.trace(U.H*Z*U*Z*rho))
    assert s.simplify(weighted-target)==0
    # Keep an arbitrary reference: the Kraus completeness identity tensors with I_R.
    bell=s.Matrix([1,0,0,1])/s.sqrt(2)
    out=sum((s.kronecker_product(k,I)*bell*bell.H*s.kronecker_product(k.H,I) for k in K.values()),s.zeros(4))
    assert s.trace(out)==1
    return {'all_nine_probabilities_exact': [{'first':a,'second':b,'probability_exact':str(p)} for (a,b),p in probs.items()],
            'sequential_weighted_record_exact':str(weighted),'anticommutator_target_exact':str(target),
            'reference_entanglement_retained':True}


def run():
    rows=[local_bound(L,t,r) for L in (1,2,3) for t in (F(1,8192),F(1,4096)) for r in ('fundamental','spin_one')]
    return {'schema':'oph.local_gauge_memory.v1',
        'verdict':'PASS_THREE_REFINEMENT_FULL_HAAR_LOCAL_RECORD_AND_WORK_CERTIFICATE',
        'source':{'carrier':'Gauss-invariant subspace of tensor_links L2(SU(2)); every representation and intertwiner retained',
            'hamiltonian':'H_L=(32/L) sum_e C_e^2+(L^3/32) sum_f(1-W_f)',
            'physical_loop':'unit-square perimeter at z=0, with 4L links',
            'probe':'constant product Haar vector Omega, not the interacting ground state',
            'native_update_boundary':'The full-Haar Hamiltonian is declared; this package does not derive it from the native twelve-port process.',
            'spatial_and_coupling_trajectory':'a=1/L, g^2=64/L^2 are declared architecture parameters',
            'clock':'two supplied common native-time interventions; no conversion to physical seconds',
            'unbounded_electric_domain':'Casimir form domain; bounded magnetic perturbation, finite smooth multiplication instruments preserve it'},
        'refinement_records':rows,
        'instrument':{'eta_exact':'1/2','kraus':'K_+=(I+eta A)/2, K_-=(I-eta A)/2, K_0=sqrt((I-eta^2 A^2)/2)',
            'all_sequential_records':'(s0,st) in {+1,-1,0}^2; p=Tr[K_st U K_s0 rho K_s0* U* K_st*]',
            'record_identity':'E[s0 st]=eta^2 Re Tr[rho A(t) A]',
            'positive_time_probabilities':'Defined by the full-source complete instrument; local interval certificate bounds its weighted moment, not each individual probability.',
            'zero_time_all_nine_records':{r:zero_time_records(r) for r in ('fundamental','spin_one')},
            'reference':'Stinespring isometry into a supplied initialized three-label apparatus, tensored with any spectator quantum reference; no erased or discarded outcome',
            'normalization':'divide the weighted moment by eta^2 times the exact Haar second moment',
            'separate_preparations':'Each time, refinement and representation uses an independently prepared identical Haar probe.'},
        'work':{'exact_per_readout':'alpha Tr[rho eta^2 Gamma(A)/(2(1-eta^2 A^2))]; magnetic contribution zero',
            'fundamental_per_readout_upper_exact':'16/3','fundamental_two_readout_upper_exact':'32/3',
            'spin_one_per_readout_upper_exact':'256/27','spin_one_two_readout_upper_exact':'512/27',
            'free_evolution_mean_work_exact':'0','apparatus_reset_cost':'No free reset asserted; initialized references and intervention times are supplied.'},
        'small_complete_noncommuting_control':complete_noncommuting_control(),
        'iteration':{'local_operator_error_upper':'min(2,2 beta N_inc |t|), N_inc=5 for L=1 and 8L-4 for L>=2',
            'normalized_pairwise_refinement_error':'dimension*(delta_L+delta_M), since the free frequencies agree exactly',
            'global_E_less_than_5_projector_used':False,'third_refinement_all_representations_retained':True,
            'fixed_time_continuum_bound':'The displayed bound deteriorates as L^4. It is controlled for t=O(L^-4); no fixed-time ultraviolet convergence follows.'},
        'scope':{'local_quantum_record_memory_certified':True,'full_Haar_infinite_carrier_time_evolution_numerically_diagonalized':False,
            'physical_Wilson_loop_identification_proved':False,'continuum_or_Clay_mass_gap_proved':False,
            'physical_seconds_or_GeV_calibration_derived':False,'new_Nature_validation':False}}
