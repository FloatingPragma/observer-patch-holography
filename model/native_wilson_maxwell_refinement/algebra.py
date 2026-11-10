"""Exact source-local cochain identities for the common comparison algebra."""
import sympy as s
from model.native_wilson_maxwell import field


def calculate():
    c=field.cube();D,C,G,PT=[c[key] for key in ('D','C','G','PT')]
    for e,(tail,head,axis) in enumerate(c['edges']):
        assert list(D[e,:])==[-int(x==tail)+int(x==head) for x in range(8)]
    assert list(C[0,:])==[1,-1,0,0,0,-1,0,0,1,0,0,0]
    PL=D*G*D.T
    assert PT+PL==s.eye(12)
    assert C*PL==s.zeros(6,12) and C*PT==C
    assert D.T*PT==s.zeros(8,12)
    assert PL*PL==PL and PT*PL==s.zeros(12)
    assert PL.T==PL and PT.T==PT
    reference=s.eye(8)-s.eye(8)[:,0]*s.ones(1,8)
    ell=D*G*reference
    assert ell.T*ell==reference.T*G*reference
    assert D.T*ell==reference
    # Conserved charge implies reference has zero commutator with an edge hop.
    assert reference*D.T==D.T
    assert ell*D.T==PL
    assert ell*D.T+PT==s.eye(12)
    # Arbitrary physical electric/magnetic Weyl strengths can use these norms.
    V,nu=field.frame()
    electric=s.simplify(V*s.diag(*(w/2 for w in nu))*V.T)
    magnetic=s.simplify(C*c['Q_vacuum_covariance']*C.T)
    electric_diag=[s.simplify(electric[i,i]) for i in range(12)]
    magnetic_diag=[s.simplify(magnetic[i,i]) for i in range(6)]
    assert len(set(electric_diag))==1 and len(set(magnetic_diag))==1
    assert s.Rational(1,2)-electric_diag[0]>0 and 2-magnetic_diag[0]>0
    edge,face=0,3
    assert C[face,edge]==0
    q=s.Matrix([s.simplify((C*V)[face,j]/s.sqrt(2*nu[j])) for j in range(5)])
    t=s.Matrix([s.simplify(V[edge,j]*s.sqrt(nu[j]/2)) for j in range(5)])
    assert s.simplify((q.T*t)[0])==0
    boundary=(q*t.T+t*q.T).applyfunc(s.simplify)
    assert boundary!=s.zeros(5)
    assert s.simplify(s.trace(boundary*boundary)-2*magnetic_diag[0]*electric_diag[0])==0
    return {'edges':12,'vertices':8,'faces':6,
        'oriented_incidence':'D_e,tail=-1; D_e,head=+1; edges positive coordinate direction',
        'oriented_face_boundary':'+(x,i)+(x+e_i,j)-(x+e_j,i)-(x,j), i<j',
        'exact_transverse_plus_longitudinal_identity':True,
        'exact_dressed_electric_hop_commutator_coefficient':'PT+D G D^T=I_12',
        'electric_hop_relation':'[E_e,F_f]=g delta_ef F_f',
        'magnetic_electric_relation':'[B_f,E_e]=i C_fe',
        'charge_hop_relation':'[rho_x,F_e]=D_ex F_e',
        'electric_charge_relation':'[E_e,rho_x]=0',
        'electric_vacuum_variance':str(electric_diag[0]),
        'magnetic_vacuum_variance':str(magnetic_diag[0]),
        'nonincident_order_test':{'electric_edge_index':edge,'magnetic_face_index':face,
            'exact_incidence':0,'unbounded_controlled_Weyl_tests_commute':True,
            'finite_linear_generator_commutator_boundary_supported':True,
            'finite_generator_commutator_norm':'K sqrt(var_B var_E)',
            'finite_generator_commutator_norm_coefficient_squared':str(s.simplify(magnetic_diag[0]*electric_diag[0])),
            'operator_norm_local_commutator_convergence_claimed':False},
        'physical_field_energy':'m/2 sum_e E_e^2 + m/2 sum_f B_f^2',
        'exact_Gauss_and_electric_energy_decomposition':True,
        'locality_scope':'exact reducing native occupation code of the unbounded comparison; retained full-source complement is not assigned this algebra',
        'full_native_or_uniform_spatial_Lieb_Robinson_cone_proved':False}
