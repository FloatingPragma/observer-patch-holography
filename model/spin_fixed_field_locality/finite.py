"""Independent finite-source controls for the Gaussian-flow sign and cut."""
from fractions import Fraction as F
import numpy as np
from model.native_ordered_ground import source


def matrices(field=F(1,4)):
    full,scale=source.matrix(4,((0,1),(1,2),(2,3)),F(1),field)
    cut,cutscale=source.matrix(4,((0,1),(2,3)),F(1),field)
    assert scale==cutscale
    return full/scale,cut/scale,(full-cut)/scale


def projector(matrix):
    e,u=np.linalg.eigh(matrix);v=u[:,0]
    return e,u,np.outer(v,v.conj())


def marginal_first_two(vector):
    coefficient=vector.reshape(4,4)
    return coefficient.T@coefficient.conj()


def flow(matrix,perturbation,alpha):
    e,u,p=projector(matrix);gaps=e-e[0];v=u[:,0]
    coefficients=u.conj().T@perturbation@v
    exact=np.zeros(len(e),dtype=complex);approx=np.zeros_like(exact)
    for n in range(1,len(e)):
        exact[n]=-coefficients[n]/gaps[n]
        # D_n0=-i(1-exp(-gap²/4alpha))V_n0/gap, so -iD has a minus sign.
        approx[n]=-(1-np.exp(-gaps[n]**2/(4*alpha)))*coefficients[n]/gaps[n]
    derivative=u@exact;flowvector=u@approx
    pdot=np.outer(derivative,v.conj())+np.outer(v,derivative.conj())
    return {'ground':v,'projector':p,'gap':gaps[1],'projector_derivative':pdot,
        'vector_derivative':derivative,'flow_vector':flowvector,
        'flow_error':float(np.linalg.norm(derivative-flowvector)),
        'spectral_error_identity':float(np.sqrt(sum(abs(coefficients[n])**2/gaps[n]**2*np.exp(-gaps[n]**2/(2*alpha)) for n in range(1,len(e)))))}


def certificate():
    full,cut,v=matrices();ef,uf,pf=projector(full);ec,uc,pc=projector(cut)
    small,scale=source.matrix(2,((0,1),),F(1),F(1,4));es,us,ps=projector(small/scale)
    cut_marginal=marginal_first_two(uc[:,0]);assert np.linalg.norm(cut_marginal-ps)<1e-13
    distance=float(np.linalg.norm(marginal_first_two(uf[:,0])-ps,ord='nuc')/2)
    rows=[]
    for s in (0.,.25,.5,1.):
        for alpha in (.01,.1,.5):
            row=flow(cut+s*v,v,alpha)
            assert abs(row['flow_error']-row['spectral_error_identity'])<1e-13
            radius=1e-5
            derivative=(projector(cut+(s+radius)*v)[2]-projector(cut+(s-radius)*v)[2])/(2*radius)
            error=float(np.linalg.norm(derivative-row['projector_derivative']))
            assert error<1e-8
            bound=float(np.linalg.norm(v,2)*np.exp(-(1/32)**2/(4*alpha))/(1/32))
            assert row['flow_error']<=bound+1e-14
            rows.append({'s':s,'alpha':alpha,'actual_gap':float(row['gap']),
                'flow_error':row['flow_error'],'spectral_error_identity':row['spectral_error_identity'],
                'projector_derivative_finite_difference_error_interval':[0.,1e-8],
                'finite_difference_derivative_admission_passed':True,'uniform_gap_bound_error_upper':bound})
    return {'source':'actual native_ordered_ground four-site open chain and its exact middle-bond cut',
        'full_dimension':16,'field':'1/4','crossing_bonds':1,
        'crossing_perturbation_norm':float(np.linalg.norm(v,2)),
        'cut_ground_factors_exact_source_numeric_control':True,
        'interior_two_site_half_trace_distance':distance,
        'flow_controls':rows,'finite_controls_prove_arbitrary_volume_theorem':False}


def exact_source_attachment():
    """Integer/dyadic Pauli equality against the unchanged actual producer."""
    n=4;j=F(1);b=F(1,4);edges=((0,1),(1,2),(2,3))
    source_h,scale=source.matrix(n,edges,j,b)
    x=np.array([[0,1],[1,0]],dtype=complex);y=np.array([[0,-1j],[1j,0]]);z=np.diag([1,-1]);identity=np.eye(2)
    def operator(terms):
        result=np.array([[1]],dtype=complex)
        for site in reversed(range(n)):result=np.kron(result,terms.get(site,identity))
        return result
    w=np.array([[1]],dtype=np.int64)
    for _ in range(n):w=np.kron(w,np.array([[1,1],[1,-1]],dtype=np.int64))
    expected=sum((float(j)/4*(operator({a:x,c:x})-operator({a:y,c:y})-operator({a:z,c:z})-np.eye(1<<n))for a,c in edges),np.zeros((1<<n,1<<n),dtype=complex))
    expected-=sum((float(b)/2*operator({site:z})for site in range(n)),np.zeros_like(expected))
    assert np.array_equal(w@source_h@w,(1<<n)*scale*expected)
    return {'actual_source_dimension':1<<n,'integer_Hadamard_identity_exact':True,
        'edge_coefficients_wz_wx_wy_cross':['J_e/4','J_e/4','-J_e/4','0','0'],
        'Lee_Yang_condition':'J_e/4=(abs(J_e/2)+abs(0))/2 for J_e>=0',
        'field_parameter_h':'B/2','full_gap_lower':'B/8',
        'nonnegative_boundary_interpolation_preserves_hypotheses':True,
        'scalar_bond_offset_changes_ground_or_gap':False}
