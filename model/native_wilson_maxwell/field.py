"""Finite open-cube cochains and the shared cutoff transverse field."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from model.affine_gauss_source.model import cutoff_control,green,simple
from model.affine_gauss_mobile.model import total_occupations

@lru_cache(None)
def cochains(n=1):
    """n open cells per axis; D is head minus tail, faces use i then j."""
    if type(n) is not int or n<1:raise ValueError('positive open-box cell count required')
    vertices=tuple(product(range(n+1),repeat=3));index={x:i for i,x in enumerate(vertices)}
    edges=[];edge_index={}
    for x in vertices:
        for axis in range(3):
            if x[axis]==n:continue
            y=tuple(a+int(j==axis) for j,a in enumerate(x))
            edge_index[x,axis]=len(edges);edges.append((index[x],index[y],axis))
    D=s.zeros(len(edges),len(vertices))
    for e,(tail,head,axis) in enumerate(edges):D[e,tail]=-1;D[e,head]=1
    faces=[];rows=[]
    for x in vertices:
        for i in range(3):
            for j in range(i+1,3):
                if x[i]==n or x[j]==n:continue
                xi=tuple(a+int(k==i) for k,a in enumerate(x))
                xj=tuple(a+int(k==j) for k,a in enumerate(x))
                faces.append((x,i,j));row=[0]*len(edges)
                for coordinate,value in (((x,i),1),((xi,j),1),((xj,i),-1),((x,j),-1)):
                    row[edge_index[coordinate]]=value
                rows.append(row)
    C=s.Matrix(rows)
    assert D.shape==(3*n*(n+1)**2,(n+1)**3)
    assert C.shape==(3*n*n*(n+1),D.rows) and C*D==s.zeros(C.rows,D.cols)
    return vertices,tuple(edges),tuple(faces),s.ImmutableMatrix(D),s.ImmutableMatrix(C)

@lru_cache(None)
def volume_boundary(n=1):
    vertices,edges,faces,D,C=cochains(n);index={face:i for i,face in enumerate(faces)}
    volumes=tuple(product(range(n),repeat=3));boundary=s.zeros(len(volumes),len(faces))
    for row,x in enumerate(volumes):
        for axis,(i,j),sign in ((0,(1,2),1),(1,(0,2),-1),(2,(0,1),1)):
            shifted=tuple(a+int(k==axis) for k,a in enumerate(x))
            boundary[row,index[(shifted,i,j)]]=sign;boundary[row,index[(x,i,j)]]=-sign
    assert boundary*C==s.zeros(len(volumes),len(edges))
    return s.ImmutableMatrix(boundary)

@lru_cache(None)
def cube():
    vertices,edges,faces,D,C=cochains(1);L=D.T*D;K=C.T*C
    inverse=green.invert([[F(value+1) for value in row] for row in L.tolist()])
    G=s.Matrix(inverse)-s.ones(8)/64
    assert L*G==s.eye(8)-s.ones(8)/8 and G==G.T and G*s.ones(8,1)==s.zeros(8,1)
    PT=s.eye(12)-D*G*D.T;P4=(6*PT-K)/2;P6=(K-4*PT)/2
    assert D.rank()==7 and C.rank()==5 and PT.trace()==5
    assert PT==PT.T and PT*PT==PT and D.T*PT==s.zeros(8,12)
    assert P4==P4.T and P6==P6.T and P4*P4==P4 and P6*P6==P6
    assert P4*P6==s.zeros(12) and P4+P6==PT and K==4*P4+6*P6
    assert P4.trace()==3 and P6.trace()==2
    x=s.Symbol('x');assert s.expand(L.charpoly(x).as_expr()-x*(x-2)**3*(x-4)**3*(x-6))==0
    covariance=P4/4+P6/(2*s.sqrt(6))
    variance=s.Rational(1,16)+1/(12*s.sqrt(6))
    assert set(P4.diagonal())=={s.Rational(1,4)} and set(P6.diagonal())=={s.Rational(1,6)}
    assert all(simple(value-variance)==0 for value in covariance.diagonal())
    return {'vertices':vertices,'edges':edges,'faces':faces,'D':D,'C':C,'D2':volume_boundary(1),'G':s.ImmutableMatrix(G),
      'K':s.ImmutableMatrix(K),'PT':s.ImmutableMatrix(PT),'P4':s.ImmutableMatrix(P4),
      'P6':s.ImmutableMatrix(P6),'Q_vacuum_covariance':s.ImmutableMatrix(covariance),
      'Q_edge_vacuum_variance':variance}

@lru_cache(None)
def frame():
    data=cube();columns=[];frequencies=[]
    for P,frequency in ((data['P4'],s.Integer(2)),(data['P6'],s.sqrt(6))):
        local=[]
        for column in P.columnspace():
            vector=column
            for u in local:vector=vector-u*(u.T*vector)[0]
            vector=vector/s.sqrt((vector.T*vector)[0]);vector=vector.applyfunc(simple)
            assert (vector.T*vector)[0]==1
            local.append(vector)
        columns.extend(local);frequencies.extend([frequency]*len(local))
    V=s.Matrix.hstack(*columns)
    assert V.T*V==s.eye(5) and V*V.T==data['PT']
    assert V.T*data['K']*V==s.diag(*(w*w for w in frequencies))
    return s.ImmutableMatrix(V),tuple(frequencies)

@lru_cache(None)
def cutoff(cap=3):
    if type(cap) is not int or cap<1:raise ValueError('positive total photon cutoff required')
    states=tuple(total_occupations(5,cap));lookup={n:i for i,n in enumerate(states)}
    assert len(states)==math.comb(5+cap,cap)
    annihilators=[]
    for j in range(5):
        entries={}
        for i,n in enumerate(states):
            if n[j]:
                target=list(n);target[j]-=1;entries[lookup[tuple(target)],i]=s.sqrt(n[j])
        annihilators.append(s.ImmutableSparseMatrix(len(states),len(states),entries))
    V,nu=frame();dimension=len(states)
    H=s.diag(*(sum(n[j]*nu[j] for j in range(5)) for n in states))
    vacuum=sum(nu)/2;top=s.diag(*(int(sum(n)==cap) for n in states))
    correction=top*(H+2*vacuum*s.eye(dimension))/2
    bare=s.zeros(dimension)
    for i,ai in enumerate(annihilators):
        for j,aj in enumerate(annihilators):
            delta=s.eye(dimension) if i==j else s.zeros(dimension)
            assert ai*aj.H-aj.H*ai==delta-top*(delta+aj.H*ai)
    for j,a in enumerate(annihilators):
        bare+=nu[j]*(a*a.H+a.H*a)/2
        q=(a+a.H)/s.sqrt(2*nu[j]);p=s.I*s.sqrt(nu[j]/2)*(a.H-a)
        assert (H*q-q*H+s.I*p).applyfunc(simple)==s.zeros(dimension)
        assert (H*p-p*H-s.I*nu[j]**2*q).applyfunc(simple)==s.zeros(dimension)
    assert bare+correction==H+vacuum*s.eye(dimension)
    if cap==3:
        # The complete total-cap carrier is invariant under mode rotations.
        # A single linear coordinate reduces to blocks of sizes1,2,3,4.
        coordinate_bound=s.simplify((3+s.sqrt(6))*cube()['Q_edge_vacuum_variance'])
        assert s.simplify(s.Rational(9,16)-coordinate_bound)>0
        inherited=cutoff_control(4)
        assert inherited['exact_Heisenberg_equations_for_finite_q_p']
    else:coordinate_bound=None
    return {'states':states,'annihilators':tuple(annihilators),'H_no_vacuum_at_a1':s.ImmutableMatrix(H),
      'vacuum_energy_at_a1':vacuum,'positive_boundary_completion_at_a1':s.ImmutableMatrix(correction),
      'Q_edge_norm_squared_at_cap3':coordinate_bound}

@lru_cache(None)
def link_coordinates():
    V,nu=frame();cut=cutoff(3);aa=cut['annihilators'];dimension=len(cut['states'])
    return tuple(s.ImmutableMatrix(sum((V[e,j]*(aa[j]+aa[j].H)/s.sqrt(2*nu[j])
      for j in range(5)),s.zeros(dimension))) for e in range(12))

def boundary_counterexample():
    coords=link_coordinates();commutator=(coords[0]*coords[1]-coords[1]*coords[0]).applyfunc(simple)
    assert commutator!=s.zeros(commutator.rows)
    states=cutoff(3)['states']
    assert all(sum(states[i])==sum(states[j])==3 for i in range(len(states))
      for j in range(len(states)) if commutator[i,j]!=0)
    i,j=next((i,j) for i in range(len(states)) for j in range(len(states)) if commutator[i,j]!=0)
    return {'edges':[0,1],'bra_occupation':states[i],'ket_occupation':states[j],
      'exact_Q_commutator_matrix_element':commutator[i,j],
      'all_nonzero_entries_supported_on_total_cutoff_boundary':True}

def certificate():
    data=cube();cut=cutoff(3);D,C,G=data['D'],data['C'],data['G']
    reference=0;Q=s.eye(8)-s.eye(8)[:,reference]*s.ones(1,8);electric=D*G*Q
    assert D.T*electric==Q and C*electric==s.zeros(6,8)
    assert electric.T*electric==Q.T*G*Q
    return {'schema':'oph.open_cube_finite_transverse_field.v1','vertices':data['vertices'],
      'edges_tail_head_axis':data['edges'],'faces_base_axes':data['faces'],
      'D':D.tolist(),'C':C.tolist(),'G':G.tolist(),'PT':data['PT'].tolist(),
      'oriented_volume_boundary':data['D2'].tolist(),'magnetic_divergence_D2_C_exactly_zero':True,
      'squared_dimensionless_frequencies':{'4':3,'6':2},'transverse_dimension':5,
      'harmonic_one_form_dimension':0,'constant_vertex_mode':'retained gauge redundancy; G kills constants',
      'Q_vacuum_covariance':data['Q_vacuum_covariance'].tolist(),
      'edge_Q_variance':data['Q_edge_vacuum_variance'],
      'exact_edge_Q_norm_squared_at_cap3':cut['Q_edge_norm_squared_at_cap3'],
      'edge_Q_norm_strict_upper':s.Rational(3,4),'total_photon_cutoff':3,'photon_dimension':56,
      'free_H':'[sum_j nu_j/2 + dGamma(nu)]/a, nu=(2,2,2,sqrt6,sqrt6)',
      'boundary_completion':'P_N=3 [dGamma(nu)+Tr(nu)]/(2a)',
      'link_variable':'Q_e=a A_phys,e','electric_flux':'P_e=a^2 E_phys,e',
      'normalized_Gauss':'D^T Ebar_L-rho-R delta_r=0, independent of g',
      'physical_electric_flux':'g Ebar_L+P_transverse','Coulomb':'g^2 rho_ref^T G rho_ref/(2a)',
      'matter_link_phase':'exp(i g Q_e) multiplies Psi_head^dagger B_e^dagger Psi_tail',
      'full_native_density_range':[-4,4],'full_longitudinal_label_range':[-4,4],
      'full_static_compensator_range':[-32,32],
      'full_static_compensator_energy':'4 abs(R), including every retained charge sector',
      'finite_photon_ladder_commutator':'[a_i,a_j^dag]=delta_ij I-P_N=3 (delta_ij I+a_j^dag a_i)',
      'link_coordinate_boundary_counterexample':boundary_counterexample(),
      'commuting_abelian_link_coordinates_claimed_on_full_cutoff':False,
      'large_occupation_CCR_or_continuum_carrier_claimed':False,
      'native_PW119_H_unchanged_or_source_transfer_claimed':False}

if __name__=='__main__':
    import json
    value=certificate();print(json.dumps(value,indent=2,default=str))
