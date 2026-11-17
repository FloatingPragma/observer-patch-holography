"""Exact full Gauss-zero support count and adopted faithful Gibbs reference."""
from collections import defaultdict
import itertools
from pathlib import Path
import numpy as np
import sympy as sy
from scipy.linalg import expm
from ... import core, spatial
from ...joint_geometry import auxiliary_columns
from ...verify import read_result, scientific_inputs
PACKAGE=Path(__file__).resolve().parents[2]

def multiplicities(charges):
    result={0:1}
    for charge in charges:
        out=defaultdict(int)
        for q,count in result.items():out[q]+=count;out[q+charge]+=count
        result=dict(out)
    assert sum(result.values())==2**len(charges)
    return result


def gauss_rank(B, cap, counts):
    rank=0;by_flux=[]
    for flux in itertools.product(range(-cap,cap+1),repeat=B.shape[1]):
        local=[sum(int(B[v,e])*flux[e] for e in range(B.shape[1])) for v in range(B.shape[0])]
        count=1
        for q in local:count*=counts.get(q,0)
        rank+=count;by_flux.append((flux,count))
    return rank,by_flux


def actual_reference():
    charges=tuple(int(x) for x in core.native_charge_coordinates())
    counts=multiplicities(tuple(q for q in charges for _ in range(6)))
    parent=read_result(PACKAGE/'RESULTS.json')
    if parent.get('scientific_inputs_sha256') != scientific_inputs() or not parent.get('all_passed'):
        raise ValueError('Reference requires the current complete native receipt')
    acols,ncols=auxiliary_columns()
    A=max(sum(abs(c) for _,c in col) for col in acols)
    N=max(sum(abs(c) for _,c in col) for col in ncols)
    records=[]
    for graph,key in [(spatial.edge(),'source_edge'),(spatial.triangle(),'source_triangle')]:
        B,_=graph.incidence();rank,by_flux=gauss_rank(B,6,counts)
        degree=[sum(v in e for e in graph.edges) for v in range(graph.nvertices)]
        M_source=parent[key]['complete_source_column_execution']['joint_source_short_time_witness']['source_operator_norm_upper_integer']
        assert parent[key]['complete_source_column_execution']['flux_cap']==6
        # ||g_v|| <=4 hopping + stiffness*(regulator/2+degree) at Phi∈[-1,1].
        local_bounds=[sy.Rational(9,4)+sy.Rational(d,2) for d in degree]
        M=int(sy.ceiling(sum(A+N*(g+sy.Rational(5,4)*M_source) for g in local_bounds)))
        geometry_dimension=48**graph.nvertices
        fullrank=rank*geometry_dimension
        beta=sy.Rational(1,1024)
        assert 0<rank<=2**(96*graph.nvertices)*13**len(graph.edges)
        records.append({'graph':graph.name,'source_digest':parent[key]['combined_source']['matter_digest'],
                        'matter_field_gauss_rank':str(rank),'geometry_dimension':geometry_dimension,
                        'complete_physical_rank':str(fullrank),'local_source_norm_bound':M_source,
                        'capacity_auxiliary_norm_upper':str(A),'capacity_lapse_norm_upper':str(N),
                        'complete_local_capacity_norm_upper_integer':M,'beta_exact':str(beta),
                        'faithfulness_floor_exact':str(sy.exp(-2*beta*M)/fullrank),
                        'flux_assignments_counted':len(by_flux),
                        'full_Fock_enumeration':False,'physical_gauss_support_nonempty':True})
    return {'actual_native_charge_coordinates':list(charges),'local_charge_multiplicities':{str(q):n for q,n in sorted(counts.items())},
            'local_mode_count':96,'local_dimension_exact':str(sum(counts.values())),'members':records}


def hermitian_log(rho):
    values,vectors=np.linalg.eigh(rho)
    assert min(values)>0
    return (vectors*np.log(values))@vectors.conj().T,values


def fixture():
    B=np.array([[1],[-1]]);charges=(-1,1);counts=multiplicities(charges)
    rank,_=gauss_rank(B,1,counts)
    occupation_charges=[sum(charges[j] for j in range(2) if mask>>j&1) for mask in range(4)]
    labels=list(itertools.product(range(4),range(4),range(-1,2)))
    keep=np.array([e==occupation_charges[a] and -e==occupation_charges[b] for a,b,e in labels])
    assert int(keep.sum())==rank
    # The omitted-charge mutation must not pass the literal Gauss support census.
    wrong=np.array([e==occupation_charges[a] for a,b,e in labels])
    assert np.any(wrong!=keep)
    pair=np.zeros((4,4));pair[0,3]=pair[3,0]=.25
    local=np.diag([0.,1.,1.5,2.5])+pair
    H=np.kron(np.kron(local,np.eye(4)),np.eye(3))+np.kron(np.kron(np.eye(4),local),np.eye(3))
    H+=np.kron(np.eye(16),np.diag([.5,0.,.5]))
    P=np.diag(keep.astype(float));assert np.max(np.abs(H@P-P@H))<1e-10
    beta=.5; M=6
    assert np.max(np.abs(np.linalg.eigvalsh(H)))<=M
    rho=P@expm(-beta*H);rho/=np.trace(rho)
    floor=np.exp(-2*beta*M)/rank
    support=rho[np.ix_(keep,keep)]
    assert min(np.linalg.eigvalsh(support))>=floor-1e-12
    assert np.linalg.norm(H@rho-rho@H)<1e-10
    rho_R=np.trace(rho.reshape(4,12,4,12),axis1=1,axis2=3)
    partial_P=np.trace(P.reshape(4,12,4,12),axis1=1,axis2=3)
    assert np.array_equal(partial_P,np.diag(np.diag(partial_P)))
    assert np.all(np.diag(partial_P)==np.round(np.diag(partial_P)))
    log_sigma,eigen=hermitian_log(rho_R);lower=float(min(eigen))
    assert lower>=floor-1e-12
    delta=np.zeros((4,4));delta[0,3]=delta[3,0]=lower/8
    tau=rho_R+delta;log_tau,_=hermitian_log(tau)
    entropy_change=float((-np.trace(tau@log_tau)+np.trace(rho_R@log_sigma)).real)
    first_law=float(-np.trace(delta@log_sigma).real)
    relative=float(np.trace(tau@(log_tau-log_sigma)).real)
    remainder_bound=float(np.linalg.norm(delta,'fro')**2/lower)
    assert abs(entropy_change-first_law+relative)<1e-10
    assert -1e-12<=relative<=remainder_bound+1e-10
    return {'full_dimension':len(labels),'exact_physical_rank':rank,'wrong_gauss_mutation_rejected':True,
            'beta':beta,'operator_norm_upper':M,
            'full_reference':rho.tolist(),'regional_reference':rho_R.tolist(),
            'regional_support_multiplicities':np.diag(partial_P).tolist(),
            'regional_eigenvalues':eigen.tolist(),'global_floor':floor,
            'entropy_change':entropy_change,'modular_first_law':first_law,'relative_entropy_remainder':relative,
            'remainder_upper_bound':remainder_bound,'all_passed':True}



def stationarity_boundary(phase):
    z=sy.sympify(phase['local_capacity_entry_exact'])
    beta=sy.Rational(1,1024)
    derivative=sy.expand(2*beta*z)
    assert z==sy.Rational(9,32) and derivative==sy.Rational(9,16384)
    return {'beta_exact':str(beta),'offdiagonal_entry_exact':str(z),
            'fixed_diagonal_area_volume_variation':True,
            'entropy_plus_area_first_variation_exact':str(derivative),
            'unrestricted_full_Gibbs_diagonal_area_stationarity_rejected':True,
            'regional_Einstein_stationarity_established':False}
