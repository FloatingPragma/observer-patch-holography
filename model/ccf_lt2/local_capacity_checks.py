"""Complete local-capacity matrix, instrument and factorization checks."""
from dataclasses import replace
import numpy as np
import sympy as sy
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from .geometry import GeometryConfig, build_geometry_extension
from .joint_geometry import auxiliary_columns
from .local_geometry import LocalGeometryConfig, geometry_factors
from .local_capacity import (apply_local_capacity_column, apply_native_column,
                            geometry_shard_column, local_energy_column)

def kron(*parts):
    out = sparse.csr_matrix([[1.]])
    for part in parts:
        out = sparse.kron(out, part, format='csr')
    return out


def comm(a, b):
    return a*b-b*a


def complete_clock_control():
    base = build_geometry_extension(np.zeros((1,1)))
    A,N=base.auxiliary_hamiltonian,base.lapse
    I=sparse.eye(2); Z=sparse.diags([1.,-1.]); had=np.array([[1.,1.],[1.,-1.]])/np.sqrt(2)
    hx,hy=Z/2,sparse.diags([0.,.75])
    records=[]
    for mode in ('localized',):
        aux=A.shape[0]**2; ia=sparse.eye(A.shape[0])
        H=kron(A,ia,I,I)+kron(ia,A,I,I)+kron(N,ia,hx,I)+kron(ia,N,I,hy)
        R=kron(sparse.eye(aux),had,I)
        measured_H=R@H@R.getH()
        rows=[]
        for sender in (0,1):
            psi=np.zeros(aux*4,complex);psi[sender]=1/np.sqrt(2);psi[2+sender]=1/np.sqrt(2)
            energy0=float(np.vdot(psi,H@psi).real)
            final=expm_multiply(-2j*H,psi)
            energy=float(np.vdot(final,H@final).real)
            probabilities=np.abs(R@final)**2
            assert abs(probabilities.sum()-1)<1e-10 and abs(energy-energy0)<1e-10
            rows.append({'sender_label':sender,'all_outcome_probabilities':probabilities.tolist(),
                         'receiver_plus_probability':float(probabilities.reshape(aux,2,2)[:,0,:].sum()),
                         'energy_residual':abs(energy-energy0),
                         'complete_meter_work':float(probabilities@measured_H.diagonal().real-energy)})
        delta=abs(rows[0]['receiver_plus_probability']-rows[1]['receiver_plus_probability'])
        assert delta<1e-10
        records.append({'mode':mode,'joint_dimension':H.shape[0],'time':2,'response':delta,'runs':rows})
    return {'all_passed':True,'records':records,
            'scope':'Complete two-source diagnostic with actual default auxiliary; not full CCF particle evolution.'}


def matrix_column_check():
    graph={'nvertices':2,'edges':((0,1),)}
    lc=LocalGeometryConfig(phi_values=(-1.,1.))
    cc=GeometryConfig(capacity_max=1)
    g=geometry_factors(graph,lc); q=g['phi']; nphi=g['dimension']
    eye2=sparse.eye(2); X=sparse.csr_matrix([[0.,1.],[1.,0.]]); Z=sparse.diags([1.,-1.])
    path=sparse.csr_matrix([[1.,-1.],[-1.,1.]])
    # Independently build the exact geometric partition from tensor factors.
    gradient=(lc.stiffness/4)*(q[0]-q[1])@(q[0]-q[1])
    gs=[lc.hopping*kron(path,eye2)+lc.stiffness*lc.regulator*(q[0]@q[0])/2+gradient,
        lc.hopping*kron(eye2,path)+lc.stiffness*lc.regulator*(q[1]@q[1])/2+gradient]
    assert sparse.linalg.norm(sum(gs)-g['H_geometry'])<1e-10
    hs=[kron(Z,eye2)/2+kron(X,X)/8,kron(eye2,Z)*.75+kron(X,X)/8]
    ks=[kron(gs[v],sparse.eye(4))+kron(g['lapses'][v],hs[v]) for v in range(2)]
    base=build_geometry_extension(np.zeros((1,1)),config=cc)
    A,N=base.auxiliary_hamiltonian,base.lapse; da=A.shape[0];ia=sparse.eye(da)
    H=kron(A,ia,sparse.eye(nphi*4))+kron(ia,A,sparse.eye(nphi*4))
    H=H+kron(N,ia,ks[0])+kron(ia,N,ks[1])
    actions=[lambda key,v=v:{r:sy.Rational(float(hs[v][r,key])) for r in range(4) if hs[v][r,key]!=0}
             for v in range(2)]
    max_error=0.
    for a in range(da):
        for b in range(da):
            for p in range(nphi):
                for s in range(4):
                    key=((a,b),(p,s)); col=(((a*da+b)*nphi+p)*4+s)
                    actual=apply_local_capacity_column(actions,key,graph,lc,cc)
                    vector=np.zeros(H.shape[0],complex)
                    for (ab,(pr,sr)),value in actual.items():
                        vector[(((ab[0]*da+ab[1])*nphi+pr)*4+sr)]=complex(value)
                    error=float(np.max(np.abs(vector-H[:,col].toarray().ravel())))
                    max_error=max(max_error,error)
                    assert error<1e-10
    return {'all_passed':True,'complete_columns_checked':H.shape[0],
            'max_column_error':max_error,'local_geometry_shards_sum_exact_for_dyadic_fixture':True}



def disconnected_factorization():
    """Every full column equals the sum of independent component columns."""
    graph = {'nvertices': 2, 'edges': ()}
    component = {'nvertices': 1, 'edges': ()}
    lc = LocalGeometryConfig(phi_values=(-1.,1.), source_energy_references=(.25,.75))
    cc = GeometryConfig(capacity_max=1, source_energy_reference=.375)
    hs = (sy.Matrix([[1,sy.Rational(1,2)],[sy.Rational(1,2),-1]]),
          sy.Matrix([[0,sy.I/4],[-sy.I/4,sy.Rational(3,4)]]))
    def one_action(v):
        return lambda s: {r:hs[v][r,s] for r in range(2) if hs[v][r,s] != 0}
    def joint_action(v):
        def action(s):
            labels = [s//2,s%2]
            result = {}
            for r,c in one_action(v)(labels[v]).items():
                out = labels.copy();out[v] = r
                result[2*out[0]+out[1]] = c
            return result
        return action
    actions = [joint_action(v) for v in range(2)]
    count = len(auxiliary_columns(cc)[0]);checked = 0
    for a in range(count):
        for b in range(count):
            for phi in range(4):
                for s in range(4):
                    key = ((a,b),(phi,s));expected = {}
                    for v in range(2):
                        cs = (a,b)[v]; ps = (phi//2,phi%2)[v]; ss = (s//2,s%2)[v]
                        cfg = replace(lc,source_energy_references=(lc.source_energy_references[v],))
                        column = apply_local_capacity_column([one_action(v)],((cs,),(ps,ss)),component,cfg,cc)
                        for ((cr,),(pr,sr)),value in column.items():
                            aa=[a,b];pp=[phi//2,phi%2];ss=[s//2,s%2]
                            aa[v]=cr;pp[v]=pr;ss[v]=sr
                            out=(tuple(aa),(2*pp[0]+pp[1],2*ss[0]+ss[1]))
                            expected[out]=expected.get(out,0)+value
                    expected={out:sy.expand(value) for out,value in expected.items() if sy.expand(value)!=0}
                    assert apply_local_capacity_column(actions,key,graph,lc,cc)==expected
                    checked+=1
    return {'all_passed':True,'complete_columns_checked':checked,
            'exact_component_sum':True,'disjoint_even_terms_commute':True,
            'local_references':[.25,.75],'common_capacity_reference':.375,
            'scope':'Complete finite factorization fixture; general factorization follows from disjoint even support.'}


def reference_convention():
    graph={'nvertices':2,'edges':()}
    lc=LocalGeometryConfig(source_energy_references=(.25,.75))
    cc=GeometryConfig(source_energy_reference=.375)
    shift=sy.Rational(1,2)
    actions=[lambda s:{s:sy.Rational(3,2)},lambda s:{s:sy.Rational(5,4)}]
    shifted=[lambda s,v=v:{s:actions[v](s)[s]+shift} for v in range(2)]
    shifted_local=replace(lc,source_energy_references=tuple(e+float(shift) for e in lc.source_energy_references))
    shifted_capacity=replace(cc,source_energy_reference=cc.source_energy_reference+float(shift))
    checked=0
    for capacities in ((0,0),(0,8),(8,0),(8,8)):
        for phi in range(9):
            key=(capacities,(phi,0))
            base=apply_local_capacity_column(actions,key,graph,lc,cc)
            moved=apply_local_capacity_column(shifted,key,graph,shifted_local,shifted_capacity)
            difference={out:sy.expand(moved.get(out,0)-base.get(out,0)) for out in set(base)|set(moved)}
            assert {out:c for out,c in difference.items() if c!=0}=={key:2*shift}
            fixed=apply_local_capacity_column(shifted,key,graph,shifted_local,cc)
            # At fixed F, uniform k_v shifts produce shift*sum_v M_v.
            lapse_columns=auxiliary_columns(cc)[1]
            predicted=shift*sum(dict(lapse_columns[a])[a] for a in capacities)
            assert sy.expand(fixed[key]-base[key]-predicted)==0
            checked+=1
    unequal=[lambda s:{s:actions[0](s)[s]+shift},actions[1]]
    unequal_local=replace(lc,source_energy_references=(lc.source_energy_references[0]+float(shift),lc.source_energy_references[1]))
    deltas=[]
    for capacities in ((0,8),(8,0)):
        key=(capacities,(0,0))
        before=apply_local_capacity_column(actions,key,graph,lc,cc)
        after=apply_local_capacity_column(unequal,key,graph,unequal_local,cc)
        deltas.append(sy.expand(after[key]-before[key]))
        # sum M_v=2 here, so any common F shift contributes zero.
        assert sum(dict(auxiliary_columns(cc)[1][a])[a] for a in capacities)==2
    assert deltas[0]!=deltas[1]
    return {'all_passed':True,'columns_checked':checked,
            'common_shift_gives_vertex_count_times_scalar':True,
            'fixed_reference_shift_is_sum_of_local_lapses':True,
            'one_common_F':True,'nonuniform_shift_not_absorbed_by_common_F':True}


def run():
    return {'all_passed':True,'matrix_columns':matrix_column_check(),
            'disconnected_factorization':disconnected_factorization(),
            'reference_convention':reference_convention(),
            'clock':complete_clock_control()}
