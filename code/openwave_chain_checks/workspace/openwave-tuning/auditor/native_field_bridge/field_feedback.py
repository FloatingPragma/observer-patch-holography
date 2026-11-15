"""Finite local-field feedback on literal native C events; development only.

The native endpoint rule and fair coin are unchanged.  The Poisson seam
scheduler, fixed compensating reference, electric energy and its inverse
temperature are explicit new controller choices.  No held-out case is run.
"""
from __future__ import annotations
import ast
import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np
import scipy.linalg as la
import sympy as sy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NATIVE = ROOT/'reverse-engineering-reality/model/observer/native.py'
ROUTING = ROOT/'reverse-engineering-reality/model/native_spatial_refinement/source.py'
N20 = ROOT/'openwave-push/issues/N20_native_green_field/n20lib.py'
N3D = ROOT/'openwave-push/issues/N3D_cubic_routing/n3dlib.py'
LINKS = (((5, 7), (6, 4)), ((9, 11), (10, 8)), ((1, 3), (2, 0)))
PINS = (NATIVE, ROUTING, N20, N3D)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def literal_native():
    tree = ast.parse(NATIVE.read_text())
    fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == 'native')
    ns = {}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(NATIVE), 'exec'), ns)
    return ns['native']


def internal_edges():
    # Independent exact Q(phi) distance construction, phi^2=phi+1.
    z, p, m, f, nf = (0,0), (1,0), (-1,0), (0,1), (0,-1)
    vertices = ((m,f,z),(p,f,z),(m,nf,z),(p,nf,z),
                (z,m,f),(z,p,f),(z,m,nf),(z,p,nf),
                (f,z,m),(f,z,p),(nf,z,m),(nf,z,p))
    out = []
    for a in range(12):
        for b in range(a+1,12):
            d0 = d1 = 0
            for u,v in zip(vertices[a],vertices[b]):
                x,y = u[0]-v[0],u[1]-v[1]
                d0 += x*x+y*y
                d1 += 2*x*y+y*y
            if (d0,d1)==(4,0): out.append((a,b))
    assert len(out)==30
    return tuple(out)


def graph(size=2):
    coords = tuple(product(range(size),repeat=3))
    ix = {x:i for i,x in enumerate(coords)}
    edges = set()
    for cell,x in enumerate(coords):
        edges.update((12*cell+a,12*cell+b) for a,b in internal_edges())
        for axis,row in enumerate(LINKS):
            y=list(x);y[axis]=(y[axis]+1)%size
            for a,b in row:edges.add(tuple(sorted((12*cell+a,12*ix[tuple(y)]+b))))
    edges = tuple(sorted(edges)); n=12*size**3
    neighbours=[[] for _ in range(n)]
    for a,b in edges:neighbours[a].append(b);neighbours[b].append(a)
    assert len(edges)==36*size**3 and all(len(row)==6 for row in neighbours)
    return edges,np.array(neighbours,dtype=int)


def exact_gap_certificate():
    # At size2 all eight Fourier characters are +-1.  Each 12x12
    # rational block minus 1/2 is positive, after lifting the sole zero.
    cert=[]
    for signs in product((-1,1),repeat=3):
        h=6*sy.eye(12)
        for a,b in internal_edges():h[a,b]-=1;h[b,a]-=1
        for sign,row in zip(signs,LINKS):
            for a,b in row:h[a,b]-=sign;h[b,a]-=sign
        if signs==(1,1,1):h += sy.ones(12)/12
        _,d=(h-sy.eye(12)/2).LDLdecomposition(hermitian=False)
        diagonal=[d[i,i] for i in range(12)]
        assert all(v>0 for v in diagonal)
        cert.append({'bloch_signs':signs,'positive_LDL_diagonal':[str(v) for v in diagonal]})
    return {'lambda2_lower_bound':'1/2','lambda_max_upper_bound':'12',
            'proof':'Eight exact Fourier blocks; rational positive LDL; degree bound for upper edge.',
            'blocks':cert}


def local_field(charges, neighbours, rounds):
    """Parallel nearest-neighbour arithmetic, starting from a zero field.

    Columns may carry independent configured probes.  No inverse or spectral
    projector appears here.  Every round reads only neighbours and local rho.
    """
    phi=np.zeros_like(charges,dtype=float)
    for _ in range(rounds):
        phi=(6*phi+phi[neighbours].sum(axis=1)+charges)/12
    return phi


def compile_controller(neighbours, reference=0, rounds=512):
    n=len(neighbours)
    charges=np.eye(n);charges[reference,:]-=1
    fields=local_field(charges,neighbours,rounds)
    energies=np.sum(charges*fields,axis=0)/2
    return {'energies':energies,'reference':reference,'rounds':rounds,
            'construction':'One local relaxation column per admitted one-load source state; cached finite circuit outputs.'}


def rates(counter, edges, energies, kappa=1., beta=4., nu=1.):
    # Only current-record state supplied by the observer is consumed.
    if any(v not in (0,1) for v in counter) or sum(counter)!=1:
        raise ValueError('outside declared one-load controller domain')
    x=int(np.argmax(counter))
    out=np.full(len(edges),nu,dtype=float)
    for e,(a,b) in enumerate(edges):
        if x==a or x==b:
            y=b if x==a else a
            out[e]=nu*np.exp(-beta*kappa*(energies[y]-energies[x])/2)
    return out


def execute(initial, edges, controller, kappa, seed, horizon=.5, cap=4096):
    """Complete competing-clock path, with all native holds and terminal wait.

    The source state goes only to literal_native.  The scheduler receives its
    own preparation/current-record counter.  Cap exhaustion is a separate
    retained output.  Log density is for ordered times and fine seam/coin marks.
    """
    n=len(controller['energies']); source=tuple(int(i==initial) for i in range(n))
    counter=list(source); flow=[0]*len(edges); native=literal_native()
    rng=np.random.default_rng(seed); t=0.; log_density=0.; records=[]
    for attempt in range(cap):
        rr=rates(counter,edges,controller['energies'],kappa)
        total=float(sum(rr)); wait=float(rng.exponential(1/total))
        if t+wait>horizon:
            log_density-=total*(horizon-t)
            return {'initial':initial,'kappa':kappa,'seed':seed,'horizon':horizon,
                    'records':records,'terminal':{'kind':'no_further_event','duration':horizon-t,
                    'total_seam_rate':total,'survival_probability':float(np.exp(-total*(horizon-t)))},
                    'log_path_density':log_density,'final_counter':counter,'final_flow':flow,
                    'source_matches_record_reconstruction':tuple(counter)==source}
        e=int(rng.choice(len(edges),p=rr/total));coin=int(rng.integers(2))
        source,(record_coin,j)=native(source,edges[e],coin)
        a,b=edges[e];counter[a]-=j;counter[b]+=j;flow[e]+=j
        log_density+=-total*wait+np.log(rr[e]/2);t+=wait
        assert tuple(counter)==source and record_coin==coin
        records.append({'attempt':attempt,'time':t,'wait':wait,'seam':e,
                        'coin':coin,'current':j,'seam_rate':float(rr[e]),
                        'total_seam_rate':total})
    return {'initial':initial,'kappa':kappa,'seed':seed,'horizon':horizon,
            'records':records,'terminal':{'kind':'cap_exhausted','time':t},
            'log_path_density':log_density,'final_counter':counter,'final_flow':flow,
            'source_matches_record_reconstruction':tuple(counter)==source}


def independent_replay(path, edges, controller):
    n=len(controller['energies']);x=path['initial']; t=0.;logp=0.;flow=[0]*len(edges)
    for r in path['records']:
        c=[int(i==x) for i in range(n)]
        rr=rates(c,edges,controller['energies'],path['kappa'])
        a,b=edges[r['seam']];old=x
        # Independent floor/ceiling endpoint rule, not the producer formula.
        total=int(x==a)+int(x==b)
        if total:
            new_a=(total+r['coin'])//2
            x=a if new_a else b
        j=int(old==a)-int(x==a)
        assert j==r['current'] and r['seam_rate']==rr[r['seam']]
        t+=r['wait'];assert abs(t-r['time'])<1e-13
        logp-=sum(rr)*r['wait'];logp+=np.log(rr[r['seam']]/2);flow[r['seam']]+=j
    if path['terminal']['kind']=='no_further_event':
        rr=rates([int(i==x) for i in range(n)],edges,controller['energies'],path['kappa'])
        logp-=sum(rr)*(path['horizon']-t)
    assert abs(logp-path['log_path_density'])<1e-10
    assert path['final_counter']==[int(i==x) for i in range(n)] and flow==path['final_flow']
    divergence=np.zeros(n,dtype=int)
    for f,(a,b) in zip(flow,edges):divergence[a]+=f;divergence[b]-=f
    initial=np.zeros(n,dtype=int);initial[path['initial']]=1
    assert np.array_equal(initial-divergence,np.array(path['final_counter']))
    return sum(r['current']==0 for r in path['records'])


def generator(edges, energies, kappa, beta=4.,nu=1.):
    n=len(energies);q=np.zeros((n,n))
    for x in range(n):
        rr=rates([int(i==x) for i in range(n)],edges,energies,kappa,beta,nu)
        for e,(a,b) in enumerate(edges):
            if x in (a,b):
                y=b if x==a else a;q[y,x]+=rr[e]/2;q[x,x]-=rr[e]/2
    return q


def run():
    before={str(p.relative_to(ROOT)):sha(p) for p in PINS}
    edges,nb=graph();n=len(nb);lap=6*np.eye(n)
    for a,b in edges:lap[a,b]-=1;lap[b,a]-=1
    certificate=exact_gap_certificate();controller=compile_controller(nb)
    # Dense inverse is an independent audit oracle; not used by the scheduler.
    green=la.inv(lap+np.ones((n,n))/n)-np.ones((n,n))/n
    true_u=(np.diag(green)+green[0,0]-2*green[:,0])/2
    eps=2*(23/24)**controller['rounds']
    actual=float(max(abs(true_u-controller['energies'])))
    assert actual<eps+2e-13
    rows=[]
    for k in (0.,1.,2.):
        q=generator(edges,controller['energies'],k)
        pi=np.exp(-4*k*controller['energies']);pi/=sum(pi)
        db=q*pi[None,:]
        assert max(abs(db-db.T).flat)<1e-15 and max(abs(q@pi))<1e-15
        p=np.ones(n)/n;out=la.expm(.5*q)@p
        drift=float(controller['energies']@(q@p))
        formula=-sum((controller['energies'][b]-controller['energies'][a])*
                     np.sinh(2*k*(controller['energies'][b]-controller['energies'][a]))
                     for a,b in edges)/n
        assert abs(drift-formula)<1e-14 and abs(sum(out)-1)<1e-13 and min(out)>=0
        rows.append({'kappa':k,'mean_base_electric_energy_at0':float(controller['energies']@p),
                     'mean_base_electric_energy_at_half_time':float(controller['energies']@out),
                     'initial_energy_current_rate':drift,'exact_edge_formula_evaluated':float(formula),
                     'detailed_balance_residual':float(max(abs(db-db.T).flat)),
                     'full_position_law_at_half_time':out.tolist()})
    assert abs(rows[0]['initial_energy_current_rate'])<1e-14
    assert rows[2]['initial_energy_current_rate']<rows[1]['initial_energy_current_rate']<0
    paths=[];holds=events=0
    for x in (0,1,5,12,35,95):
        for k in (0.,1.,2.):
            for seed in (173,804):
                path=execute(x,edges,controller,k,seed)
                holds+=independent_replay(path,edges,controller);events+=len(path['records']);paths.append(path)
    # Explicit retained cap failure; no conclusive-only normalization.
    cap_path=execute(35,edges,controller,1.,57,horizon=100.,cap=1)
    assert cap_path['terminal']['kind']=='cap_exhausted';independent_replay(cap_path,edges,controller)
    noevent=execute(35,edges,controller,1.,57,horizon=0.,cap=1)
    assert noevent['records']==[] and noevent['terminal']['survival_probability']==1.
    # All source columns, seam choices and coins agree independently.
    native=literal_native();branches=0
    for x in range(n):
        for a,b in edges:
            for coin in (0,1):
                state=tuple(int(i==x) for i in range(n));out,(_,j)=native(state,(a,b),coin)
                target=x
                if x in (a,b):target=a if coin else b
                assert out==tuple(int(i==target) for i in range(n))
                assert j==int(x==a)-int(target==a);branches+=1
    after={str(p.relative_to(ROOT)):sha(p) for p in PINS};assert before==after
    # Rate/path perturbation bound for the exact finite arithmetic circuit.
    # Only six incident seam rates change between approximate/exact fields.
    path_bound=.5*6*np.exp(4)*np.expm1(4*eps)
    result={'schema':'oph.native_field_feedback.development.v1',
      'checks':{'independent_actual_native_all_branches':True,'local_relaxation_no_inverse_in_controller':True,
                'exact_rational_gap_certificate':True,'finite_field_error_bound':True,
                'all_current_only_scheduled_paths_replayed':True,'complete_holds_waits_and_cap_failure':True,
                'energy_bias_and_detailed_balance':True,'source_bytes_unchanged':True},
      'domain':{'ports':n,'seams':len(edges),'source':'one native load; known initial preparation',
                'reference':'fixed compensating -1 charge at port0, supplied apparatus',
                'native_branches_checked':branches,'executed_paths':len(paths),'native_events':events,
                'retained_zero_currents':holds,'official_rows_changed':[]},
      'controller':{'kappa_cases':[0,1,2],'beta':4,'base_seam_rate':1,
                    'local_relaxation_rounds':controller['rounds'],'round_clock':'compiler work before the native Poisson experiment',
                    'online_input':'initial one-load preparation and subsequent native C current records only',
                    'energy':'U_K=kappa rho^T G_K rho/2; G_K=(1/12)sum_{j<K}(I-L/12)^j on neutral vectors'},
      'field_error':{'exact_operator_norm_bound':eps,'observed_max_energy_error':actual,
                     'total_variation_complete_path_bound_time_half_kappa1_exact_arithmetic':min(1.,path_bound),
                     'float_status':'Numerical execution checked to stated residuals; exact-arithmetic truncation bound does not include floating rounding.'},
      'gap_certificate':certificate,'interventions':rows,
      'clock':{'unchanged_kappa0_attempt_rate':len(edges),'unchanged_kappa0_load_transfer_rate':3,
                'fine_path_density':'exp(-sum Lambda_i Delta_t_i) product_i(r_e_i/2), with terminal survival',
                'uniform_attempt_clock_diffusion_coefficient':'1/(16M)',
                'Poisson_per_seam_nu_diffusion_coefficient':'nu/16',
                'coherent_Hm_speed_if_separately_installed':'1/sqrt(8); not this dissipative clock'},
      'source_pins':before,'complete_paths':paths,'cap_failure_example':cap_path,
      'zero_duration_example':noevent,'auditor_sha256':sha(Path(__file__))}
    (HERE/'FIELD_FEEDBACK.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'domain':result['domain'],'field_error':result['field_error'],
                      'energy_response':[{'kappa':r['kappa'],'initial_energy_current_rate':r['initial_energy_current_rate']} for r in rows]}))


if __name__=='__main__':run()
