"""Independent exact source audit and a declared current-square extension.

No empirical hadron value enters a Hamiltonian. Comparisons are evaluated
after operators and certified intervals have been built. All inputs are frozen.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('frozen_colour_source', HERE/'inputs/0_source.py')
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def inertia(matrix):
    """Exact symmetric congruence; count strictly negative eigenvalues."""
    A = [row[:] for row in matrix]
    count = 0
    while A:
        n = len(A)
        j = next((i for i in range(n) if A[i][i]), None)
        if j is None:
            pair = next(((i,k) for i in range(n) for k in range(i) if A[i][k]), None)
            if pair is None:
                break
            i,k = pair
            old = [row[:] for row in A]
            for r in range(n):
                A[i][r] = old[i][r]+old[k][r]
                A[r][i] = old[r][i]+old[r][k]
            A[i][i] = old[i][i]+old[k][k]+2*old[i][k]
            j = i
        A[0],A[j] = A[j],A[0]
        for row in A:
            row[0],row[j] = row[j],row[0]
        pivot = A[0][0]
        count += pivot < 0
        A = [[A[i][k]-A[i][0]*A[0][k]/pivot for k in range(1,n)] for i in range(1,n)]
    return count


def extension_values(d, kind):
    """One inherited U=1 coefficient, no fit or search.

    b=N_colour/3 and q=actual colourful electric current. The amended
    interaction is sum_x(b_x-q_x)^2; flipped swaps up and down as a control.
    EM is the separate positive grounded-chain electric field at a=1/16.
    """
    vals=[]
    L=d['L']
    for occ in d['counts']:
        ns=[sum(occ[4*x:4*x+4]) for x in range(L)]
        qs=[sum((F(2,3) if z%2==0 else -F(1,3))*occ[4*x+z] for z in range(4)) for x in range(L)]
        if kind=='none':
            value=F()
        elif kind=='current_square':
            value=sum((F(n,3)-q)**2 for n,q in zip(ns,qs))
        elif kind=='flavour_reflected':
            value=sum(q*q for q in qs)
        elif kind=='electric':
            value=sum(qs[i]*qs[j]*F((min(i,j)+1)*(L-max(i,j)),L+1) for i in range(L) for j in range(L))/32
        else:
            raise ValueError(kind)
        assert value>=0
        vals.append(value)
    return vals+[F()]*(d['code_dimension']-d['dimension'])


def rational_matrix(d, vals):
    # Congruence by sqrt(prod n_i!) removes every bosonic normalization radical.
    gs=[]
    for occ in d['counts']:
        g=1
        for n in occ:
            g*=factorial(n)
        gs.append(g)
    n=d['dimension']
    A=[[F() for _ in range(n)] for _ in range(n)]
    for (i,j),z in d['entries'].items():
        z=source.mul(z,source.root(gs[i]*gs[j]))
        assert set(z)<={1}
        A[i][j]=z.get(1,F())
    for i in range(n):
        A[i][i]+=vals[i]*gs[i]
    assert all(A[i][j]==A[j][i] for i in range(n) for j in range(n))
    return A,gs


def ground_interval(A,G):
    def nbelow(x):
        B=[row[:] for row in A]
        for i,g in enumerate(G):
            B[i][i]-=x*g
        return inertia(B)
    lo,hi=F(0),F(16)
    assert nbelow(lo)==0 and nbelow(hi)>0
    while hi-lo>F(1,1<<36):
        mid=(lo+hi)/2
        if nbelow(mid)==0:
            lo=mid
        else:
            hi=mid
    assert nbelow(lo)==0 and nbelow(hi)>0
    return lo,hi


def run():
    pins=json.loads((HERE/'INPUTS.json').read_text())
    for pin in pins:
        assert hashlib.sha256((HERE/pin['snapshot']).read_bytes()).hexdigest()==pin['sha256']
    trial_checks=[]
    for L in [2,3,4,5]:
        for kappa in [F(1,2),F(1),F(2)]:
            ds=[source.data(L,kappa,w) for w in [1,2]]
            target={lab:i for i,lab in enumerate(ds[1]['labels'])}
            perm=[target[tuple(sorted(z^1 for z in lab))] for lab in ds[0]['labels']]
            for (i,j),v in ds[0]['entries'].items():
                assert v==ds[1]['entries'][perm[i],perm[j]]
            for d in ds:
                # The exact endpoint local spin-half eigenvalue is 57/16.
                for x in [0,L-1]:
                    ids=[i for i,p in enumerate(d['patterns']) if p[x]==3]
                    A,G=rational_matrix(d,extension_values(d,'none'))
                    B=[[A[i][j]-(F(57,16)*G[i] if i==j else 0) for j in ids] for i in ids]
                    assert len(ids)==2 and B[0][0]*B[1][1]==B[0][1]*B[1][0]
                    assert B[0][0]>=0 and B[1][1]>=0
                    assert all(d['span'][i]==0 for i in ids)
                    trial_checks.append([L,str(kappa),d['weakup'],x])
    ratios=json.loads((HERE/'inputs/3_results.json').read_text())['exact_free_mass_ratio_intervals']
    a=list(map(F,ratios['up']));b=list(map(F,ratios['down']))
    assert 0<a[0]<=a[1]<b[0]<=b[1]<1
    gaps={}
    eigen=[]
    for kind in ['none','electric','current_square','flavour_reflected']:
        ints=[]
        for weakup in [1,2]:
            d=source.data(2,F(1),weakup)
            vals=extension_values(d,kind)
            base=ground_interval(*rational_matrix(d,vals))
            shiftlo=weakup*(a[0]-1)+(3-weakup)*(b[0]-1)
            shifthi=weakup*(a[1]-1)+(3-weakup)*(b[1]-1)
            E=[base[0]+shiftlo,base[1]+shifthi]
            ints.append(E)
            eigen.append(dict(kind=kind,weakup=weakup,complete_sector_dimension=d['dimension'],energy_exact=list(map(str,E))))
        gaps[kind]=[ints[0][0]-ints[1][1],ints[0][1]-ints[1][0]]
    assert 0<gaps['none'][0]<gaps['none'][1]<1
    assert gaps['electric'][1]<gaps['none'][0]
    assert gaps['current_square'][0]>1
    assert gaps['flavour_reflected'][1]<0
    # Full one-cell colourful occupation check. u,d and their conjugates each
    # have 3 colours times 2 spins. B is signed quark number, so charge
    # conjugation changes both B and Q, preserving the squared current.
    local_controls=0
    for u in range(7):
        for d in range(7):
            for au in range(7):
                for ad in range(7):
                    B=F(u+d-au-ad,3)
                    Q=F(2*u-d-2*au+ad,3)
                    P=(B-Q)**2
                    Bc=F(au+ad-u-d,3)
                    Qc=F(2*au-ad-2*u+d,3)
                    assert P>=0 and P==(Bc-Qc)**2
                    local_controls+=1
    unit=F('0.00051099895069')
    proton=F('0.93827208943')
    neutron=F('0.93956542194')
    nuclear=json.loads((HERE/'inputs/4_results.json').read_text())
    nuclear_checks=[]
    for row in nuclear['complete_ground_bounds']:
        if row['name'] not in ['helion','alpha']:
            continue
        lo,hi=map(F,row['full_ground_GeV_interval']);tar=F(row['target_GeV'])
        assert not lo<=tar<=hi
        nuclear_checks.append(dict(name=row['name'],interval_GeV=[str(lo),str(hi)],target_GeV=str(tar),
            distance_to_interval_GeV=str(max(lo-tar,tar-hi)),maximum_field_lowering_GeV=row['maximum_field_induced_ground_lowering_GeV_exact']))
    out=dict(schema='oph.obstruction_repair.hadrons.v1',input_pins_verified=True,
        exact_endpoint_trials=len(trial_checks),trial_configurations=trial_checks,
        colour_singlet_trial_bound='57/16',trial_mass_upper_GeV=str(F(57,16)*unit),
        proton_over_trial_bound=float(proton/(F(57,16)*unit)),
        existing_colour_bridge_already_has_SU3=True,common_EM_transfer_proved=False,
        energy_intervals=eigen,gap_intervals_exact={k:list(map(str,v)) for k,v in gaps.items()},
        gap_intervals_decimal={k:list(map(float,v)) for k,v in gaps.items()},
        amended_beta_available_energy_GeV=list(map(str,[(v-1)*unit for v in gaps['current_square']])),
        physical_np_gap_comparison_GeV=str(neutron-proton),
        amended_np_gap_GeV=list(map(str,[v*unit for v in gaps['current_square']])),
        negative_controls=['switch off current square: threshold closed','reflect flavour current: order reverses','add only positive electric field: gap decreases'],
        nuclear_exclusions=nuclear_checks,
        current_square_fixed_coefficient='inherited U=1; no numerical tuning',
        current_square_ancestry='New architecture chosen to test a known obstruction; no independent empirical prediction or forced OPH derivation claimed.',
        current_square_full_Fock_positive=True,current_square_electron_sector_unchanged=True,
        local_charge_conjugation_and_positivity_checks=local_controls,
        structural_beta_threshold_repaired_on_declared_amendment=True,
        physical_mass_rows_solved=False,physical_beta_decay_executed=False,halo_lifetime_solved=False,
        retained_obligations=['physical mass scale','correct neutron-proton gap','common field/source clock map','weak outgoing decay on same amended source','native full preparation and observation program','nuclear continuum and halo channels'])
    (HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['exact_endpoint_trials','proton_over_trial_bound','gap_intervals_decimal','structural_beta_threshold_repaired_on_declared_amendment','physical_mass_rows_solved']},indent=2))


if __name__=='__main__':
    run()
