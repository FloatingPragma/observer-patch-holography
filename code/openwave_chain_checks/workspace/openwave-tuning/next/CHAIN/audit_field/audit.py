"""Independent arithmetic, symbol, integer-fiber and complete-mark audit.

Producer execution and verification routines are not called. Source classes
and saved artifacts are comparison targets for independently rebuilt objects.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product, combinations
import hashlib
import importlib.util
import json
import math
import sys

import numpy as np
import scipy.linalg as la
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
import sympy as sy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CHAIN = HERE.parent
BRIDGE = ROOT / "openwave-tuning/auditor/native_field_bridge"
PINS = {
    "field_limit/MANIFEST.json": "612aedd40f13faea2b80b5c31d4afbb19cb2332ff22db5b459ac2ae936879008",
    "dynamical_field/MANIFEST.json": "2aa4fe69b85d389949d4bbaed733b5cd69aafab4760cb4a2a8cf68019c5fe2cc",
    "force/MANIFEST.json": "f4f2f6b490aef1a87eb0fd6b5c07c429d5f394baf4414e0f39b7dc309babba1e",
}
LINKS = (((5, 7), (6, 4)), ((9, 11), (10, 8)), ((1, 3), (2, 0)))
HALF = ((0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (0, 5))
SEAMS = ((0, 1), (1, 2), (2, 0))
T, RATE, BETA = F(37, 100), F(4, 5), F(361, 50)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    sys.modules[name] = value
    spec.loader.exec_module(value)
    return value


def verify_pins():
    result = {}
    for relative, expected in PINS.items():
        path = CHAIN / relative
        assert sha(path) == expected, relative
        manifest = read(path)
        for entry, digest in manifest.get("files", manifest.get("sha256", {})).items():
            target = path.parent / entry if "files" in manifest else ROOT / entry
            assert sha(target) == digest, str(target)
            result[str(target.relative_to(ROOT))] = digest
        result[str(path.relative_to(ROOT))] = expected
    for relative, expected in read(CHAIN / "field_limit/RESULT.json")["source_pins"].items():
        assert sha(ROOT / relative) == expected, relative
        result[relative] = expected
    return result


def frame_and_graph():
    # Coordinates in Q(sqrt(5)); comparison of exact squared distances.
    p = (1 + sy.sqrt(5)) / 2
    vertices = [(-1,p,0),(1,p,0),(-1,-p,0),(1,-p,0),
                (0,-1,p),(0,1,p),(0,-1,-p),(0,1,-p),
                (p,0,-1),(p,0,1),(-p,0,-1),(-p,0,1)]
    internal = [(a,b) for a,b in combinations(range(12),2)
                if sy.simplify(sum((vertices[a][i]-vertices[b][i])**2 for i in range(3))) == 4]
    assert len(internal) == 30
    cells = list(product(range(2), repeat=3)); edges = set()
    for ci, cell in enumerate(cells):
        edges.update((12*ci+a,12*ci+b) for a,b in internal)
        for axis, links in enumerate(LINKS):
            other = list(cell); other[axis] ^= 1
            cj = cells.index(tuple(other))
            edges.update(tuple(sorted((12*ci+a,12*cj+b))) for a,b in links)
    edges = sorted(edges)
    incidence = np.zeros((96,len(edges)),int)
    for col,(a,b) in enumerate(edges): incidence[a,col]=1; incidence[b,col]=-1
    lap = incidence @ incidence.T
    neighbours = [list(np.flatnonzero(lap[x] == -1)) for x in range(96)]
    assert len(edges)==288 and all(len(row)==6 for row in neighbours)
    return internal, edges, incidence, lap, neighbours


def symbol(internal, momentum):
    matrix = 6*np.eye(12,dtype=complex)
    for a,b in internal: matrix[a,b]-=1; matrix[b,a]-=1
    for angle,links in zip(momentum, LINKS):
        for a,b in links:
            matrix[a,b]-=np.exp(1j*angle); matrix[b,a]-=np.exp(-1j*angle)
    return matrix


def gap_certificates(internal):
    pivots = []
    for signs in product((-1,1),repeat=3):
        h=6*sy.eye(12)
        for a,b in internal: h[a,b]-=1; h[b,a]-=1
        for sign,links in zip(signs,LINKS):
            for a,b in links: h[a,b]-=sign; h[b,a]-=sign
        if signs==(1,1,1): h+=sy.ones(12)/12
        h-=sy.eye(12)/2
        # Exact symmetric Schur elimination, independent of producer LDL.
        local=[]
        while h.rows:
            pivot=h[0,0]; assert pivot>0
            local.append(str(pivot))
            h=h[1:,1:]-h[1:,0]*h[0,1:]/pivot
        pivots.append({"signs":list(signs),"positive_pivots":local})
    zero=sy.Matrix(np.rint(symbol(internal,(0,0,0)).real).astype(int))
    assert zero.eigenvals()=={sy.Integer(0):1,sy.Integer(4):3,sy.Integer(6):2,sy.Integer(8):6}
    return pivots


def fixed_solve(z, rho, nb, rounds, bits=40):
    scale=2**bits
    for _ in range(rounds):
        nxt=[]
        for x,row in enumerate(nb):
            numerator=6*z[x]+sum(z[y] for y in row)+scale*int(rho[x])
            quotient,remainder=divmod(numerator,12)
            nxt.append(quotient+int(remainder>=6))
        z=nxt
    return z


def upward(value, places):
    scale=10**places
    return F((value.numerator*scale+value.denominator-1)//value.denominator,scale)


def arithmetic_audit(edges,incidence,lap,nb):
    prod=module("audited_static_field",CHAIN/"field_limit/controlled_field.py")
    for numerator in range(-1200,1201):
        q,r=divmod(numerator,12); rounded=q+int(r>=6)
        assert rounded==prod.round_ratio(numerator)
        assert abs(F(rounded)-F(numerator,12))<=F(1,2)
    assert F(12*96,4)<17**2
    theta=F(23,24); alpha=theta**32; floor=F(408,2**40)
    green=la.pinvh(lap.astype(float))
    source=read(BRIDGE/"RECIPROCAL_FEEDBACK.json")
    paths=[p for p in source["complete_paths"] if p["initial"]==[1,35] and p["kappa"] in (0.,1.) and p["seed"]==804]+[source["cap_failure"]]
    stored=read(CHAIN/"field_limit/RESULT.json")["local_tracking"]
    rows=[]; signs_detected=holds=0
    for path,recorded in zip(paths,stored):
        positions=list(path["initial"]); rho=np.zeros(96,int)
        rho[positions[0]]+=1;rho[positions[1]]-=1
        z=[0]*96; bound=F(2); flow=np.zeros(288,int)
        initial=rho.copy(); snapshots=[]; maxratio=0.; max_residual=0.
        for index,record in enumerate([None]+path["records"]):
            moved=0
            if record is not None:
                side,e,coin=record["copy"],record["seam"],record["coin"]
                a,b=edges[e]; old=positions[side]
                dest=(a if coin else b) if old in (a,b) else old
                current=int(old==a)-int(dest==a)
                assert current==record["current"]==record["observer_current"]
                positions[side]=dest; sign=1-2*side
                flow[e]+=sign*current; rho=initial-incidence@flow
                direct=np.zeros(96,int); direct[positions[0]]+=1;direct[positions[1]]-=1
                assert np.array_equal(rho,direct)
                holds+=current==0;signs_detected+=current!=0
                moved=abs(current)
            z=fixed_solve(z,rho,nb,32)
            bound=upward(alpha*(bound+moved)+floor*(1-alpha),18)
            residual=[2**40*int(rho[x])-6*z[x]+sum(z[y] for y in nb[x]) for x in range(96)]
            assert sum(residual)==0
            certificate=F(2*sum(v*v for v in residual),2**80)
            assert certificate==prod.residual_squared_certificate(z,rho,nb)
            error=np.array(z)/2**40-green@rho
            actual=float(np.linalg.norm(incidence.T@error))
            assert actual<=float(bound)+1e-12 and actual**2<=float(certificate)+1e-23
            assert abs(float(rho@error)/2)<=actual+1e-12
            maxratio=max(maxratio,actual/float(bound));max_residual=max(max_residual,float(certificate))
            if moved or index in (0,len(path["records"])):
                snapshots.append({"after_record_count":index,"tracking_bound_exact_upper":str(bound),
                                  "residual_squared_bound_exact_upper":str(upward(certificate,24))})
        assert [np.eye(96,dtype=int)[p].tolist() for p in positions]==path["final_counters"]
        assert rho.tolist()==recorded["final_charge"]
        for expected,actual in zip(recorded["snapshot_certificates"],snapshots):
            assert all(expected[key]==value for key,value in actual.items())
        assert len(snapshots)==len(recorded["snapshot_certificates"])
        rows.append({"records":len(path["records"]),"cap_kind":path["terminal"]["kind"],
                     "maximum_error_to_bound":maxratio,"snapshot_count":len(snapshots)})
    assert sum(row["records"] for row in rows)==539 and holds>0 and signs_detected>0
    # Uniform branch generator, in source time, from all incident fair coins.
    generator=np.zeros((96,96))
    for x in range(96):
        for a,b in edges:
            for coin in (0,1):
                y=(a if coin else b) if x in (a,b) else x
                generator[y,x]+=.5;generator[x,x]-=.5
    assert np.array_equal(generator,-lap/2)
    caps=[]
    for row in read(BRIDGE/"CAPPED_LAW.json")["count_controls"]:
        if row["copies"]==2 and row["kappa"]==0:
            mean=576*row["horizon"]
            complete=math.exp(-mean)*sum(mean**j/math.factorial(j) for j in range(row["cap"]))
            assert abs(complete-row["success_probability"])<1e-13
            assert 1-complete>.01
            caps.append({"cap":row["cap"],"discarded_probability_if_flag_deleted":1-complete})
    return {"paths":rows,"records":539,"holds":holds,"moving_records":signs_detected,
            "rounding_residues_tested":2401,"cap_mutations":caps}, stored[0]["final_charge"], prod


def pi_enclosure():
    # Independent tighter alternating Machin interval.
    intervals=[]
    for inverse,count in ((5,40),(239,12)):
        terms=[F((-1)**j,(2*j+1)*inverse**(2*j+1)) for j in range(count)]
        partial=sum(terms,F());next_term=F((-1)**count,(2*count+1)*inverse**(2*count+1))
        intervals.append((min(partial,partial+next_term),max(partial,partial+next_term)))
    return 16*intervals[0][0]-4*intervals[1][1],16*intervals[0][1]-4*intervals[1][0]


def static_audit(internal,charge,prod):
    lo,hi=pi_enclosure();unit=np.ones(12)/np.sqrt(12);qproj=np.eye(12)-np.outer(unit,unit)
    modes=[m for m in product((-1,0,1),repeat=3) if any(m)]
    saved=read(CHAIN/"field_limit/RESULT.json")["static_refinement"]
    rows=[]
    for size,recorded in zip((16,32,64),saved):
        sums=[F(),F(),F()];maxratio=0.;total_diagnostic=0.
        for mode in modes:
            r=2*hi/size;q=(2*lo/size)**2*sum(v*v for v in mode)/8
            rem=F(3,8)*r**4; beta=4-r; defect=r*r/(4-r-r*r/2)**2
            resolvent=rem/(q*(q-rem))+defect/(q-rem)+defect/beta
            cross=r/(beta*(q-rem))
            bound=(resolvent/6+F(3,5)*cross)/size**2
            averaged=resolvent/(6*size**2)
            pbound,pavg,pres=prod.mode_bounds(size,mode)
            assert bound<=pbound and averaged<=pavg and resolvent<=pres
            weight=F(1,2**sum(v!=0 for v in mode));sums[0]+=weight*pbound;sums[1]+=weight*pavg
            sums[2]+=weight*pbound*2*prod.pi_interval()[1]*sum(abs(v) for v in mode)
            xi=2*np.pi*np.array(mode);k=xi/size;matrix=symbol(internal,k)
            values,vectors=la.eigh(matrix);inv=la.inv(matrix)
            actualq=float(k@k)/8;actualr=float(max(abs(k)))
            assert abs(values[0]-actualq)<=3*actualr**4/8+1e-13
            assert min(values[1:])>=4-actualr-1e-13
            optical_weight=sum(abs(unit@vectors[:,j])**2 for j in range(1,12))
            assert optical_weight<=float(defect)+1e-13
            assert abs(unit@inv@unit-1/actualq)<=float(resolvent)+1e-11
            assert la.norm(unit@inv@qproj)<=float(cross)+1e-11
            rho=np.zeros(12,complex)
            for port,value in enumerate(charge):
                if value:
                    cell,internal_port=divmod(port,12)
                    rho[internal_port]+=value*np.exp(-1j*xi@np.array(np.unravel_index(cell,(2,2,2)))/2)
            finite=unit@inv@rho/(np.sqrt(12)*size**2);target=2*rho.sum()/(3*float(xi@xi))
            error=abs(finite-target)
            assert error<=float(bound)+1e-13
            maxratio=max(maxratio,error/float(bound));total_diagnostic+=float(weight)*error
        keys=("all_position_supnorm_bound_exact_upper","uniform_internal_supnorm_bound_exact_upper","all_position_electric_bound_exact_upper")
        assert all(F(recorded[key])==upward(value,15) for key,value in zip(keys,sums))
        k=2*np.pi*np.array((1,0,1))/size;full=symbol(internal,k)
        reduced=symbol(internal,(k[0],0,0));target=1/np.pi**2
        _,_,allowance=prod.mode_bounds(size,(1,0,1))
        mutations={"remove_z_routing":abs(unit@la.solve(reduced,unit)/size**2-target)>float(allowance)/size**2,
                   "replace_acoustic_with_optical":abs(1/(size**2*la.eigvalsh(full)[1])-target)>float(allowance)/size**2}
        assert all(mutations.values())
        rows.append({"N":size,"all_twelve_modes_checked_at":len(modes),"localized_potential_bound":str(sums[0]),
                     "averaged_potential_bound":str(sums[1]),"electric_vector_bound":str(sums[2]),
                     "maximum_mode_error_to_bound":maxratio,"mutations":{k:bool(v) for k,v in mutations.items()}})
    return rows


def car_move(mask,src,dst):
    if not(mask>>src&1) or mask>>dst&1:return None
    between=sum(mask>>j&1 for j in range(min(src,dst)+1,max(src,dst)))
    return mask^(1<<src)^(1<<dst),(-1)**between


def integer_charts():
    incidence=sy.zeros(6,6)
    for h,(a,b) in enumerate(HALF):incidence[a,h]=1;incidence[b,h]=-1
    # Solve five independent rows directly, with F3 as free coordinate.
    tree=[0,1,2,4,5]; inverse=incidence[:5,tree].inv()
    offsets={}
    for mask in range(64):
        rhs=sy.Matrix([int(z<3)-(mask>>z&1)-int(z==0)*(3-mask.bit_count()) for z in range(6)])
        f=sy.zeros(6,1);values=inverse*rhs[:5,0]
        for h,value in zip(tree,values):f[h]=value
        assert incidence*f==rhs
        offsets[mask]=np.array(f,dtype=int).ravel()
    null=sy.zeros(6,1);null[3]=1
    for h,value in zip(tree,-inverse*incidence[:5,3]):null[h]=value
    assert incidence*null==sy.zeros(6,1)
    return offsets,np.array(null,dtype=int).ravel(),np.array(incidence,dtype=int)


def independent_family(modulus,offsets,cycle):
    ks=range(-modulus//2,(modulus-1)//2+1) if modulus%2==0 else range(-(modulus//2),(modulus-1)//2+1)
    states=[(mask,k) for mask in range(64) for k in ks];index={state:i for i,state in enumerate(states)}
    center=lambda k:(k+modulus//2)%modulus-modulus//2
    n=len(states);records=np.array([offsets[mask]+k*cycle for mask,k in states]);r=center(records)
    electric=.5*.37*(r*r).sum(axis=1);mass=np.array([.23*sum((1 if z<3 else -1)*(mask>>z&1) for z in range(6)) for mask,k in states])
    rows=[];cols=[];vals=[];hop_edges=[];groups={};stepmax=0
    for col,(mask,k) in enumerate(states):
        f=records[col]
        for h,(a,b) in enumerate(HALF):
            for src,dst,delta in ((a,b,1),(b,a,-1)):
                moved=car_move(mask,src,dst)
                if moved is None:continue
                target_f=f.copy();target_f[h]+=delta
                target=(moved[0],center(int(target_f[3])))
                rows.append(index[target]);cols.append(col);vals.append((-1)**h*moved[1]);hop_edges.append((index[target],col))
                exact=offsets[moved[0]]+int(target_f[3])*cycle
                assert np.array_equal(exact,target_f)
        for sign in (-1,1):
            rows.append(index[mask,center(k+sign)]);cols.append(col);vals.append(-.61)
        for e,(a,b) in enumerate(SEAMS):
            na,nb=mask>>a&1,mask>>b&1
            for coin in (0,1):
                outa=(na+nb+coin)//2; current=na-outa
                moved=(mask,1) if current==0 else car_move(mask,a if current>0 else b,b if current>0 else a)
                tf=f.copy();tf[2*e]+=current;tf[2*e+1]-=current
                exact=offsets[moved[0]]+int(tf[3])*cycle
                assert np.array_equal(exact,tf)
                stepmax=max(stepmax,int(abs(tf-f).max()))
                target=index[moved[0],center(int(tf[3]))]
                for law,value in (("C",current),("E",na-nb)):
                    key=e,law,coin,value
                    groups.setdefault(key,[]).append((target,col,moved[1]))
    h=sparse.csr_matrix((vals,(rows,cols)),shape=(n,n))+sparse.diags(mass+electric+2*.61)
    operators={}
    for key,entries in sorted(groups.items()):
        rr,cc,vv=zip(*entries)
        # Two law choices, two coins and three seams.
        operators[key]=sparse.csr_matrix((np.array(vv)*(.5/np.sqrt(3)),(rr,cc)),shape=(n,n))
    assert stepmax==1 and len(operators)==30
    assert max(np.bincount([a for a,b in hop_edges],minlength=n))<=6
    assert max(np.bincount([b for a,b in hop_edges],minlength=n))<=6
    return states,index,records,electric,h,operators


def rational_tail(x,d):
    if d<=0 or x>=d+1:return F(1)
    return min(F(1),2*x**d/F(math.factorial(d))/(1-x/F(d+1)))


def dynamical_audit():
    sys.path.insert(0,str(ROOT/"openwave-tuning/next/C10"))
    from joint_source import JointNativeSource
    actual=JointNativeSource()
    prod=module("audited_dynamic_field",CHAIN/"dynamical_field/dynamical_field.py")
    offsets,cycle,incidence=integer_charts()
    assert np.array_equal(cycle,-np.array([1,-1,1,-1,1,-1]))
    source_labels=[(mask,int((f[3]+3)%6-3)) for mask,f in actual.states]
    six=independent_family(6,offsets,cycle);states,index,records,electric,h,ops=six
    order=[index[state] for state in source_labels]
    h_actual=h.toarray()[np.ix_(order,order)]
    herr=float(abs(h_actual-actual.hamiltonian).max());assert herr<2e-15
    # Rebuild every actual M6 Hamiltonian entry, then use source-order
    # diagonalization independently of producer's number-block reference.
    evals,evecs=la.eigh(h_actual);psi=evecs[:,0];rho=np.outer(psi,psi)
    assert evals[1]-evals[0]>.001
    source_rho,energy,multiplicity=actual.reference()
    assert multiplicity==1 and abs(source_rho-rho).max()<2e-13
    arr=np.load(CHAIN/"dynamical_field/dynamical_field_states.npz")
    clock=np.load(CHAIN/"clock/common_clock_states.npz")
    source_columns=arr["actual_preparation_columns"]
    reservoir=[];branches=[];max_preparation_error=0.
    # Rebuild fair reservoir swap from seven-orbital CAR actions.
    for initial in (0,1):
        for final in (0,1):
            op=np.zeros((384,384))
            for col,(mask,f) in enumerate(actual.states):
                full=mask+(initial<<6);n0=mask&1
                if n0==initial:
                    target=full;phase=1-n0-initial
                else:
                    target,phase=car_move(full,0 if n0 else 6,6 if n0 else 0)
                if target>>6==final:
                    op[actual.source.index[(target&63,f)],col]=phase/math.sqrt(2)
            prescribed=actual.source.reservoir_instrument()[2*initial+final][1]
            assert np.array_equal(op,prescribed)
            branch=op@rho@op.T;branches.append(branch);reservoir.append(op@psi)
            key=f"reservoir_{initial}_{final}_input"
            max_preparation_error=max(max_preparation_error,float(abs(branch-clock[key]).max()))
            assert abs(branch-np.outer(source_columns[:,2*initial+final],source_columns[:,2*initial+final].conj())).max()<2e-13
    columns=np.column_stack(reservoir)
    expected=sum(branches);assert abs(expected-clock["input_state"]).max()<2e-13
    assert max_preparation_error<2e-13 and abs(np.linalg.norm(columns)**2-1)<1e-13
    clock_contract=read(CHAIN/"clock/CLOCK_CONTRACT.json")
    contract=read(CHAIN/"dynamical_field/DYNAMICAL_FIELD_CONTRACT.json")
    assert sha(CHAIN/"clock/CLOCK_CONTRACT.json")==read(CHAIN/"dynamical_field/MANIFEST.json")["sha256"]["openwave-tuning/next/CHAIN/clock/CLOCK_CONTRACT.json"]
    digest=hashlib.sha256(np.asarray(actual.hamiltonian,dtype="<f8").tobytes()).hexdigest()
    assert digest==contract["source"]["H_sha256_le_float64"]==clock_contract["source_H_array_sha256_le_float64"]
    assert (clock_contract["numerical_member"]["duration"],clock_contract["numerical_member"]["native_attempt_rate"],clock_contract["numerical_member"]["uniformization_rate"],clock_contract["numerical_member"]["law"])==(.37,.8,1.6,"MIX")
    max_kerror=0.
    for seam in actual.field.operators_by_law["MIX"]:
        for meta,op in seam:
            key=tuple(meta[k] for k in ("seam","law","coin","record_value"))
            max_kerror=max(max_kerror,float(abs(ops[key].toarray()[np.ix_(order,order)]-op/math.sqrt(3)).max()))
    assert max_kerror<1e-15
    for mask,k in states:
        f=offsets[mask]+k*cycle
        assert tuple(np.array([mask>>j&1 for j in range(6)])+incidence@f+np.array([3-mask.bit_count(),0,0,0,0,0]))==(1,1,1,0,0,0)
        assert tuple(f)==prod.integer_record(mask,k)
    source_f=records[order];radius=int(abs(source_f).max());assert radius==5
    periodic=actual.field.electric_values
    energy_defect=.185*(source_f*source_f).sum(axis=1)-periodic
    assert abs(abs(energy_defect).max()-13.32)<1e-13
    probability=np.sum(abs(columns)**2,axis=1);defect_weight=float(probability[abs(energy_defect)>1e-12].sum())
    assert defect_weight>0
    # The spatial chart is read only as a comparison map, with its own tree.
    spatial=np.array([prod.lift_record(list(f),3-mask.bit_count())["lifted_records"] for mask,f in actual.states])
    assert np.all((source_f-spatial)%6==0) and np.all((source_f-spatial)@incidence.T==0)
    x=BETA*T;z=RATE*T;y=x+z
    assert y==F(14837,5000)
    # Exact binomial convolution for all coefficients used by finite controls.
    for degree in range(41):
        lhs=sum((z**n/F(math.factorial(n))*x**(degree-n)/F(math.factorial(degree-n)) for n in range(degree+1)),F())
        assert lhs==y**degree/F(math.factorial(degree))
    results=read(CHAIN/"dynamical_field/RESULT.json");members=[]
    for saved in results["members"]:
        modulus=saved["modulus"];d=(modulus-1)//2-radius+1
        full=rational_tail(y,d);terminal=rational_tail(y,d-1);passive=rational_tail(x,d)
        for key,value in (("complete_marked_clock_TV_theorem_bound",full),("marked_clock_plus_terminal_MIX_TV_bound",terminal),("passive_TV_theorem_bound",passive)):
            assert F(saved[key]["exact_rational_upper"])==value
        members.append({"modulus":modulus,"full_clock_TV_bound":str(full),"terminal_MIX_TV_bound":str(terminal)})
    assert rational_tail(y,16)<F(4185417,10**12) and rational_tail(y,-2)==1
    diagnostics=[];lost_holds=0.;complete_keys=None
    for modulus in (6,42):
        fam=six if modulus==6 else independent_family(modulus,offsets,cycle)
        states,index,records,electric,h,ops=fam
        initial=np.zeros((len(states),4),complex)
        for row,state in enumerate(source_labels):initial[index[state]]=columns[row]
        producer_family=prod.Family(modulus)
        assert producer_family.states==tuple(states)
        hdifference=(h-producer_family.h).tocoo()
        assembly_error=float(abs(hdifference.data).max()) if hdifference.nnz else 0.
        # Independent summation groups the diagonal terms differently.
        # Off-diagonal entries agree exactly; this is a floating diagnostic.
        assert np.all(hdifference.row==hdifference.col) and assembly_error<1e-12
        if modulus==6:
            assert np.array_equal(producer_family.h.toarray()[np.ix_(order,order)],actual.hamiltonian)
        for key,op in producer_family.instruments():
            difference=op-ops[key]
            assert difference.nnz==0 or abs(difference.data).max()<1e-15
        complete=sum(op.T@op for op in ops.values())-sparse.eye(len(states))
        assert complete.nnz==0 or abs(complete.data).max()<1e-14
        final=expm_multiply(-.37j*h,initial)
        reference=arr[f"M{modulus}_passive_final_columns"]
        phase=np.vdot(reference,final)/abs(np.vdot(reference,final))
        error=float(np.linalg.norm(final-phase*reference));assert error<1e-11
        marks=list(ops);complete_keys=marks
        probs=[float(np.linalg.norm(ops[key]@final)**2) for key in marks]
        saved=next(r for r in results["members"] if r["modulus"]==modulus)
        assert max(abs(p-row["probability"]) for p,row in zip(probs,saved["complete_terminal_MIX_records"]))<2e-12
        assert abs(sum(probs)-1)<2e-11
        if modulus==6:
            lost_holds=sum(p for p,(e,law,coin,value) in zip(probs,marks) if (value if law=="C" else ((value-(2*coin-1)*abs(value))//2))==0)
            changed=expm_multiply(-.37j*(h-sparse.diags(electric)),initial)
            removal_error=float(np.linalg.norm(changed-final));assert removal_error>.1
        # Retain all 30 native marks and four reservoir columns in one event.
        middle=expm_multiply(-.17j*h,initial)
        marked=np.concatenate([op@middle for op in ops.values()],axis=1)
        continued=expm_multiply(-.20j*h,marked)
        archived=arr[f"M{modulus}_one_event_final_columns"]
        event_error=float(np.linalg.norm(continued-phase*archived));assert event_error<2e-11
        assert continued.shape[1]==120 and abs(np.linalg.norm(continued)**2-1)<2e-11
        diagnostics.append({"M":modulus,"independent_diagonal_assembly_roundoff":assembly_error,
                            "passive_archive_error":error,"one_event_archive_error":event_error,
                            "terminal_marks":len(marks),"one_event_columns":continued.shape[1]})
    assert lost_holds>0.1
    mutants={}
    for name,hm,prep,time in (("electric_removed",actual.hamiltonian-actual.field.electric,expected,T),
                               ("reservoir_replaced_by_reference",actual.hamiltonian,rho,T),
                               ("clock_doubled",actual.hamiltonian,expected,2*T)):
        try:prod.admit_actual_input(hm,prep,time,actual);mutants[name]=False
        except ValueError:mutants[name]=True
    assert all(mutants.values())
    return {"M6_H_difference":herr,"M6_all_MIX_difference":max_kerror,"H_sha256":digest,
            "preparation_branch_max_error":max_preparation_error,"reservoir_probabilities":[float(np.linalg.norm(v)**2) for v in columns.T],
            "initial_record_radius_all384":radius,"M6_energy_chart_max_defect":float(abs(energy_defect).max()),
            "M6_preparation_energy_defect_weight":defect_weight,"M6_preparation_mean_energy_defect":float(probability@energy_defect),
            "bounded_offdiagonal_norm":"361/50","coefficient_convolutions_exact":41,"members":members,
            "executed_diagnostics":diagnostics,"drop_electric_vector_distance":removal_error,
            "terminal_probability_lost_by_deleting_holds":lost_holds,"admission_mutations":mutants}


def force_adapter(nb,lap):
    contract=read(CHAIN/"force/FORCE_CONTRACT.json")
    transfer=read(CHAIN/"force/RESULT.json")["fixedpoint_variant_transfer"]
    theta=F(23,24);alpha=theta**512;eta=2*alpha+F(408,2**40)*(1-alpha)
    delta=eta+2*alpha;ell=4*delta;tv=F(96,1000)*ell/(1-ell)
    assert ell<1
    # Rho-pair energy is resistance/2. The resistance metric implies
    # a one-edge change at most 1/2; two ideal-stencil errors add 2 epsilon.
    # e<11/4 follows by bounding its factorial tail after order three.
    e_upper=sum((F(1,math.factorial(j)) for j in range(4)),F())+F(1,24)/(1-F(1,5))
    assert e_upper<F(11,4)
    assert F(11,4)**2/(1-16*alpha)<8
    for key,value in (("zero_start_field_error_upper",eta),("unit_energy_difference_to_ideal_frozen_upper",delta),("log_seam_rate_ratio_absolute_upper",ell),("complete_marked_capped_path_TV_upper",tv)):
        assert F(transfer[key])>=value
    lower=F(1,1000)-F(660,2*1000**2);assert lower==F(67,100000)
    score=2*F(transfer["complete_marked_capped_path_TV_upper"])/lower
    assert F(contract["rational_force_score_transfer"]["frozen_score_error_upper"])==score
    green=la.pinvh(lap.astype(float));maximum=0.
    for row in transfer["compiled_pairs"]:
        x,y=row["state"];rho=np.zeros(96,int);rho[x]+=1;rho[y]-=1
        z=fixed_solve([0]*96,rho,nb,512);energy=F(z[x]-z[y],2**41)
        assert energy==F(row["exact_energy"])
        residual=2**40*rho-np.array([6*z[x]-sum(z[y] for y in nb[x]) for x in range(96)],dtype=object)
        cert=F(2*sum(int(r)**2 for r in residual),2**80)
        assert cert==F(row["residual_field_error_squared_upper"])
        measured=abs(float(energy)-float(rho@green@rho/2));maximum=max(maximum,measured)
        assert measured<float(eta)+1e-13
    assert len(transfer["compiled_pairs"])==22
    assert contract["boundaries"]["static_field_refinement_implies_dynamic_force_refinement"] is False
    return {"exact_eta":str(eta),"exact_delta":str(delta),"exact_log_rate_ratio_bound":str(ell),
            "complete_marked_capped_path_TV_upper":transfer["complete_marked_capped_path_TV_upper"],
            "rational_score_error_upper":str(score),"incident_rate_upper":8,
            "total_attempt_rate_upper":660,"compiled_pair_count":22,"maximum_Green_energy_error_diagnostic":maximum,
            "comparison_domain":"ideal exact-real exponential rates and ideal frozen local stencil; floating library error excluded"}


def main():
    before=verify_pins()
    internal,edges,incidence,lap,nb=frame_and_graph()
    pivots=gap_certificates(internal)
    arithmetic,charge,prod=arithmetic_audit(edges,incidence,lap,nb)
    static=static_audit(internal,charge,prod)
    print("Independent fixed-point, 539-record and full-symbol checks passed.",flush=True)
    dynamical=dynamical_audit()
    print("Independent actual M6, full preparation and M42 marked-event checks passed.",flush=True)
    force=force_adapter(nb,lap)
    after=verify_pins();assert before==after
    report={"schema":"oph.chain.independent_field_audit.v1","independent_review":True,
            "producer_execute_or_verify_called":False,"source_artifacts_unchanged":True,
            "frozen_manifests":PINS,"checked_file_sha256":before,
            "side_two_exact_positive_pivots":pivots,"fixed_point_and_records":arithmetic,
            "static_full_symbol":static,"dynamical_integer_field":dynamical,"force_adapter":force,
            "official_rows_changed":[]}
    (HERE/"AUDIT_RESULT.json").write_text(json.dumps(report,indent=2)+"\n")
    claims=[
        ("fixed_point_tracking","PASS","Exact integer rounding, energy-norm contraction, every recorded move and residual certificate independently reconstructed."),
        ("static_smoothed_3D_Coulomb","PASS","All twelve Bloch branches and 26 Fejer modes retained; localized O(h), internally averaged O(h^2), electric vector bounds on the unit torus."),
        ("actual_C10_M6_baseline","PASS","Actual H, complete MIX matrices, global reference, all four reservoir outcomes and clock pins independently compared."),
        ("interacting_integer_register_limit","PASS","Exact Gauss fiber, bounded Dyson support and all-count convolution give the complete marked-law bound; M42 TV below 4.185417e-6."),
        ("extra_terminal_MIX","PASS","One extra support increment is charged; M42 TV below 2.286972e-5."),
        ("classical_fixedpoint_energy_rate_adapter","PASS","Exact energy and logarithmic-rate allowances imply complete capped-path TV; exact-real evaluation domain retained."),
        ("M6_small_full_clock_limit_error","FAIL","The support theorem yields TV<=1 at actual M6; passive numerical discrepancy is not a full-clock certificate."),
        ("M6_integer_chart_equals_periodic_electric_energy","FAIL","The maximal electric-energy discrepancy is 13.32 and the full actual preparation has positive defect weight."),
        ("classical_force_to_C10_join","FAIL","No map identifies classical reciprocal scheduler, C10 Hamiltonian, full preparation and clock; labels do not supply a transfer."),
        ("static_charge_lift_is_dynamic_refinement","FAIL","The same physical coordinate charge lift does not intertwine a coarse edge hop or clock with the fine source."),
        ("Maxwell_mechanical_force_or_physical_seconds","FAIL","No common spatial quantum refinement, transverse electromagnetic modes, mechanical inertia, calibrated force units or native coherent clock compiler is supplied."),
    ]
    contract={"schema":"oph.chain.field_audit_contract.v1","independent_review":True,
              "scope":"Independent audit of frozen field-limit, C10 integer-field and classical arithmetic force-adapter claims.",
              "audited_manifests":PINS,"claims":[{"claim":name,"status":status,"finding":finding} for name,status,finding in claims],
              "supported_claims_pass":True,"excluded_physical_claims_are_failed_transfers":True,
              "exact_remaining_physical_gap":"An executable common source, preparation and clock transfer connecting the classical cubic Coulomb/feedback instrument to C10, together with a spatial dynamical quantum field refinement with transverse electromagnetic observables and a calibrated mechanical force/inertia dictionary. The proved C10 limit changes flux-register capacity on a fixed one-cycle triangle; it supplies no spatial Maxwell limit. Coherent H services have no supplied native compiler or physical-second identification.",
              "numeric_certification_boundary":"Analytic rational certificates are independent of floating propagation diagnostics. General floating-library exponential/eigensolver error is unproved. Ground preparation follows the declared source prescription and is compared numerically.",
              "audit_artifacts":{"audit.py":sha(HERE/"audit.py"),"PROOF.md":sha(HERE/"PROOF.md"),"AUDIT_RESULT.json":sha(HERE/"AUDIT_RESULT.json")},
              "producer_mutations":[],"official_rows_changed":[]}
    (HERE/"AUDIT_CONTRACT.json").write_text(json.dumps(contract,indent=2)+"\n")
    print(json.dumps({"supported_claims_pass":True,"records":539,"M42_bound":dynamical["members"][-1]["full_clock_TV_bound"],"contract_sha256":sha(HERE/"AUDIT_CONTRACT.json")},indent=2))


if __name__=="__main__":main()
