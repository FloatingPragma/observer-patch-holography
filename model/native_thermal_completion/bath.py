"""Complete reversible proposal records, equilibrium and finite cooling laws."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from model.native_abundance_epochs import model as native
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
QUANTA=(F(1),F(3,4),F(1,2),F(1,4))


@lru_cache(maxsize=1)
def source():
    states=native.STATES;index={x:i for i,x in enumerate(states)};kernel=native.kernel()
    energy=[native.energy(x) for x in states]
    directed={(index[x],index[y]):p for x,row in kernel.items() for y,p in row.items() if x!=y}
    proposal={key:max(directed.get(key,F(0)),directed.get(key[::-1],F(0)))
              for key in set(directed)|{key[::-1] for key in directed}}
    sums=[sum((p for (x,_),p in proposal.items() if x==i),F(0)) for i in range(len(states))]
    uniformization=2*max(sums)
    graph=csr_matrix((np.ones(len(proposal)),tuple(zip(*proposal))),shape=(len(states),len(states)))
    count,labels=connected_components(graph,directed=False)
    components=[np.where(labels==k)[0].tolist() for k in range(count)]
    if count!=26 or any(sum(states[x][:3])!=sum(states[y][:3]) for x,y in proposal):
        raise ArithmeticError('thermal support or conserved load changed')
    if any((energy[y]-energy[x])%2 for x,y in proposal):raise ArithmeticError('energy quantum is not integral')
    conjugation=[index[tuple(5-z for z in state[:3])+((-state[3])%6,)] for state in states]
    if any(proposal.get((conjugation[x],conjugation[y]))!=p for (x,y),p in proposal.items()):
        raise ArithmeticError('proposal is not charge-conjugation equivariant')
    return dict(states=states,index=index,energy=energy,proposal=proposal,row_sums=sums,
        uniformization=uniformization,components=components,conjugation=conjugation,
        reverse_edges_added=sum(key not in directed for key in proposal),native_directed_edges=len(directed))


@lru_cache(maxsize=8)
def instrument(q):
    q=F(q)
    if not 0<q<=1:raise ValueError('0<q=exp(-2 beta)<=1 required')
    d=source();n=len(d['states']);lam=d['uniformization'];records=[];rows=[{} for _ in range(n)]
    for (x,y),base in sorted(d['proposal'].items()):
        delta=(d['energy'][y]-d['energy'][x])//2
        accept=q**max(0,delta);weight=base/lam
        records.append((x,y,y,'accepted',weight*accept))
        records.append((x,y,x,'rejected',weight*(1-accept)))
        rows[x][y]=rows[x].get(y,F(0))+weight*accept
        rows[x][x]=rows[x].get(x,F(0))+weight*(1-accept)
    for x in range(n):
        idle=1-d['row_sums'][x]/lam
        records.append((x,x,x,'idle',idle));rows[x][x]=rows[x].get(x,F(0))+idle
    if any(sum(row.values())!=1 or min(row.values())<0 for row in rows):
        raise ArithmeticError('complete bath row lost original probability')
    for (x,y) in d['proposal']:
        delta=(d['energy'][y]-d['energy'][x])//2
        if rows[x][y]/rows[y][x]!=q**delta:raise ArithmeticError('exact local detailed balance failed')
    conj=d['conjugation']
    if any(rows[conj[x]].get(conj[y],F(0))!=p for x,row in enumerate(rows) for y,p in row.items()):
        raise ArithmeticError('bath violates native conjugation')
    entries=[(x,y,float(p)) for x,row in enumerate(rows) for y,p in row.items()]
    matrix=csr_matrix(([p for _,_,p in entries],([x for x,_,_ in entries],[y for _,y,_ in entries])),shape=(n,n))
    return rows,records,matrix


def component_weights():
    d=source();p=native.initial_distribution()
    return [sum((p.get(d['states'][i],F(0)) for i in comp),F(0)) for comp in d['components']]


def equilibrium(q):
    d=source();q=F(q);out=[F(0)]*len(d['states']);weights=component_weights();parts=[]
    for component,probability in zip(d['components'],weights,strict=True):
        floor=min(d['energy'][i] for i in component)
        w=[q**((d['energy'][i]-floor)//2) for i in component];z=sum(w,F(0))
        for i,a in zip(component,w,strict=True):out[i]=probability*a/z
        parts.append({'states':len(component),'total_load':sum(d['states'][component[0]][:3]),
            'energy_floor':floor,'partition_sum_relative':str(z),'preparation_component_weight':str(probability)})
    rows,_,_=instrument(q);advanced=[F(0)]*len(out)
    for x,row in enumerate(rows):
        for y,p in row.items():advanced[y]+=out[x]*p
    if advanced!=out or sum(out)!=1:raise ArithmeticError('exact equilibrium is not stationary')
    return out,parts


def entropy(p):
    positive=p[p>0]
    return float(-np.dot(positive,np.log(positive)))


def shares(p):
    result=np.zeros(6)
    for state,value in zip(source()['states'],p,strict=True):
        result[state[3]]+=value/2;result[-state[3]%6]+=value/2
    return result.tolist()


def cooling(steps):
    if type(steps) is not int or steps<1:raise ValueError('positive steps per temperature required')
    d=source();initial,_=equilibrium(F(1));p=np.array([float(x) for x in initial]);energy=np.array(d['energy'])
    initial_error=sum((abs(F.from_float(float(x))-x) for x in initial),F(0))
    error=initial_error
    plan=(F(3,4),F(1,2),F(1,4),F(1,2),F(3,4),F(1));rows=[];total_sigma=0.;total_heat=0.;worst=0.
    for q in plan:
        exact,_,matrix=instrument(q);target,_=equilibrium(q);gibbs=np.array([float(x) for x in target]);beta=-math.log(float(q))/2
        coefficient_error=max(sum((abs(F.from_float(float(x))-x) for x in row.values()),F(0)) for row in exact)
        indegree=max(Counter(y for row in exact for y,p in row.items() if p).values())
        gamma=F(2*indegree,2**53-2*indegree)
        delta=coefficient_error+gamma*(1+coefficient_error)
        smallest_matrix_value=min(matrix.data[matrix.data>0])
        start_energy=float(p@energy);start_entropy=entropy(p);stage_sigma=0.;smallest=math.inf
        for _ in range(steps):
            if float(min(p[p>0]))*smallest_matrix_value<=np.finfo(float).tiny:
                raise ArithmeticError('relative-roundoff domain includes underflow')
            following=np.asarray(p@matrix).ravel();heat=float((following-p)@energy)
            sigma=entropy(following)-entropy(p)-beta*heat
            if sigma< -2e-12:raise ArithmeticError('finite bath entropy production became negative')
            smallest=min(smallest,sigma);stage_sigma+=sigma;worst=max(worst,abs(float(following.sum())-1))
            p=following
        # Positive sparse multiplication: at most two rounded operations per incoming edge.
        exact_error=(1+error)*(1+delta)**steps-1
        scaled=exact_error*2**80
        error=F((scaled.numerator+scaled.denominator-1)//scaled.denominator,2**80)
        heat=float(p@energy)-start_energy
        rows.append({'q_exp_minus_two_beta':str(q),'beta_in_declared_cost_units':beta,'steps':steps,
            'expected_cost':float(p@energy),'cost_heat':heat,'entropy_change':entropy(p)-start_entropy,
            'entropy_production':stage_sigma,'minimum_step_entropy_production':smallest,
            'distance_to_component_Gibbs_TV':float(np.sum(abs(p-gibbs))/2),
            'state_probability_L1_roundoff_bound':str(error),
            'state_probability_L1_roundoff_bound_numeric':float(error),
            'roundoff_domain':'IEEE binary64 products/sums, no overflow or underflow; no renormalization',
            'face_class_shares':shares(p),'complete_state_probabilities':p.tolist()})
        total_sigma+=stage_sigma;total_heat+=heat
    if worst>2e-12 or np.max(abs(p-p[np.array(d['conjugation'])]))>2e-12:
        raise ArithmeticError('complete law or conjugation was lost')
    return {'steps_per_temperature':steps,'total_steps':6*steps,'initial_state':'componentwise q=1 Gibbs law with the original native preparation component weights',
        'initial_expected_cost':float(np.array([float(x) for x in initial])@energy),
        'total_heat':total_heat,'total_entropy_production':total_sigma,
        'maximum_probability_mass_error':worst,'stages':rows,'state_probabilities_renormalized':False}


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in {'.py','.md'}}
    paths.update(['model/native_abundance_epochs/PROOF.md','model/defect_integer_lift/federation.json',
                  'Lean/Dynamics/ChargeConjugation.lean'])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def certificate():
    pins=bindings();d=source();temperatures=[]
    for q in QUANTA:
        _,records,_=instrument(q);pi,parts=equilibrium(q)
        encoded=[[x,proposal,target,label,str(p)] for x,proposal,target,label,p in records]
        temperatures.append({'q':str(q),'complete_proposal_records':len(records),
            'zero_probability_records_retained':sum(p==0 for *_,p in records),
            'record_sha256':hashlib.sha256(json.dumps(encoded,separators=(',',':')).encode()).hexdigest(),
            'detailed_balance_exact':True,'stationarity_exact':True,'conjugation_exact':True,
            'components':parts,'equilibrium_face_class_shares':[str(sum((p*(int(state[3]==k)+int((-state[3])%6==k))/2
                 for state,p in zip(d['states'],pi,strict=True)),F(0))) for k in range(6)]})
    result={'schema':'oph.native_thermal_completion.v1','source_bindings':pins,'states':len(d['states']),
        'state_order':[list(x) for x in d['states']],
        'complete_proposal_edges':[[x,y,str(weight),d['energy'][y]-d['energy'][x]] for (x,y),weight in sorted(d['proposal'].items())],
        'complete_idle_weights':[str(1-x/d['uniformization']) for x in d['row_sums']],
        'native_directed_nonidentity_edges':d['native_directed_edges'],'reverse_edges_added':d['reverse_edges_added'],
        'weak_components':len(d['components']),'uniformization_rate_exact':str(d['uniformization']),
        'thermal_energy':'E=5 sum(load^2)+2[C(k)+C(-k)]','q_definition':'q=exp(-2 beta)',
        'proposal_rule':'s_xy=max(K_xy,K_yx); accepted=s_xy/Lambda * min(1,q^((E_y-E_x)/2)); rejected and idle records retained',
        'temperature_rows':temperatures,'cooling_experiments':[cooling(n) for n in (32,256)],
        'scope':{'law_is_explicit_bath_extension':True,'unchanged_native_repair_claimed_thermal':False,
                 'component_weights_selected_by_Gibbs_law':False,'physical_energy_or_cosmic_clock_calibrated':False,
                 'source_particles_or_cosmological_relic_abundances_identified':False,
                 'all_proposal_rejections_and_idle_records_retained':True}}
    if bindings()!=pins:raise ArithmeticError('thermal completion source changed during execution')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
    result=certificate();text=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n';path=HERE/'results.json'
    if a.write:path.write_text(text)
    if a.check and path.read_text()!=text:raise ArithmeticError('thermal completion receipt changed')
    print('Exact reversible bath, complete records and finite cooling laws verified.')
