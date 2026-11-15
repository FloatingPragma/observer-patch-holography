"""Opposing review of the uncapped field-feedback law; no acceptance run."""
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh, expm
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
from model.native_spatial_refinement import source as routing


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def run():
    paths=[HERE/'field_feedback.py',HERE/'FIELD_FEEDBACK.json',
           ROOT/'reverse-engineering-reality/model/native_spatial_refinement/source.py']
    pins={str(p.relative_to(ROOT)):sha(p) for p in paths}
    receipt=json.loads(paths[1].read_text())
    g=routing.graph(2)
    edges=g.seams
    n=g.ports
    lap=np.zeros((n,n),dtype=int)
    for a,b in edges:
        lap[a,a]+=1;lap[b,b]+=1;lap[a,b]-=1;lap[b,a]-=1
    # Derive each exact block from the full actual graph, independently of
    # the producer's link-symbol formula and LDL decomposition.
    minors=[]
    columns=[]
    for frequency in np.ndindex(2,2,2):
        chars=np.array([(-1)**sum(a*b for a,b in zip(frequency,x))
                        for x in g.cell_coordinates],dtype=int)
        basis=np.kron(chars[:,None],np.eye(12,dtype=int))
        raw=basis.T@lap@basis
        assert np.all(raw%8==0)
        block=raw//8
        assert np.array_equal(lap@basis,basis@block)
        positive=12*s.Matrix(block)-6*s.eye(12)
        if frequency==(0,0,0): positive+=s.ones(12)
        dets=[positive[:k,:k].det(method='domain-ge') for k in range(1,13)]
        assert all(d>0 for d in dets)
        minors.append({'frequency':frequency,'positive_leading_minors':list(map(str,dets))})
        columns.append(basis)
    orth=np.concatenate(columns,axis=1)
    assert np.array_equal(orth.T@orth,8*np.eye(n,dtype=int))
    vals,vecs=eigh(lap)
    f=np.zeros(n);f[1:]=-np.expm1(512*np.log1p(-np.minimum(vals[1:],12)/12))/vals[1:]
    field=(vecs*f)@vecs.T
    energies=(np.diag(field)+field[0,0]-2*field[:,0])/2
    max_law=max_balance=max_drift=0.
    for row in receipt['interventions']:
        k=row['kappa'];q=np.zeros((n,n))
        for a,b in edges:
            diff=energies[b]-energies[a]
            ab=math.exp(-2*k*diff)/2;ba=math.exp(2*k*diff)/2
            q[b,a]+=ab;q[a,a]-=ab;q[a,b]+=ba;q[b,b]-=ba
        stationary=np.exp(-4*k*energies);stationary/=sum(stationary)
        balance=q*stationary[None,:]
        max_balance=max(max_balance,float(np.max(abs(balance-balance.T))))
        output=expm(q/2)@np.full(n,1/n)
        max_law=max(max_law,float(np.max(abs(output-row['full_position_law_at_half_time']))))
        drift=float(energies@(q@np.full(n,1/n)))
        max_drift=max(max_drift,abs(drift-row['initial_energy_current_rate']))
        assert abs(sum(output)-1)<1e-12 and min(output)>=0
    assert max_law<1e-12 and max_balance<1e-14 and max_drift<1e-12
    max_log=0.;events=holds=0
    for path in receipt['complete_paths']+[receipt['cap_failure_example'],receipt['zero_duration_example']]:
        x=path['initial'];t=0.;logp=0.;k=path['kappa']
        for event in path['records']:
            def rate(a,b):
                if x not in (a,b): return 1.
                y=b if x==a else a
                return math.exp(-2*k*(energies[y]-energies[x]))
            rates=[rate(a,b) for a,b in edges];total=sum(rates)
            assert abs(total-event['total_seam_rate'])<1e-10
            a,b=edges[event['seam']]
            nxt=(a if event['coin'] else b) if x in (a,b) else x
            current=int(x==a)-int(nxt==a)
            assert current==event['current']
            logp-=total*event['wait'];logp+=math.log(rates[event['seam']]/2)
            x=nxt;t+=event['wait'];events+=1;holds+=current==0
        if path['terminal']['kind']=='no_further_event':
            total=len(edges)
            for a,b in edges:
                if x in (a,b):
                    y=b if x==a else a
                    total+=math.exp(-2*k*(energies[y]-energies[x]))-1
            logp-=total*(path['horizon']-t)
        max_log=max(max_log,abs(logp-path['log_path_density']))
        assert path['final_counter']==[int(i==x) for i in range(n)]
    assert max_log<1e-9
    assert pins=={str(p.relative_to(ROOT)):sha(p) for p in paths}
    result={'schema':'oph.native_field_feedback.opposing_uncapped.v1',
      'checks':{'actual_cubic_routing_used':True,'complete_exact_Fourier_decomposition':True,
                'gap_half_by_independent_Sylvester_minors':True,'complete_uncapped_output_laws_match':True,
                'detailed_balance_and_energy_response_match':True,'fine_native_holds_waits_and_path_densities_match':True,
                'source_pins_unchanged':True},
      'exact_blocks':minors,'position_law_max_error':max_law,'energy_drift_max_error':max_drift,
      'detailed_balance_max_error':max_balance,'fine_path_log_density_max_error':max_log,
      'replayed_events_including_cap_control':events,'retained_holds':holds,
      'source_sha256':pins,'reviewer_sha256':sha(Path(__file__)),
      'scope':'Finite, installed energy-dependent scheduler on actual native C process. Full96-position law here is uncapped. Finite-cap probability is a separate count-augmented instrument obligation; its omission cannot be normalized away. No coupling-strength, Maxwell, Lorentz or orbital criterion is awarded.',
      'official_rows_changed':[]}
    (HERE/'ROOT_UNCAPPED_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('checks','position_law_max_error','fine_path_log_density_max_error')}))


if __name__=='__main__':run()
