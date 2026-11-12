"""Independent compressed-PW string and constrained-energy review."""
from collections import deque
from itertools import product
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

import sympy as s

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=ROOT/'model/openwave_reference_binding'


def run():
    receipt=json.loads((PACKAGE/'results.json').read_text())
    for path,expected in receipt['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms
    from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current,label_index
    from model.coherent_defect_wave.physical_attachment.chiral import representation
    edges=[(term.path[0],term.left_cell,term.right_cell) for term in kinetic_terms(Configuration())]
    adjacency={x:[] for x in range(8)}
    for edge,a,b in edges:
        adjacency[a].append((edge,a,b,1));adjacency[b].append((edge,b,a,0))
    distances={0:0};todo=deque([0])
    while todo:
        vertex=todo.popleft()
        for _,_,target,_ in adjacency[vertex]:
            if target not in distances:distances[target]=distances[vertex]+1;todo.append(target)
    paths={0:[()]}
    for level in (1,2,3):
        for target in [x for x in distances if distances[x]==level]:
            paths[target]=[path+((edge,start,target,forward),)
                           for start in distances if distances[start]==level-1
                           for edge,_,end,forward in adjacency[start] if end==target
                           for path in paths[start]]
    for row in receipt['admissible_native_strings']['rows']:
        actual={tuple(tuple(edge) for edge in path) for path in row['shortest_paths']}
        assert actual==set(paths[row['cell']])

    # Compress the actual119-state currents to each complete path-irrep
    # tensor factor, then apply independent Kronecker Gauss operators.
    gauges=representation()[2];labels=label_index();checks=0
    for family,indices in (('ec',(9,)),('L',(10,11))):
        d=len(indices)
        for length in range(4):
            for orientations in product((0,1),repeat=length):
                factors=[d]+[d*d]*length+[d]
                dimension=s.prod(factors)
                state=s.MutableSparseMatrix(dimension,1,{})
                def linear(address):
                    result=0
                    for coordinate,size in zip(address,factors):result=result*size+coordinate
                    return result
                for word in product(range(d),repeat=length+1):
                    address=[word[0]]+[word[j]*d+word[j+1] for j in range(length)]+[word[-1]]
                    state[linear(address)]=d**(-s.Rational(length+1,2))
                assert (state.H*state)[0]==1
                rowlists=[]
                for forward in orientations:
                    rowlists.append([labels[(1,indices[a],indices[b])] if forward
                                     else labels[(-1,indices[b],indices[a])]
                                     for a,b in product(range(d),repeat=2)])
                def apply_factor(matrix,slot):
                    # Apply the tensor factor directly to the finitely many
                    # nonzero contracted-string entries, without constructing
                    # a mostly zero 256-by256 Kronecker operator.
                    stride=s.prod(factors[slot+1:]);size=factors[slot]
                    result=s.MutableSparseMatrix(dimension,1,{})
                    for (address,_),amplitude in state.todok().items():
                        coordinate=(address//stride)%size
                        for target in range(size):
                            value=matrix[target,coordinate]
                            if value:result[address+(target-coordinate)*stride]+=value*amplitude
                    return result
                for current,t in enumerate(gauges):
                    native=t.extract(indices,indices)
                    link_actions=[]
                    for forward,rows in zip(orientations,rowlists):
                        pair=[]
                        for endpoint in ((0,1) if forward else (1,0)):
                            full=endpoint_current(current,endpoint)
                            selected=full.extract(rows,rows)
                            outside=[r for r in range(119) if r not in rows]
                            assert full.extract(outside,rows)==s.zeros(len(outside),len(rows))
                            pair.append(selected)
                        link_actions.append(pair)
                    for node in range(length+1):
                        action=s.MutableSparseMatrix(dimension,1,{})
                        if node==0:action+=apply_factor(native,0)
                        if node==length:action+=apply_factor(-native.T,len(factors)-1)
                        if node>0:action+=apply_factor(link_actions[node-1][1],node)
                        if node<length:action+=apply_factor(link_actions[node][0],node+1)
                        assert action==s.zeros(dimension,1)
                        checks+=1
                    if current==11:
                        wrong=action+2*apply_factor(native.T,len(factors)-1)
                        assert wrong!=s.zeros(dimension,1)

    # Rebuild the actual weak endpoint from the local highest-weight
    # basis; infer its representation through the positive Gram inverse.
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import target
    from model.coherent_defect_wave.chiral_joint.local_reference import physical
    from model.coherent_defect_wave.chiral_joint.local_reference.reduction import quark_module,lepton_operators
    from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
    qw,lw,sw=physical.labels();q=quark_module();leptons=lepton_operators();scalar=higgs_memory()[1]
    endpoint_checks=[]
    for family,indices,number,hyper in (('ec',(9,),4,-2),('L',(10,11),2,1)):
        block=target(family,0)
        keys=tuple((i,mask,h) for i,(_,qy) in enumerate(qw)
                   for mask,(n,_,ly,spin) in enumerate(lw) if (n,spin)==(number,0)
                   for h,(_,sy) in enumerate(sw) if qy+ly+sy==hyper)
        positions={key:i for i,key in enumerate(keys)}
        high=s.zeros(len(keys),block['basis'].cols)
        for row,key in enumerate(block['keys']):high[positions[key],:]=block['basis'][row,:]
        currents=[physical.tensor_map(keys,keys,quark=q['currents'][j],lepton=leptons['currents'][j],scalar=scalar[8+j]) for j in range(4)]
        basis=high if family=='ec' else s.Matrix.hstack(-(currents[0]-s.I*currents[1])*high,high)
        metric=physical.metric_on(keys);gram=basis.H*metric*basis
        assert gram==s.kronecker_product(s.eye(len(indices)),block['gram'])
        left=gram.inv()*basis.H*metric
        for j,current in enumerate(currents):
            induced=left*current*basis
            expected=s.kronecker_product(-gauges[8+j].extract(indices,indices).T,s.eye(high.cols))
            assert induced==expected and current*basis==basis*induced
        induced_h=left*physical.hamiltonian_on(keys)*basis
        assert induced_h==s.kronecker_product(s.eye(len(indices)),block['H'])
        endpoint_checks.append({'family':family,'exact_endpoint_weight_rows':len(keys),'multiplicity':high.cols})

    source=json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json').read_text())
    local=source['exact_local_channels'];x=s.symbols('x')
    for family in ('ec','L'):
        result=next(row for row in receipt['constrained_minima'] if row['family']==family)
        gap=s.Rational(local[family]['full_compensated_local_gap_above_neutral_ground_lower'])
        trial=s.Rational(result['actual_native_local_trial_gap_upper'])
        assert gap>0 and trial>=gap
        for row in result['rows']:
            ell=row['distance'];bounds=list(map(s.Rational,row['binding_cost_relative_to_origin_interval']))
            if family=='L' or ell==0:assert bounds==[ell,ell]
            else:assert bounds==[ell+gap-trial,ell+trial-gap] and bounds[0]>0
        if family=='L':
            blocks=local[family]['complete_compensated_local_maximal_module_blocks']
            minimum=min(blocks,key=lambda row:s.Rational(row['ground_interval'][0]))
            poly=s.Poly(s.sympify(minimum['complete_characteristic_polynomial']),x)
            visible=s.Poly(s.sympify(local[family]['complete_target_annihilator']),x)
            interval=list(map(s.Rational,local[family]['complete_local_energy_records'][0]['energy_interval']))
            assert s.gcd(poly,visible).count_roots(*interval)==1
            assert interval[1]<s.Rational(local[family]['every_other_color_admissible_quark_module_lower'])
        else:assert not local[family]['lowest_compensated_local_energy_exactly_identified']
    # A same-sector compressed minimum is not a stationary source eigenvalue.
    assert not receipt['scope']['Pi_y_invariant_under_actual_H']
    assert not receipt['scope']['reference_binding_removed_by_energy_zero_subtraction']
    assert not receipt['scope']['screened_cloud_minimum_or_string_breaking_computed']
    for row in receipt['screened_cloud_counterexample']['records']:
        assert row['single_unscreened_charge_projectors_all_annihilate_state']
        assert row['joint_Gauss_constraints_preserved']
    return {'schema':'oph.independent_reference_binding_review.v1','all_pass':True,
            'all_actual_shortest_stored_paths_independently_reconstructed':True,
            'independent_compressed119_link_tensor_Gauss_checks':checks,
            'actual_local_endpoint_embeddings':endpoint_checks,
            'L_minimum_trial_root_matches_complete_local_minimum':True,
            'ec_visible_pole_not_promoted_to_global_compensated_ground':True,
            'Pi_y_T_Pi_y_reason':'Every hopping monomial has a nonarrival matter endpoint; its nontrivial fermion operator has zero gauge-singlet compression.',
            'scope':'Constrained expectation minima on complete Pi_y; screened complement separate, no invariant charged band or Coulomb interpretation.',
            'reviewed_sha256':{str((PACKAGE/name).relative_to(ROOT)):hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest()
                               for name in ('strings.py','binding.py','PROOF.md','results.json')}}


if __name__=='__main__':
    result=run();target=Path(__file__).resolve().parents[1]/'reference_binding_review.json'
    target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: independent PW119 strings, native endpoints and constrained binding minima')
