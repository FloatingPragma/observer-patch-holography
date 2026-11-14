"""Particle/current/field records from the actual common autonomous API.

No Hamiltonian or force array is supplied by this consumer. Each call uses
the producer's complete sparse source map with every amplitude retained.
"""
from collections import defaultdict
from pathlib import Path
import hashlib
import gzip
import json
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'vector_field'))
import autonomous as producer


def joint_law(state):
    """Every terminal position/event/field Born probability, including zeros."""
    return [{'position':x,'all288_Z16_field_values_hex':field.hex(),
             'all12_event_probabilities':(abs(amplitude)**2).tolist()}
            for (x,field),amplitude in sorted(state.items())]


def write_joint_laws(filename, records):
    payload=json.dumps({'schema':'oph.local12.complete_particle_field_laws.v1',
                       'records':records,
                       'omitted_field_position_blocks':'exactly zero in this ideal finite source law',
                       'native_complements':'Complete source/reset/padding/pointer laws and their unconditional TV error are supplied by the native composition certificate; these ideal probabilities are not success-renormalized.'},
                      separators=(',',':')).encode()
    encoded=gzip.compress(payload,mtime=0)
    (HERE/filename).write_bytes(encoded)
    return {'path':str((HERE/filename).relative_to(ROOT)),
            'sha256':hashlib.sha256(encoded).hexdigest(),
            'uncompressed_sha256':hashlib.sha256(payload).hexdigest()}


def measure(state,charge):
    """Complete coarse law, with unnormalized branches admitted unchanged."""
    graph,B,A,labels=producer.data()
    positions=np.zeros(96);events=np.zeros(12);field_energy=0.;total=0.
    for (x,F),amplitude in state.items():
        probability=abs(amplitude)**2;weight=float(probability.sum())
        field=np.frombuffer(F,dtype=np.uint8).astype(np.int64)
        n=np.zeros(96,dtype=int);n[x]=charge
        reference=np.zeros(96,dtype=int);reference[0]=charge
        assert not np.any((n+B@field-reference)%producer.MODULUS)
        positions[x]+=weight;events+=probability;field_energy+=weight*producer.electric_energy(F);total+=weight
    return {'unconditional_norm_squared':total,'all96_position_probabilities':positions.tolist(),
            'all12_event_probabilities':events.tolist(),'electric_energy_moment':field_energy,
            'retained_source_field_blocks':len(state),'renormalized':False}


def current(state,charge):
    graph,B,A,labels=producer.data()
    neighbours,edges=producer.mb.common.geometry_tables()
    mean=np.zeros(288)
    for (x,F),amplitude in state.items():
        reflected=amplitude.sum()/6-amplitude
        for d in range(6):
            y=int(neighbours[x,d]);edge=int(edges[x,d]);a,b=graph.seams[edge]
            mean[edge]+=charge*(int(x==a)-int(y==a))*abs(reflected[2*d+1])**2
    return mean


def run():
    paths=[Path(__file__),Path(producer.__file__),Path(producer.mb.__file__),
           Path(producer.mode_field.__file__),HERE.parent/'manybody/MANIFEST.json',
           HERE.parent/'clock/MANIFEST.json',HERE/'PARTICLE_RECORDS.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    graph,B,A,labels=producer.data()
    experiments=[];finals={};current_cases=[];complete_laws=[]
    for charge,alpha,beta in ((1,0.,0.),(1,.17,.41),(-1,.17,.41)):
        initial=producer.initial(charge=charge)
        state=initial;rows=[]
        for n in range(4):
            complete_laws.append({'charge':charge,'alpha':alpha,'beta':beta,'step':n,
                                  'joint_law':joint_law(state)})
            row=measure(state,charge)
            overlap=producer.phase_correlator(initial,state)
            row.update(step=n,return_correlation=[float(overlap.real),float(overlap.imag)])
            assert abs(row['unconditional_norm_squared']-1)<2e-10
            rows.append(row)
            if n<3:
                J=current(state,charge)
                next_state=producer.step(state,alpha,beta,charge=charge)
                nxt=measure(next_state,charge)
                difference=charge*(np.asarray(nxt['all96_position_probabilities'])-np.asarray(row['all96_position_probabilities']))+B@J
                error=float(np.max(abs(difference)))
                assert error<2e-11
                current_cases.append({'charge':charge,'alpha':alpha,'beta':beta,'step':n,
                                      'all288_native_mean_currents':J.tolist(),'continuity_error':error})
                state=next_state
        finals[charge,alpha,beta]=state
        experiments.append({'charge':charge,'alpha':alpha,'beta':beta,
                            'other207_face_couplings':0.,'records':rows})
    base=np.asarray(experiments[0]['records'][-1]['all96_position_probabilities'])
    coupled=np.asarray(experiments[1]['records'][-1]['all96_position_probabilities'])
    hole=np.asarray(experiments[2]['records'][-1]['all96_position_probabilities'])
    response=float(abs(base-coupled).sum()/2)
    conjugation=float(abs(coupled-hole).sum()/2)
    assert response>1e-5 and conjugation<2e-10
    # An actual coarse position intervention at step1, followed by the SAME
    # augmented source. Its branches keep their original probabilities.
    first=producer.step(producer.initial(),.17,.41)
    branches=defaultdict(dict)
    for (x,F),value in first.items():branches[x][x,F]=value
    total=np.zeros(96);norm=0.;joint=[]
    for x,branch in sorted(branches.items()):
        for _ in range(2):branch=producer.step(branch,.17,.41)
        record=measure(branch,1)
        norm+=record['unconditional_norm_squared'];total+=record['all96_position_probabilities']
        joint.append({'first_position':x,**record})
        complete_laws.append({'alpha':.17,'beta':.41,'first_step':1,'first_position':x,
                              'last_step':3,'joint_law':joint_law(branch)})
    assert abs(norm-1)<2e-10
    backaction=float(abs(total-coupled).sum()/2)
    # Compact held-event dark mode never enters a moving slot, even under
    # this autonomous field update. Only its field state evolves.
    dark=producer.initial();key=next(iter(dark));dark[key][:]=0
    dark[key][0]=1/np.sqrt(2);dark[key][2]=-1/np.sqrt(2)
    for _ in range(3):dark=producer.step(dark,.17,.41)
    dark_record=measure(dark,1)
    assert abs(dark_record['all96_position_probabilities'][0]-1)<2e-10
    result={'schema':'oph.local12.particle_coupled_response.v1','candidate_id':producer.MODEL_ID,
        'checks':{'actual_autonomous_provider_used':True,'all_complete_source_outcomes_retained':True,
                  'native_mean_current_continuity':True,'interaction_changes_particle_position_records':True,
                  'same_source_opposite_current_hole':True,'coarse_position_backaction_not_renormalized':True,
                  'held_dark_mode_survives_field_coupling':True,'source_bytes_unchanged':True},
        'experiments':experiments,'native_current_cases':current_cases,
        'interaction_position_TV':response,'particle_hole_position_TV':conjugation,
        'two_time_position':{'first_step':1,'last_step':3,'joint_records':joint,
                             'complete_probability':norm,'measurement_backaction_TV':backaction},
        'held_dark_record':dark_record,
        'complete_terminal_and_two_time_laws':write_joint_laws('COUPLED_JOINT_LAWS.json.gz',complete_laws),
        'apparatus':'PARTICLE_RECORDS native_joint_position_full_Z16_field is a finite compute/read/uncompute native program for these one-particle/one-hole outputs. Readout source/reset errors add to the autonomous producer implementation error. No full source-error number is inferred from sparse numerical residuals.',
        'clock':'Each complete V is one declared macro step at t=n*tau0*h. Actual augmented pulse duration/error and field initialization are the producer contract; bare-W clock resource counts are not transferred.',
        'scope':'One finite face-coupled source, complete numerical quantum laws, separate particle/hole references. This is measured interaction response and continuity, not a derived two-charge force law, vacuum pole, observed mass or homogeneous Maxwell limit.',
        'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    (HERE/'COUPLED_RESPONSE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'interaction_position_TV':response,
                      'particle_hole_TV':conjugation,'measurement_backaction_TV':backaction}))


if __name__=='__main__':run()
