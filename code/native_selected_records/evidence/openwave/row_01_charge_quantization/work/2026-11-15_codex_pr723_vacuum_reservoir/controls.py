"""Feasible source-state sensitivity, retained timeout and vacuum boundary."""
import json
from pathlib import Path
from model import Source, graph, history, tree
from reader import observe
from verify import record_reading

HERE = Path(__file__).resolve().parent


def main():
    edges = graph(json.loads((HERE/'seams.json').read_text()), 2)
    states = {'zero':[0]*24, 'single':[1]+[0]*23,
              'double':[2]+[0]*23, 'negative_double':[-2]+[0]*23,
              'filled_first_cell':[2]*12+[0]*12}
    trials = []
    for i,(name,initial) in enumerate(states.items()):
        body,probe,outcome = history(initial,'native',72312900+i,edges,0,64)
        reading = observe(24,edges,body,probe,64,0,0)
        independent = record_reading(24,edges,body,probe,64,0,0)
        assert reading == independent and reading['initial_from_records'] == initial
        trials.append({'name':name,'initial_audit_only':initial,'body':body,
                       'probe':probe,'outcome':outcome,'reading':reading})
    source = Source(states['double'],'native')
    first_edge = tree(24,edges)[0]
    timeout_probe = [source.event(first_edge,0,probe=True) for _ in range(2)]
    timeout = observe(24,edges,[],timeout_probe,2,0,0)
    assert timeout == record_reading(24,edges,[],timeout_probe,2,0,0) == {'outcome':'inconclusive'}
    boundary = []
    for background in (0,2):
        initial = [background]*24
        body,probe,outcome = history(initial,'native',72312999,edges,64,64)
        boundary.append({'background_audit_only':background,'body':body,'probe':probe})
    assert boundary[0]['body'] == boundary[1]['body']
    assert boundary[0]['probe'] == boundary[1]['probe']
    return {'feasible_detector_trials':trials,
            'forced_timeout_probe':{'records':timeout_probe,'outcome':timeout},
            'uniform_background_boundary':boundary,
            'boundary_reason':'All-two over all24 ports violates the declared initial vacuum reference; it remains fixed and record-indistinguishable from zero without that reference.',
            'passed':True,'same_group':True}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
