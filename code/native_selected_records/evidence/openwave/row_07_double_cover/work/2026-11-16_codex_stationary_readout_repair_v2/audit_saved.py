"""Read-only semantic checks on saved laws, in addition to numerical replay.

Probability nonnegativity and exact control availability are scientific
constraints, not fields to which a generic absolute roundoff tolerance applies.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import confirm

ROOT=Path(__file__).resolve().parent


def require(value,message):
    if not value:raise RuntimeError(message)


def probabilities(values,name):
    require(isinstance(values,list)and values,name+' missing probability law')
    require(all(type(v)is float and math.isfinite(v)and v>=0 for v in values),name+' negative or invalid probability')
    require(abs(sum(values)-1)<2e-10,name+' probability sum')


def case(record):
    require(type(record['selected_clock_available'])is bool,'Clock availability must be boolean')
    for stage in record['target_first_failures']:
        if stage['kind']=='declared_filter_hold':continue
        probabilities(stage['probabilities'],'target flag')
        require(stage['first_failure_mass']>=0,'Negative target first-failure mass')
    for marker in record['marker_records']:
        require(marker['first_failure_mass']>=0,'Negative marker first-failure mass')
        probabilities([float(r['probability'])for r in marker['histories']],'marker history')
    require(abs(record['raw_complete_marker_target_mass']-1)<1e-8,'Stopped-path total mass')
    returned=record['returned']
    probabilities(returned['positive_reencode_probabilities'],'returned counter')
    probabilities(returned['path_probabilities_plus_minus'],'physical path')
    visibility=returned['path_probabilities_plus_minus'][0]-returned['path_probabilities_plus_minus'][1]
    require(abs(visibility-record['visibility'])<1e-12,'Reader visibility disagrees with both-port law')
    require(record['selected_path']==('minus'if record['turns']%2 else'plus'),'Law-dependent continuation setting')
    confirm.exact_path_semantics(record)
    for branch in returned['conditional_clocks']:
        require(type(branch['available'])is bool,'Conditional clock availability must be boolean')
        for setting in('ramsey_stride1','ramsey_stride8'):
            for row in branch[setting]:
                for key in('x_plus','y_plus'):
                    require(-2e-12<=row[key]<=1+2e-12,'Ramsey floating center outside allowed roundoff enclosure')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=Path)
    args=parser.parse_args()
    if args.case:
        opener=gzip.open if args.case.suffix=='.gz'else open
        with opener(args.case,'rt')as stream:case(json.load(stream))
        print(json.dumps({'case_semantics':True}));return
    folder=ROOT/'CONFIRMATION';confirm.verify_freeze();confirm.verify_custody(folder)
    with gzip.open(folder/'PREPARATION.json.gz','rt')as stream:preparation=json.load(stream)
    for group in preparation.values():
        for source in group.values():
            for stage in source['stages']:probabilities(stage['probabilities'],'preparation flag')
            require(abs(source['complete_outcome_mass']-1)<1e-9,'Preparation stopped-path mass')
    with gzip.open(folder/'SUMMARY.json.gz','rt')as stream:summary=json.load(stream)
    confirm.scientific_gates(preparation,summary['rows'])
    count=0
    for item in json.loads((folder/'INDEX.json').read_text())['files']:
        with gzip.open(folder/item['path'],'rt')as stream:record=json.load(stream)
        if 'target_first_failures'in record:case(record);count+=1
        if 'source'in record and'counter_probabilities'in record['source']:
            probabilities(record['source']['counter_probabilities'],'initial counter')
            probabilities(record['source']['decode_probabilities'],'initial decoder')
    require(count==256,'Saved case coverage')
    print(json.dumps({'saved_semantics_pass':True,'cases':count,'read_only':True}))


if __name__=='__main__':main()
