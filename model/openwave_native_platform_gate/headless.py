"""Fresh native instrument computation, with an independent endpoint oracle."""
from collections import defaultdict
from fractions import Fraction as F
import json
from model.openwave_native_process.source import unread_column,raw_local_branches

def independent_column(left,right):
    result=defaultdict(F)
    for law in ('C','E'):
        for coin in (0,1):
            outputs=[];labels=[]
            for x in (left,right):
                loads=(x&1,(x>>1)&1)
                total=sum(loads)
                left_output=total//2+coin*(total%2)
                y=(left_output,total-left_output);j=loads[0]-left_output
                outputs.append(y[0]+2*y[1]);labels.append(j if law=='C'else loads[0]-loads[1])
            if labels[0]==labels[1]:result[tuple(outputs)]+=F(1,4)
    return dict(result)

def run():
    rows=[]
    for i in range(4):
        for j in range(4):
            actual=unread_column(i,j,(0,1));expected=independent_column(i,j)
            assert actual==expected,(i,j,actual,expected)
            rows.append({'input_matrix_unit':[i,j],
                'output':[[*key,str(v)]for key,v in sorted(actual.items())]})
    branch_maps=raw_local_branches()
    errors=[]
    for b in branch_maps:
        for x,y in b.pairs:
            if b.weight:errors.append(abs(x.bit_count()-y.bit_count()))
    assert max(errors)==0
    altered=[(1,0)] # Deliberate destruction of an occupied endpoint.
    bad=max(abs(x.bit_count()-y.bit_count())for x,y in altered)
    assert bad==1
    return {'schema':'oph.native_headless.v1','verdict':'PASS_LITERAL_NATIVE_CHANNEL_AND_INDEPENDENT_ENDPOINT_ORACLE',
        'complete_raw_slots':len(branch_maps),'all_endpoint_matrix_units':rows,
        'maximum_load_conservation_error':max(errors),'altered_rule_negative_control_error':bad,
        'empirical_values_used':[],'stored_result_substitution':False,
        'criterion_observable_identified_with_Nature':False}

if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
