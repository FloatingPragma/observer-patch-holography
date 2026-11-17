"""Post-run exact bookkeeping of already-certified two-source success costs."""
from fractions import Fraction as F
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def endpoints(x):
    def exact(v):return (-1 if v[0] else 1)*F(int(v[1]))*F(2)**int(v[2])
    return tuple(exact(x[k]) for k in ('lo_dyadic','hi_dyadic'))

def main():
    records=json.loads((HERE/'INTERVALS.json').read_text())['rows'];rows=[];error=F(2*2*2**22+1,2**52)
    for law in ('native','swap'):
        for k in (0,6):
            for turns in (1,2):
                for direction in (-1,1):
                    group={r['charge_sign']:r for r in records if (r['clock_law'],r['k'],r['turns'],r['direction'])==(law,k,turns,direction)}
                    if set(group)!={-1,1}:raise ValueError('cost case coverage')
                    for a in (-1,1):
                        for b in (-1,1):
                            x=endpoints(group[a]['unconditional_probability_lower']);y=endpoints(group[b]['unconditional_probability_lower']);lo,hi=x[0]*y[0],x[1]*y[1]
                            if lo<=F(9,10000) or error/lo>=F('0.00002'):raise ValueError('full pair cost/error gate')
                            rows.append({'clock_law':law,'k':k,'turns':turns,'direction':direction,'charge_signs':[a,b],'pair_unconditional_probability_enclosure':[str(lo),str(hi)],'pair_perturbed_acceptance_lower':str(lo-error),'conditional_binary_error_upper':str(error/lo)})
    return {'passed':True,'verification_only':True,'pair_cost_rows':rows,'count':len(rows),'entry_bound_per_source':'1816/2^64','pair_cost_rule':'Product of both complete selected-clock probabilities, each including admission, TWO calibrations, quantum preparation, counter/decoder, frame/marker, re-encoding and fixed path flag.','joint_initial_error_counted_once':True}

if __name__=='__main__':print(json.dumps(main(),indent=2))
