"""Exact phase-verifier boundary tests; diagnostic, not an alternative source law."""
import argparse
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
from exact_phase import certify_source,certify_five_selector,require
from circuits import reversal,TRI

HERE=Path(__file__).resolve().parent


def record_tangent(t):
    n=int((t*(2**40)).to_integral_value(rounding=ROUND_HALF_EVEN));q=F(n,2**40)
    return {'bits':40,'numerator':n,'real':str((1-q*q)/(1+q*q)),'imag':str(2*q/(1+q*q))}


def compile_known(counts):
    c=sum(counts[0]);tr=sum(counts[i][i]for i in range(3))-c
    sq=sum(counts[i][j]*counts[j][i]for i in range(3)for j in range(3))-c*c
    d=2*sq-tr*tr
    with localcontext()as ctx:
        ctx.prec=320
        dec=lambda x:Decimal(x.numerator)/Decimal(x.denominator)
        s=(abs(dec(tr))+dec(d).sqrt())/2
        boost=s<Decimal('0.5');b=3*s-4*s**3 if boost else s
        if b==Decimal('0.5'):
            return {'bits':40,'endpoint':'minus_one','real':'-1','imag':'0'},boost
        return record_tangent((1/(4*b*b-1)).sqrt()),boost


def isotropic(s):
    return [[s*int(i==j)+(1-s)/3 for j in range(3)]for i in range(3)]


def calculate():
    rows=[]
    for m in (16,32):
        for slot in range(3*m):
            vertex,k=divmod(slot,m)
            edges=[e for e in TRI if vertex in e]*(m-k)+[e for e in TRI if (vertex+1)%3 in e]*k
            for law in ('native','swap'):
                _,counts,_=reversal(edges,4*m,law,1)
                record,boost=compile_known(counts)
                rows.append({'M':m,'slot':slot,'law':law,'phase':record,'preboost':boost,
                             'certificate':certify_source(counts,record,boost)})
    # Includes a deliberately tiny nonzero distance from half. The exact
    # certificate, unlike a finite Decimal equality test, must distinguish it.
    boundary=[]
    for label,s in [('half',F(1,2)),('above_half',F(1,2)+F(1,10**120)),
                    ('below_half',F(1,2)-F(1,10**120)),('lower_domain',F(1,5))]:
        counts=isotropic(s);record,boost=compile_known(counts)
        boundary.append({'label':label,'norm':str(s),'record':record,'preboost':boost,
                         'certificate':certify_source(counts,record,boost)})
    rejected=[]
    def rejects(label,callback):
        try:callback()
        except AssertionError as error:
            rejected.append({'label':label,'error':str(error)});return
        raise AssertionError('Accepted malformed phase certificate '+label)
    endpoint={'bits':40,'endpoint':'minus_one','real':'-1','imag':'0'}
    rejects('false_endpoint_above_half',lambda:certify_source(isotropic(F(1,2)+F(1,10**120)),endpoint,False))
    rejects('false_endpoint_below_half',lambda:certify_source(isotropic(F(1,2)-F(1,10**120)),endpoint,False))
    counts=isotropic(F(3,5));record,boost=compile_known(counts)
    wrong=dict(record);n=record['numerator']+1;t=F(n,2**40)
    wrong.update(numerator=n,real=str((1-t*t)/(1+t*t)),imag=str(2*t/(1+t*t)))
    rejects('adjacent_wrong_rounding_cell',lambda:certify_source(counts,wrong,boost))
    rejects('wrong_preboost_branch',lambda:certify_source(counts,record,True))
    rejects('nonunitary_phase',lambda:certify_source(counts,dict(record,imag='0'),boost))
    with localcontext()as ctx:
        ctx.prec=320
        selector=record_tangent((Decimal(2)/3).sqrt())
    selector_certificate=certify_five_selector(selector)
    return {'passed':True,'verification_only':True,'external_G6':False,
            'source_norm_cases':rows,'source_norm_case_count':len(rows),
            'boundary_cases':boundary,'rejected_mutations':rejected,
            'five_selector':{'record':selector,'certificate':selector_certificate},
            'interpretation':'Exact rounding integrity checks. No source-law sensitivity or new emergence credit.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record-to',type=Path)
    args=parser.parse_args()
    if args.record_to:require(not args.record_to.exists(),'First-record destination exists')
    start=time.perf_counter();result=calculate()
    if args.record_to:
        with args.record_to.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:require(result==json.loads((HERE/'PHASE_TESTS.json').read_text()),'Saved exact phase tests differ')
    print(json.dumps({'passed':True,'source_cases':result['source_norm_case_count'],
                     'rejections':len(result['rejected_mutations']),'seconds':time.perf_counter()-start}))
