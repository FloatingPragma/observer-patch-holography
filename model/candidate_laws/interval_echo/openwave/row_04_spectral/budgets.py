"""Verify the source ledger at small budgets without underflow in a successful branch."""

import argparse
from decimal import Decimal as D
import json
from pathlib import Path

from tools.strict_json import load
from .experiment import HERE
from .forks import summary
from .local import fuel_probability
from .source import decimal_branches


def run():
    scopes = summary()
    cases = []
    for size,edges in ((20,630),(80,2520),(320,10080)):
        for depth in (1,2):
            for cap in (1,8,16,32,96):
                for arm,p in (("matched",.5),("bare_128",fuel_probability(128,matched=False))):
                    good,fail = decimal_branches(edges,depth,cap,-1,p)
                    cases.append({"cells":size,"prefix":depth,"cap":cap,"arm":arm,
                                  "accepted":format(good,".15E"),
                                  "failures":{key:format(val,".15E") for key,val in fail.items()}})
    # An exact witness against replacing a small survival by 1 - large failure.
    value,_ = decimal_branches(630,1,1,-1,.5)
    if not abs(value/D(2)**(-635)-1)<D("1e-24"):
        raise ArithmeticError("Single-attempt survival lost its 635-bit factor")
    return {"scope":"Post-run source-denominator and complete fork controls", "forks":scopes,"budgets":cases}


def verify(expected,actual,key=""):
    if type(actual) is not type(expected):
        raise ValueError("Unexpected ledger type")
    if isinstance(actual,dict):
        if actual.keys()!=expected.keys():
            raise ValueError("Ledger key mismatch")
        for k,v in actual.items():
            verify(expected[k],v,k)
    elif isinstance(actual,list):
        for x,y in zip(expected,actual,strict=True):
            verify(x,y,key)
    elif isinstance(actual,float):
        if not abs(actual-expected)<2e-10:
            raise ValueError("Different source control")
    elif key in {"accepted","writer_before_return","commissioning_hop","last_write","route_to_parent",
                 "failed_reset","no_accepted_word"}:
        a,b = D(actual),D(expected)
        if not (a.is_finite() and b.is_finite() and a>=0 and b>=0
                and (a==b or (min(a,b)>0 and abs(a-b)/max(a,b)<D("2e-11")))):
            raise ValueError("Lost or changed finite branch")
    elif actual!=expected:
        raise ValueError("Exact budget ledger differs")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--record",type=Path)
    args=ap.parse_args()
    record=json.loads(json.dumps(run(),allow_nan=False))
    if args.record:
        with args.record.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(record,f,separators=(",",":"),allow_nan=False)
            f.write("\n")
    else:
        verify(load(HERE/"budgets.json"),record)
    print("Verified",len(record["budgets"]),"finite-budget ledgers and all eligible three-word preparations.")
