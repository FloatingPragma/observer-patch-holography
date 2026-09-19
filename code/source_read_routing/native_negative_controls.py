"""Run coherent malformed histories through the independent native checker."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from check_control import rows_for
from run import command, native_path
from verify import load

HERE = Path(__file__).resolve().parent


def check(binary, work):
    folder = HERE/"controls"
    packet = load(folder/"q3_baseline.json")
    original = np.asarray(list(rows_for(packet,folder)),dtype="<u8")
    work.mkdir(parents=True,exist_ok=True)
    args = command(binary,native_path(folder/"q3.input"),"baseline",native_path(work/"logical.bin"))
    positive = subprocess.run(args,input=original.tobytes(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if positive.returncode:
        raise ValueError("native positive control failed: "+positive.stderr.decode())
    cases = {}
    first = lambda op: int(np.flatnonzero(original[:,0] == op)[0])
    for name,op,column in [("mean arithmetic",5,7),("foreign owner",3,4),
                           ("writer forgery",5,5),("unsupported route",5,2),
                           ("capture law",6,0),("receiver reset payload",4,1)]:
        data=original.copy(); data[first(op),column] += np.uint64(1); cases[name]=data
    removed=first(6)+1
    data=np.delete(original,removed,axis=0)
    for column in (5,6):
        mask=(data[:,column]>removed)&(data[:,column]<np.uint64(2**64-1))
        data[mask,column] -= np.uint64(1)
    cases["coherently renumbered missing cleanup"]=data
    data=original.copy()
    committed=first(9)
    j=next(i for i in range(committed+1,len(data)) if data[i,0]==8)
    old=first(2); data[j,2]=data[old,3]; data[j,6]=old
    cases["stale immutable version"]=data
    cases["truncated final commit"]=original[:-1]
    cases["extra final event"]=np.concatenate((original,original[-1:]))
    for name,data in cases.items():
        result=subprocess.run(args,input=data.tobytes(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if result.returncode==0:
            raise ValueError("native verifier accepted "+name)
    print(f"Native verifier: positive replay and {len(cases)} semantic mutations checked.",flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--work",type=Path)
    args=parser.parse_args()
    if args.work:
        check(args.binary,args.work)
    else:
        with tempfile.TemporaryDirectory(prefix="oph-routing-negative-") as directory:
            check(args.binary,Path(directory))
