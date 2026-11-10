"""Optionally generate an untrusted full-sector trial, without changing receipts."""
import argparse
from pathlib import Path
import json
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from .operator import ActualOperator,algebraic_frame


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('choose a new trial output path')
    operator=ActualOperator(algebraic_frame())
    matrix=LinearOperator((70928,70928),matvec=operator.apply,dtype=float)
    initial=np.zeros((13,5456));initial[:,0]=1/np.sqrt(13)
    values,vectors=eigsh(matrix,k=1,which='SA',v0=initial.ravel(),tol=1e-12,maxiter=1000)
    with args.output.open('xb') as stream:np.save(stream,vectors[:,0],allow_pickle=False)
    print(json.dumps({'output':str(args.output),'energy_diagnostic':float(values[0]),
        'certified':False,'next_requirement':'actual complete-source residual verification'}))


if __name__=='__main__':main()
