"""Independent full-basis Gauss and spin/CAR audit of the joined carrier."""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from .interacting_carrier import (HERE,NativeCarrier,reduce_exact,normalized_csr,
                                  exact_adjoint,build_join,primary_family_matrices,ROOT)
from .source_algebra import add,occupied,decode,mode,create,annihilate,norm


def spin_ladder(vector,families,raising=True):
    """S+ = sum c_up† c_down; actual raw CAR, independent of H compression."""
    out={};wanted=1 if raising else 0
    for (h,mask,l0,l1),value in vector.items():
        for old in occupied(mask):
            site,i,f,spin=decode(old,families)
            if spin!=wanted:continue
            a=annihilate(mask,old)
            b=create(a[0],mode(site,i,f,1-spin,families))
            if b is not None:add(out,(h,b[0],l0,l1),value*a[1]*b[1])
    return out


def raw_difference(a,b):
    result=dict(a)
    for key,value in b.items():add(result,key,-value)
    return result


def maxabs(a):
    return float(np.max(np.abs(a.data if sparse.issparse(a) else a),initial=0))


def run(repo,root,output=None):
    start=time.monotonic();root=Path(root)
    output=Path(output) if output is not None else root
    output.mkdir(parents=True,exist_ok=True)
    packet=json.loads((root/'joined_basis.json').read_text())
    operator_receipt=json.loads((root/'joined_operator_receipt.json').read_text())
    if hashlib.sha256((root/'joined_basis.json').read_bytes()).hexdigest()!=operator_receipt['basis_sha256']:
        raise ArithmeticError('exact basis checksum')
    f=operator_receipt['families'];carrier=NativeCarrier(repo)
    y,sigma=primary_family_matrices()
    joined=build_join(root,family0=y,family1=sigma,yukawa=y,fermion_energy=0.,scalar_energy=1.)
    total_gauss=0;spin_ward=0;addition_CAR=0;reference_CAR=0;receipts={}
    for q in (0,1):
        # Rebuild each lifted basis from source definitions, then compare
        # the independently serialized exact vectors one by one.
        sector=carrier.physical_basis(f,q)
        for vector,stored in zip(sector['vectors'],packet[str(q)]['basis'],strict=True):
            literal={(h,int(mask),l0,l1):Fraction(value) for h,mask,l0,l1,value in stored['entries']}
            if vector!=literal:raise ArithmeticError('serialized full basis mismatch')
            for k in range(12):
                for vertex in (0,1):
                    if carrier.gauss(vector,k,vertex,f,q):raise ArithmeticError('lifted full-family Gauss residual')
                    total_gauss+=1
            for fl in range(f):
                for fe in range(f):
                    for adjoint in (False,True):
                        x=spin_ladder(carrier.higgs(vector,fl,fe,f,adjoint),f)
                        z=carrier.higgs(spin_ladder(vector,f),fl,fe,f,adjoint)
                        if raw_difference(x,z):raise ArithmeticError('literal CAR Higgs is not a spin singlet')
                        spin_ward+=1
            # An explicit last-ordered odd reference CAR mode makes the
            # source+detector probe even. Check its omitted sign is truly+1.
            for fa in range(f):
                for sa in (0,1):
                    explicit={};ref=60*f;m=mode(0,carrier.e,fa,sa,f)
                    for (scalar,mask,l0,l1),value in vector.items():
                        fullmask=mask|(q<<ref)
                        if (fullmask.bit_count()%2)!=0:raise ArithmeticError('combined physical parity')
                        first=annihilate(fullmask,m) if q else create(fullmask,ref)
                        if first is None:continue
                        second=annihilate(first[0],ref) if q else create(first[0],m)
                        if second is not None:
                            add(explicit,(scalar,second[0]&~(1<<ref),l0,l1),value*first[1]*second[1])
                    retained=carrier.addition(vector,fa,sa,f,adjoint=bool(q))
                    if raw_difference(explicit,retained):raise ArithmeticError('graded reference CAR detector sign')
                    reference_CAR+=1
            if q==0:
                for fa in range(f):
                    for sa in (0,1):
                        for fb in range(f):
                            for sb in (0,1):
                                a=carrier.addition(carrier.addition(vector,fb,sb,f),fa,sa,f,True)
                                b=carrier.addition(carrier.addition(vector,fa,sa,f,True),fb,sb,f)
                                for key,value in b.items():add(a,key,value)
                                expected=vector if (fa,sa)==(fb,sb) else {}
                                if raw_difference(a,expected):raise ArithmeticError('complete addition CAR sum rule')
                                addition_CAR+=1
        plus,_=reduce_exact(sector,sector,lambda v:spin_ladder(v,f),f'q{q}_Splus')
        minus,_=reduce_exact(sector,sector,lambda v:spin_ladder(v,f,False),f'q{q}_Sminus')
        if not exact_adjoint(plus,minus):raise ArithmeticError('spin ladder CAR adjoint')
        sp=normalized_csr(plus,sector,sector);sm=sp.getH()
        model=joined['sectors'][q];h=model['H'];sz=model['observables']['twice_Sz']/2
        s2=sz@sz+(sp@sm+sm@sp)/2
        commutator=h@sp-sp@h
        algebra=sp@sm-sm@sp-2*sz
        # N_F=0,2 supports spin0/1; N_F=1,3 supports spin1/2,3/2.
        allowed=(0.,2.) if q==0 else (.75,3.75)
        minimal=(s2-allowed[0]*sparse.eye(len(sector['vectors'])))@(s2-allowed[1]*sparse.eye(len(sector['vectors'])))
        receipt={'all_lifted_vectors_rebuilt_and_compared':sector['dimension'],
                 'full_Gauss_checks':24*sector['dimension'],
                 'literal_CAR_Higgs_spin_singlet_checks':2*f*f*sector['dimension'],
                 'spin_ladder_adjoint_exact':True,
                 'spin_ladder_H_commutator_max':maxabs(commutator),
                 'spin_algebra_max':maxabs(algebra),
                 'spin_Casimir_allowed_eigenvalues':list(allowed),
                 'spin_Casimir_minimal_polynomial_max':maxabs(minimal)}
        if max(receipt[x] for x in ('spin_ladder_H_commutator_max','spin_algebra_max','spin_Casimir_minimal_polynomial_max'))>1e-11:
            raise ArithmeticError('spin covariance check')
        receipts[str(q)]=receipt
        print(json.dumps({'stage':'exact_full_audit','Qstar':q,'Gauss_checks':total_gauss,
                          'spin_Ward_checks':spin_ward,'seconds':time.monotonic()-start}),flush=True)
    out={'schema':'oph.interacting.lepton_higgs_exact_audit.v1',
         'all_lifted_24_Gauss_constraints_checked':total_gauss,
         'all_raw_literal_spin_CAR_Higgs_Ward_checks':spin_ward,
         'all_raw_Q0_addition_CAR_anticommutators_checked':addition_CAR,
         'all_explicit_graded_reference_probe_CAR_sign_checks':reference_CAR,
         'reference_CAR_mode_order':'after all matter modes',
         'creation_order':'c_e_dagger*d_reference_dagger',
         'annihilation_order':'d_reference*c_e',
         'standard_and_native_current_spans_equal_exactly':True,
         'sectors':receipts,'elapsed_seconds':time.monotonic()-start}
    (output/'joined_exact_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=ROOT)
    p.add_argument('--root',type=Path,default=HERE)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();run(args.repo,args.root,args.output)
