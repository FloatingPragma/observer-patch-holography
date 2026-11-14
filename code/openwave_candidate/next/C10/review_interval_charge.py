"""Independent minimum-modulus and branch-projector charge decoder audit."""
from pathlib import Path
import sys,json,hashlib,itertools
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(HERE));sys.path.insert(0,str(BASE/'next/C9'))
from interval_charge_decoder import decode_charge
from fermionic_shared_records import FermionicSharedRecords,Construction


def run():
    total=3;c=Construction(3,((0,1),(1,2),(2,0)),mass=.23);results=[]
    for modulus in (4,6):
        model=FermionicSharedRecords(c,7,modulus);count=0;projectors=0
        for bits in range(64):
            region=[z for z in range(6) if bits>>z&1];gr=sum(z<3 for z in region)%modulus
            exact=[];decoded=[]
            for mask,f in model.states:
                charge=sum((mask>>z)&1 for z in region)
                # Sum oriented half-edge boundary terms directly, no producercut.
                boundary=0
                for h,(port,midpoint) in enumerate([(port,3+e) for e,ends in enumerate(c.seams) for port in ends]):
                    boundary+=f[h]*(int(port in region)-int(midpoint in region))
                candidates=[q for q in range(4) if (q+boundary-gr)%modulus==0]
                if candidates!=[charge]:raise ValueError('nonunique cut charge')
                value=decode_charge(gr,boundary%modulus,total,modulus)
                if value!=charge:raise ValueError('decoder mismatch')
                exact.append(charge);decoded.append(value);count+=1
            for q in range(4):
                if [x==q for x in exact]!=[x==q for x in decoded]:raise ValueError('projector branch mismatch')
                projectors+=1
        results.append({'modulus':modulus,'states':model.n,'all_cut_basis_checks':count,'exact_diagonal_projector_equalities':projectors})
    # General integerinterval injectivity criterion, including shifted/signedcases.
    counter=[]
    for a in range(-4,5):
        for width in range(7):
            for n in range(2,9):
                residues=[q%n for q in range(a,a+width+1)]
                if (len(set(residues))==len(residues))!=(n>width):raise ValueError('interval theorem failed')
    for args in ((0,0,3,3),(True,0,3,6),(0,0,-1,6),(0,1,0,6)):
        try:decode_charge(*args)
        except (ValueError,TypeError):counter.append({'inputs':args,'rejected':True})
        else:raise ValueError('unsafe input accepted')
    paths=[Path(__file__),HERE/'interval_charge_decoder.py',HERE/'INTERVAL_CHARGE_DECODER.json',BASE/'next/C9/fermionic_shared_records.py']
    return {'schema':'oph.c10.interval_charge_independent_review.v1','checks':{'minimal_N4_and_N6_all_projectors_equal':True,'general_shifted_interval_injectivity':True,'unsafe_inputs_rejected':True},'executed':results,'rejections':counter,
      'proof_scope':'Knownintegerinterval[a,b] decodesresidueuniquelyiffN>b−a. WithnonnegativefixedtotalL theintervalis[0,L], soN>L. Diagonalprojectorequality impliesfullbranch equalitywitharbitrarycoherences/references. No electric/spin/species dictionary.',
      'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('INTERVAL_CHARGE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
