"""Reuse modular Gauss readout with the actual nonnegative occupation interval.

Q_R lies in [0,L], so N>L suffices to recover it from G_R-cut_F modulo N.
This is narrower than decoding an arbitrary signed cut in [-L,L].
"""
from pathlib import Path
import sys,json,hashlib,operator
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'next/C9'))
from fermionic_shared_records import FermionicSharedRecords,Construction


def decode_charge(gauss_residue,cut_residue,total,modulus):
    args=(gauss_residue,cut_residue,total,modulus)
    if any(isinstance(x,bool) for x in args):raise ValueError('integer residues and bounds required')
    gauss,cut,total,modulus=map(operator.index,args)
    if total<0 or modulus<=total or modulus<2:raise ValueError('modulus must exceed the width of the admitted charge interval')
    value=(gauss-cut)%modulus
    if value>total:raise ValueError('residue is outside the admitted nonnegative charge interval')
    return value


def cut_record(model,state,region):
    _,records=state;value=0;c=model.construction
    for e,(a,b) in enumerate(c.seams):
        t=c.ports+e
        value+=(int(a in region)-int(t in region))*records[2*e]
        value+=(int(b in region)-int(t in region))*records[2*e+1]
    return value%model.modulus


def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    model=FermionicSharedRecords(c,7,6);g=model.gauss(model.states[0]);checked=0
    for bits in range(64):
        region={i for i in range(6) if bits>>i&1};gr=sum(g[i] for i in region)%6
        for state in model.states:
            true=sum((state[0]>>i)&1 for i in region)
            decoded=decode_charge(gr,cut_record(model,state,region),3,6)
            if decoded!=true:raise ValueError('boundary apparatus differs from full occupation charge')
            checked+=1
    aliased=FermionicSharedRecords(c,7,3);region={0,1,2};pairs={}
    for state in aliased.states:
        q=sum((state[0]>>i)&1 for i in region);r=cut_record(aliased,state,region)
        pairs.setdefault(r,set()).add(q)
    witness=next(({'cut_residue':r,'charges':sorted(q)} for r,q in pairs.items() if 0 in q and 3 in q),None)
    if witness is None:raise ValueError('expected zero/three alias absent')
    rejected=False
    try:decode_charge(0,0,3,3)
    except ValueError:rejected=True
    files=[Path(__file__),BASE/'next/C9/fermionic_shared_records.py',BASE/'candidate/oph_openwave/construction.py']
    return {'schema':'oph.c10.interval_charge_decoder.v1','checks':{'all_complete_regional_charge_projectors_equal':checked==120*64,'Z3_zero_three_alias_retained':witness is not None,'unsafe_interval_decoder_rejected':rejected},
      'reuse':'Existing modular Gauss/cut identity and complete regional occupation PVM; reduce residue on known [0,L] interval instead of arbitrary [-L,L] cut interval.',
      'proof':'If two charges q,q_prime in [0,L] have the same residue and N>L, their difference is divisible by N and has absolute value at most L<N, hence it is zero. Therefore Q_R=(G_R-cut_F) mod N whenever this representative lies in [0,L]. Equality on every basis state gives equality of every diagonal charge-projector branch, for any density and external reference.',
      'parameters':{'fixed_nonnegative_total':3,'successful_modulus':6,'alias_control_modulus':3,'regions':64,'states':120},
      'projector_basis_comparisons':checked,'alias_control':witness,
      'scope':'Fixed nonnegative occupation domain and known regional Gauss residue. No arbitrary signed charge, vacuum-subtracted quantum regional charge, all-N C10mod3 boundary recovery or physical electric-charge dictionary. N>2L remains a valid sufficient condition for the separate generic signed-cut decoder; it is not necessary for this known charge interval.',
      'full_criterion_met':False,'official_rows_changed':[],
      'source_sha256':{str(p.resolve().relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}

if __name__=='__main__':
    r=run();Path(__file__).with_name('INTERVAL_CHARGE_DECODER.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
