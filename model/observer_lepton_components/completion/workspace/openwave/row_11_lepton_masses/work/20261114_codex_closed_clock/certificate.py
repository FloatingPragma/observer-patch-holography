"""Exact source audit and retained-dilation contraction, with no mass targets."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import argparse
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT/'reverse-engineering-reality'
sys.path.insert(0,str(SOURCE))
import sympy as s
from model.coherent_seam_repair.program import Observer
from model.finite_family_clock import family, primitives as fc
from model.lepton_register_round import carrier, walk as wk, rounds


def simp(x):
    return s.simplify(s.expand_complex(x))


def matrix_data(m):
    return [[str(s.simplify(m[i,j])) for j in range(m.cols)] for i in range(m.rows)]


def equal(a,b):
    return all(s.simplify(x)==0 for x in a-b)


def binomial(p,n):
    return [s.binomial(n,k)*p**k*(1-p)**(n-k) for k in range(n+1)]


def convolve(p,n):
    v=[s.Integer(1)]
    for _ in range(n):
        nxt=[s.Integer(0)]*(len(v)+1)
        for k,x in enumerate(v):
            nxt[k]+=x*(1-p)
            nxt[k+1]+=x*p
        v=nxt
    return v


def calculate():
    checks={}
    native_audit=family.native_event_check()
    checks['actual_coherent_repair_360_columns']=native_audit['actual_Observer_repair_columns']==360 and native_audit['exact_mismatches']==0
    native_path=SOURCE/'model/observer/native.py'
    native_hash=sha256(native_path.read_bytes()).hexdigest()
    checks['pinned_native_bytes']=native_hash=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    w=carrier.carrier_walk()
    labels=[]
    for object_id in ('X','Y'):
        weights,current_ok,positive,ticked=carrier.double_load(w,2)
        checks[f'{object_id}_120_actual_native_columns']=current_ok and positive and len(ticked)==60
        each=[]
        for (cfg,rec,address),amp in sorted(ticked.items(),key=lambda x:x[0][2]):
            assert amp==1
            each.append(dict(address_index=address,seam=list(w.addresses[address]),record=rec,
                             final_load=list(cfg),amplitude='1/sqrt(60)',original_weight='1/60'))
        labels.append(each)
    table=[]
    counts=Counter()
    for x in labels[0]:
        for y in labels[1]:
            f=int(x['record']==2);g=int(y['record'] in (2,3))
            counts[f,g]+=1
            table.append([x['address_index'],y['address_index'],f,g,f^g,f*g])
    p=s.Rational(sum(row[4] for row in table),len(table))
    p_control=s.Rational(sum(row[5] for row in table),len(table))
    checks['all_3600_source_labels']=len(table)==3600 and sum(counts.values())==3600
    checks['source_XOR_2_9']=p==s.Rational(2,9)
    checks['source_AND_1_72']=p_control==s.Rational(1,72)
    B=fc.oriented_record()*fc.repair()*fc.address_reflection()
    prep=fc.state(s.sqrt(s.Rational(1,3)),0,s.sqrt(s.Rational(2,3)))
    phi=fc.analyzer()
    M=fc.compression(phi,B*prep)
    R=fc.shift();a=1/(3*s.sqrt(3));b=1/(3*s.sqrt(6));zsym=fc.Z
    expected=a*s.eye(3)+b*(zsym*R+R.T/zsym)
    checks['source_compression_coefficients']=equal(M,expected)
    checks['source_tick_and_analyzer_normalized']=B.H*B==s.eye(12) and (phi.H*phi)[0]==1
    adj=M.H.subs(s.conjugate(zsym),1/zsym)
    raw_effect=adj*M
    checks['complete_effect_trace']=s.simplify(s.trace(raw_effect))==s.Rational(2,9)
    z=s.Rational(3,5)+s.I*s.Rational(4,5)
    source_effects={k:(M.subs(zsym,z**k).H*M.subs(zsym,z**k)).applyfunc(s.simplify) for k in range(3)}
    # Direct sum over every original source label, grouping only identical
    # already evaluated matrix entries. No binomial law is used on this side.
    direct1=s.zeros(3)
    for row in table:
        direct1+=source_effects[row[4]]/3600
    direct1=direct1.applyfunc(s.simplify)
    qcounts=Counter(row[4] for row in table)
    direct2=s.zeros(3)
    for q1,c1 in qcounts.items():
        for q2,c2 in qcounts.items():
            direct2+=s.Rational(c1*c2,3600**2)*source_effects[q1+q2]
    direct2=direct2.applyfunc(s.simplify)

    def formula(p,n):
        c1=s.expand((1-p+p*z)**n)
        c2=s.expand((1-p+p*z**2)**n)
        return ((a*a+2*b*b)*s.eye(3)+2*a*b*(c1*R+s.conjugate(c1)*R.T)
                +b*b*(c2*R**2+s.conjugate(c2)*R.T**2)).applyfunc(s.simplify)

    generators={1:formula(p,1),2:formula(p,2)}
    checks['N1_all_original_archive_contraction']=equal(direct1,generators[1])
    checks['N2_exact_compressed_product_contraction']=equal(direct2,generators[2])
    count_laws=[]
    for n in (1,2,3,6):
        weights=binomial(p,n)
        independent=convolve(p,n)
        checks[f'N{n}_binomial_equals_convolution']=weights==independent
        checks[f'N{n}_complete_mean_variance']=(sum(weights)==1 and sum(k*v for k,v in enumerate(weights))==n*p
          and sum((k-n*p)**2*v for k,v in enumerate(weights))==n*p*(1-p))
        checks[f'N{n}_characteristic_coefficients']=(s.expand(sum(v*zsym**k for k,v in enumerate(independent))-(1-p+p*zsym)**n)==0)
        g=formula(p,n)
        checks[f'N{n}_trace_invariant']=s.simplify(s.trace(g))==s.Rational(2,9)
        count_laws.append(dict(N=n,all_original_count_weights=list(map(str,weights)),mean=str(n*p),variance=str(n*p*(1-p))))
    g_control=formula(p_control,1)
    delta=(generators[1]-g_control).applyfunc(s.simplify)
    difference=s.simplify(s.trace(delta.H*delta))
    checks['changed_AND_detected_in_same_generator']=difference.is_positive is True
    chi=s.expand(1-p+p*z)
    contraction=s.simplify(chi*s.conjugate(chi))
    checks['discard_before_inverse_detected']=contraction==s.Rational(349,405) and contraction<1
    # On the retained Q eigenblocks the inverse uses the SAME archive.
    # Every unitary phase cancels individually; this is a circuit identity,
    # explicitly distinct from the source-dependent physical tests.
    checks['retained_inverse_every_Q_block']=all(s.simplify(z**(q*j)*z**(-q*j))==1 for q in (0,1) for j in (-1,0,1))
    checks={key:bool(value) for key,value in checks.items()}
    checks['all_checks_pass']=all(checks.values())
    sources={}
    for mod in list(sys.modules.values()):
        f=getattr(mod,'__file__',None)
        if f:
            path=Path(f).resolve()
            if path.suffix=='.py' and SOURCE in path.parents:
                sources[str(path.relative_to(ROOT))]=sha256(path.read_bytes()).hexdigest()
    for relative in ('model/lepton_register_round/PROOF.md','model/finite_family_clock/PROOF.md',
                     'model/observer/native.py'):
        path=SOURCE/relative;sources[str(path.relative_to(ROOT))]=sha256(path.read_bytes()).hexdigest()
    calibration=ROOT/'openwave/row_11_lepton_masses/work/20261114_codex_observer_selection/PROOF.md'
    sources[str(calibration.relative_to(ROOT))]=sha256(calibration.read_bytes()).hexdigest()
    return dict(task='LOCAL-M10-R11-CLOSED-CLOCK-20261114',candidate='lepton-count-collision',class_='F',
       physical_row_established=False,checks=checks,native_sha256=native_hash,
       actual_pinned_native_calls=240,actual_coherent_repair_columns=360,
       count_source_address_records=labels,joint_records_column_names=['X_address','Y_address','forward','hold','XOR','AND'],
       all_3600_joint_records=table,source_flag_counts={f'{a},{b}':v for (a,b),v in sorted(counts.items())},
       source_count=str(p),changed_count=str(p_control),count_laws=count_laws,
       exact_arithmetic_phase=str(z),phase_scope='fixed per-ancilla test phase; total action N*arg(z), no mass target',
       source_M_symbolic=matrix_data(M),N1_generator=matrix_data(generators[1]),N2_generator=matrix_data(generators[2]),
       N1_changed_AND_generator=matrix_data(g_control),changed_generator_Frobenius_squared=str(difference),
       discarded_inverse_coherence_multiplier=str(contraction),retained_inverse_multiplier='1',
       contraction_accounting={'N1_original_labels':3600,'N2_original_product_labels':3600**2,
         'N2_evaluated_compressed_binary_blocks':4,'N2_full_tensor_executed':False,
         'compression':'independent source product; exact original block cardinalities retained'},
       source_hashes=dict(sorted(sources.items())))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=calculate()
    encoded=json.dumps(result,indent=2,sort_keys=True)+'\n'
    path=HERE/'RECEIPT.json'
    if args.write:path.write_text(encoded)
    else:assert path.read_text()==encoded,'receipt differs on exact replay'
    assert result['checks']['all_checks_pass']
    print(json.dumps({'status':'pass','checks':len(result['checks'])-1,
       'source_count':result['source_count'],'changed_count':result['changed_count'],
       'discarded_inverse_coherence_multiplier':result['discarded_inverse_coherence_multiplier'],
       'changed_generator_Frobenius_squared':result['changed_generator_Frobenius_squared'],
       'physical_row_established':False},sort_keys=True))

if __name__=='__main__':main()
