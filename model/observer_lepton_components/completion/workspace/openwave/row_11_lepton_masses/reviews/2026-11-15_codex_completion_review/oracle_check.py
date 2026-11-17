#!/usr/bin/env python3
"""Exact independent oracle reconstruction, plus exhaustive reader input classes.

Primary arithmetic/source-factor checks import no producer implementation.
Only the final reader adapter imports outcome.py/headline.py, as explicitly
reported. All values tested there have first been independently reconstructed.
"""
from fractions import Fraction as F
from itertools import product, combinations_with_replacement
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time


# Independently multiply as complex pairs over the real quadratic field.
class Number(tuple):
    def __new__(cls, values=(0, 0, 0, 0)):
        return tuple.__new__(cls, tuple(map(F, values)))

    @staticmethod
    def lift(value):
        return value if isinstance(value, Number) else Number((value, 0, 0, 0))

    def __add__(self, other):
        other = Number.lift(other)
        return Number(a+b for a, b in zip(self, other))

    __radd__ = __add__

    def __neg__(self):
        return Number(-a for a in self)

    def __sub__(self, other):
        return self + -Number.lift(other)

    def __rsub__(self, other):
        return Number.lift(other) + -self

    def __mul__(self, other):
        other = Number.lift(other)
        a,b,c,d = self
        e,f,g,h = other
        return Number((a*e+2*b*f-c*g-2*d*h,
                       a*f+b*e-c*h-d*g,
                       a*g+2*b*h+c*e+2*d*f,
                       a*h+b*g+c*f+d*e))

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        return Number(x/F(divisor) for x in self)

    def conjugate(self):
        a,b,c,d = self
        return Number((a,b,-c,-d))

    def __bool__(self):
        return any(x != 0 for x in self)


ZERO = Number()
ONE = Number((1,0,0,0))
Z = Number((F(3,5),0,F(4,5),0))


def sparse(rows):
    return {int(i): Number(v) for i,v in rows if any(F(x) for x in v)}


def dot(a,b):
    return sum((x.conjugate()*b.get(i,ZERO) for i,x in a.items()), ZERO)


def accumulate(out, i, value):
    value = out.get(i,ZERO)+value
    if value:
        out[i] = value
    elif i in out:
        del out[i]


def linear(columns, vector):
    out = {}
    for i,a in vector.items():
        for j,u in columns[i].items():
            accumulate(out,j,u*a)
    return out


def adjoint(columns):
    rows = [{} for _ in columns]
    for i,column in enumerate(columns):
        for j,u in column.items():
            rows[j][i] = u.conjugate()
    return rows


def family_checks(receipt):
    data = receipt['retained_family_query']
    B = [[F(x) for x in row] for row in data['B12']]
    assert len(B) == 12 and all(len(row) == 12 for row in B)
    assert all(sum(B[k][i]*B[k][j] for k in range(12)) == (i == j)
               for i in range(12) for j in range(12))
    prep, phi = sparse(data['prep_numerator']), sparse(data['phi_numerator'])
    assert dot(prep,prep) == Number.lift(36) and dot(phi,phi) == Number.lift(24)
    preps = [{4*i+f:a for i,a in prep.items()} for f in range(4)]
    all_y = {}
    effects = {}
    inverse_columns = projector_columns = 0
    for b in (0,1):
        # Explicit tensor matrix from the definition, rather than the
        # producer's sequence of factor()/shift()/phase() vector routines.
        U = []
        for jc,je,f in product(range(12),range(12),range(4)):
            column = {}
            for ic in range(12):
                if not B[ic][jc]:
                    continue
                signed_record = (1 if ic%2 == 0 else -1) if (ic//2)%2 else 0
                target_f = (f+signed_record)%3 if f < 3 else 3
                phase = Z if b and signed_record == 1 else Z.conjugate() if b and signed_record == -1 else ONE
                for ie in range(12):
                    if B[ie][je]:
                        column[4*(12*ic+ie)+target_f] = phase*(B[ic][jc]*B[ie][je])
            U.append(column)
        Ud = adjoint(U)
        for i,col in enumerate(U):
            assert linear(Ud,col) == {i:ONE}
            inverse_columns += 1
        vectors = [linear(Ud,{4*i+f:a for i,a in phi.items()}) for f in range(4)]
        assert vectors == [sparse(v) for v in data['projector_rank4_numerators'][str(b)]]
        assert all(dot(vectors[i],vectors[j]) == Number.lift(24*(i==j))
                   for i in range(4) for j in range(4))

        def projection(v):
            out = {}
            for source in vectors:
                coefficient = dot(source,v)/24
                for j,a in source.items():
                    if coefficient:
                        accumulate(out,j,a*coefficient)
            return out

        # Independently apply U^* Pi U to every source column, retaining the
        # whole576 domain. No 3x3 contraction is used to define the inverse.
        for i,col in enumerate(U):
            pi_col = {}
            for f in range(4):
                coefficient = sum((a.conjugate()*col.get(4*obj+f,ZERO) for obj,a in phi.items()),ZERO)/24
                for obj,a in phi.items():
                    if coefficient:
                        pi_col[4*obj+f] = a*coefficient
            assert linear(Ud,pi_col) == projection({i:ONE})
            projector_columns += 1
        ys = [projection(v) for v in preps]
        assert ys == [sparse(v) for v in data['query0_family_factors'][str(b)]]
        all_y[b] = ys
        A = [[dot(preps[i],ys[j])/36 for j in range(4)] for i in range(4)]
        effects[b] = A
        assert all(dot(ys[i],ys[j])/36 == A[i][j] for i in range(4) for j in range(4))
        for f in range(4):
            complement = dict(preps[f])
            for i,a in ys[f].items():
                accumulate(complement,i,-a)
            assert not projection(complement)
            assert dot(ys[f],complement) == ZERO
            assert dot(ys[f],ys[f])+dot(complement,complement) == Number.lift(36)
        null = sparse(data['null_plants'][str(b)])
        assert null and projection(null) == {}
        z = Z if b else ONE
        off = Number((0,1,0,0))*z/27+z.conjugate()*z.conjugate()/54
        for j in range(3):
            assert A[j][j] == Number.lift(F(2,27))
            assert A[(j+1)%3][j] == off and A[(j-1)%3][j] == off.conjugate()
        # W decoding restores processor3 on every prepared physical input:
        # these logical source factors never populate processor3 from f<3.
        assert all(i%4 < 3 for v in ys[:3] for i in v)
    rays = {}
    for name,p in (('native_xor',F(2,9)),('changed_AND',F(1,72)),('count_HOLD',F(0))):
        A = [[(1-p)*effects[0][i][j]+p*effects[1][i][j] for j in range(3)] for i in range(3)]
        stored = data['effects_and_controls'][name]
        assert A == [[Number(x) for x in row] for row in stored['compression']]
        ray = sum((x for row in A for x in row),ZERO)/3
        assert ray == Number(stored['query0_invariant_family_ray'])
        assert ONE-ray == Number(stored['query1_invariant_family_ray'])
        rays[name] = list(map(str,ray))
        assert sum((A[i][i] for i in range(3)),ZERO) == Number.lift(F(2,9))
    assert len({tuple(v) for v in rays.values()}) == 3
    return prep,all_y,{'full_inverse_columns':inverse_columns,'full_projector_columns':projector_columns,
                       'rank_per_count_block':4,'complement_dimension_per_block':572,
                       'control_ray_probabilities':rays}


def count_checks(receipt):
    data = receipt['retained_count_source']
    C = [[F(x) for x in row] for row in data['local24_completion']]
    R = [[F(x) for x in row] for row in data['R60']]
    assert all(C[i][j] == C[j][i] for i in range(24) for j in range(24))
    assert all(sum(C[k][i]*C[k][j] for k in range(24)) == (i==j)
               for i in range(24) for j in range(24))
    assert R == [[F(1,30)-(i==j) for j in range(60)] for i in range(60)]
    assert all(sum(R[k][i]*R[k][j] for k in range(60)) == (i==j)
               for i in range(60) for j in range(60))
    configurations=[]
    for i,j in combinations_with_replacement(range(12),2):
        cfg=tuple((k==i)+(k==j) for k in range(12))
        configurations.append(cfg)
    assert len(set(configurations)) == 78
    addresses = data['addresses']
    assert len(addresses) == len({tuple(a) for a in addresses}) == 60
    local = ((0,0),(1,0),(0,1),(1,1),(2,0),(0,2))
    allowed = set(configurations)

    def repair(cfg,record,address):
        a,b=address
        column = 4*local.index((cfg[a],cfg[b]))+record
        out={}
        for row in range(24):
            if C[row][column]:
                source,raw=divmod(row,4)
                dst=list(cfg)
                dst[a],dst[b]=local[source]
                dst=tuple(dst)
                assert dst in allowed
                out[(dst,raw)] = C[row][column]
        return out

    n=0
    for cfg,raw,address in product(configurations,range(4),addresses):
        forward=repair(cfg,raw,address)
        back={}
        for (next_cfg,next_raw),coefficient in forward.items():
            for key,value in repair(next_cfg,next_raw,address).items():
                back[key]=back.get(key,F(0))+coefficient*value
        back={key:value for key,value in back.items() if value}
        assert back == {(cfg,raw):F(1)}
        n+=1
    assert n == data['full_count_source_dimension_per_cell'] == 18720
    blank=tuple(2*(i==2) for i in range(12))
    records=[]
    for address in addresses:
        out=repair(blank,0,address)
        assert len(out)==1 and list(out.values())==[1]
        records.append(next(iter(out))[1])
    flags=([int(r==2) for r in records],[int(r in (2,3)) for r in records])
    assert list(flags)==data['count_flags']
    for side,flag in enumerate(flags):
        matrix=[[sum(R[i][k]*flag[k]*R[k][j] for k in range(60)) for i in range(60)] for j in range(60)]
        assert matrix==[[F(x) for x in row] for row in data['one_cell_projector_columns'][side]]
    for ix,iy in product(range(60),repeat=2):
        fx,fy=flags[0][ix],flags[1][iy]
        assert data['joint_labels'][60*ix+iy]==[ix,iy,records[ix],records[iy],fx,fy,fx^fy,fx*fy]
    for name,p in (('xor',F(2,9)),('and',F(1,72))):
        stored=data['retained_count_vectors'][name]
        q=[]
        for ix,iy in product(range(60),repeat=2):
            x=F(sum(flags[0]),30)-flags[0][ix]
            y=F(sum(flags[1]),30)-flags[1][iy]
            q.append(x+y-2*x*y if name=='xor' else x*y)
        assert q==list(map(F,stored['q1_numerator']))
        assert [1-x for x in q]==list(map(F,stored['q0_numerator']))
        assert sum(q)/3600==sum(x*x for x in q)/3600==p
        assert sum(x*(1-x) for x in q)==0
    assert F(7,9)**2+F(2,9)**2==F(53,81)
    return {'complete_count_repair_columns':n,'count_original_labels':3600,
            'retained_count_return':'1','committed_count_return':'53/81',
            'committed_count_nonreturn':'28/81','full_tensor_formula_retained':True}


def arm_checks(receipt):
    data=receipt['whole_arm_native_query']
    states=data['input_states'];work=data['working_states']
    assert len(states)==60 and len(work)==78
    assert len({json.dumps(s,sort_keys=True) for s in work})==78
    total=queries=0
    for law in ('C','E'):
        pairs={}
        for source,flavors,destination,raw,sign in data['signed_pairs'][law]:
            a=(source,tuple(flavors),0);b=(destination,tuple(flavors),raw)
            assert a not in pairs and b not in pairs and sign in (-1,1)
            pairs[a]=(b,sign);pairs[b]=(a,sign)
        def Q(state):
            return pairs.get(state,(state,1))
        for source,fa,fb,raw,query in product(range(78),range(6),range(6),range(16),(0,1)):
            key=(source,(fa,fb),raw)
            mid,sign=Q(key);back,sign2=Q(mid)
            assert back==key and sign*sign2==1
            total+=1
        for row in data['query_rows'][law]:
            source,flavors,query,stored_mid,stored_copy,stored_out,stored_sign=row
            mid,sign=Q((source,tuple(flavors),0))
            bit=query ^ int(mid[2]%4 in (1,2))
            back,sign2=Q(mid)
            out=[back[0],back[2],bit^1]
            assert [mid[0],mid[2],query]==stored_mid
            assert [mid[0],mid[2],bit]==stored_copy
            assert out==stored_out==[source,0,query ^ int(states[source]['tagged_port']!=7)]
            assert sign*sign2==stored_sign==1
            queries+=1
        for row in data['HOLD_query_rows'][law]:
            source,flavors,query,out,sign=row
            assert out==[source,0,query^1] and sign==1
    assert total==179712 and queries==8640
    transfer={tuple(r['input']):tuple(r['output']) for r in data['neutral_transfer']}
    assert len(transfer)==24 and len(set(transfer.values()))==24
    assert all(transfer[v]==k for k,v in transfer.items())
    assert all(transfer[(1,f,3)]==(1,0,f) and transfer[(0,f,3)]==(0,f,3) for f in range(3))
    assert len(data['invalid_witnesses'])==5
    assert all(w['actual_query']=='invalid_identity' for w in data['invalid_witnesses'])
    overflow=next(w for w in data['invalid_witnesses'] if w['name']=='outward_auxiliary_boundary')
    assert overflow['hypothetical_native']['overflow']
    assert overflow['hypothetical_native']['native_current']==1
    assert overflow['hypothetical_native']['realized_signed_current']==0
    return {'signed_completion_columns':total,'whole_arm_query_columns':queries,
            'W_columns':24,'invalid_witnesses_retained':5,
            'invalid_domain': 'Identity invalid outcome outside admitted chart; not a full-Fock native-unitary emergence claim.'}


def reader_checks(receipt,packet,prep,ys):
    # Exact complete equivalence classes for all scalar values the reader uses.
    count=receipt['retained_count_source']['retained_count_vectors']['xor']
    count_reps={}
    for j,(a,b) in enumerate(zip(count['q0_numerator'],count['q1_numerator'])):
        count_reps.setdefault((F(a),F(b)),divmod(j,60))
    family_reps={}
    for obj,fi,fo in product(range(144),range(3),range(3)):
        initial=prep.get(obj,ZERO) if fi==fo else ZERO
        key=(initial,ys[0][fi].get(4*obj+fo,ZERO),ys[1][fi].get(4*obj+fo,ZERO))
        family_reps.setdefault(key,(obj,fi,fo))
    arms={}
    for i,state in enumerate(receipt['whole_arm_native_query']['input_states']):
        arms.setdefault(state['tagged_port']==7,i)
    assert set(arms)=={False,True}
    # Import exclusively for the adapter comparison; main scientific
    # calculations above did not use these arithmetic or producer routines.
    sys.path.insert(0,str(packet))
    spec=importlib.util.spec_from_file_location('reviewed_outcome',packet/'outcome.py')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls=0
    for (h0,h1),(cx,cy) in count_reps.items():
        for (initial,y0,y1),(obj,fi,fo) in family_reps.items():
            projected=h0*y0+h1*y1
            for reference,source in arms.items():
                for qi,qo in product((0,1),repeat=2):
                    if reference:
                        amp=initial/360 if qo==1-qi else ZERO
                    else:
                        blocks=((projected,initial-projected),(initial-projected,-projected))
                        amp=blocks[qi][qo]/360
                    got=module.read(receipt,source,cx,cy,obj,fi,fo,qi,qo)
                    assert Number(got['amplitude_coefficients'])==amp
                    assert Number(got['original_probability_coefficients'])==amp.conjugate()*amp
                    assert got['physical_source_base_unchanged']==source and got['physical_family_output']==fo
                    assert got['neutral_processor_output']==3 and got['current_raw_output']==got['chi_work_output']==0
                    calls+=1
    assert calls==len(count_reps)*len(family_reps)*2*4
    return {'count_equivalence_classes':len(count_reps),'family_equivalence_classes':len(family_reps),
            'reader_equivalence_cases':calls,'original_coordinate_domain':60*3600*144*3*3*2*2,
            'coverage_argument':'All3600 count labels and all1296 family-coordinate tuples are mapped to their exact scalar input classes; both source arms and all four query-bit pairs are exercised. Source index metadata is checked on each representative and the remaining source behavior is the independently checked arm table.',
            'producer_import_scope':'Only final outcome adapter imports headline.K; primary checks above import no producer functions.'}


def main():
    assert not sys.flags.optimize, 'assertions must be enabled'
    ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,required=True);args=ap.parse_args()
    start=time.perf_counter()
    receipt=json.loads((args.packet/'RECEIPT.json').read_text())
    arm=arm_checks(receipt)
    count=count_checks(receipt)
    prep,ys,family=family_checks(receipt)
    reader=reader_checks(receipt,args.packet,prep,ys)
    result={'schema':'row11.retained_oracle_completion.v1','all_pass':True,
            'arm':arm,'count':count,'family':family,'reader':reader,
            'new_native_calls':0,'source_dependency':'Original frozen native execution is checked separately by isolated byte replay; independent calculations use its bound source matrices and full records.',
            'scope':'Complete declared logical construction and prepared-output reader; no coherent-lift emergence, physical clock, phase selection or mass ratios.',
            'seconds':time.perf_counter()-start}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
