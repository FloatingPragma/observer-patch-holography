"""An actual native-link / CAR-Higgs join and charged-addition probe.

The exact reducing control has no colored fermions, the source charged-lepton
Higgs vertex at vertex0, both full native60 links and their color loop, and
the same two-site Gauss law. Quark pair couplings are explicitly zero.
All physical states in the declared Q*=0 and Q*=1 sectors are retained.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations,product
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import sympy as s
from scipy import sparse

from .native_inputs import native_data,manifest,assert_zero,ROOT
from .source_algebra import (add,rational,norm,dot,primitive,standard_generators,
                            rank_mod_prime,create,annihilate,occupied,mode,decode,links)

HERE=Path(__file__).resolve().parent


class NativeCarrier:
    def __init__(self,repo=ROOT):
        self.repo=Path(repo).resolve();self.source=manifest(self.repo)
        self.data=native_data(self.repo)
        self.anti,self.generators,self.scalar_generators=standard_generators(self.data)
        self.records,self.reverse,self.link_blocks=links(self.data['blocks'])
        channel=self.data['channels']['charged_lepton']
        if len(channel['left'])!=2 or len(channel['right'])!=1:raise ArithmeticError('native lepton support')
        self.left=tuple(channel['left']);self.e=channel['right'][0]
        self.leptons=set((*self.left,self.e))
        self.partner={i:(r,tuple(b['indices'])) for r,b in enumerate(self.data['blocks']) for i in b['indices']}
        self.colored={r for r,b in enumerate(self.data['blocks']) if b['color_casimir']!=0}
        self.columns=[{i:[(j,int(g[j,i])) for j in range(15) if g[j,i]] for i in range(15)} for g in self.generators]
        self.rows=[{i:[(j,int(g[i,j])) for j in range(15) if g[i,j]] for i in range(15)} for g in self.generators]
        self.scalar_columns=[{i:[(j,int(h[j,i])) for j in range(2) if h[j,i]] for i in range(2)} for h in self.scalar_generators]
        self.gamma=[]
        for h,component in enumerate(channel['components']):
            for i in self.left:
                value=s.conjugate(component[i,self.e])
                if value:
                    if not value.is_Integer:raise ArithmeticError('literal real integral Higgs tensor expected')
                    self.gamma.append((h,i,int(value)))
        self.proto_cache={};self.prototype_receipts=[]

    def gauss(self,vector,k,vertex,families,reference):
        out={};cols=self.columns[k];rows=self.rows[k]
        for (h,mask,l0,l1),coefficient in vector.items():
            for old in occupied(mask):
                site,i,f,spin=decode(old,families)
                if site!=vertex:continue
                removed,sign1=annihilate(mask,old)
                for j,amount in cols[i]:
                    created=create(removed,mode(site,j,f,spin,families))
                    if created is not None:
                        new_mask,sign2=created;add(out,(h,new_mask,l0,l1),coefficient*amount*sign1*sign2)
            if vertex==0:
                if reference:add(out,(h,mask,l0,l1),-coefficient*int(self.generators[k][self.e,self.e]))
                if h>=0:
                    for new,amount in self.scalar_columns[k][h]:add(out,(new,mask,l0,l1),coefficient*amount)
            for edge,lid in enumerate((l0,l1)):
                if not lid:continue
                r,a,b=self.records[lid]
                if vertex==edge:
                    for j,amount in rows[a]:
                        new=self.reverse[r,j,b]
                        dest=(h,mask,new,l1) if edge==0 else (h,mask,l0,new)
                        add(out,dest,-coefficient*amount)
                else:
                    for j,amount in cols[b]:
                        new=self.reverse[r,a,j]
                        dest=(h,mask,new,l1) if edge==0 else (h,mask,l0,new)
                        add(out,dest,coefficient*amount)
        return out

    def cartan(self,key,k,vertex,reference):
        h,mask,l0,l1=key;g=self.generators[k];weight=0
        for m in occupied(mask):
            site,i,_,_=decode(m,1)
            if site==vertex:weight+=int(g[i,i])
        if vertex==0:
            if reference:weight-=int(g[self.e,self.e])
            if h>=0:weight+=int(self.scalar_generators[k][h,h])
        for edge,lid in enumerate((l0,l1)):
            if lid:
                _,a,b=self.records[lid];weight+=-int(g[a,a]) if vertex==edge else int(g[b,b])
        return weight

    def prototype(self,reference,e0,e1,lsite,r0,r1):
        """Entire exact Gauss kernel for one immutable occupation/link block."""
        cache_key=(reference,e0,e1,lsite,r0,r1)
        if cache_key in self.proto_cache:return self.proto_cache[cache_key]
        mask_e=sum(1<<mode(site,self.e,0,spin,1) for site,count in ((0,e0),(1,e1)) for spin in range(count))
        raw=[]
        local=[(-1,mask_e)] if lsite<0 else [(h,mask_e|(1<<mode(lsite,i,0,0,1))) for h in range(2) for i in self.left]
        raw_size=len(local)*len(self.link_blocks[r0])*len(self.link_blocks[r1])
        for (h,mask),l0,l1 in product(local,self.link_blocks[r0],self.link_blocks[r1]):
            key=(h,mask,l0,l1)
            if all(self.cartan(key,k,v,reference)==0 for k in (6,7,10,11) for v in (0,1)):raw.append(key)
        if not raw:
            self.proto_cache[cache_key]=[]
            return []
        equations={}
        for col,key in enumerate(raw):
            for k in range(12):
                for vertex in (0,1):
                    for dest,value in self.gauss({key:Fraction(1)},k,vertex,1,reference).items():
                        row_key=(k,vertex,dest)
                        row=equations.setdefault(row_key,{})
                        add(row,col,value)
        row_list=list(equations.values())
        matrix=s.SparseMatrix(len(row_list),len(raw),{(i,j):s.Rational(v.numerator,v.denominator) for i,row in enumerate(row_list) for j,v in row.items()})
        null=matrix.nullspace()
        result=[]
        for vector in null:
            work={key:rational(vector[i]) for i,key in enumerate(raw) if vector[i]}
            for previous in result:
                coefficient=dot(work,previous)/norm(previous)
                for key,value in previous.items():add(work,key,-coefficient*value)
            work=primitive(work)
            for k in range(12):
                for vertex in (0,1):
                    if self.gauss(work,k,vertex,1,reference):raise ArithmeticError('prototype Gauss residual')
            result.append(work)
        certificate=rank_mod_prime(row_list,len(raw),len(raw)-len(null))
        self.prototype_receipts.append({'reference_ec_dual':reference,'e_at0':e0,'e_at1':e1,
            'L_site':lsite,'link0_block':r0,'link1_block':r1,'raw_internal_dimension':raw_size,
            'Cartan_neutral_dimension':len(raw),'exact_Gauss_dimension':len(result),
            'all_24_native_constraints_checked_on_each_kernel_vector':True,
            'rank_certificate':certificate,
            'raw_basis':[{'norm':str(norm(v)),'entries':[[h,str(mask),l0,l1,str(c)] for (h,mask,l0,l1),c in sorted(v.items())]} for v in result]})
        self.proto_cache[cache_key]=result
        return result

    def physical_basis(self,families,reference):
        if type(families) is not int or families < 1:
            raise ValueError('a positive integer number of families is required')
        if reference not in (0,1):raise ValueError('Q*=0 or1 only')
        basis=[];labels=[]
        charge={-1:0}|{r:int(b['q6Y']) for r,b in enumerate(self.data['blocks'])}
        for nh in (0,1):
            ne=reference+nh
            for e0 in range(ne+1):
                e1=ne-e0
                for lsite in ((-1,) if nh==0 else (0,1)):
                    m0=6*e0-3*int(lsite==0)-3*nh-6*reference
                    for r0,r1 in product(charge,repeat=2):
                        if charge[r0]-charge[r1]!=m0:continue
                        prototypes=self.prototype(reference,e0,e1,lsite,r0,r1)
                        for occ0 in combinations(range(2*families),e0):
                            for occ1 in combinations(range(2*families),e1):
                                emask=0
                                for site,occ in ((0,occ0),(1,occ1)):
                                    for fs in occ:
                                        f,spin=divmod(fs,2);emask|=1<<mode(site,self.e,f,spin,families)
                                for lfs in ((-1,) if lsite<0 else range(2*families)):
                                    for kernel,proto in enumerate(prototypes):
                                        vector={}
                                        for (h,mask,l0,l1),value in proto.items():
                                            full_mask=emask
                                            if lsite>=0:
                                                li=[i for old in occupied(mask) for site,i,_,_ in (decode(old,1),) if i in self.left]
                                                if len(li)!=1:raise ArithmeticError('one L prototype')
                                                f,spin=divmod(lfs,2);full_mask|=1<<mode(lsite,li[0],f,spin,families)
                                            vector[h,full_mask,l0,l1]=value
                                        twice_sz={sum(1 if decode(m,families)[3]==0 else -1 for m in occupied(key[1])) for key in vector}
                                        if len(twice_sz)!=1:raise ArithmeticError('fixed spin weight')
                                        basis.append(vector)
                                        labels.append({'N_H':nh,'N_F':reference+2*nh,'Qstar':reference,
                                            'e_at0':list(occ0),'e_at1':list(occ1),'L_site':lsite,'L_family_spin':lfs,
                                            'link0_block':r0,'link1_block':r1,'kernel':kernel,'twice_Sz':twice_sz.pop(),
                                            'raw_norm_squared':str(norm(vector))})
        expected=6+48*families**2 if reference==0 else 16*families+64*families**3-20*families**2
        if len(basis)!=expected:raise ArithmeticError(('complete Gauss count',reference,families,len(basis),expected))
        inverse=defaultdict(list)
        for j,vector in enumerate(basis):
            for key,value in vector.items():inverse[key].append((j,value))
        gram={}
        for values in inverse.values():
            for i,a in values:
                for j,b in values:add(gram,(i,j),a*b)
        for (i,j),value in gram.items():
            if value!=(norm(basis[i]) if i==j else 0):raise ArithmeticError('exact full Gram matrix')
        return {'vectors':basis,'labels':labels,'inverse':inverse,'norms':[norm(v) for v in basis],
                'families':families,'reference':reference,'dimension':len(basis)}

    def hopping(self,vector,edge,fa,fb,families,adjoint=False):
        out={};x,y=edge,1-edge
        for (h,mask,l0,l1),value in vector.items():
            lid=(l0,l1)[edge]
            for old in occupied(mask):
                site,i,f,spin=decode(old,families)
                if i not in self.leptons:raise ArithmeticError('colored fermion leakage')
                if not adjoint and lid==0 and site==y and f==fb:
                    r,ids=self.partner[i]
                    removed,sign1=annihilate(mask,old)
                    for j in ids:
                        created=create(removed,mode(x,j,fa,spin,families))
                        if created is not None:
                            newmask,sign2=created;link=self.reverse[r,j,i]
                            key=(h,newmask,link,l1) if edge==0 else (h,newmask,l0,link)
                            add(out,key,value*sign1*sign2)
                elif adjoint and lid and site==x and f==fa:
                    r,a,b=self.records[lid]
                    if i!=a:continue
                    removed,sign1=annihilate(mask,old)
                    created=create(removed,mode(y,b,fb,spin,families))
                    if created is not None:
                        newmask,sign2=created
                        key=(h,newmask,0,l1) if edge==0 else (h,newmask,l0,0)
                        add(out,key,value*sign1*sign2)
        return out

    def higgs(self,vector,fl,fe,families,adjoint=False):
        out={}
        for (h,mask,l0,l1),value in vector.items():
            if not adjoint and h<0:
                for scalar,i,gamma in self.gamma:
                    for sl,se,epsilon in ((0,1,1),(1,0,-1)):
                        right=mode(0,self.e,fe,se,families);left=mode(0,i,fl,sl,families)
                        one=create(mask,right)
                        if one is None:continue
                        two=create(one[0],left)
                        if two is not None:add(out,(scalar,two[0],l0,l1),value*gamma*epsilon*one[1]*two[1])
            elif adjoint and h>=0:
                for scalar,i,gamma in self.gamma:
                    if scalar!=h:continue
                    for sl,se,epsilon in ((0,1,1),(1,0,-1)):
                        right=mode(0,self.e,fe,se,families);left=mode(0,i,fl,sl,families)
                        one=annihilate(mask,left)
                        if one is None:continue
                        two=annihilate(one[0],right)
                        if two is not None:add(out,(-1,two[0],l0,l1),value*gamma*epsilon*one[1]*two[1])
        return out

    def loop(self,vector,adjoint=False):
        out={}
        for (h,mask,l0,l1),value in vector.items():
            if not adjoint and l0==l1==0:
                for r in self.colored:
                    for a,b in product(self.data['blocks'][r]['indices'],repeat=2):
                        add(out,(h,mask,self.reverse[r,a,b],self.reverse[r,b,a]),value)
            elif adjoint and l0 and l1:
                r,a,b=self.records[l0];z,c,d=self.records[l1]
                if r==z and r in self.colored and a==d and b==c:add(out,(h,mask,0,0),value)
        return out

    def addition(self,vector,f,spin,families,adjoint=False):
        out={};m=mode(0,self.e,f,spin,families)
        for (h,mask,l0,l1),value in vector.items():
            result=annihilate(mask,m) if adjoint else create(mask,m)
            if result is not None:add(out,(h,result[0],l0,l1),value*result[1])
        return out


def reduce_exact(domain,codomain,apply,label):
    numerators={};leakage_checks=0
    for col,vector in enumerate(domain['vectors']):
        action=apply(vector);dots={}
        for key,value in action.items():
            for row,coefficient in codomain['inverse'].get(key,()):add(dots,row,value*coefficient)
        residual=dict(action)
        for row,numerator in dots.items():
            numerators[row,col]=numerator
            factor=numerator/codomain['norms'][row]
            for key,value in codomain['vectors'][row].items():add(residual,key,-factor*value)
        if residual:raise ArithmeticError(f'{label}: physical leakage in column{col}, {len(residual)} entries')
        leakage_checks+=1
    return numerators,leakage_checks


def normalized_csr(numerators,domain,codomain):
    rows=[];cols=[];values=[]
    for (i,j),value in numerators.items():
        rows.append(i);cols.append(j)
        values.append(float(value)/math.sqrt(float(codomain['norms'][i]*domain['norms'][j])))
    return sparse.csr_matrix((values,(rows,cols)),shape=(codomain['dimension'],domain['dimension']),dtype=complex)


def exact_adjoint(a,b):
    return a=={(j,i):value for (i,j),value in b.items()}


def diagonal_operators(carrier,sector):
    out={name:[] for name in ('N_H','N_F','N_at1','N_e_total','N_L_total',
                              'matter_scalar_charge6','reference_charge6',
                              'electric_native_trace','electric_color','twice_Sz')}
    for label,vector in zip(sector['labels'],sector['vectors'],strict=True):
        counts=set();energies=set();colors=set()
        for h,mask,l0,l1 in vector:
            counts.add(sum(decode(m,sector['families'])[0]==1 for m in occupied(mask)))
            ne=sum(decode(m,sector['families'])[1]==carrier.e for m in occupied(mask))
            nl=sum(decode(m,sector['families'])[1] in carrier.left for m in occupied(mask))
            nh=int(h>=0)
            if nl!=nh or ne!=sector['reference']+nh:
                raise ArithmeticError('complete lepton-sector charge identity')
            if 6*ne-3*nl-3*nh!=6*sector['reference']:
                raise ArithmeticError('matter plus scalar charge must equal6Qstar')
            energy=color=Fraction(0)
            for lid in (l0,l1):
                if lid:
                    r,a,b=carrier.records[lid]
                    energy+=rational(carrier.data['casimir'][a,a])
                    color+=rational(carrier.data['operators']['color_C2'][a,a])
            energies.add(energy);colors.add(color)
        if len(counts)!=1 or len(energies)!=1 or len(colors)!=1:raise ArithmeticError('diagonal observable is not scalar on invariant basis')
        out['N_H'].append(label['N_H']);out['N_F'].append(label['N_F'])
        out['N_e_total'].append(sector['reference']+label['N_H'])
        out['N_L_total'].append(label['N_H'])
        out['matter_scalar_charge6'].append(6*sector['reference'])
        out['reference_charge6'].append(-6*sector['reference'])
        out['N_at1'].append(counts.pop());out['electric_native_trace'].append(energies.pop())
        out['electric_color'].append(colors.pop());out['twice_Sz'].append(label['twice_Sz'])
    return {name:sparse.diags([float(x) for x in values],format='csr',dtype=complex) for name,values in out.items()},out


def produce_operators(repo,out,families=3):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    carrier=NativeCarrier(repo);sectors={q:carrier.physical_basis(families,q) for q in (0,1)}
    print(json.dumps({'stage':'exact_bases','families':families,'dimensions':{q:x['dimension'] for q,x in sectors.items()},'seconds':time.monotonic()-start}),flush=True)
    matrices={};exact={};leakage=0;adjoints=0
    for q,sector in sectors.items():
        diag,diag_exact=diagonal_operators(carrier,sector)
        matrices.update({f'q{q}_{k}':v for k,v in diag.items()})
        exact[f'q{q}_diagonal']={k:[str(x) for x in v] for k,v in diag_exact.items()}
        for edge in (0,1):
            for fa,fb in product(range(families),repeat=2):
                name=f'q{q}_D{edge}_{fa}_{fb}'
                a,n=reduce_exact(sector,sector,lambda v:carrier.hopping(v,edge,fa,fb,families),name);leakage+=n
                b,n=reduce_exact(sector,sector,lambda v:carrier.hopping(v,edge,fa,fb,families,True),name+'_adjoint');leakage+=n
                if not exact_adjoint(a,b):raise ArithmeticError('CAR/link adjoint '+name)
                adjoints+=1;matrices[name]=normalized_csr(a,sector,sector)
                exact[name]=[[i,j,str(value)] for (i,j),value in sorted(a.items())]
        for fl,fe in product(range(families),repeat=2):
            name=f'q{q}_W_{fl}_{fe}'
            a,n=reduce_exact(sector,sector,lambda v:carrier.higgs(v,fl,fe,families),name);leakage+=n
            b,n=reduce_exact(sector,sector,lambda v:carrier.higgs(v,fl,fe,families,True),name+'_adjoint');leakage+=n
            if not exact_adjoint(a,b):raise ArithmeticError('CAR/Higgs adjoint '+name)
            adjoints+=1;matrices[name]=normalized_csr(a,sector,sector)
            exact[name]=[[i,j,str(value)] for (i,j),value in sorted(a.items())]
        name=f'q{q}_loop'
        a,n=reduce_exact(sector,sector,carrier.loop,name);leakage+=n
        b,n=reduce_exact(sector,sector,lambda v:carrier.loop(v,True),name+'_adjoint');leakage+=n
        if not exact_adjoint(a,b):raise ArithmeticError('loop adjoint')
        adjoints+=1;matrices[name]=normalized_csr(a,sector,sector)
        exact[name]=[[i,j,str(value)] for (i,j),value in sorted(a.items())]
        print(json.dumps({'stage':'operators','Qstar':q,'leakage_columns_checked':leakage,'seconds':time.monotonic()-start}),flush=True)
    for f,spin in product(range(families),(0,1)):
        name=f'addition_{f}_{spin}'
        a,n=reduce_exact(sectors[0],sectors[1],lambda v:carrier.addition(v,f,spin,families),name);leakage+=n
        # The annihilator can leave Q*=0 only by removing a nonexistent reference;
        # every nonzero result lies in the full physical sibling by Gauss covariance.
        b,n=reduce_exact(sectors[1],sectors[0],lambda v:carrier.addition(v,f,spin,families,True),name+'_adjoint');leakage+=n
        if not exact_adjoint(a,b):raise ArithmeticError('gauge-invariant addition adjoint')
        adjoints+=1;matrices[name]=normalized_csr(a,sectors[0],sectors[1])
        exact[name]=[[i,j,str(value)] for (i,j),value in sorted(a.items())]
    arrays={}
    for name,matrix in matrices.items():
        matrix.eliminate_zeros();arrays[name+'_data']=matrix.data;arrays[name+'_indices']=matrix.indices
        arrays[name+'_indptr']=matrix.indptr;arrays[name+'_shape']=np.asarray(matrix.shape,dtype=np.int64)
    np.savez_compressed(out/'joined_operators.npz',**arrays)
    basis_packet={str(q):{'dimension':sector['dimension'],'labels':sector['labels'],
                         'raw_norms':[str(x) for x in sector['norms']],
                         'basis':[{'entries':[[h,str(mask),l0,l1,str(c)] for (h,mask,l0,l1),c in sorted(v.items())]} for v in sector['vectors']]}
                  for q,sector in sectors.items()}
    (out/'joined_basis.json').write_text(json.dumps(basis_packet,separators=(',',':'))+'\n')
    exact_packet={'normalization':'matrix_ij = numerator_ij / sqrt(raw_norm_i*raw_norm_j)',
                  'families':families,'operators':exact}
    (out/'joined_exact_operators.json').write_text(json.dumps(exact_packet,separators=(',',':'))+'\n')
    receipt={'schema':'oph.interacting.lepton_higgs_join.v1','source':carrier.source,'families':families,
             'complete_Gauss_dimensions':{str(q):x['dimension'] for q,x in sectors.items()},
             'dimension_formulas':{'q0':'6+48*f^2','q1':'16*f+64*f^3-20*f^2'},
             'ambient_global_charge_compatible_dimensions':{'q0':(1+64*families**2)*60**2,
                'q1':(4*families+2*math.comb(4*families,2)*8*families)*60**2},
             'all_24_native_Gauss_constraints_on_each_prototype':True,
             'prototype_Gauss_checks':24*sum(r['exact_Gauss_dimension'] for r in carrier.prototype_receipts),
             'exact_prototype_rank_exhaustion':carrier.prototype_receipts,
             'every_lifted_basis_Gram_checked_exactly':True,'complete_operator_leakage_columns_checked':leakage,
             'independent_literal_adjoint_checks':adjoints,'native_Higgs_tensor_Ward_checks':carrier.data['higgs_ward_checks'],
             'operator_names':sorted(matrices),'numeric_packet_sha256':hashlib.sha256((out/'joined_operators.npz').read_bytes()).hexdigest(),
             'basis_sha256':hashlib.sha256((out/'joined_basis.json').read_bytes()).hexdigest(),
             'exact_operator_sha256':hashlib.sha256((out/'joined_exact_operators.json').read_bytes()).hexdigest(),
             'scope':{'up_and_down_pair_vertices_disabled':True,'colored_fermion_occupation_exactly_zero':True,
                 'all_native60_link_states_allowed':True,'colored_gauge_loop_active':True,
                 'Higgs_at_vertex0_only':True,'scalar_occupancy_at_most_one':True,
                 'static_reference_ec_dual_present_only_in_Q1':True,'full_Gauss_space_in_each_declared_reducing_sector':True,
                 'reference_vacuum_even_and_ec_dual_odd_fermion_parity_declared':True,
                 'charged_addition_even_on_combined_matter_reference_carrier':True,
                 'full_three_channel_or_continuum_theory_claimed':False,'physical_particle_masses_claimed':False},
             }
    (out/'joined_operator_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'stage':'complete','families':families,'dimensions':receipt['complete_Gauss_dimensions'],
                      'exact_leakage_checks':leakage,'seconds':time.monotonic()-start}),flush=True)
    return receipt


def load_operators(root):
    root=Path(root);receipt=json.loads((root/'joined_operator_receipt.json').read_text())
    for name,key in (('joined_operators.npz','numeric_packet_sha256'),
                     ('joined_basis.json','basis_sha256'),
                     ('joined_exact_operators.json','exact_operator_sha256')):
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=receipt[key]:
            raise ArithmeticError('operator packet checksum: '+name)
    out={}
    with np.load(root/'joined_operators.npz',allow_pickle=False) as data:
        for name in receipt['operator_names']:
            out[name]=sparse.csr_matrix((data[name+'_data'].copy(),data[name+'_indices'].copy(),data[name+'_indptr'].copy()),shape=tuple(data[name+'_shape']))
            if not np.all(np.isfinite(out[name].data)):
                raise ArithmeticError('nonfinite operator packet: '+name)
    return out,receipt


def primary_family_matrices(path=None):
    """Evaluate the declared exact family reference and ordered-history control."""
    path=Path(path) if path is not None else HERE.parent/'dependencies'/'exact_family_operators.json'
    packet=json.loads(path.read_text())
    values={key:np.asarray(s.Matrix(packet[key]).applyfunc(s.sympify).evalf(40),dtype=complex)
            for key in ('Y','D','J')}
    if any(x.shape!=(3,3) or not np.all(np.isfinite(x)) for x in values.values()):
        raise ValueError('three finite exact family matrices are required')
    angle=math.log(5/4)
    return values['Y'],values['Y']+(1-math.cos(angle))*values['D']+math.sin(angle)*values['J']


def build_join(root=HERE,*,family0=None,family1=None,yukawa=None,
               electric=1.,magnetic=1.,hopping=1.,higgs_coupling=1.,
               fermion_energy=0.,scalar_energy=1.,reference_energy=0.,
               phase0=0.,phase1=0.,higgs_phase=0.):
    """Complete sparse H0,H1 and all charged-addition/record probes.

    Families are the produced packet's dimension. Hopping kernels must be
    Hermitian; the literal Higgs family matrix may be any complex matrix.
    The external reference energy is explicit and shifts all Q1 additions.
    """
    matrices,receipt=load_operators(root);f=receipt['families']
    fs=[np.eye(f) if x is None else np.asarray(x,complex) for x in (family0,family1,yukawa)]
    for k,x in enumerate(fs):
        if x.shape!=(f,f) or not np.all(np.isfinite(x)):raise ValueError('finite equal family matrices required')
        if k<2 and np.linalg.norm(x-x.conj().T)>1e-12:raise ValueError('hopping kernels must be Hermitian')
    controls=[electric,magnetic,hopping,higgs_coupling,fermion_energy,scalar_energy,reference_energy,phase0,phase1,higgs_phase]
    if not all(np.isscalar(x) and np.isreal(x) and np.isfinite(x) for x in controls):
        raise ValueError('finite real scalar controls required for a Hermitian Hamiltonian')
    models={}
    for q in (0,1):
        n=receipt['complete_Gauss_dimensions'][str(q)];zero=sparse.csr_matrix((n,n),dtype=complex)
        ds=[]
        for edge in (0,1):
            d=zero.copy()
            for fa,fb in product(range(f),repeat=2):d+=fs[edge][fa,fb]*matrices[f'q{q}_D{edge}_{fa}_{fb}']
            ds.append(d)
        w=zero.copy()
        for fl,fe in product(range(f),repeat=2):w+=fs[2][fl,fe]*matrices[f'q{q}_W_{fl}_{fe}']
        field=electric*matrices[f'q{q}_electric_native_trace']-magnetic*(matrices[f'q{q}_loop']+matrices[f'q{q}_loop'].getH())
        bare=fermion_energy*matrices[f'q{q}_N_F']+scalar_energy*matrices[f'q{q}_N_H']
        if q:bare+=reference_energy*sparse.eye(n,format='csr')
        c0=np.exp(-1j*phase0)*ds[0];c1=np.exp(-1j*phase1)*ds[1]
        cw=np.exp(-1j*higgs_phase)*w
        hop=-hopping*(c0+c0.getH()+c1+c1.getH())
        higgs=higgs_coupling*(cw+cw.getH());h=(field+bare+hop+higgs).tocsr()
        if not np.all(np.isfinite(h.data)):
            raise ValueError('configured Hamiltonian has nonfinite entries')
        observables={name:matrices[f'q{q}_{name}'] for name in ('N_H','N_F','N_at1','twice_Sz',
                     'N_e_total','N_L_total','matter_scalar_charge6','reference_charge6')}
        observables.update(J0=1j*hopping*(c0-c0.getH()),J1=1j*hopping*(c1-c1.getH()),
                           JH=1j*higgs_coupling*(cw.getH()-cw))
        models[q]={'H':h,'H_field':field,'H_bare':bare,'H_hop':hop,'H_Higgs':higgs,
                   'W':w,'D0':ds[0],'D1':ds[1],'observables':observables,
                   'phase_second':{'edge0':hopping*(c0+c0.getH()),'edge1':hopping*(c1+c1.getH()),
                                   'Higgs':-higgs},
                   'dimension':n}
    return {'sectors':models,'addition':{(f0,spin):matrices[f'addition_{f0}_{spin}'] for f0,spin in product(range(f),(0,1))},
            'source':receipt['source'],'scope':receipt['scope'],
            'family_matrices':{'edge0':fs[0],'edge1':fs[1],'Higgs':fs[2]},
            'parameters':dict(electric=electric,magnetic=magnetic,hopping=hopping,higgs_coupling=higgs_coupling,
                fermion_energy=fermion_energy,scalar_energy=scalar_energy,reference_energy=reference_energy,
                phase0=phase0,phase1=phase1,higgs_phase=higgs_phase),
            'reference_grading':{'reference_CAR_mode_order':'after all matter modes',
                'reference_mode_index':60*f,'Q0_reference_parity':0,'Q1_reference_parity':1,
                'complete_carrier_parity':'even',
                'probe_order':'Bdagger=c_e_dagger*d_reference_dagger; reference creates first',
                'adjoint_order':'B=d_reference*c_e; matter annihilates first',
                'reference_CAR_sign_on_both_physical_maps':1,
                'reference_grading_is_declared':True},
            'source_selected_scalar_parameters':False,'reference_energy_offset_is_declared':True,
            'common_gauge_potential_cancels_when_reference_is_shifted_consistently':True}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=ROOT)
    p.add_argument('--out',type=Path,default=HERE);p.add_argument('--families',type=int,default=3)
    a=p.parse_args();produce_operators(a.repo,a.out,a.families)
