"""Sparse local C/E instrument on a modular path ball, with explicit rejections.

The finite regulator is CPTP, with rejected native jumps separately retained.
An a priori one-interval trace-distance bound concerns its exact mathematical
map versus the complete modular sector, from the initial basis state only.
"""
from collections import deque
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from .construction import Construction
from .geometry import flat_ring
from .modular_native import LocalCoarseModular


def require(ok,message):
    if not ok:raise ValueError(message)


def upper_sqrt(value):
    scale=2**32;top=value.numerator*scale*scale;bottom=value.denominator
    root=math.isqrt(top//bottom)
    if root*root*bottom<top:root+=1
    return Fraction(root,scale)


def upper_exp(x):
    require(x>=0,'nonnegative exponent required')
    terms=max(32,math.ceil(x)+2);term=Fraction(1);total=term
    for k in range(1,terms+1):term=term*x/k;total+=term
    following=term*x/(terms+1);ratio=x/(terms+2)
    require(ratio<1,'invalid exponential tail ratio')
    return total+following/(1-ratio)


def trace_bound(construction,radius,time,probability):
    """Certified rational majorant, no finite reference or quadrature input."""
    require(type(radius) is int and radius>=0,'nonnegative radius required')
    require(math.isfinite(time) and math.isfinite(probability) and 0<=probability<=1,'invalid clock')
    degree=[0]*construction.ports
    for a,b in construction.seams:degree[a]+=1;degree[b]+=1
    mass=Fraction.from_float(float(construction.mass))
    norm=upper_sqrt(mass*mass+2*max(degree,default=0))
    x=abs(Fraction.from_float(float(time)))*norm
    exponential=upper_exp(x)
    def tail(order):return exponential*x**order/math.factorial(order)
    coherent=min(Fraction(1),2*tail(radius+1))
    # All local jumps span <=2 half-hops. The projected coherent state's
    # outside-(R-2) norm is bounded by the order-(R-1) exponential tail.
    outer=Fraction(1) if radius<2 else min(Fraction(1),tail(radius-1))
    result=min(Fraction(1),coherent+2*Fraction.from_float(float(probability))*outer)
    rounded=math.nextafter(float(result),math.inf) if result not in (0,1) else float(result)
    return {'trace_distance_upper_bound':rounded,
            'coherent_part_upper':math.nextafter(float(coherent),math.inf) if coherent else 0.,
            'commit_outer_norm_upper':math.nextafter(float(outer),math.inf) if outer else 0.,
            'method':'exact rational sqrt/exponential majorants, outward float conversion; trace distance = half trace norm',
            'scope':'One interval from initial basis state, complete modular sector at same modulus; exact mathematical map, no certified floating evolution error and no multi-step bound.'}


def outward(value,upper=True):
    value=Fraction(value)
    if value==0:return 0.
    number=float(value)
    return math.nextafter(number,math.inf if upper else -math.inf)


def _rational_taylor_action(matrix,mass,max_degree,time,tolerance=Fraction(1,10**14)):
    """Internal helper: matrix must be the generated projected incidence H_m, with matching mass and maximum port degree. Exact rational polynomial and tail."""
    require(math.isfinite(time),'finite time required')
    degree=Fraction.from_float(float(mass))**2+2*max_degree
    x=abs(Fraction.from_float(float(time)))*upper_sqrt(degree)
    exponential=upper_exp(x);order=0;tail=exponential*x
    while tail>tolerance:
        order+=1
        require(order<=128,'certified Taylor budget exceeded')
        tail=exponential*x**(order+1)/math.factorial(order+1)
    size=matrix.shape[0];real=[Fraction(0)]*size;imag=[Fraction(0)]*size
    real[0]=Fraction(1);total_real=real.copy();total_imag=imag.copy()
    coo=matrix.tocoo()
    entries=[(int(i),int(j),Fraction.from_float(float(v))) for i,j,v in zip(coo.row,coo.col,coo.data)]
    duration=Fraction.from_float(float(time))
    for step in range(1,order+1):
        hr=[Fraction(0)]*size;hi=[Fraction(0)]*size
        for i,j,value in entries:hr[i]+=value*real[j];hi[i]+=value*imag[j]
        real=[duration*value/step for value in hi]
        imag=[-duration*value/step for value in hr]
        total_real=[a+b for a,b in zip(total_real,real)]
        total_imag=[a+b for a,b in zip(total_imag,imag)]
    return total_real,total_imag,tail,order


def branch_metadata(key):
    seam,law,coin,value,status=key
    proposed=value if law=='C' else (value-(2*coin-1)*(abs(value)%2))//2
    return {'seam':seam,'law':law,'coin':coin,'proposed_record_value':value,
            'proposed_current':proposed,'applied_current':proposed if status=='native' else 0,
            'native_record':[coin,value] if status=='native' else None,'status':status}


class SparseCoarseWindow:
    def __init__(self,construction,modulus=6,radius=4,max_states=2000):
        require(type(modulus) is int and modulus>=2,'integer modulus>=2 required')
        require(type(radius) is int and radius>=0,'nonnegative integer radius required')
        require(bool(construction.seams),'event instrument requires seams')
        self.construction=construction;self.modulus=modulus;self.radius=radius
        initial=(0,(0,)*(2*len(construction.seams)))
        self.states=[initial];self.index={initial:0};self.distances=[0];queue=deque([initial])
        while queue:
            state=queue.popleft();distance=self.distances[self.index[state]]
            if distance==radius:continue
            for target,_ in self.half_hops(state):
                if target not in self.index:
                    require(len(self.states)<max_states,'path ball exceeds explicit state budget')
                    self.index[target]=len(self.states);self.states.append(target);self.distances.append(distance+1);queue.append(target)
        rows=[];cols=[];data=[]
        for col,state in enumerate(self.states):
            rows.append(col);cols.append(col);data.append(construction.mass*(1 if state[0]<construction.ports else -1))
            for target,sign in self.half_hops(state):
                if target in self.index:rows.append(self.index[target]);cols.append(col);data.append(sign)
        self.h=csr_matrix((data,(rows,cols)),shape=(len(self.states),)*2)
        require((self.h-self.h.T).nnz==0,'projected generator not Hermitian')
        groups={}
        for seam in range(len(construction.seams)):
            for law in ('C','E'):
                for coin in (0,1):
                    for source,state in enumerate(self.states):
                        target,current,contrast=self.native_target(state,seam,coin)
                        record=current if law=='C' else contrast
                        status='native' if target in self.index else 'boundary_rejected'
                        destination=self.index[target] if status=='native' else source
                        groups.setdefault((seam,law,coin,record,status),[]).append((source,destination))
        self.maps=[]
        for key,pairs in sorted(groups.items()):
            src,dst=map(np.array,zip(*pairs))
            require(len(set(dst))==len(dst),'noninjective outcome map; regulator must not merge branches')
            self.maps.append((key,src,dst))
        # Every event branch contributes (1/(4M)) times a partial permutation.
        counts=np.zeros(len(self.states),int)
        for _,src,_ in self.maps:counts[src]+=1
        require(np.all(counts==4*len(construction.seams)),'Kraus completeness failed')

    def half_hops(self,state):
        pos,f=state;c=self.construction
        for edge,(a,b) in enumerate(c.seams):
            for side,port in enumerate((a,b)):
                if pos==port:target,delta=c.ports+edge,1
                elif pos==c.ports+edge:target,delta=port,-1
                else:continue
                record=list(f);record[2*edge+side]=(record[2*edge+side]+delta)%self.modulus
                yield (target,tuple(record)),1 if side==0 else -1

    def native_target(self,state,seam,coin):
        pos,f=state;c=self.construction;a,b=c.seams[seam]
        contrast=int(pos==a)-int(pos==b)
        current=(contrast-(2*coin-1)*(abs(contrast)%2))//2
        if not current:return state,current,contrast
        record=list(f);record[2*seam]=(record[2*seam]+current)%self.modulus
        record[2*seam+1]=(record[2*seam+1]-current)%self.modulus
        return (b if current==1 else a,tuple(record)),current,contrast

    def initial(self):
        psi=np.zeros(len(self.states),complex);psi[0]=1;return psi

    def interval(self,time,probability,density=False):
        require(math.isfinite(time) and math.isfinite(probability) and 0<=probability<=1,'invalid clock')
        psi=expm_multiply(-1j*time*self.h,self.initial())
        probabilities=(1-probability)*abs(psi)**2
        rho=(1-probability)*np.outer(psi,psi.conj()) if density else None
        ledger=[{'status':'no_commit','probability':1-probability}]
        scale=math.sqrt(probability/(4*len(self.construction.seams)))
        for key,src,dst in self.maps:
            branch=np.zeros(len(psi),complex);branch[dst]=scale*psi[src]
            weight=float(np.vdot(branch,branch).real)
            probabilities+=abs(branch)**2
            if density:rho+=np.outer(branch,branch.conj())
            ledger.append({**branch_metadata(key),'probability':weight})
        require(abs(probabilities.sum()-1)<1e-12,'finite instrument trace failed')
        positions=np.zeros(self.construction.ports+len(self.construction.seams))
        for i,state in enumerate(self.states):positions[state[0]]+=probabilities[i]
        omitted={}
        for source,state in enumerate(self.states):
            for target,sign in self.half_hops(state):
                if target not in self.index:omitted[target]=omitted.get(target,0j)+sign*psi[source]
        shell=float(sum(abs(psi[i])**2 for i,distance in enumerate(self.distances) if distance>self.radius-2))
        return {'position_probabilities':positions.tolist(),
                'coherent_outer_shell_probability':shell,
                'omitted_coherent_derivative_norm':float(math.sqrt(sum(abs(value)**2 for value in omitted.values()))),
                'boundary_diagnostic_scope':'Outer-shell population and omitted-H derivative at interval end, not a finite-time leakage probability. Coherent projection error is bounded separately.',
                'boundary_rejection_probability':sum(r['probability'] for r in ledger if r['status']=='boundary_rejected'),
                'ledger_probability_sum':sum(r['probability'] for r in ledger),
                'ledger':ledger,'bound':trace_bound(self.construction,self.radius,time,probability)},rho

    def certified_interval(self,time,probability):
        """Rational polynomial CP output with explicit trace defect and intervals.

        No normalization is applied to the polynomial. For an effect E, its
        error is bounded by (epsilon*(2+epsilon)+abs(norm_defect))/2.
        The extra trace term is necessary for a subnormalized approximation.
        """
        require(math.isfinite(probability) and 0<=probability<=1,'invalid commit probability')
        degree=[0]*self.construction.ports
        for a,b in self.construction.seams:degree[a]+=1;degree[b]+=1
        real,imag,epsilon,order=_rational_taylor_action(self.h,self.construction.mass,max(degree),time)
        weights=[a*a+b*b for a,b in zip(real,imag)];norm=sum(weights)
        defect=norm-1
        numerical_effect=(epsilon*(2+epsilon)+abs(defect))/2
        regulator=trace_bound(self.construction,self.radius,time,probability)
        error=Fraction.from_float(regulator['trace_distance_upper_bound'])+numerical_effect
        p=Fraction.from_float(float(probability));factor=p/(4*len(self.construction.seams))
        positions=[Fraction(0)]*(self.construction.ports+len(self.construction.seams))
        for i,state in enumerate(self.states):positions[state[0]]+=(1-p)*weights[i]
        ledger=[{'status':'no_commit','probability':(1-p)*norm}]
        for key,src,dst in self.maps:
            weight=factor*sum(weights[int(i)] for i in src)
            for source,target in zip(src,dst):positions[self.states[int(target)][0]]+=factor*weights[int(source)]
            ledger.append({**branch_metadata(key),'probability':weight})
        require(sum(positions)==norm,'exact polynomial trace accounting failed')
        require(sum(row['probability'] for row in ledger)==norm,'exact polynomial ledger accounting failed')
        def effect_record(q):
            nearest=float(q);rounding=abs(Fraction.from_float(nearest)-q)
            return {'polynomial_probability':nearest,'conversion_abs_error_upper':outward(rounding),
                    'full_modular_probability_interval':[outward(max(Fraction(0),q-error),False),
                                                          outward(min(Fraction(1),q+error),True)]}
        output_ledger=[]
        for row in ledger:
            copy={key:value for key,value in row.items() if key!='probability'}
            copy.update(effect_record(row['probability']));output_ledger.append(copy)
        return {'taylor_order':order,'discarded_state_norm_tail_upper':outward(epsilon),
                'polynomial_norm_defect':float(defect),'polynomial_norm_defect_abs_upper':outward(abs(defect)),
                'numerical_effect_error_upper':outward(numerical_effect),
                'combined_effect_error_upper':outward(error),
                'regulator_bound':regulator,'positions':[effect_record(q) for q in positions],
                'ledger':output_ledger,
                'scope':'Exact rational polynomial with retained norm defect; certified complete-modular event intervals from initial basis for one interval. No polynomial renormalization or normalized-CPTP claim.'}


def run(include_ledger=False):
    """Small default assessment; --full-ledger retains complete certified output."""
    triangle=Construction(3,((0,1),(1,2),(2,0)),.23)
    full=LocalCoarseModular(triangle,6)
    exact,_=full.step(full.initial_density(),.37,.29)
    rows=[]
    for radius in (2,4):
        window=SparseCoarseWindow(triangle,6,radius)
        result,rho=window.interval(.37,.29,True)
        embedding=np.zeros((full.n,len(window.states)))
        for i,state in enumerate(window.states):embedding[full.index[state],i]=1
        error=float(np.abs(np.linalg.eigvalsh(embedding@rho@embedding.T-exact)).sum()/2)
        require(error<=result['bound']['trace_distance_upper_bound']+1e-12,'one-interval bound violated')
        rows.append({'radius':radius,'states':len(window.states),'actual_trace_distance':error,
                     'bound':result['bound'],'rejection_probability':result['boundary_rejection_probability']})
    carrier=SparseCoarseWindow(flat_ring(1,1,.23),6,4)
    actual,_=carrier.interval(.1,.29);certificate=carrier.certified_interval(.1,.29)
    accepted=np.zeros(len(carrier.states),int)
    for key,source,_ in carrier.maps:
        if key[-1]=='native':accepted[source]+=1
    null_key,null_source,null_target=next(item for item in carrier.maps if item[0][1]=='C' and item[0][3]==0 and item[0][-1]=='native' and len(item[1])>=2)
    coherent=np.zeros(len(carrier.states));coherent[null_source[:2]]=1/np.sqrt(2)
    branch=np.zeros(len(carrier.states));branch[null_target]=coherent[null_source]
    null_coherence=float(abs(branch[null_target[0]]*branch[null_target[1]]))
    ring={'geometry':'flat_ring(1,1)','states':len(carrier.states),'sparse_h_nonzeros':carrier.h.nnz,
          'source_port_probability':actual['position_probabilities'][0],
          'boundary_rejection_probability':actual['boundary_rejection_probability'],
          'coherent_outer_shell_probability':actual['coherent_outer_shell_probability'],
          'omitted_coherent_derivative_norm':actual['omitted_coherent_derivative_norm'],
          'regulator_trace_distance_bound':actual['bound']['trace_distance_upper_bound'],
          'certificate':{key:value for key,value in certificate.items() if key not in ('positions','ledger','regulator_bound')},
          'source_port_certified_probability':certificate['positions'][0],
          'native_outcome_branch_count':sum(key[-1]=='native' for key,_,_ in carrier.maps),
          'rejected_outcome_branch_count':sum(key[-1]=='boundary_rejected' for key,_,_ in carrier.maps)}
    if include_ledger:ring.update({'complete_window_readout':actual,'complete_certificate':certificate})
    package=Path(__file__).resolve().parent
    return {'scope':'Sparse path-window local-coarse C/E instrument on canonical ring; certified one initial-basis interval against complete same-modulus sector including classical outcomes. No multi-step, continuum or stable-particle theorem.',
            'checks':{'triangle_exact_full_bound':all(row['actual_trace_distance']<=row['bound']['trace_distance_upper_bound'] for row in rows),
                      'finite_trace':bool(abs(sum(actual['position_probabilities'])-1)<1e-12),
                      'ledger_complete':bool(abs(actual['ledger_probability_sum']-1)<1e-12),
                      'rational_tail':bool(certificate['discarded_state_norm_tail_upper']<1e-14),
                      'actual_ring_error_bound':bool(certificate['combined_effect_error_upper']<.006)},
            'mutations':{'dropping_rejection_branches_loses_trace':bool(np.any(accepted<4*len(carrier.construction.seams))),
                         'global_monitoring_changes_null_coherence':null_coherence>.4,
                         'zero_projection_error_rejected':any(row['actual_trace_distance']>1e-6 for row in rows)},
            'full_criterion_met':False,'triangle_full_comparison':rows,'actual_ring':ring,
            'supplied':{'mass':.23,'modulus':6,'radius':4,'time':.1,'commit_probability':.29,
                        'preparation':'port0 with zero modular half-seam records',
                        'instrument':'local-coarse C/E, uniform seam, fair law and coin, transit=elsewhere',
                        'window_boundary':'separate retained native-jump rejection and projected Hermitian H'},
            'source_sha256':{name:hashlib.sha256((package/name).read_bytes()).hexdigest() for name in
                             ('sparse_native.py','modular_native.py','construction.py','geometry.py')}}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-ledger',action='store_true')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=run(args.full_ledger)
    encoded=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(encoded)
    print(encoded,end='')
