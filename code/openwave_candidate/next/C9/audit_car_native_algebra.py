"""Independent finite combinatorial proof controls, no producer import."""
import itertools,json,cmath,hashlib
from pathlib import Path


def phase(occupation,a,b):
    return (-1)**sum(occupation[min(a,b)+1:max(a,b)])


def gauss(occupation,records,seams):
    ports=3;g=list(occupation)
    for e,(a,b) in enumerate(seams):
        g[a]+=records[2*e];g[b]+=records[2*e+1]
        g[ports+e]-=records[2*e]+records[2*e+1]
    return g


def run():
    seams=((0,1),(1,2),(2,0));modulus=3;checks=0;wrap=None;signs=set()
    for total in (1,3):
        words=[u for u in itertools.product((0,1),repeat=6) if sum(u)==total]
        for seam,(a,b) in enumerate(seams):
            groups={}
            for u in words:
                for x,y in itertools.product(range(modulus),repeat=2):
                    f=[0]*6;f[2*seam:2*seam+2]=[x,y];d=u[a]-u[b]
                    for coin in (0,1):
                        j=(d-(2*coin-1)*abs(d))//2
                        out=list(u);out[a]-=j;out[b]+=j
                        g=f.copy();g[2*seam]=(g[2*seam]+j)%modulus;g[2*seam+1]=(g[2*seam+1]-j)%modulus
                        before=gauss(u,f,seams);after=gauss(out,g,seams)
                        if any((r-s)%modulus for r,s in zip(before,after)):raise ValueError('Gauss residue changed')
                        sign=phase(u,a,b) if j else 1;signs.add(sign)
                        if total==1 and sign!=1:raise ValueError('one-load phase mismatch')
                        for law in ('C','E'):
                            key=(law,coin,j if law=='C' else d)
                            target=(tuple(out),tuple(g))
                            group=groups.setdefault(key,{})
                            if target in group and group[target]!=(u,tuple(f)):raise ValueError('noninjective outcome')
                            group[target]=(u,tuple(f));checks+=1
                        if before!=after and wrap is None:
                            angle=.2
                            wrap={'before_integer_gauss':before,'after_integer_gauss':after,
                                  'continuous_representative_phase_difference':abs(cmath.exp(1j*angle*before[a])-cmath.exp(1j*angle*after[a])),
                                  'discrete_phase_difference':abs(cmath.exp(2j*cmath.pi*before[a]/modulus)-cmath.exp(2j*cmath.pi*after[a]/modulus))}
            # Every basis input participates once for each law/coin; each amplitude1/2.
            source_counts={}
            for group in groups.values():
                for source in group.values():source_counts[source]=source_counts.get(source,0)+1
            if set(source_counts.values())!={4}:raise ValueError('incomplete instrument')
    if signs!={-1,1}:raise ValueError('fermionic negative sign control not exercised')
    return {'schema':'oph.c9.car_native_algebra_review.v1','checked_outcome_arrows':checks,
        'checks':{'each_group_partial_bijection':True,'four_quarter_weights_complete':True,
                  'fixed_number_and_Pauli_preserved':True,'modular_Gauss_preserved':True,
                  'one_load_JW_phase_is_identity':True,'nontrivial_fermion_sign_exercised':True,
                  'continuous_integer_representative_gauge_rejected':wrap['continuous_representative_phase_difference']>.1},
        'wrap_control':wrap,'scope':'Exact local algebra for arbitrary unmodified spectator bits/records; no physical species, sea preservation or continuumU1 claim.',
        'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(Path(__file__).resolve().relative_to(Path(__file__).resolve().parents[2])):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}

if __name__=='__main__':
    r=run();Path(__file__).with_name('CAR_NATIVE_ALGEBRA_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
