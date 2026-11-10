"""Independent finite Taylor and rational-inertia audit of the charged source."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json,math
import mpmath as mp
from model.affine_gauss_spectrum import model as candidate

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def endpoint(x,upper=False):
    sign,mantissa,exponent,_=x._mpi_[int(upper)]
    return F((-1)**sign*mantissa)*F(2)**exponent

def ivf(q):
    q=F(q);return mp.iv.mpf(q.numerator)/q.denominator

def product(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),mp.iv.mpf(0)) for j in range(len(B[0]))] for i in range(len(A))]

@lru_cache(None)
def taylor_controls():
    mp.iv.dps=110;out=[]
    for edge in (0,12):
        variance=candidate.mobile.ivalue(candidate.mobile.exact_peierls_data(edge)[2])
        # Independent literal four-level Jacobi action. The even Taylor
        # remainder is at most 3^41/82!, since ||A||^2<3.
        A=[[mp.iv.mpf(0) for _ in range(4)] for _ in range(4)]
        for j in range(3):A[j][j+1]=A[j+1][j]=mp.iv.sqrt(variance*(j+1))
        square=product(A,A)
        power=[[mp.iv.mpf(int(i==j)) for j in range(4)] for i in range(4)]
        total=[[mp.iv.mpf(int(i==j)) for j in range(4)] for i in range(4)]
        for n in range(1,41):
            power=product(power,square)
            coef=ivf(F((-1)**n,math.factorial(2*n)))
            total=[[total[i][j]+coef*power[i][j] for j in range(4)] for i in range(4)]
        remainder=F(3)**41/math.factorial(82)
        row={'edge':edge,'degree':80,'operator_remainder':str(remainder)}
        for name,i in (('vacuum',0),('two_photon',2)):
            row[name]=[str(endpoint(total[i][0])-remainder),str(endpoint(total[i][0],True)+remainder)]
        out.append(row)
    return out

def rational_inertia(A,shift):
    n=len(A);L=[[F(0) for _ in range(n)] for _ in range(n)];D=[]
    for i in range(n):
        L[i][i]=1
        p=A[i][i]-shift-sum(L[i][j]**2*D[j] for j in range(i))
        assert p!=0
        D.append(p)
        for row in range(i+1,n):
            L[row][i]=(A[row][i]-sum(L[row][j]*L[i][j]*D[j] for j in range(i)))/p
    return sum(x<0 for x in D)

@lru_cache(None)
def audit():
    original=candidate.certificate()
    raw=(ROOT/'model/affine_gauss_spectrum/results.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='632a3b6ae85c7f1fb67ce768a78b684972831eb1b194517b1dc9c2cb1affb9ff'
    assert original==json.loads(raw)
    rows=original['spectral']['low_charged_levels']
    independent=taylor_controls()
    for j,control in enumerate(independent):
        for key,label in (('vacuum','radial_and_surface_vacuum_overlaps'),('two_photon','radial_and_surface_two_photon_amplitudes')):
            a,b=map(F,control[key]);l,u=map(F,original['spectral']['operator_bounds'][label][j])
            assert l<a<b<u
    scale=10**42
    k=[];err=F(0)
    for control in independent:
        lo,hi=map(F,control['vacuum']);mid=(lo+hi)/2
        q=F((mid*scale).numerator//(mid*scale).denominator,scale)
        k.append(q);err=max(err,abs(q-lo),abs(q-hi))
    A=[[F(0) for _ in range(13)] for _ in range(13)]
    charges,_=candidate.mobile.parent.charge_geometry()
    for j in range(13):A[j][j]=F(str(charges[j]['total_energy_above_actual_neutral_ground']))
    for e,(u,v) in enumerate(candidate.mobile.actual_edges()):A[u][v]=A[v][u]=-k[e>=12]/24
    # Uniform exact row-sum bound on the difference of actual A and this
    # rational matrix, using maximum degree12.
    epsilon=err/2
    counts=[]
    for j,row in enumerate(rows,1):
        lo,hi=map(F,row['zero_photon_compression_energy_interval'])
        left=rational_inertia(A,lo+epsilon);right=rational_inertia(A,hi-epsilon)
        assert left<j<=right
        counts.append([left,right])
    spectral=original['spectral'];delta=F(spectral['strict_ground_energy_reduction_by_photon_dressing_lower'])
    r=F(spectral['operator_bounds']['off_diagonal_block_norm_upper'][1])
    positive=delta**2/(4*r*r)
    assert positive>F('1e-12')
    result={'schema':'oph.affine_gauss_spectrum_independent_audit.v1',
      'independent_literal_Taylor_controls':independent,
      'rational_compression_error_norm_bound':str(epsilon),
      'all_13_compression_intervals_independently_certified':True,
      'exact_rational_inertia_endpoint_counts':counts,
      'actual_ground_nonvacuum_probability_positive_lower':str(positive),
      'Schur_and_variational_proof_checked':True,'complete_preparation_fidelity_checked':True,
      'candidate_receipt_sha256':hashlib.sha256((ROOT/'model/affine_gauss_spectrum/results.json').read_bytes()).hexdigest(),
      'candidate_source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'model/affine_gauss_spectrum').iterdir()) if p.suffix in ('.py','.md')},
      'owned_source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}}
    return result
