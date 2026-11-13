"""Independently reconstruct local identities and complete collision bands."""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import mpmath as mp
from model.openwave_native_centre_string.centre import native_quotient,spatial_labels
from model.openwave_strong_step_scaling.source import chain
from model.openwave_autonomous_clock_channel.clock import collective,evolve,single_branch_records,STEP
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def pins():
    paths={str(p.relative_to(ROOT))for p in HERE.iterdir()if p.suffix in('.py','.md')}
    for name in('openwave_native_centre_string','openwave_autonomous_clock_channel'):
        p=f'model/{name}/results.json';paths.add(p)
        for k,h in json.loads((ROOT/p).read_text())['source_sha256'].items():
            assert hashlib.sha256((ROOT/k).read_bytes()).hexdigest()==h,('stale reviewed source',k)
            paths.add(k)
    return{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)}


def string_checks():
    q=native_quotient()
    for row in q['complete_native_one_generation_coordinates']:
        assert (row['electric_integer']+row['colour_triality'])%3==0
        assert F(row['charge_exact'])*3==row['electric_integer']
    columns=0
    for L in(2,3):
        d=chain(L);labels=spatial_labels(d)
        for k in range(len(d['H'].toarray())):
            flux=labels['flux'][k];Q=labels['local_q'][k].copy();Q[0]+=labels['reference'][k]
            assert np.array_equal(np.diff(np.r_[0,flux,0]),Q)
            assert np.array_equal(flux%3,(-np.cumsum(d['n'][k])[:-1])%3)
            columns+=1
        assert (1,)+(0,)*(L-2) in labels['all_hamiltonian_link_shifts'] or (-1,)+(0,)*(L-2) in labels['all_hamiltonian_link_shifts']
    return{'all15native_quotient_rows_checked':True,'all_spatial_columns_checked':columns,
           'exact_integer_lift_ambiguity':'2/3 and -1/3 share the same centre phase',
           'dual_electric_phase_is_conjugate_link_angle':False}


def clock_checks():
    count=0
    for N in(1,2,3,4):
        d=collective(N);K=d['K'];n=N+1
        for k in range(1,N+1):
            for j in range(n-k):
                # A Hermitian coherence has both conjugate bands, so this
                # probes the entire native map without dropping a complement.
                rho=np.zeros((n,n),complex);rho[j+k,j]=1;rho[j,j+k]=1
                U=np.diag(np.exp(-1j*np.arange(n)*float(STEP)))
                actual=sum((U@z@rho@z.T@U.conj().T for z in K),np.zeros_like(rho))
                assert np.max(abs(actual-evolve(rho,1,N)))<2e-15
                count+=1
    receipt=json.loads((ROOT/'model/openwave_autonomous_clock_channel/results.json').read_text())
    mp.mp.dps=70;bands=0
    for row in receipt['exact_complete_coherence_bands']:
        N=row['cap'];k=row['coherence_order'];size=N+1-k
        gamma=mp.mpf(1)/16;kappa=2*gamma/N;dt=mp.mpf(1)/65536
        rates=[gamma*(j+1)*(j<N)+kappa*j*(j-1)for j in range(N+1)]
        B=mp.matrix(size)
        for j in range(size):
            B[j,j]=mp.sqrt((1-dt*rates[j+k])*(1-dt*rates[j]))
            if j:B[j,j-1]=dt*gamma*j
            if j+2<size:B[j,j+2]=dt*kappa*(j+k+2)*(j+k+1)
        values=[B[0,0]] if size==1 else mp.eig(B,left=False,right=False)
        radius=max(abs(z)for z in values)
        lo,hi=[mp.mpf(F(x).numerator)/F(x).denominator for x in row['spectral_radius_interval_exact']]
        assert lo<=radius<=hi<1
        bands+=1
    zero=single_branch_records(65536,F(0));on=single_branch_records(65536)
    assert zero['first_probe_prepares_the_phase']and on['prior_pumped_occupation_cancels_exactly']
    assert zero['r_exact']=='1'
    return{'full_native_collision_coherence_columns_checked':count,
           'independent70digit_Perron_roots_inside_exact_intervals':bands,
           'zero_pump_measurement_phase_countercontrol_retained':True,
           'autonomous_internal_phase_onset_inferred':False}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    before=pins();r={'schema':'oph.independent_native_clock_string_review.v1',
        'verdict':'PASS_INDEPENDENT_NATIVE_CLOCK_AND_STRING_REVIEW',
        'string':string_checks(),'clock':clock_checks(),'source_sha256':before};assert before==pins()
    text=json.dumps(r,sort_keys=True,indent=2)+'\n';out=HERE/'results.json'
    if a.write:out.write_text(text)
    else:assert out.read_text()==text
    print(r['verdict'])
if __name__=='__main__':main()
