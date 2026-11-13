"""High-precision direct Fourier reconstruction of every portable probe row."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def high(value):
    f=F.from_float(float(value));return mp.mpf(f.numerator)/f.denominator

def pins():
    from model.openwave_boosted_records.__main__ import pins as parent_pins
    target=ROOT/'model/openwave_boosted_records/results.json';saved=json.loads(target.read_text());out=parent_pins()
    assert saved['source_sha256']==out
    for f in(target,HERE/'__init__.py',HERE/'__main__.py'):
        out[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest()
    return dict(sorted(out.items()))

def calculate():
    target=ROOT/'model/openwave_boosted_records/results.json';saved=json.loads(target.read_text())
    packet=json.loads((ROOT/'model/openwave_boosted_records/source_packet.json').read_text());bindings=pins()
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in bindings.items())
    mp.mp.dps=65;T=np.array(packet['native_positive_frame_real'])+1j*np.array(packet['native_positive_frame_imag']);Q=np.array(packet['charges'])
    assert np.linalg.norm(T.conj().T@T-np.eye(48))<3e-15
    assert np.linalg.norm(np.r_[Q,-Q][:,None]*T-T*Q[None,:])==0
    comparisons=[]
    for row in saved['portable_Taichi_execution']['rows']:
        case=next(c for c in packet['cases']if c['name']==row['case']);n=row['side'];j=case['mode'];m=high(packet['masses_electron_units'][j]);a=8*mp.pi/(3*n);time=high(case['time']);norm=mp.sqrt(sum(mp.mpf(w)**2 for w in case['amplitudes_real']));psi=[];ref=[];mean=mp.mpf(0)
        modes=[]
        for k,w in zip(case['momenta_in_grid_units'],case['amplitudes_real']):
            p=mp.mpf(3)*k/4;E=mp.sqrt(m*m+(2*mp.sin(a*p/2)/a)**2);Ec=mp.sqrt(m*m+p*p);amp=mp.mpf(w)/norm;modes.append((k,amp,E,Ec));mean+=amp*amp*E
        for x in range(n):
            psi.append(sum(amp*mp.exp(-1j*time*E+2j*mp.pi*k*x/n)for k,amp,E,Ec in modes)/mp.sqrt(n))
            ref.append(sum(amp*mp.exp(-1j*time*Ec+2j*mp.pi*k*x/n)for k,amp,E,Ec in modes)/mp.sqrt(n))
        nativeweights=[abs(mp.mpc(high(z.real),high(z.imag)))**2 for z in T[:,j]]
        prob=np.array([[float(mp.mpf('0.5')*abs(z)**2*w)for w in nativeweights]for z in psi]);cont=np.array([[float(mp.mpf('0.5')*abs(z)**2*w)for w in nativeweights]for z in ref]);record=np.array(row['complete_position_native_probe_probabilities'])
        delta=float(np.max(abs(prob-record)));tv=float(np.sum(abs(prob-cont))/2)
        assert delta<3e-12 and abs(prob.sum()+row['preparation_failure_probability']-1)<3e-14
        assert abs(tv-row['same_input_continuum_complete_record_TV_diagnostic'])<3e-12
        assert tv<=float(F(row['continuum_record_TV_upper_exact']))+3e-12
        assert abs(float(mean)-row['energy_mean_electron_units_diagnostic'])<2e-10
        assert row['success_preparation_work_electron_units_diagnostic']==-row['absorbing_probe_source_energy_change_electron_units_diagnostic']
        comparisons.append({'case':row['case'],'side':n,'all_native_position_effects':96*n,'maximum_probability_difference_65_digits':delta,'independent_complete_TV_diagnostic':tv})
    assert len(comparisons)==18
    for row in saved['boosted_energy_and_velocity_checks']:
        n=row['side'];a=8*mp.pi/(3*n);p=mp.mpf(3)/4;E=mp.sqrt(1+(2*mp.sin(a*p/2)/a)**2);v=mp.sin(a*p)/(a*E)
        defect=abs(1-mp.mpf(5)/4*(E-mp.mpf(3)/5*p))
        assert defect<=mp.mpf(str(float(F(row['boost_energy_shell_defect_upper_exact']))))
        assert abs(v-mp.mpf(3)/5)<=mp.mpf(str(float(F(row['group_velocity_error_upper_exact']))))
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in bindings.items())
    return {'verdict':'PASS_INDEPENDENT_COMPLETE_PORTABLE_BOOST_RECORD_REVIEW','source_sha256':bindings,
      'direct_Fourier_precision_decimal_digits':65,'complete_cases':comparisons,'independent_boost_defect_checks':4,
      'native_probe_isometry_and_charge_sectors_verified':True,
      'review_scope':'Independent direct Fourier reconstruction of all18 portable probability tables and all retained native rows, same-input continuum law, work account and exact stated boost bounds. This reviews actual reducing-sector records, not a full interacting Fock evolution or a finite-lattice Lorentz representation.'}

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();r=calculate();out=HERE/'results.json'
    if a.write:out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    else:assert json.loads(out.read_text())==r
    print(r['verdict'])
if __name__=='__main__':main()
