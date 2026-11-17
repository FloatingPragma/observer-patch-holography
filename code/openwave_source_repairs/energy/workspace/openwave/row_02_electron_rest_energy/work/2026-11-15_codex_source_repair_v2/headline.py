#!/usr/bin/env python3
"""Small offline reconstruction; reads bounded arrays, never invokes a producer.

This regenerates the prepared-source/observed-energy headline. The separate
independent review checks every native signed-pair column and error derivation.
"""
from pathlib import Path
from fractions import Fraction as F
import json,math,time
import numpy as np
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def main():
    t=time.monotonic();coupled=read(HERE/'coupled_energy/RECEIPT.json')['scientific']
    joined=read(HERE/'routed_join/RECEIPT.json')['scientific']
    sector=np.load(HERE/'poisson_energy/source_sector_maps.npz');source=np.load(HERE/'routed_join/maps.npz')
    all_records=[]
    for handoff in joined['handoffs']:
        index=handoff['trial_index'];p,n=divmod(handoff['actual_routed_pair_code'],256)
        x=handoff['actual_source'];assert x[p]==1 and x[n]==-1 and sum(abs(v) for v in x)==2
        flags=[];fail=[]
        for name,offset,address,prefix in [('first_cell_certified',0,p,'first'),('second_cell',12,n,'second')]:
            data=np.load(HERE/name/'maps.npz');q=data['certified_Q_proposal'];d=data['dark_integer_basis'][offset:offset+12]
            P=d@np.linalg.inv(d.T@d)@d.T;rho=np.zeros((12,12),complex);rho[address-offset,address-offset]=1
            assert np.linalg.norm(P@rho)==0
            A=np.eye(12)-q.conj().T@q
            for _ in range(512):rho=q@rho@q.conj().T+np.trace(A@rho)*np.eye(12)/12
            assert np.max(abs(rho-source['trial_%02d_%s'%(index,prefix)]))<1e-12
            assert np.trace(P@rho).real>1-1e-7
            flags.append(np.array([np.trace(e@rho).real for e in sector[prefix+'_effects']]))
            fail.append(float(np.trace(sector[prefix+'_processor_failure']@rho).real))
        laws={'single':flags[0],'pair':np.array([flags[0][a]*flags[1][b] for b in (0,1) for a in (0,1)])}
        fringes={}
        for duration,table in coupled['native_flag_controlled_words'].items():
            fringes[duration]={}
            for name,meta in table.items():
                pair=[]
                for setting in ('X','Y'):
                    law=laws[name]@np.asarray(meta['complete_per_input_flag_output_centres'][setting])
                    expected=handoff['complete_energy_output_centres'][duration][name][setting]
                    assert np.max(abs(law-expected['diagnostic_probability_centres']))<1e-12
                    a,b=(law[j] for j in meta['success_outputs']);assert a+b>.99
                    pair.append(float((a-b)/(a+b)))
                fringes[duration][name]=complex(*pair)
        all_records.append(fringes)
    assert len(all_records)==12
    intervals={};observed=all_records[0]
    for name in ('single','pair'):
        epsilon=joined['energy_intervals'][name]['complete_error_used']
        radius=2*math.sqrt(2)*epsilon/(1-2*epsilon);phase=radius/(1-radius)
        z=observed['1/10'][name];e=-10*math.atan2(z.imag,z.real);short=[e-10*phase,e+10*phase]
        z=observed['4'][name];base=-math.atan2(z.imag,z.real)/4
        candidate=[base+k*math.pi/2 for k in range(32) if short[0]-phase/4<=base+k*math.pi/2<=short[1]+phase/4]
        assert len(candidate)==1
        intervals[name]=[candidate[0]-phase/4,candidate[0]+phase/4]
        assert max(abs(np.asarray(intervals[name])-joined['energy_intervals'][name]['long_interval']))<1e-12
    ratio=[intervals['pair'][0]/intervals['single'][1],intervals['pair'][1]/intervals['single'][0]]
    assert ratio[0]<2<ratio[1] and ratio[1]-ratio[0]<.001
    old=coupled['fixed_count_control_records'];a=complex(*old['single']['complex_accepted_fringe']);b=complex(*old['pair']['complex_accepted_fringe'])
    mismatch=abs(b-a*a);assert mismatch>.05
    for record in coupled['negative_source_records'].values():assert record['energy_success_probability']<1e-12
    assert F(joined['complete_join_error_upper'])<F(1,1000)
    print(json.dumps({'headline_pass':True,'actual_source_handoffs':len(all_records),
        'energy_intervals':intervals,'pair_single_ratio_interval':ratio,
        'pair_keV_at_single511_anchor_interval':[511*v for v in ratio],
        'fixed_count_fringe_mismatch':mismatch,'complete_output_error_upper':float(F(joined['complete_join_error_upper'])),
        'fresh_native_calls':0,'scope':'Reconstructed certified-source preparation and retained native detector records',
        'seconds':time.monotonic()-t},sort_keys=True))
if __name__=='__main__':main()
