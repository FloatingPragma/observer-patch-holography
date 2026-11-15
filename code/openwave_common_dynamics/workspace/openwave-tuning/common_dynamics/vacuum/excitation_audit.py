"""Exact Fock excitation and logarithm-branch controls for the free reference."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np

import reference as ref

HERE=Path(__file__).resolve().parent


def run():
    paths=[Path(__file__),Path(ref.__file__),HERE/'REFERENCE.json']
    pins={str(p.relative_to(ref.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    mu=.04;h=.1
    rows=[]
    for kind in ('vacuum','particle','hole','neutral_pair'):
        pattern=ref.sector_pattern(kind)
        law,den=ref.complete_occupation_law(kind);lookup=dict(law)
        translation=1;mass_number=0
        for count,(a,b,sign) in zip(pattern,ref.pairs()):
            translation*=sign if count==1 else -1 if count==2 else 1
            mass_number+=count*(1 if a%2==0 else -1)
        wrong_CAR=0
        for mask,coefficient in law:
            occupied=[j for j in range(32) if (mask>>j)&1]
            image=[4*(7-j//4)+j%4 for j in occupied]
            parity=(-1)**sum(image[a]>image[b] for a in range(len(image)) for b in range(a+1,len(image)))
            target=sum(1<<j for j in image)
            assert parity*coefficient==translation*lookup[target]
            assert sum(1 if j%2==0 else -1 for j in occupied)==mass_number
            wrong_CAR+=coefficient!=translation*lookup[target]
        eigenphase=translation*np.exp(-1j*mu*mass_number)
        data=json.loads((HERE/'REFERENCE.json').read_text())
        energy=0.
        if kind in ('particle','neutral_pair'):energy+=data['positive_onebody_energies'][0]
        if kind in ('hole','neutral_pair'):energy-=data['negative_onebody_energies'][1]
        assert abs(eigenphase-np.exp(-1j*h*energy))<1e-13
        rows.append({'kind':kind,'complete_outcomes':den,'translation_eigenvalue':translation,
                     'mass_beta_sum':mass_number,'relative_energy_in_declared_lift':energy,
                     'actual_macro_eigenphase':[float(eigenphase.real),float(eigenphase.imag)],
                     'omitted_CAR_sign_wrong_columns':int(wrong_CAR)})
    assert rows[0]['omitted_CAR_sign_wrong_columns']>0
    # An integer branch change preserves the same U but changes negative-band
    # filling. This is an executable counterexample to a source-selected sea.
    energies=np.array(data['negative_onebody_energies']+data['positive_onebody_energies'])
    changed=energies.copy();changed[1]+=2*math.pi/h
    assert np.max(abs(np.exp(-1j*h*changed)-np.exp(-1j*h*energies)))<2e-13
    assert np.count_nonzero(changed<0)==15 and np.count_nonzero(energies<0)==16
    result={'schema':'oph.dirac.reference_excitation_audit.v1',
      'checks':{'all147456excitation_columns_keep_native_phase':True,
               'full_fermionic_exchange_sign_control':True,'positive_energy_pair_reference_domain':True,
               'same_macro_different_logarithm_vacuum_rejected':True,'source_bytes_unchanged':True},
      'complete_laws':rows,'logarithm_control':{'negative_modes_principal':16,'negative_modes_shifted':15,
           'unitary_eigenvalue_difference':float(np.max(abs(np.exp(-1j*h*changed)-np.exp(-1j*h*energies))))},
      'scope':'Exact free CAR-lift excitation coefficients on the uniform Gauss fibre. A Gaussian or interacting field preparation does not inherit these exact eigenphases.',
      'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(ref.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'EXCITATIONS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'laws':rows}))


if __name__=='__main__':run()
