"""Independent exact arithmetic, character counting and record checks."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import numpy as np
import sympy as s
from model.openwave_axial_state.source import sector, parameters, weights
from model.openwave_baryon_weak_formfactors.native import currents
from model.openwave_native_force.force import source_matrices
from model.lepton_mass_compiler.analyzer import rational_effect_intervals
from model.particle_common_contract.provenance import import_closure
from model.openwave_larmor_dynamics.certificate import equivalent

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PACKAGES=('openwave_axial_state','openwave_hadron_spin_constraints',
          'openwave_nuclear_distance','openwave_neutrino_calibration',
          'openwave_joint_native_clocks')

def pins():
    seeds={str(p.relative_to(ROOT)) for package in (*PACKAGES,HERE.name)
           for p in (ROOT/'model'/package).iterdir()
           if p.suffix in ('.py','.md','.tex') or p.name in ('results.json','inputs.json','oscillation_inputs.json')}
    seeds.discard(str((HERE/'results.json').relative_to(ROOT)))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(import_closure(seeds))}

def rational_action(matrix,v):
    a=matrix.tocsr()
    return [sum(F(int(a.data[k]))*v[a.indices[k]]
                for k in range(a.indptr[i],a.indptr[i+1])) for i in range(a.shape[0])]

def run():
    before=pins()
    receipts={p:json.loads((ROOT/'model'/p/'results.json').read_text()) for p in PACKAGES}
    for receipt in receipts.values():
        for p,h in receipt['source_sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    ground=receipts[PACKAGES[0]]['current_ground_certificate']
    d=sector('baryon');v=list(map(int,ground['integer_native_occupation_trial']))
    norm=sum(x*x for x in v);assert norm==int(ground['integer_trial_norm_squared'])
    hv=[F(0) for _ in v]
    for k,c in weights(parameters()).items():
        hv=[a+c*b for a,b in zip(hv,rational_action(d['parts'][k],v))]
    energy=sum(a*b for a,b in zip(v,hv))/norm
    variance=sum((a-energy*b)**2 for a,b in zip(hv,v))/norm
    assert energy==F(ground['trial_energy_rational'])
    assert variance==F(ground['trial_residual_squared_rational'])
    native=currents()['full']
    vector=sum(native[x][0][0] for x in range(2))
    axial=sum(native[x][3][0] for x in range(2))
    vv=np.asarray(v,float);native_value=np.vdot(vector@vv,axial@vv).real/norm
    lo,hi=map(F,ground['native_ground_axial_interval_rational'])
    assert float(lo)-1e-12<native_value<float(hi)+1e-12
    records=ground['complete_trial_original_records']
    assert sum(map(F,records['success_native_occupation_probabilities']))==F(1,2024)
    assert sum(map(F,records['failure_native_occupation_probabilities']))==F(2023,2024)
    # Independent exact operator proof of the fixed-D shortcut, before reduction.
    B=d['Binteger'];identity=d['parts']['C48']+32*d['parts']['B9']-96*d['parts']['N']
    assert (identity@B).nnz==0
    spin=receipts[PACKAGES[1]]['common_spin_obstruction']
    required=spin['rho_required_global_spin_interval']
    assert float(spin['Delta_fit_global_spin'])<float(required[0])
    # Independent SU(3) Weyl-character count on all 24 choose 6 occupations.
    modes=[(cw[0],cw[1],sz,iz) for cell in range(2)
           for sz,iz,cw in product((-1,1),(-1,1),((1,0),(0,1),(-1,-1)))]
    counts=Counter();zero=Counter()
    for occupied in combinations(modes,6):
        weight=tuple(map(sum,zip(*occupied)))
        counts[weight[:2]]+=1
        if weight[2:]==(0,0):zero[weight[:2]]+=1
    roots=((1,-1),(2,1),(1,2));full=fixed=0
    for bits in product((0,1),repeat=3):
        shift=tuple(sum(b*r[j] for b,r in zip(bits,roots)) for j in range(2))
        sign=(-1)**sum(bits);full+=sign*counts[shift];fixed+=sign*zero[shift]
    assert (full,fixed)==(2520,240)
    nuclear=receipts[PACKAGES[2]]['nuclear']
    census=nuclear['all_spin_isospin_parity_channels']
    assert sum(c['fixed_Ms0_Mi0_dimension'] for c in census)==fixed
    assert sum(c['fixed_Ms0_Mi0_dimension']*(2*c['spin']+1)*(2*c['isospin']+1) for c in census)==full
    # Independent family rotation and exact mass-splitting arithmetic.
    neutral=receipts[PACKAGES[3]]['calibrated_neutral_catalogue'];inp=neutral['calibration_inputs']
    sq=[F(inp['declared_lightest_mass_eV'])**2]
    sq.extend(sq[0]+F(inp[k]) for k in ('delta_m21_squared_eV2','delta_m31_squared_eV2'))
    for row,x in zip(neutral['rows'],sq):
        l,u=map(F,row['mass_eV_interval']);assert l*l<=x<=u*u
    U=np.asarray(neutral['PMNS_real'])+1j*np.asarray(neutral['PMNS_imag'])
    assert np.linalg.norm(U.conj().T@U-np.eye(3))<1e-14
    for key,i,j in (('sin2_theta12',0,1),('sin2_theta23',1,2)):
        assert abs(abs(U[i,j])**2/(1-float(inp['sin2_theta13']))-float(inp[key]))<1e-14
    assert abs(abs(U[0,2])**2-float(inp['sin2_theta13']))<1e-14
    tau,e,mu=[tuple(map(F,b)) for b in rational_effect_intervals()['effect_intervals_k0_k1_k2']]
    rigid=(mu[1]**2-e[0]**2)/(tau[0]**2-e[1]**2)
    assert rigid<F(inp['solar_3sigma_eV2'][0])/F(inp['normal_atmospheric_3sigma_eV2'][1])
    eta=F(5,10**10)**2/(4*F(49,100)**2)
    assert (F(41,10**10)-F(1,10**14))**2>64*eta
    # Independently reconstruct native spin-chain clocks with exact Pauli matrices.
    I=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1)
    majorana=[s.kronecker_product(*([Z]*j+[end]+[I]*(2-j))) for j in range(3) for end in (X,Y)]
    reconstructed=3*s.eye(8)
    for j in range(3):reconstructed+=s.I*majorana[2*j]*majorana[2*j+1]/2
    for j in range(2):reconstructed+=s.I*majorana[2*j+1]*majorana[2*j+2]/7
    assert reconstructed==source_matrices()[2]
    coefficient=F(1,2)*(F(121,100)*8)**2*(F(13,4)**2+F(61,50)**2)
    assert coefficient<600 and F(600)/F(4,5)==750
    clocks=receipts[PACKAGES[4]]
    assert 31<float(clocks['event16_relative_response_difference_ppm_interval'][0])
    assert not clocks['physical_clock_G_or_equivalence_principle_violation_identified']
    assert pins()==before
    return {'schema':'oph.independent_source_constraints_review.v1','verdict':'PASS_INDEPENDENT_SOURCE_CONSTRAINTS_REVIEW',
            'source_sha256':before,'exact_native_trial_energy':str(energy),'exact_native_trial_variance':str(variance),
            'native_axial_diagnostic':native_value,'full_original_trial_records_exactly_normalized':True,
            'fixed_D_native_operator_identity_checked':True,'independent_complete_six_quark_character_count':full,
            'independent_zero_spin_weak_weight_count':fixed,'calibrated_neutral_mass_brackets_and_rotation_checked':True,
            'rigid_neutral_hierarchy_excluded':True,'uniform_neutral_complete_record_bound_checked':True,
            'independent_exact_native_clock_reconstruction':True,'weak_probe_remainder_coefficient':str(coefficient),
            'scope':'Independent checks of stated finite algebra, calibration ancestry and selected bounds. Does not certify a common physical continuum or convert calibrated inputs into predictions.'}

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    value=run();path=HERE/'results.json'
    if a.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),value)
    print(value['verdict'])

if __name__=='__main__':main()
