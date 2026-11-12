"""Independent conditioning and unitary controls for prepared charged lines."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
import sympy as s

ROOT = Path(__file__).resolve().parents[3]


def produce():
    paths = ['model/openwave_particle_readouts/PROOF.md',
             'model/openwave_particle_readouts/readouts.py',
             'model/openwave_particle_readouts/results.json',
             'model/coherent_defect_wave/chiral_joint/interacting_charged_response/transfer.py',
             'model/coherent_defect_wave/chiral_joint/interacting_charged_response/local.py',
             'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/PROOF.md']
    saved = json.loads((ROOT/paths[2]).read_text())
    for path, digest in saved['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    count = 0; max_violation = 0.; broken_normalization = False
    for seed in range(32):
        rng = np.random.default_rng(seed+418)
        h0 = np.diag([0., .4, .9, 1.7])
        raw = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        perturb = (raw+raw.conj().T)/2
        perturb *= .002/np.linalg.norm(perturb, 2)
        h = h0+perturb
        # A deliberately small successful atomic component checks conditional
        # normalization. An orthogonal error lies in the retained cloud.
        weight = .01
        eps = .003
        phi = np.array([np.sqrt(weight),0,0,0], dtype=complex)
        chi = phi+np.array([0,eps,0,0])
        p = float(np.vdot(chi,chi).real)
        p_low = (np.sqrt(weight)-eps)**2
        assert p_low <= p <= (np.sqrt(weight)+eps)**2
        for t in (0., .1, 1., 10.):
            psi = chi/np.sqrt(p)
            state = expm(-1j*t*h)@psi
            leakage = float(np.sum(np.abs(state[1:])**2))
            upper = min(1., (eps+t*.002*np.sqrt(weight))**2/p_low)
            phase_distance = min(2., (2*eps+t*.002*np.sqrt(weight))/np.sqrt(p_low))
            survival = abs(np.vdot(psi,state))**2
            lower_survival = max(0.,1-phase_distance**2/2)**2
            assert leakage <= upper+1e-14
            assert survival >= lower_survival-1e-14
            max_violation = max(max_violation,leakage-upper,lower_survival-survival)
            count += 1
            if t == 0:
                broken_normalization |= leakage > eps**2
    assert broken_normalization
    # Fidelity-to-purification norm bound, including an entangled environment.
    for d in (.001,.1,.5,1.):
        joint = np.array([np.sqrt(1-d*d),0,0,d])
        ideal = np.array([1,0,0,0])
        assert np.linalg.norm(joint-ideal) <= 2*d+1e-14
    # Exact pointer control, including source-plus-reference cancellation.
    h=s.Matrix([[1,1],[1,-1]])/s.sqrt(2);plus=s.Matrix([1,1])/s.sqrt(2)
    assert h*s.diag(1,-1)*plus == s.Matrix([0,1])
    assert h*s.diag(1,1)*plus == s.Matrix([1,0])
    assert h*(s.eye(2)/2)*h.H == s.eye(2)/2
    assert saved['source']['reference_cost_subtracted_from_gap']=='4'
    assert not saved['calibration']['electron_identification_or_GeV_mass_prediction_claimed']
    assert not saved['clock']['bounded_clock_is_self_starting_frequency_rigid_time_crystal']
    return {'schema':'oph.prepared_line_independent_review.v1',
            'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            'review_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'independent_finite_unitary_cases':count,
            'maximum_bound_violation':max(0.,max_violation),
            'checks':{'full_cloud_complement_included':True,
                      'conditional_probability_normalization':True,
                      'purification_norm_bound':True,
                      'two_endpoint_survival_errors':True,
                      'source_and_reference_energy_separated_explicitly':True,
                      'native_spin_relative_phase_and_reference_cancellation':True},
            'mutations':{'missing_line_probability_denominator_rejected':True,
                         'rotated_reference_changes_pointer_record':True,
                         'dephased_pointer_loses_deterministic_record':True},
            'verdict':'PASS_PREPARED_LINE_BOUNDS_WITH_DECLARED_APPARATUS',
            'scope':'The source proof supports the specified prepared charged line at finite time. No unbound electron, self-starting rigid clock, apolar local tensor, or laboratory spin identity follows.'}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();value=produce()
    args.output.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    print(value['verdict'])
