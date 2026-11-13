"""Portable tests executed inside the pinned OpenWave package.

Imports only the copied original native Taichi engine, not RER modules.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import taichi as ti
from .engine import NativeKernel, matrices_from_packet, decode


def run(packet, test):
    ti.init(arch=ti.cpu, default_fp=ti.f64, cpu_max_num_threads=4, offline_cache=False)
    mats = matrices_from_packet(packet)
    frames = [decode(x) for x in packet['native_charged_frames']]
    masses = packet['native_mass_coefficients']
    rows = []
    side = 6
    kernel = NativeKernel(side, mats)
    for family, frame in enumerate(frames):
        momentum = np.array([np.sin(2*np.pi/side), 0., 0.]) if test == 'mass_shell' else np.zeros(3)
        symbol = mats[0]+sum(p*g for p,g in zip(momentum,mats[1:]))
        block = frame.conj().T@symbol@frame
        values, vectors = np.linalg.eigh(block)
        expected = np.sqrt(masses[family]**2+momentum@momentum)
        value, state = values[-1], frame@vectors[:,-1]
        residual = float(np.linalg.norm(symbol@state-value*state))
        assert abs(value-expected) < 3e-14 and residual < 3e-14
        x = np.indices((side,side,side))[0]
        wave = np.exp(2j*np.pi*x/side) if test == 'mass_shell' else np.ones_like(x)
        initial = wave[...,None]*state/(side**1.5*np.sqrt(2))
        time = 2. if test == 'stability' else .75
        steps = 80
        out = kernel.evolve(initial,time/steps,steps)
        target = initial*np.exp(-1j*value*time)
        error = float(np.linalg.norm(out-target))
        integration_bound = time**5*packet['symbol_norm_bound']**5/(120*steps**4)
        assert error < integration_bound+2e-12
        probabilities = abs(out)**2
        total = float(probabilities.sum()+.5)
        assert abs(total-1) < 2*integration_bound+2e-12
        # A broken propagator must fail the same full-vector gate.
        wrong = initial*np.exp(-1j*(value+.03)*time)
        wrong_error = float(np.linalg.norm(out-wrong))
        assert wrong_error > 100*integration_bound
        # Every retained coordinate and the original failure event are checked.
        exact_probabilities = abs(target)**2
        variation = float(np.sum(abs(probabilities-exact_probabilities))/2)
        assert variation < 2*integration_bound+2e-12
        rows.append({'family_index':family, 'side':side,
            'full_native_Nambu_coordinates':180*side**3,
            'source_momentum':momentum.tolist(),
            'rest_coefficient_from_source':float(masses[family]),
            'full_symbol_energy':float(value),
            'full_symbol_eigenvector_residual_diagnostic':residual,
            'event_parameter':time,'RK4_steps':steps,
            'exact_arithmetic_integration_error_bound':integration_bound,
            'full_state_error_diagnostic':error,
            'all_coordinate_record_TV_diagnostic':variation,
            'complete_probability_including_failure_diagnostic':total,
            'preparation_failure_probability':.5,
            'changed_frequency_negative_control_error':wrong_error})
    return {'test':test,'verdict':'PASS_FINITE_NATIVE_PLATFORM_TEST','runs':rows,
        'official_admission_or_cell':False,
        'roundoff_certified_by_analytic_truncation_bound':False}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--test',choices=('stability','rest_spectrum','mass_shell'),required=True)
    a=p.parse_args()
    packet=json.loads(a.source.read_text())
    result=run(packet,a.test)
    import openwave
    result.update({'upstream_version':openwave.__version__,
        'taichi_version':list(ti.__version__),
        'source_packet_sha256':hashlib.sha256(a.source.read_bytes()).hexdigest()})
    a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['verdict'],a.test)


if __name__=='__main__':main()
