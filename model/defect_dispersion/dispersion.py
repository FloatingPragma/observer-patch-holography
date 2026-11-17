"""Periodic geometry controls; crystal band energies are not physical masses."""
from __future__ import annotations
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy import linalg as la
from model.compact_defect_pair.spectral import SIGMA
from model.compact_defect_pair.topology import direction

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def bloch(side, mass, theta, *, uniform=False):
    if type(side) is not int or side < 6 or side % 2:
        raise ValueError('an even supercell side at least six is required')
    if not np.isfinite(mass) or mass < 8:
        raise ValueError('finite mass at least eight is required')
    t = np.asarray(theta)
    if t.shape != (3,) or not np.isrealobj(t) or not np.isfinite(t).all():
        raise ValueError('three finite real Bloch phases are required')
    points = tuple(product(range(-side//2, side//2), repeat=3))
    address = {p: i for i, p in enumerate(points)}
    h = np.zeros((2*len(points),)*2, complex)
    for p, i in address.items():
        n = np.array((0, 0, 1) if uniform else direction(p), float)
        length = la.norm(n)
        if length:
            n /= length
        h[2*i:2*i+2, 2*i:2*i+2] = mass*sum(n[j]*SIGMA[j] for j in range(3))
        for axis in range(3):
            q = list(p)
            q[axis] += 1
            phase = 1
            if q[axis] == side//2:
                q[axis] = -side//2
                phase = np.exp(1j*t[axis])
            j = address[tuple(q)]
            block = .5j*phase*SIGMA[axis]
            h[2*i:2*i+2, 2*j:2*j+2] = block
            h[2*j:2*j+2, 2*i:2*i+2] = block.conj().T
    return h


def homogeneous_spectrum(side, mass, theta):
    q = np.asarray(list(product(range(side), repeat=3)))
    s = np.sin((2*np.pi*q + theta)/side)
    energy = np.sqrt(s[:, 0]**2+s[:, 1]**2+(mass-s[:, 2])**2)
    return np.sort(np.r_[-energy, energy])


def two_cell_fourier_control(side=6, mass=8, axis=2):
    """All rows of the actual two-cell periodic source, without truncation."""
    h0 = bloch(side, mass, [0, 0, 0])
    theta = np.zeros(3); theta[axis] = np.pi
    hp = bloch(side, mass, theta)
    # Construct every actual edge on a two-cell periodic graph independently.
    points = tuple(product(range(-side//2, side//2), repeat=3))
    address = {p: i for i, p in enumerate(points)}
    n = len(h0)
    full = np.zeros((2*n, 2*n), complex)
    for cell in range(2):
        for p, i in address.items():
            v = np.array(direction(p), float)
            if la.norm(v):
                v /= la.norm(v)
            at = cell*n+2*i
            full[at:at+2, at:at+2] = mass*sum(v[j]*SIGMA[j] for j in range(3))
            for j in range(3):
                q = list(p); q[j] += 1; other = cell
                if q[j] == side//2:
                    q[j] = -side//2
                    if j == axis:
                        other = 1-cell
                to = other*n+2*address[tuple(q)]
                block = .5j*SIGMA[j]
                full[at:at+2, to:to+2] = block
                full[to:to+2, at:at+2] = block.conj().T
    eye = np.eye(len(h0))
    f = np.block([[eye, eye], [eye, -eye]])/np.sqrt(2)
    diagonal = la.block_diag(h0, hp)
    error = la.norm(full@f-f@diagonal)
    if error > 1e-10:
        raise ArithmeticError('complete two-cell Fourier identity failed')
    return {'one_body_dimension': len(full), 'all_row_intertwiner_error': float(error),
            'crystal_momenta': [0, 'pi/side'], 'discarded_rows': 0}


@lru_cache(maxsize=1)
def calculate():
    cases = []
    # pi/2 is an exact momentum of a four-supercell periodic ring.
    phases = [(None, 0)] + [(a, t) for a in range(3) for t in [-np.pi/2, np.pi/2]]
    for side, mass in [(6, 8), (6, 64), (8, 8), (8, 64)]:
        rows = []
        for axis, phase in phases:
            theta = np.zeros(3)
            if axis is not None:
                theta[axis] = phase
            h = bloch(side, mass, theta)
            e, w = la.eigh(h, check_finite=True)
            band = abs(e) < 3.000000001
            residual = la.norm(h@w-w*e, 'fro')
            if sum(band) != 4 or min(abs(e[~band])) < mass-3-1e-9 or residual > 1e-8:
                raise ArithmeticError('complete band or spectral residual failed')
            rows.append({'axis': axis, 'phase': float(phase),
                'complete_energies': e.tolist(), 'low_energies': e[band].tolist(),
                'complement_count': int(sum(~band)),
                'minimum_complement_abs_energy': float(min(abs(e[~band]))),
                'all_eigenvector_residual_Frobenius': float(residual),
                'uniform_eigenmode_add_remove_probability': 1/len(e),
                'nonzero_probe_outcomes': len(e), 'zero_probe_outcomes': len(e)})
        e0 = np.array(rows[0]['low_energies'])[2:]
        for row in rows[1:]:
            e = np.array(row['low_energies'])[2:]
            row['positive_branch_E2_secant'] = ((e**2-e0**2)/(row['phase']/side)**2).tolist()
            row['positive_branch_energy_change'] = (e-e0).tolist()
        cases.append({'side': side, 'texture_mass': mass,
            'one_body_dimension': 2*side**3, 'rows': rows,
            'isolated_band_dimension_exact': 4,
            'isolated_band_abs_energy_upper_exact': 3,
            'complement_abs_energy_lower_exact': mass-3,
            'full_finite_ring_records_computed': False})
    control = la.eigvalsh(bloch(6, 8, [.2, -.3, .4], uniform=True))
    control_error = float(max(abs(control-homogeneous_spectrum(6, 8, [.2, -.3, .4]))))
    if control_error > 1e-10:
        raise ArithmeticError('homogeneous full-spectrum control failed')
    return {'schema': 'oph.defect_dispersion.v1', 'cases': cases,
        'two_cell_Fourier': two_cell_fourier_control(),
        'homogeneous_complete_spectrum_error': control_error,
        'readout': {'state': 'filled negative-energy sea of the quadratic periodic member',
            'clock': 'hbar_model=1; same hopping coefficient 1/2 on every bond',
            'instrument': 'all eigenmode additions and removals, each divided by sqrt(2*side^3)',
            'energy_transfer': 'abs(lambda); addition momentum theta/side, hole momentum -theta/side',
            'completeness': 'sum_j (c_j c_j^dagger+c_j^dagger c_j)/(2*side^3)=I',
            'reference': 'dual charge and opposite momentum carried by each probe branch'},
        'scope': {'geometry_is_open_compact_member': False,
            'all_interacting_common_source_spectra_computed': False,
            'selected_Bloch_blocks_are_complete': True,
            'all_crystal_momenta_sampled': False,
            'relativistic_mass_assignment_established': False,
            'GeV_calibration': None,
            'numerical_secants_are_certified_curvature_intervals': False,
            'claim': 'anisotropic nearly flat defect bands; no shared relativistic rest interpretation from Gamma energies'},
        'source_bindings': bindings()}


def bindings():
    names = ['model/defect_dispersion/dispersion.py', 'model/defect_dispersion/PROOF.md',
             'model/compact_defect_pair/spectral.py', 'model/compact_defect_pair/topology.py',
             'model/lattice_texture_selection/shells.py']
    return {n: hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}


def compare(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(compare(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(compare(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return type(b) in (float, int) and np.isfinite(b) and abs(a-b) <= 2e-10*max(1, abs(a))
    return a == b


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    result = calculate(); path = HERE/'results.json'
    if args.check:
        if not compare(result, json.loads(path.read_text())):
            raise SystemExit('periodic dispersion receipt differs')
    else:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('four complete-band geometry controls; physical mass assignment unproved')

if __name__ == '__main__':
    main()
