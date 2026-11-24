"""All normalized null contractions from diagonal radar and matter records."""
import argparse
import json
from pathlib import Path

import numpy as np

from . import measure


def axial_residual(observation, coupling):
    curvature = observation['temporal'] + observation['spatial']
    matter = observation['energy'] + observation['pressure']
    return curvature - coupling * matter


def run():
    bank = json.loads((measure.HERE / 'results.json').read_text())
    diagonals = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1]]) / np.sqrt(3)
    directions = np.vstack([np.eye(3), -np.eye(3), diagonals])
    fixed_coupling = measure.classical.COUPLING
    surveys = []
    for survey in bank['surveys']:
        measurements = {
            name: [measure.read(tape, bank['recovered_mass2']) for tape in tapes]
            for name, tapes in survey['records'].items()
        }
        coupling, _ = measure.inferred(measurements['ground'][0], measurements['x'][0])
        held = [('x', j, item) for j, item in enumerate(measurements['x']) if j]
        held.extend((name, j, item) for name in ('y', 'z', 'even_x')
                    for j, item in enumerate(measurements[name]))
        rows = []
        for name, appointment, item in held:
            axes = axial_residual(item, coupling)
            rows.append(dict(source=name, appointment=appointment, axial=axes.tolist(),
                             nine_rays=(directions**2 @ axes).tolist(),
                             all_sky_min=float(axes.min()), all_sky_max=float(axes.max())))
        largest = max(max(abs(row['all_sky_min']), abs(row['all_sky_max'])) for row in rows)
        np.testing.assert_allclose(largest, survey['max_held_axial_null_error'], atol=1e-12)
        fixed_error = max(float(np.max(np.abs(axial_residual(item, fixed_coupling))))
                          for _, _, item in held)
        control_error = max(float(np.max(np.abs(axial_residual(item, coupling))))
                            for item in measurements['decoupled_y'])
        surveys.append(dict(spacing=survey['spacing'], recovered_coupling=coupling,
                            held_all_sky_max_abs=largest, held=rows,
                            held_fixed_hamiltonian_coupling_max_abs=fixed_error,
                            decoupled_all_sky_max_abs=control_error))
    return dict(null_spatial_unit_directions=directions.tolist(),
                fixed_hamiltonian_coupling=fixed_coupling, surveys=surveys)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    out = run()
    dest = Path(__file__).with_name('null_sky.json')
    if parser.parse_args().write:
        dest.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    else:
        stored = json.loads(dest.read_text())
        np.testing.assert_equal(stored['null_spatial_unit_directions'], out['null_spatial_unit_directions'])
        np.testing.assert_equal(stored['fixed_hamiltonian_coupling'], out['fixed_hamiltonian_coupling'])
        for old, new in zip(stored['surveys'], out['surveys'], strict=True):
            for key in ('spacing', 'held_all_sky_max_abs', 'decoupled_all_sky_max_abs',
                        'held_fixed_hamiltonian_coupling_max_abs', 'recovered_coupling'):
                np.testing.assert_allclose(old[key], new[key], atol=2e-13, rtol=0)
            for a, b in zip(old['held'], new['held'], strict=True):
                assert (a['source'], a['appointment']) == (b['source'], b['appointment'])
                for key in ('nine_rays', 'axial', 'all_sky_min', 'all_sky_max'):
                    np.testing.assert_allclose(a[key], b[key], atol=2e-13, rtol=0)
    for survey in out['surveys']:
        print(survey['spacing'], survey['recovered_coupling'], survey['held_all_sky_max_abs'],
              survey['held_fixed_hamiltonian_coupling_max_abs'], survey['decoupled_all_sky_max_abs'])
