"""Generate numerical proof coordinates; the certificate checks them exactly."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh
from .certificate import source_digest
from model.openwave_coupled_haar_source.source import matrices, hamiltonian


def main():
    records = {}
    for shared, g2, cap in ((True,2,4),(True,2,6),(True,2,8),(True,4,6),(True,8,6),(False,2,6)):
        basis, _ = matrices(cap)
        source = hamiltonian(cap, g2, shared=shared)
        energies, vectors = eigh(np.array(source.evalf(18), float))
        D = 2**40
        key = f'{"shared" if shared else "independent"}_g2_{g2}_cap{cap}'
        records[key] = {'cap': cap, 'g_squared_exact': str(g2), 'spacing_exact': '1', 'shared_edge': shared,
            'exact_source_matrix_sha256': source_digest(source), 'ordered_spin_network_labels': [list(x) for x in basis],
            'ground_denominator': 2**48,
            'ground_integer_coordinates': [str(round(float(x)*2**48)) for x in vectors[:, 0]],
            'low_trial_integer_coordinates': [[str(round(float(x)*2**48)) for x in vectors[:, j]] for j in range(3)],
            'complete_finite_frame': {'dyadic_denominator': D,
                'dyadic_eigenvalues': [str(round(float(x)*D)) for x in energies],
                'dyadic_eigenframe': [[str(round(float(x)*D)) for x in row] for row in vectors]}}
        print(key, len(basis), flush=True)
    Path(__file__).with_name('numerical_witnesses.json').write_text(json.dumps(records, sort_keys=True, indent=2)+'\n')

if __name__ == '__main__':
    main()
