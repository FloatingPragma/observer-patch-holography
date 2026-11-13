"""Untrusted integer-frame proposal; the producer checks every finite bound."""
import json
from fractions import Fraction as F
import numpy as np
from scipy.optimize import brentq
from .schur import HERE, blocks, shifted
from model.openwave_coupled_haar_spectrum.certificate import source_digest


def main():
    result = {}
    for which, window in (('coarse', (95., 98.)), ('refined', (55., 59.))):
        source, grouped, norm2, count = blocks(which)
        h = np.array(source['H'], float)
        arrays = [(float(d), np.array(__import__('sympy').SparseMatrix(h.shape[0], h.shape[1], terms), float)) for d, terms in grouped.items()]
        def eigen(z):
            matrix = h.copy()
            for d, gram in arrays:
                matrix -= gram/(d-z)
            return np.linalg.eigvalsh(matrix)[1]-z
        center = brentq(eigen, *window, xtol=1e-12)
        floor, magnetic = float(source['floor']), float(source['faces']*source['beta'])
        squared_gram = sum((gram/(d-center)**2 for d, gram in arrays), np.zeros_like(h))
        err = magnetic*np.max(np.sum(abs(squared_gram), axis=1))/(1-magnetic/(floor+magnetic-center))
        radius = 1.02*err + 2e-8
        packets = {}
        for label, value in (('lower', center-radius), ('upper', center+radius)):
            z = F(int(np.floor(value*10**10)) if label == 'lower' else int(np.ceil(value*10**10)), 10**10)
            exact, gram, error = shifted(which, z)
            values, frame = np.linalg.eigh(np.array(exact, float)); D = 2**40
            packets[label] = {'energy_exact': str(z), 'exact_matrix_sha256': source_digest(exact),
                              'complete_frame': {'dyadic_denominator': D,
                                'dyadic_eigenframe': [[str(round(float(v)*D)) for v in row] for row in frame],
                                'dyadic_eigenvalues': [str(round(float(v)*D)) for v in values]}}
        result[which] = packets
    (HERE/'numerical_witnesses.json').write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')


if __name__ == '__main__':
    main()
