"""Source-derived portable matrices, all native coordinates retained."""
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from model.openwave_lepton_mass_shell.spectrum import native_frames
from model.lepton_mass_compiler.analyzer import rational_effect_intervals
from model.common_source.chiral_candidate.kinetic import bdg_symbol

HERE = Path(__file__).resolve().parent


def source_packet():
    data = native_frames()
    intervals = rational_effect_intervals()['effect_intervals_k0_k1_k2']
    masses = np.array([float(sum(map(F, b), F(0))/2) for b in intervals])
    assert all(F(b[1]) < 1 for b in intervals)
    D = np.array(data['D']).astype(complex)
    K = np.kron(D, np.diag(masses))
    K = K+K.T
    M = bdg_symbol((0, 0, 0), K, speed=1)
    matrices = [M]+[bdg_symbol(p, K, speed=1)-M for p in np.eye(3)]
    frames = [np.array(x).astype(complex) for x in data['frames']]
    sparse = []
    for matrix in matrices:
        sparse.append([[int(i), int(j), float(matrix[i,j].real), float(matrix[i,j].imag)]
                       for i,j in zip(*np.nonzero(matrix))])
    for i, a in enumerate(matrices[1:]):
        assert np.linalg.norm(a@M+M@a) < 1e-15
        for j, b in enumerate(matrices[1:]):
            assert np.linalg.norm(a@b+b@a-(2*np.eye(180) if i==j else 0)) < 1e-15
    assert np.linalg.norm(M, 2) < 1
    return {'schema':'oph.native_openwave_port.v1',
        'native_H0_and_three_kinetic_sparse_matrices':sparse,
        'native_charged_frames':[np.stack((x.real,x.imag),axis=-1).tolist() for x in frames],
        'native_mass_coefficients':masses.tolist(), 'rational_mass_intervals':intervals,
        'symbol_norm_bound':2, 'native_CAR_modes_per_site':90, 'Nambu_rows_per_site':180,
        'declared_source':'Periodic kinetic-plus-analyzer-pairing member, phase2/9, speed1, spacing1',
        'source_parameter_is_physical_seconds':False,
        'massless_internal_coefficients_retained':39,
        'port_is_full_doubled_source_not_only_three_Dirac_frames':True,
        'physical_species_or_continuum_identification':False}


if __name__ == '__main__':
    p = HERE/'source_packet.json'
    p.write_text(json.dumps(source_packet(), indent=2, sort_keys=True)+'\n')
    print(p)
