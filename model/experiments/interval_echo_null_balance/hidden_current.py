"""A pair of opposite stationary seam currents with identical passive exterior histories."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import compile_wave
from model.executable_model import cell, federation

HERE = Path(__file__).resolve().parent


def exact_modes(car, gateways):
    adjacency = sp.Matrix(car.adjacency().toarray().astype(int))
    local_laplacian = 5*sp.eye(12) - adjacency
    kernel = (local_laplacian - 6*sp.eye(12)).col_join(sp.eye(12)[gateways, :]).nullspace()
    u, v = sp.GramSchmidt(kernel[:2], True)
    for w in (u, v):
        assert local_laplacian*w == 6*w and all(w[i] == 0 for i in gateways)
    return u, v


def seam_laws(car, u, v):
    rows = []
    for a, b in car.seams:
        # For (u+i v)/sqrt(2), the continuity current from a to b is Im(psi_a* psi_b).
        j = sp.simplify((u[a]*v[b]-v[a]*u[b])/2)
        rows.append(dict(ports=[int(a),int(b)], number_current=str(j), kinetic_current=str(3*j)))
    return rows


def read_current(car, psi, ports):
    a, b = ports
    probabilities = []
    for sign in (1, -1):
        target = np.zeros(12, complex)
        target[[a,b]] = np.array([1, sign*1j])/np.sqrt(2)
        analyzer = compile_wave(car, target)
        probabilities.append(float(abs(np.vdot(analyzer,psi))**2))
    return dict(yes_probabilities=probabilities, kinetic_current=1.5*(probabilities[0]-probabilities[1]))


def run():
    local = cell()
    boards = []
    orientation = None
    for level in range(3):
        board = federation(level)
        gateways = np.flatnonzero(board.degree[:12] == 6).tolist()
        u, v = exact_modes(local, gateways)
        laws = seam_laws(local, u, v)
        if orientation is None:
            orientation = max(laws, key=lambda row: abs(float(sp.sympify(row['kinetic_current']))))['ports']
        exact = next(sp.sympify(row['kinetic_current']) for row in laws if row['ports']==orientation)
        u, v = [np.asarray(w, dtype=float).ravel() for w in (u,v)]
        h = candidate.generator(board)
        offset = 6*candidate.KAPPA
        reads = []
        for name, wave in [('positive',(u+1j*v)/np.sqrt(2)),('negative',(u-1j*v)/np.sqrt(2)),
                           ('standing',u.astype(complex))]:
            state = compile_wave(local, wave)
            full = np.zeros(board.n, complex)
            full[:12] = state
            np.testing.assert_allclose(h @ full,(offset+3)*full,rtol=0,atol=1e-12)
            read = read_current(local,state,orientation)
            a,b = orientation
            # Direct anticommutator of the generator above the uniform orbital with the
            # continuity current i[h_ab,P_b], tested independently of the two analyzers.
            jpsi = np.zeros_like(full)
            jpsi[a] = -.5j*full[b]
            jpsi[b] = .5j*full[a]
            direct = np.vdot(h @ full - offset*full, jpsi).real
            np.testing.assert_allclose(read['kinetic_current'],direct,atol=1e-13)
            expected = {'positive':float(exact),'negative':-float(exact),'standing':0.}[name]
            np.testing.assert_allclose(direct,expected,atol=1e-13)
            reads.append(dict(preparation=name, **read, operator_current=float(direct)))
        boards.append(dict(cells=board.cells,gateways=gateways,selected_ports=orientation,
                           kinetic_current_magnitude=str(abs(exact)), all_internal_seams=laws,
                           all_null_error_floor_per_unit_coupling=str(2*abs(exact)),
                           quadrature_reads=reads))
    return dict(boards=boards)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    result = run()
    path = HERE/'hidden_current.json'
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=2,allow_nan=False)+'\n')
    else:
        from .regional import verify
        verify(json.loads(path.read_text()), result)
    print([(b['cells'], b['selected_ports'], b['kinetic_current_magnitude']) for b in result['boards']])
