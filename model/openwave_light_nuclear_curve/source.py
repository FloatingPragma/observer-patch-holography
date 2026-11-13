"""The native three- and four-baryon quartet compressions."""
import json
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import sys
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients

RER = Path(__file__).resolve().parents[2]
PARENT = RER / 'model/openwave_common_hyperfine_spectrum/inputs.json'
RAW = json.loads(PARENT.read_text())
PARAMETERS = {k: F(v) for k, v in RAW['declared_rational_parameters'].items()}
UNIT = F(RAW['electron_unit_GeV'])


def lift(a, x, count=3):
    result = np.ones((1, 1), complex)
    for i in range(count):
        result = np.kron(result, a if i == x else np.eye(4))
    return result


@lru_cache(None)
def projected_source(cells=3, three_body=F(0), vacuum_compatible=False):
    """Literal compressed current algebra, not a replacement for full Fock H."""
    assert cells in (3, 4)
    p = PARAMETERS
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    S = [[lift(np.kron(a / 2, np.eye(2)), x, cells) for a in pauli] for x in range(cells)]
    I = [[lift(np.kron(np.eye(2), a / 2), x, cells) for a in pauli] for x in range(cells)]
    identity = np.eye(4**cells)
    local = lambda a: sum(x @ x for x in a)
    total = lambda a: [sum(a[x][k] for x in range(cells)) for k in range(3)]
    TS, TI = total(S), total(I)
    S2, I2 = local(TS), local(TI)
    scalar = 3 * cells * p['mbar'] + cells * p['A'] + F(3*cells, 4) * (p['Js'] + p['Jw'])
    # The positive graph Laplacian has only its occupation diagonal in this P.
    scalar += (3 * cells if vacuum_compatible else 6 * (cells-1)) * p['t']
    terms=[(scalar,identity),(-p['delta'],TI[2])]
    for x in range(cells-1):
        terms.append((p['spin_gradient'],local([S[x][k]-S[x+1][k] for k in range(3)])))
        terms.append((p['eta'],local([I[x][k]+I[x+1][k] for k in range(3)])))
    if vacuum_compatible:
        terms.append((-F(3,4)*(cells-2)*(p['spin_gradient']+p['eta']),identity))
        # All local P cells are colour singlets with exactly three quarks.
        dU,dm=degree_coefficients()
        terms.append((3*(cells-2)*dm,identity))
    for x in range(cells-2):
        window = [sum(S[j][k] for j in range(x, x+3)) for k in range(3)]
        terms.append((three_body,local(window)-7/4*identity))
    h=sum(float(c)*a for c,a in terms)
    assert np.max(abs(h-h.conj().T)) < 1e-12
    for g in (S2, I2, TS[2], TI[2]):
        assert np.max(abs(h@g-g@h)) < 2e-10
    return {'H': h.real, 'S2': S2.real, 'I2': I2.real,
            'Sz': TS[2].real, 'Iz': TI[2].real, 'dimension': len(h),
            'rational_terms':tuple((F(c),a.real) for c,a in terms)}


def levels(cells=3, three_body=F(0), vacuum_compatible=False):
    d = projected_source(cells, three_body, vacuum_compatible)
    rows = []
    for twice_s in range(cells % 2, cells+1, 2):
        for twice_i in range(cells % 2, cells+1, 2):
            s, i = twice_s/2, twice_i/2
            iz = -.5 if cells % 2 else 0.
            z = np.eye(d['dimension'])
            selector = sum((g-target*z)@(g-target*z) for g,target in
                           ((d['S2'],s*(s+1)),(d['I2'],i*(i+1)),(d['Sz'],s),(d['Iz'],iz)))
            w, v = eigh(selector)
            B = v[:, w<1e-10]
            if B.shape[1]:
                energy = eigh(B.T@d['H']@B,eigvals_only=True)
                rows.append({'spin':s,'isospin':i,'Iz':iz,'multiplicity_at_fixed_projections':len(energy),
                             'all_projected_energies_GeV':(energy*float(UNIT)).tolist()})
    return rows

