"""Complete native quark-removal target blocks and metric CAR images."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source
from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from model.coherent_defect_wave.chiral_joint.local_reference.algebra import particle_hole
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.coherent_defect_wave.chiral_joint.local_reference.physical import columns
from model.coherent_defect_wave.chiral_joint.local_reference.reduction import (
    coordinate_kernel, lepton_operators,
)
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory

from .quark import target_module


# Original removal in the actual representative component (0, 12, 22, 23).
# u^c and d^c are particle-hole transformed; Q is not.
CHANNELS = {'uc': (True, 0, 0, 4), 'dc': (True, 1, 0, -2),
            'Q': (False, 3, 1, -1)}


@lru_cache(maxsize=None)
def weight_basis(create, two_spin, two_weak, six_hyper):
    q = target_module(create)
    lepton = lepton_operators()
    scalar = higgs_memory()[1]
    labels = tuple((mask.bit_count(), int(2*lepton['spin'][2][mask, mask]),
                    int(2*lepton['currents'][2][mask, mask]),
                    int(6*lepton['currents'][3][mask, mask])) for mask in range(64))
    return tuple((i, mask, h) for i, (w, y) in enumerate(q['weights'])
                 for mask, (n, spin, lw, ly) in enumerate(labels)
                 if (n, spin) == (3, two_spin)
                 for h in range(5)
                 if w+lw+2*scalar[10][h, h] == two_weak
                 and y+ly+6*scalar[11][h, h] == six_hyper)


def tensor_map(keys, target, quark=None, lepton=None, scalar=None):
    positions = {key: i for i, key in enumerate(target)}
    entries = defaultdict(lambda: s.Integer(0))
    for axis, matrix in enumerate((quark, lepton, scalar)):
        if matrix is None:
            continue
        cols = columns(matrix)
        for j, key in enumerate(keys):
            for destination, value in cols.get(key[axis], ()):
                changed = list(key)
                changed[axis] = destination
                assert tuple(changed) in positions, ('omitted tensor image', key, changed)
                entries[positions[tuple(changed)], j] += value
    return s.ImmutableSparseMatrix(len(target), len(keys), entries)


def metric_on(keys, create):
    cols = columns(target_module(create)['metric'])
    positions = {key: i for i, key in enumerate(keys)}
    entries = {}
    for j, (a, mask, h) in enumerate(keys):
        for destination, value in cols.get(a, ()):
            assert (destination, mask, h) in positions
            entries[positions[destination, mask, h], j] = value
    return s.ImmutableSparseMatrix(len(keys), len(keys), entries)


def hamiltonian_on(keys, create):
    q = target_module(create)
    lepton = lepton_operators()
    qc = tuple(columns(t) for t in q['transitions'])
    lc = tuple(columns(t) for t in lepton['transitions'])
    positions = {key: i for i, key in enumerate(keys)}
    entries = defaultdict(lambda: s.Integer(0))
    for j, (a, mask, h) in enumerate(keys):
        if h:
            entries[j, j] += 1
        branches = tuple((r, r+1) for r in range(4)) if h == 0 else (((h+1) % 4, 0),)
        for r, target_h in branches:
            for destination, value in qc[r].get(a, ()):
                assert (destination, mask, target_h) in positions
                entries[positions[destination, mask, target_h], j] += value
            for destination, value in lc[r].get(mask, ()):
                assert (a, destination, target_h) in positions
                entries[positions[a, destination, target_h], j] += value
    return s.ImmutableSparseMatrix(len(keys), len(keys), entries)


@lru_cache(maxsize=None)
def target(family):
    create, mode, weak, hyper = CHANNELS[family]
    q = target_module(create)
    lepton = lepton_operators()
    scalar = higgs_memory()[1]
    keys = weight_basis(create, 1, weak, hyper)
    weak_map = tensor_map(keys, weight_basis(create, 1, weak+2, hyper),
                          quark=q['currents'][0]+s.I*q['currents'][1],
                          lepton=lepton['currents'][0]+s.I*lepton['currents'][1],
                          scalar=scalar[8]+s.I*scalar[9])
    spin_map = tensor_map(keys, weight_basis(create, 3, weak, hyper),
                          lepton=lepton['spin'][0]+s.I*lepton['spin'][1])
    basis, coordinates = coordinate_kernel(weak_map.col_join(spin_map))
    expected = (len(keys)-len(weight_basis(create, 1, weak+2, hyper))
                -len(weight_basis(create, 3, weak, hyper))
                +len(weight_basis(create, 3, weak+2, hyper)))
    assert basis.cols == expected
    full_metric = metric_on(keys, create)
    h = hamiltonian_on(keys, create)
    assert full_metric*h == h.H*full_metric
    image = h*basis
    reduced_h = s.ImmutableSparseMatrix(image.extract(coordinates, range(basis.cols)))
    assert basis*reduced_h == image
    gram = s.ImmutableSparseMatrix(basis.H*full_metric*basis)
    assert gram*reduced_h == reduced_h.H*gram
    return dict(keys=keys, basis=basis, coordinates=coordinates, H=reduced_h, gram=gram)


@lru_cache(maxsize=None)
def transfer(family):
    create, mode, weak, hyper = CHANNELS[family]
    q = target_module(create)
    sb, tb = source(1), target(family)
    native_mode = q['native_component'][mode]
    assert native_mode % 15 in FAMILIES[family]
    assert particle_hole(car.monomial(((native_mode, False),))) == car.monomial(((native_mode, create),))
    positions = {key: i for i, key in enumerate(tb['keys'])}
    entries = {}
    cols = columns(q['maps'][mode])
    for j, (a, mask, h) in enumerate(sb['keys']):
        for destination, value in cols.get(a, ()):
            assert (destination, mask, h) in positions
            entries[positions[destination, mask, h], j] = value
    raw = s.ImmutableSparseMatrix(len(tb['keys']), len(sb['keys']), entries)
    image = raw*sb['basis']
    matrix = s.ImmutableSparseMatrix(image.extract(tb['coordinates'], range(sb['basis'].cols)))
    assert tb['basis']*matrix == image
    adjoint = s.ImmutableSparseMatrix(sb['gram'].inv()*matrix.H*tb['gram'])
    assert sb['gram']*adjoint == matrix.H*tb['gram']
    return dict(H=tb['H'], gram=tb['gram'], map=matrix, adjoint=adjoint,
                full_weight_dimension=len(tb['keys']),
                complete_selected_target_dimension=tb['basis'].cols,
                native_removed_mode=native_mode)


if __name__ == '__main__':
    energy = s.Symbol('E')
    for family in CHANNELS:
        data = transfer(family)
        print(family, data['native_removed_mode'], data['full_weight_dimension'],
              data['H'].rows, s.factor(data['H'].charpoly(energy).as_expr()), flush=True)
