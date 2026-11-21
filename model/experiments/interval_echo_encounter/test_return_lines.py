"""Independent checks of the distance, moment and sharp-read arguments."""
import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse import csgraph

from model.candidate_laws.contrast_memory_spherical_v1.sectors import MultisetSector, PairSector
from model.executable_model.carrier import cell, federation
from . import return_lines as lines


def test_distance_witness_on_full_pair_graph():
    board = federation(0)
    pairs = PairSector(board)
    adjacency = pairs.exchange().copy()
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    adjacency.data[:] = 1
    home = pairs.index(0,1)
    dist = csgraph.shortest_path(adjacency,directed=False,indices=home,unweighted=True)
    bound = lines.line_bound(board,(0,1))
    assert dist[pairs.index(*bound['target'])]>=bound["moves_at_least"]
    assert bound["moves_at_least"]==27
    # At each exact-distance shell the first nonzero moment has the same sign.
    v = np.zeros(pairs.P,dtype=np.int64)
    v[home] = 1
    twice_h = (2*pairs.hamiltonian(8)).astype(np.int64)
    for j in range(1,9):
        v = twice_h @ v
        np.testing.assert_array_equal(v[dist>j],0)
        assert np.min((-1)**j*v[dist==j])>0


def test_finite_reference_capacity():
    assert lines.counter(17)==list(range(18))
    with pytest.raises(ValueError):
        lines.counter(32)


def test_finite_read_error_budget():
    cfg = lines.json.loads((lines.HERE/"return_protocol.json").read_text())
    name, joint = next(lines.acquire_laws(cfg))
    assert name == "native"
    admit = joint[:, :2].sum(axis=1)
    assert np.min(admit) > .000060
    delta = 2*16*2**-32
    assert 2*delta/(min(admit)-delta)<.00025
    assert 17*.00025 < .15795626-.01309490


def test_local_survival_against_full_qutrit_generator_and_energy_zero():
    board = cell()
    full = MultisetSector(board,2)
    hard = PairSector(board)
    e = hard.index(0,1)
    sharp = np.zeros(len(full.states))
    sharp[full.index[(1,1)+(0,)*10]] = 1
    H = full.hamiltonian(8)
    phase = expm(-1j*.21*H)
    U = expm(-1j*.21*hard.hamiltonian(8).toarray())
    for n in (1,2,3):
        p = abs(sharp @ np.linalg.matrix_power(phase,n) @ sharp)**2
        assert abs(p-abs(np.linalg.matrix_power(U,n)[e,e])**2)<1e-12
        shifted = expm(-1j*.21*n*(H+13*np.eye(len(full.states))))
        assert abs(p-abs(sharp @ shifted @ sharp)**2)<1e-12


def test_frequency_count_and_mass_outside_levels():
    energies = np.array([0,.51,1.63,2.75,4.23])
    weights = np.array([.54,.19,.15,.10,.02])
    shots = 10**12
    for n in range(1,6):
        p = weights[:n]/sum(weights[:n])
        amplitudes = np.exp(-1j*np.outer(np.arange(17),energies[:n])) @ p
        back = np.rint((shots//2)*abs(amplitudes)**2).astype(np.int64)
        counts = np.stack((back,shots//2-back,np.zeros(17,dtype=np.int64),
                           np.full(17,shots//2)),axis=1)
        measured = lines.read(counts,shots,19.)
        assert measured["energies_at_least"]<=n
        for d,bound in measured["outside_weight_lower"].items():
            assert bound<=sum(p[int(d):])+1e-14
        if n<=2:
            assert measured["first_three_excess_lower"]<=0
        if n==5:
            assert measured["positive_gaps_at_least"]>1
            assert measured["outside_weight_lower"]["2"]>0


@pytest.mark.parametrize("mutation",["missing","unsigned","oversized","nonrepeat","noherald","negative"])
def test_record_reader_refuses_invalid_evidence(mutation):
    counts = np.tile(np.array([100,0,0,100],np.int64),(6,1))
    if mutation=="missing":
        counts = counts[:,:3]
    elif mutation=="unsigned":
        counts = counts.astype(np.uint64)
    elif mutation=="oversized":
        counts[2,0] += 1
    elif mutation=="nonrepeat":
        counts[0] += [-1,1,0,0]
    elif mutation=="noherald":
        counts[3] = [0,0,0,200]
    else:
        counts[-1,1:3] = [-3,3]
    with pytest.raises(ValueError):
        lines.read(counts,200,15)
