"""Capacity-independent strong-source stability; no occupation cutoff."""
import json
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

RER = Path(__file__).resolve().parents[2]


@lru_cache(None)
def inputs():
    root = RER / 'model/openwave_common_hyperfine_spectrum'
    values = json.loads((root / 'inputs.json').read_text())
    receipt = json.loads((root / 'results.json').read_text())
    p = {k: F(v) for k, v in values['declared_rational_parameters'].items()}
    kappa = F(receipt['common_length_attachment']['kappa_required_interval'][1])
    return p, kappa, F(values['electron_unit_GeV'])


@lru_cache(None)
def declared_kappa():
    receipt=json.loads((RER/'model/openwave_common_hyperfine_spectrum/results.json').read_text())
    return F(receipt['common_length_attachment']['kappa'])


@lru_cache(None)
def uniform_bound(vacuum_compatible=False):
    p, kappa, unit = inputs()
    c = max(F(0), kappa / 4 - p['U'])
    h = p['pair_triplet']
    mass = p['mbar'] - abs(p['delta']) / 2
    alpha, beta = h / 4 - c / 3, p['A'] / 9 - h / 4
    gamma = mass - h / 2 - c
    one = mass + p['A'] / 9 - F(4, 3) * c
    if vacuum_compatible:
        removed = p['spin_gradient'] + p['eta']
        alpha -= removed / 4
        gamma -= removed / 2 + p['t']
        one -= F(3, 4) * removed + p['t']
    assert alpha > 0 and beta > 0 and one > 0 and 2 * alpha + gamma > 0
    return {'alpha': alpha, 'beta': beta, 'gamma': gamma,
            'N1_floor': one, 'N_ge_2_per_particle_floor': 2 * alpha + gamma,
            'uniform_floor': min(one, 2 * alpha + gamma), 'unit': unit,
            'kappa_upper': kappa, 'parameters': p}


def direct_count_bound(nq, na, vacuum_compatible=False):
    """The exact polynomial lower bound, valid for every nonnegative count."""
    b = uniform_bound(vacuum_compatible)
    n, d = nq + na, nq - na
    return (b['alpha'] * n * n + b['beta'] * d * d
            + b['gamma'] * n + b['parameters']['pair_triplet'] * abs(d) / 2)


def joint_stability(three_body=F(0)):
    """C completion and occupancy-gated nuclear term on the same full carrier.

    Each triple projector requires N_x=3 on all its cells. A cell belongs to
    at most three triples, so their sum is bounded by the number of N_x=3
    cells. Global min(Nq,Na)<=N/2 supplies the separate C-completion bound.
    """
    z=uniform_bound(True);p=z['parameters'];c=F(3,4)*(p['Jw']+p['eta'])
    from .parameters import degree_coefficients
    dU,dm=degree_coefficients()
    # Discard the positive number shift. The negative colour shift uses the
    # same capacity-independent Casimir bound C²<=N(N+3)/3.
    z=dict(z);z['alpha']-=max(F(0),-dU)/3;z['gamma']-=max(F(0),-dU)
    z['N1_floor']-=F(4,3)*max(F(0),-dU)
    candidates={
        'one':z['N1_floor']-c/2,
        'two':2*z['alpha']+z['gamma']-c/2,
        'three':3*z['alpha']+z['gamma']+z['beta']/3+p['pair_triplet']/6-c/2-three_body/3,
        'four_or_more':4*z['alpha']+z['gamma']-c/2,
    }
    floor=min(candidates.values());assert floor>0
    for nq in range(49):
        for na in range(49):
            n=nq+na
            actual=direct_count_bound(nq,na,True)-max(F(0),-dU)*n*(n+3)/3-c*n/2-(three_body if n==3 else 0)
            assert actual>=floor*n
    return {'uniform_floor':floor,'cases':candidates,'C_cost_per_pair':c,
            'three_body_coefficient':three_body,'domain':'Every local occupation, arbitrary family capacity and every open-chain length.'}


if __name__ == '__main__':
    for changed in (False, True):
        b = uniform_bound(changed)
        for nq in range(49):
            for na in range(49):
                assert direct_count_bound(nq, na, changed) >= b['uniform_floor'] * (nq + na)
        print('vacuum_compatible', changed, 'floor', str(b['uniform_floor']), float(b['uniform_floor']))
