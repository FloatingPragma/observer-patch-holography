"""Full Cartan carrier with external fundamental/antifundamental probes."""
from itertools import combinations, product
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent


def fermion_hop(word, dest, source):
    if not (word >> source) & 1:
        return None, 0
    sign = (-1)**((word & ((1 << source)-1)).bit_count())
    word ^= 1 << source
    if (word >> dest) & 1:
        return None, 0
    sign *= (-1)**((word & ((1 << dest)-1)).bit_count())
    return word | (1 << dest), sign


def current(state, sites, has_q, has_a, a, b):
    word, q, anti = state
    out = {}
    for x in sites:
        w, sign = fermion_hop(word, 3*x+a, 3*x+b)
        if w is not None:
            z = w, q, anti
            out[z] = out.get(z, 0)+sign
    if has_q and q == b:
        z = word, a, anti
        out[z] = out.get(z, 0)+1
    if has_a and anti == a:
        z = word, q, b
        out[z] = out.get(z, 0)-1
    return {z: c for z, c in out.items() if c}


def carrier(k):
    states = []
    for q, anti in product(range(3), repeat=2):
        counts = [k-int(q == c)+int(anti == c) for c in range(3)]
        if min(counts) < 0 or max(counts) > 4:
            continue
        choices = [list(combinations(range(4), n)) for n in counts]
        for xs3 in product(*choices):
            word = sum(1 << (3*x+c) for c, xs in enumerate(xs3) for x in xs)
            states.append((word, q, anti))
    return states


def casimir6(states, sites, has_q, has_a):
    ix = {s: i for i, s in enumerate(states)}
    rows, cols, values = [], [], []
    mask = sum(7 << (3*x) for x in sites)
    for j, state in enumerate(states):
        n = (state[0] & mask).bit_count()+int(has_q)-int(has_a)
        column = {j: -n*n}
        for a, b in product(range(3), repeat=2):
            for middle, c1 in current(state, sites, has_q, has_a, b, a).items():
                for target, c2 in current(middle, sites, has_q, has_a, a, b).items():
                    i = ix[target]
                    column[i] = column.get(i, 0)+3*c1*c2
        for i, value in column.items():
            if value:
                rows.append(i); cols.append(j); values.append(value)
    result = sp.csr_matrix((values, (rows, cols)), shape=(len(states), len(states)), dtype=np.int64)
    assert (result-result.T).nnz == 0
    return result


def base_operators(states):
    ix = {s: i for i, s in enumerate(states)}
    rows, cols, values = [], [], []
    for j, (word, q, anti) in enumerate(states):
        for x, c in product(range(3), range(3)):
            for u, v in [(x, x+1), (x+1, x)]:
                out, sign = fermion_hop(word, 3*u+c, 3*v+c)
                if out is not None:
                    rows.append(ix[(out, q, anti)]); cols.append(j); values.append(-sign)
    kinetic = sp.csr_matrix((values, (rows, cols)), shape=(len(states), len(states)), dtype=np.int64)
    mass = sp.diags([sum((-1)**x*((word >> (3*x)) & 7).bit_count() for x in range(4))
                     for word, _, _ in states], dtype=np.int64).tocsr()
    colour = casimir6(states, range(4), True, True)
    assert (kinetic-kinetic.T).nnz == 0
    for op in [kinetic, mass]:
        assert (colour@op-op@colour).nnz == 0
    return kinetic, mass, colour


def main():
    paths = [HERE/'STATIC_COLOUR_PLAN.md', HERE/'hadrons/RAW_COEFFICIENT_INPUTS.json',
             HERE/'hadrons/STAGGERED_CHAIN4_RESULT.json', HERE/'hadrons/CHAIN4_MANIFEST.json']
    pins = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    inp = json.loads(paths[1].read_text())
    reference = json.loads(paths[2].read_text())
    mass = float(F(inp['432']['declared_light_midpoint']))
    # Every nontrivial SU3 representation has C2>=4/3. Even this weaker
    # bound suffices: H>=-6/4-6m; penalty100 gives >131.8. A physical
    # localized singlet trial has energy<=8+6m for our whole parameter set.
    assert 100*4/3-1.5-6*mass > 8+6*mass
    output = []
    for k in [0, 2]:
        states = carrier(k)
        kinetic, staggered, colour6 = base_operators(states)
        print('carrier', 3*k, len(states), flush=True)
        for d in range(4):
            electric6 = sum((casimir6(states, range(cut), 0 < cut, d < cut)
                             for cut in [1, 2, 3]), sp.csr_matrix(colour6.shape, dtype=np.int64))
            assert (colour6@electric6-electric6@colour6).nnz == 0
            for kap in [.5, 1., 2.]:
                H = kinetic.astype(float)/8+mass*staggered+kap*electric6/6+100*colour6/6
                if len(states) < 5:
                    energies, vectors = np.linalg.eigh(H.toarray())
                    energies, vectors = energies[:3], vectors[:, :3]
                else:
                    v0 = np.random.default_rng(50701).normal(size=len(states))
                    energies, vectors = sla.eigsh(H, k=3, which='SA', tol=1e-12, maxiter=100000, v0=v0)
                    order = np.argsort(energies)
                    energies, vectors = energies[order], vectors[:, order]
                residuals = np.linalg.norm(H@vectors-vectors*energies, axis=0)
                colour = np.sum(vectors*(colour6@vectors), axis=0)/6
                assert residuals[0] < 1e-8 and abs(colour[0]) < 1e-9
                if k == 0:
                    assert abs(energies[0]-4*kap*d/3) < 1e-10
                    e0 = 0.
                else:
                    ref = next(r for r in reference['results'] if float(F(r['kappa'])) == kap)
                    vac = next(r for r in ref['complete_sector_records'] if r['N'] == 6)
                    e0 = sum(float(F(z)) for z in vac['ground_interval_exact'])/2
                    if d == 0:
                        # Full coincident probe carrier includes adjoint pair.
                        # The separately prepared singlet-pair reference is E0.
                        assert energies[0] <= e0+1e-8
                row = dict(fermion_number=3*k, dimension=len(states), kappa=kap,
                           separation=d, low_ritz_values=energies.tolist(), residuals=residuals.tolist(),
                           total_colour_casimirs=colour.tolist(), reference_vacuum=e0,
                           potential_from_original_vacuum=float(energies[0]-e0), bare_cost=4*kap*d/3,
                           rigorous_ground_interval=False)
                output.append(row)
                print(3*k, d, kap, row['potential_from_original_vacuum'], flush=True)
    assert pins == {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result = dict(schema='oph.obstruction_repair.static_colour_probe.v1',
                  complete_Cartan_carrier=True, exact_colour_commutators=True,
                  nonsinglet_penalty_exclusion_proved=True, fixed_fermion_number=6,
                  source_changed_by_declared_static_probes=True, rows=output,
                  infinitely_separated_probe_limit=False, physical_string_tension_solved=False,
                  common_spatial_clock_transfer=False, native_probe_program_executed=False, pins=pins)
    (HERE/'STATIC_COLOUR_RESULT.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
