"""Complete projected triple-CAR hopping; no producer code imported."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt, lcm
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy import sparse
import sympy as sy

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent/'hadrons/shared_auxiliary/LOCAL_OPERATORS.json'


def annihilate(mode, word):
    if not (word >> mode) & 1:
        return None, 0
    return word ^ (1 << mode), (-1 if (word & ((1 << mode)-1)).bit_count() % 2 else 1)


def create(mode, word):
    if (word >> mode) & 1:
        return None, 0
    return word | (1 << mode), (-1 if (word & ((1 << mode)-1)).bit_count() % 2 else 1)


def raw_onebody(i, j):
    rows, cols, values = [], [], []
    for word in range(4096):
        middle, a = annihilate(j, word)
        if middle is None:
            continue
        out, b = create(i, middle)
        if out is not None:
            rows.append(out); cols.append(word); values.append(a*b)
    return sparse.csr_matrix((values, (rows, cols)), shape=(4096, 4096), dtype=np.int64)


def rational_scaled(A, denominator):
    rows, cols, vals = [], [], []
    for (i, j), value in sy.SparseMatrix(A).todok().items():
        v = value*denominator
        assert v.q == 1
        rows.append(int(i)); cols.append(int(j)); vals.append(int(v))
    return sparse.csr_matrix((vals, (rows, cols)), shape=A.shape, dtype=np.int64)


def sqbounds(q):
    scale = 1 << 100
    n = isqrt(q.numerator*scale*scale//q.denominator)
    a, b = F(n, scale), F(n+1, scale)
    assert a*a <= q <= b*b
    return a, a if a*a == q else b


def exact_control():
    # Complete symmetric SU4 baryon chain, including both colour-flux states.
    # Even and odd reflection blocks have diagonals (0,Delta -/+ 2t).
    t = F(1, 8)
    out = []
    for kap in (F(1, 2), F(1), F(2)):
        Delta = 4*kap/3
        eigen = []
        for middle in (Delta-2*t, Delta+2*t):
            sl, sh = sqbounds(middle*middle+12*t*t)
            eigen.append(((middle-sh)/2, (middle-sl)/2))
        exact_hop = ((eigen[1][0]-eigen[0][1])/2, (eigen[1][1]-eigen[0][0])/2)
        third = 6*t**3/Delta**2
        assert exact_hop[0] > 0
        out.append({'kappa': str(kap), 'Delta': str(Delta),
                    'exact_lowest_even_odd_intervals': [[str(v) for v in pair] for pair in eigen],
                    'exact_half_splitting_interval': list(map(str, exact_hop)),
                    'third_order_hopping': str(third),
                    'exact_to_third_order_interval': [float(x/third) for x in exact_hop],
                    'control_mass': 0, 'control_auxiliary_coupling': 0,
                    'actual_shared_vacuum_identification': False})
    return out


def main():
    source = json.loads(INPUT.read_text())['sectors']
    rows, cols, vals, counts = [], [], [], []
    offset = 0
    for sector in source:
        block = np.array(sector['R'], dtype=np.int64)
        for i, j in zip(*np.nonzero(block)):
            rows.append(sector['words'][i]); cols.append(offset+j); vals.append(block[i, j])
        counts.extend([sector['N']]*sector['dimension'])
        offset += sector['dimension']
    assert offset == 92
    R = sparse.csr_matrix((vals, (rows, cols)), shape=(4096, 92), dtype=np.int64)
    G = (R.T@R).tocsr()
    Gsym = sy.SparseMatrix(G.toarray().tolist())
    Gi = sy.SparseMatrix(Gsym.inv(method='DM'))
    assert Gsym*Gi == sy.eye(92)
    # Independently check that the consumed frame is indeed colour-singlet.
    colour_checks = 0
    for a in range(3):
        for b in range(3):
            if a != b:
                generator = sum((raw_onebody(3*f+a, 3*f+b) for f in range(4)),
                                sparse.csr_matrix((4096, 4096), dtype=np.int64))
                diff = generator@R; diff.eliminate_zeros()
                assert diff.nnz == 0
                colour_checks += 1
    for a, b in ((0, 1), (1, 2)):
        difference = np.array([sum((w >> (3*f+a) & 1)-(w >> (3*f+b) & 1) for f in range(4))
                               for w in range(4096)], dtype=np.int64)
        diff = sparse.diags(difference)@R; diff.eliminate_zeros()
        assert diff.nnz == 0
        colour_checks += 1
    print('exact local frame and eight colour generators passed', flush=True)
    templates, triples = {}, []
    nonzero = 0
    for modes in combinations(range(12), 3):
        rr, cc, vv = [], [], []
        for word in range(4096):
            state, sign = word, 1
            for mode in modes:
                state, step_sign = annihilate(mode, state)
                if state is None:
                    break
                sign *= step_sign
            if state is not None:
                rr.append(state); cc.append(word); vv.append(sign)
        raw = sparse.csr_matrix((vv, (rr, cc)), shape=(4096, 4096), dtype=np.int64)
        form = (R.T@raw@R).tocsr(); form.eliminate_zeros()
        if form.nnz == 0:
            assert len(set(m % 3 for m in modes)) != 3
            triples.append({'modes': modes, 'projected_zero': True})
            continue
        assert len(set(m % 3 for m in modes)) == 3
        nonzero += 1
        key = tuple(sorted(m//3 for m in modes))
        if key not in templates:
            templates[key] = {'form': form, 'multiplicity': 0, 'representative': modes}
        reference = templates[key]['form']
        minus = form-reference; minus.eliminate_zeros()
        plus = form+reference; plus.eliminate_zeros()
        assert minus.nnz == 0 or plus.nnz == 0
        sign = 1 if minus.nnz == 0 else -1
        templates[key]['multiplicity'] += 1
        triples.append({'modes': modes, 'projected_zero': False,
                        'symmetric_flavours': key, 'relative_sign': sign})
    assert len(triples) == 220 and nonzero == 64 and len(templates) == 20
    print('all220 triples checked;64 survive in20 exact flavour templates', flush=True)
    denominators = []
    for info in templates.values():
        form = sy.SparseMatrix(info['form'].toarray().tolist())
        info['lower'] = Gi*form
        info['raise'] = Gi*form.T
        for name in ('lower', 'raise'):
            denominators.extend(int(v.q) for v in sy.SparseMatrix(info[name]).todok().values())
    denom = lcm(*denominators)
    parity = sparse.diags([(-1)**n for n in counts], dtype=np.int64, format='csr')
    transfer = sparse.csr_matrix((8464, 8464), dtype=np.int64)
    local_exports = []
    for key, info in templates.items():
        low = rational_scaled(info['lower'], denom)
        high = rational_scaled(info['raise'], denom)
        assert np.max(np.abs(low.data), initial=0) < 100000
        assert np.max(np.abs(high.data), initial=0) < 100000
        transfer += 6*info['multiplicity']*(sparse.kron(high@parity, low, format='csr')
                                           + sparse.kron(parity@low, high, format='csr'))
        local_exports.append({'flavours': key, 'representative_modes': info['representative'],
            'multiplicity': info['multiplicity'], 'operator_denominator': denom,
            'lower_integer_matrix': low.toarray().tolist(), 'raise_integer_matrix': high.toarray().tolist()})
    transfer.eliminate_zeros()
    gram = sparse.kron(G, G, format='csr')
    form = gram@transfer
    defect = form-form.T; defect.eliminate_zeros()
    assert defect.nnz == 0
    total = np.array([a+b for a in counts for b in counts], dtype=np.int64)
    left = np.repeat(np.array(counts, dtype=np.int64), 92)
    coo = transfer.tocoo()
    assert np.all(total[coo.row] == total[coo.col])
    assert np.all(np.abs(left[coo.row]-left[coo.col]) == 3)
    # A finite current operator is retained; it changes baryon placement by one.
    numerator_velocity = sparse.csr_matrix(
        (coo.data*((left[coo.row]-left[coo.col])//3), (coo.row, coo.col)), shape=transfer.shape)
    numerator_velocity.eliminate_zeros()
    assert numerator_velocity.nnz == transfer.nnz
    assert np.array_equal(np.sort(np.abs(numerator_velocity.data)), np.sort(np.abs(transfer.data)))
    generators = {}
    for key in ('twice_Sz', 'twice_Iz', 'Splus', 'Iplus'):
        block_form = sy.diag(*[sy.Matrix(r[key]) for r in source])
        local_action = Gi*block_form
        gd = lcm(*(int(v.q) for v in sy.SparseMatrix(local_action).todok().values()))
        integer = rational_scaled(local_action, gd)
        global_generator = sparse.kron(integer, sparse.eye(92, dtype=np.int64), format='csr')+sparse.kron(sparse.eye(92, dtype=np.int64), integer, format='csr')
        commutator = global_generator@transfer-transfer@global_generator
        commutator.eliminate_zeros()
        assert commutator.nnz == 0
        generators[key] = {'exact_commutator_zero': True, 'generator_denominator': gd}
    sectors = []
    for N in range(0, 25, 3):
        ix = np.flatnonzero(total == N)
        block = transfer[ix, :][:, ix]
        bound = F(int(np.max(np.asarray(abs(block).sum(axis=1)), initial=0)), denom**2)
        if N in (3, 21):
            check = block@block-36*denom**4*sparse.eye(len(ix), dtype=np.int64)
            check.eliminate_zeros()
            assert check.nnz == 0
        sectors.append({'N': N, 'dimension': len(ix), 'nonzero_transfer_entries': block.nnz,
                        'rigorous_physical_norm_upper': str(bound),
                        'exact_norm_six': N in (3, 21)})
    assert sum(r['dimension'] for r in sectors) == 8464
    sparse.save_npz(HERE/'TRIPLE_TRANSFER_NUMERATOR.npz', transfer)
    sparse.save_npz(HERE/'PRODUCT_GRAM.npz', gram)
    (HERE/'LOCAL_TRIPLES.json').write_text(json.dumps({'denominator': denom,
        'all220_triples': triples, 'templates': local_exports}, indent=2)+'\n')
    result = {'schema': 'oph.baryon_mobility.complete_projection.v1',
        'source_pins': {str(INPUT.relative_to(HERE.parent)): hashlib.sha256(INPUT.read_bytes()).hexdigest()},
        'local_singlet_dimension': 92, 'complete_product_dimension': 8464,
        'independent_colour_generator_checks': colour_checks, 'all_triples': 220,
        'nonzero_projected_triples': 64, 'symmetric_flavour_templates': 20,
        'transfer_numerator_denominator': denom**2,
        'dimensionless_operator': 'P K^3 P, K=sum_mu(cL_mu^dagger cR_mu+h.c.)',
        'third_order_Hamiltonian_coefficient': '-t^3/(4*kappa/3)^2',
        'fermion_parity_retained': True, 'exact_Gram_self_adjoint': True,
        'number_conserved_and_left_number_changes_by_three': True,
        'symmetry_checks': generators, 'number_sectors': sectors,
        'massless_zero_auxiliary_control': exact_control(),
        'full_H_eff3_claim': False, 'full625_spectral_transfer': False,
        'physical_halo_or_lifetime_claim': False,
        'all_other_channels_retained_in_full_source_obligation': True}
    (HERE/'RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'denominator': denom**2, 'sectors': sectors,
                      'control': result['massless_zero_auxiliary_control']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
