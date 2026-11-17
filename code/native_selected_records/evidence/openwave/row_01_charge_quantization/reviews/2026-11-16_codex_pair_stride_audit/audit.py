"""Same-group independent pair-readout audit; imports no producer modules."""
from pathlib import Path
from itertools import product
import gzip
import hashlib
import json
import math
import random
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROW = HERE.parent.parent
SOURCE = ROW / 'work/2026-11-16_codex_common_refinement_family'
PAIR = ROW / 'work/2026-11-16_codex_pair_stride_readout'
ETA, DELTA, STRIDE = .01, 1e-5, 8
SCALE, OLD = 511 * math.sqrt(5) / 2, 1024.4508487023415


def unpack(x):
    return np.array(x['real']) + 1j*np.array(x['imag']) if isinstance(x, dict) else np.array(x)


def close(x, y, tolerance=2e-10):
    assert np.shape(x) == np.shape(y), (np.shape(x), np.shape(y))
    residual = float(np.max(np.abs(np.asarray(x)-np.asarray(y)))) if np.size(x) else 0.
    assert residual < tolerance, residual
    return residual


def principal(record, eta, delta):
    for axis in ('x', 'y'):
        p, q = record[axis+'_plus'], record[axis+'_minus']
        assert math.isfinite(p) and math.isfinite(q) and 0 <= min(p,q) <= max(p,q) <= 1
        assert abs(p+q-1) < 1e-12
    r = 2*record['x_plus']-1
    lo = max(-1., min(1., (r-eta-2*delta)/(1-eta)))
    hi = max(-1., min(1., (r+eta+2*delta)/(1-eta)))
    return math.acos(hi), math.acos(lo)


def intervals(first, last, eta=ETA, delta=DELTA, stride=STRIDE):
    """Independent enumeration via monotone half-waves, not signed windings."""
    i0, i1 = principal(first, eta, delta), principal(last, eta, delta)
    result = []
    for segment in range(stride):
        if segment % 2 == 0:
            left, right = ((segment*math.pi+x)/stride for x in i1)
        else:
            left, right = (((segment+1)*math.pi-x)/stride for x in reversed(i1))
        left, right = max(left, i0[0]), min(right, i0[1])
        if left <= right:
            if result and left <= result[-1][1]+1e-14:
                result[-1][1] = max(right, result[-1][1])
            else:
                result.append([left, right])
    return result


def contains(parts, value):
    return any(a-1e-10 <= value <= b+1e-10 for a,b in parts)


def contained(parts, outer):
    return all(any(c-1e-10 <= a <= b <= d+1e-10 for c,d in outer) for a,b in parts)


def probabilities(z):
    return {axis+'_'+side: float((1+sign*part)/2)
            for axis,part in (('x',z.real),('y',z.imag))
            for side,sign in (('plus',1),('minus',-1))}


def check_records(a, b):
    return max(close(a[key], b[key]) for key in a)


def shifted(record, dx, dy):
    result = dict(record)
    for axis, error in (('x',dx), ('y',dy)):
        result[axis+'_plus'] += error
        result[axis+'_minus'] -= error
    return result


def clock(k):
    """Literal floor/ceil endpoint averaging, independent expression of the law."""
    reversal = np.zeros((33,33)); columns = 0
    for position in range(3):
        reversal[11*position+10,11*position+10] = 1
        for edge,(u,v) in enumerate(((0,1),(0,1),(0,2),(0,2),(1,2))):
            for coin in (0,1):
                q = [int(position==j) for j in range(3)]
                total, old = q[u]+q[v], q[u]
                q[u] = (total+coin)//2
                q[v] = total-q[u]
                current = old-q[u]
                target = q.index(1)
                next_coin = coin ^ int(current != 0)
                reversal[11*target+2*edge+next_coin,11*position+2*edge+coin] = 1
                columns += 1
    epsilon = 2.**-k
    choice = np.r_[np.full(10,epsilon/math.sqrt(10)), math.sqrt(1-epsilon*epsilon)]
    encoding = np.kron(np.eye(3),choice[:,None])
    operator = reversal@(2*encoding@encoding.T-np.eye(33))
    close(operator.T@operator,np.eye(33)); close(reversal@reversal,np.eye(33))
    return operator, columns


def run():
    begin = time.monotonic()
    data = json.loads((PAIR/'OUTPUT.json').read_text())
    bindings = json.loads((HERE/'INPUTS.json').read_text())
    for path,digest in bindings['files'].items():
        assert hashlib.sha256((ROW/path).read_bytes()).hexdigest() == digest, path
    assert len(data['rows']) == 112 and data['passed']
    counts = {'source_columns':0, 'source_cases':0, 'pair_moments':0,
              'corner_cases':0, 'interior_perturbations':0, 'controls':0,
              'full_chart_plants':0, 'false_principal_branch_rejections':0}
    rng = random.Random(72316008)
    extrema = {'pair_off_band_max':0., 'minimum_density_eigenvalue':1.,
               'minimum_positive_f_weight':1., 'maximum_negative_f_weight':0.,
               'maximum_probability_residual':0., 'maximum_certificate_residual':0.}
    final_bounds, corner_bounds, rows = [], [], []
    by_case = {(r['k'], *(r['case'][x] for x in ('M','turns','direction','law'))):r for r in data['rows']}
    controls = {(tuple(c['case'][x] for x in ('M','turns','direction','law'))):c
                for c in data['controls'] if 'case' in c}
    for k in range(7):
        base = json.loads(gzip.decompress((SOURCE/f'k{k}/EVIDENCE.json.gz').read_bytes()))['base_execution']
        W, n = clock(k); counts['source_columns'] += n
        close(W, unpack(base['clock']['W']))
        step = np.linalg.matrix_power(W,2**k)
        later = np.linalg.matrix_power(step,STRIDE)
        eigenvalues, eigenvectors = np.linalg.eig(W.astype(complex))
        theta = 2*math.asin(2.**-k/math.sqrt(5))
        desired = []
        for sign in (1,-1):
            j = int(np.argmin(abs(eigenvalues-np.exp(1j*sign*theta))))
            assert abs(eigenvalues[j]-np.exp(1j*sign*theta)) < 1e-11
            desired.append(eigenvectors[:,j]/np.linalg.norm(eigenvectors[:,j]))
        frequency = 2*theta*2**k
        for case in base['cases']:
            identity = tuple(case[x] for x in ('M','turns','direction','law'))
            filed = by_case[(k,*identity)]
            rho = unpack(case['clock_density'])
            close(rho,rho.conj().T); close(np.trace(rho),1)
            mineig = float(np.linalg.eigvalsh(rho).min())
            assert mineig >= -1e-10
            extrema['minimum_density_eigenvalue'] = min(extrema['minimum_density_eigenvalue'],mineig)
            weights = [float(np.vdot(v,rho@v).real) for v in desired]
            off = 1-sum(w*w for w in weights)
            assert 0 <= off <= ETA
            extrema['pair_off_band_max'] = max(extrema['pair_off_band_max'],off)
            extrema['minimum_positive_f_weight'] = min(extrema['minimum_positive_f_weight'],weights[0])
            extrema['maximum_negative_f_weight'] = max(extrema['maximum_negative_f_weight'],weights[1])
            cert = filed['numerical_source_certificate']
            residual = close([off,weights[0]**2,weights[1]**2],
                             [cert['eta'],cert['positive_pair_weight'],cert['negative_pair_weight']])
            extrema['maximum_certificate_residual'] = max(extrema['maximum_certificate_residual'],residual)
            # Product trace identity is checked against producer's full 1089x1089 calculation;
            # actual two-unit native admissibility is separately manifest-bound below.
            moments = {1:complex(np.trace(rho@step)**2), STRIDE:complex(np.trace(rho@later)**2)}
            records = {s:probabilities(z) for s,z in moments.items()}
            for s in (1,STRIDE):
                residual = check_records(records[s],filed['pair_records'][str(s)])
                extrema['maximum_probability_residual'] = max(extrema['maximum_probability_residual'],residual)
                counts['pair_moments'] += 1
            for delta,key in ((0,'nominal_reader'),(DELTA,'bounded_reader'),(2*DELTA,'worst_corner_envelope_reader')):
                out = intervals(records[1],records[STRIDE],delta=delta)
                close(out,filed[key]['frequency_components'])
                assert len(out)==1 and contains(out,frequency)
            envelope = intervals(records[1],records[STRIDE],delta=2*DELTA)
            for signs,corner in zip(product((-1,1),repeat=4),filed['corner_perturbations']):
                first = shifted(records[1],DELTA*signs[0],DELTA*signs[1])
                last = shifted(records[STRIDE],DELTA*signs[2],DELTA*signs[3])
                check_records(first,corner['first_record']); check_records(last,corner['stride8_record'])
                out = intervals(first,last)
                close(out,corner['reader']['frequency_components'])
                assert contains(out,frequency) and contained(out,envelope)
                counts['corner_cases'] += 1
            for _ in range(32):
                errors = [rng.uniform(-DELTA,DELTA) for _ in range(4)]
                out = intervals(shifted(records[1],*errors[:2]),shifted(records[STRIDE],*errors[2:]))
                assert contains(out,frequency) and contained(out,envelope)
                counts['interior_perturbations'] += 1
            false = [[x/STRIDE for x in principal(records[STRIDE],ETA,DELTA)]]
            assert not contains(false,frequency)
            counts['false_principal_branch_rejections'] += 1
            if k == 6:
                bounds = intervals(records[1],records[STRIDE])[0]
                final_bounds.append([SCALE*x for x in bounds]); corner_bounds.append([SCALE*x for x in envelope[0]])
                assert SCALE*envelope[0][1] < OLD
                control = controls[identity]
                shift = (OLD-1022)/SCALE
                rec = {s:probabilities(moments[s]*np.exp(1j*s*shift)) for s in moments}
                for s in rec: check_records(rec[s],control['records'][str(s)])
                out = intervals(rec[1],rec[STRIDE])
                close(out,control['reader']['frequency_components'])
                assert 1-weights[0]**2 <= ETA and contains(out,frequency+shift) and not contains(out,frequency)
                counts['controls'] += 1
            rows.append({'k':k,'case':identity,'pair_off_band':off,
                         'positive_f_weight':weights[0],'negative_f_weight':weights[1]})
            counts['source_cases'] += 1
        if k == 6:
            higher_theta = 2*math.asin(2.**-k*math.sqrt(.3))
            j = int(np.argmin(abs(eigenvalues-np.exp(1j*higher_theta))))
            v = eigenvectors[:,j]/np.linalg.norm(eigenvectors[:,j])
            rec = {s:probabilities(complex(np.vdot(v,U@v)**2)) for s,U in ((1,step),(8,later))}
            control = data['controls'][-1]
            for s in rec: check_records(rec[s],control['records'][str(s)])
            out = intervals(rec[1],rec[8]); close(out,control['reader']['frequency_components'])
            assert contains(out,2*higher_theta*2**k) and not contains(out,frequency)
            counts['controls'] += 1
    # Entire unsigned chart including all turning points, two slopes and aliases.
    for index in range(4097):
        omega = math.pi*index/4096
        for eta in (0.,.01,.2):
            rec = {s:probabilities(complex(np.exp(1j*s*omega))) for s in (1,STRIDE)}
            result = intervals(rec[1],rec[STRIDE],eta=eta,delta=1e-12)
            assert contains(result,omega), (index,eta,result)
            counts['full_chart_plants'] += 1
    return {'passed':True,'same_group':True,'external_G6':False,'producer_imports':[],
            'counts':counts,'extrema':extrema,'rows':rows,
            'k6_bounded_energy_union':[min(x[0] for x in final_bounds),max(x[1] for x in final_bounds)],
            'k6_all_allowed_perturbations_energy_envelope':[min(x[0] for x in corner_bounds),max(x[1] for x in corner_bounds)],
            'old_pair_keV':OLD,'seconds':time.monotonic()-begin}


if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
