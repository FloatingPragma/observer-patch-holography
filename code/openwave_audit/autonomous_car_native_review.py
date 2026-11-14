"""Independent N16 source/native audit; producer modules are not imported."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACK = ROOT / 'openwave-tuning/common_dynamics/autonomous_car_native'
FIELD = ROOT / 'openwave-tuning/next/CHAIN3/native_field_connection'
EXPECTED = '1f510cbd194afa43f4bb37fc1cc45792d582aff67bee201b06504141a6f892d7'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def ceil(x): return -(-x.numerator // x.denominator)
def sqrt_bounds(x, bits=256):
    x = F(x); scale = 1 << bits
    a = math.isqrt((x.numerator * scale * scale) // x.denominator)
    return F(a, scale), F(a + 1, scale)
def up(x): return sqrt_bounds(x)[1]


def custody():
    assert sha(PACK / 'MANIFEST.json') == EXPECTED
    m = read(PACK / 'MANIFEST.json')
    pins = {str((PACK / 'MANIFEST.json').relative_to(ROOT)): EXPECTED}
    for section in ('files', 'source_pins', 'consumed_manifests'):
        for key, digest in m[section].items():
            assert key not in pins or pins[key] == digest, key
            assert sha(ROOT / key) == digest, key
            pins[key] = digest
    return pins


def budget_and_source():
    receipt = read(PACK / 'BUDGET.json'); config = receipt['config']
    coord = read(FIELD / 'COORDINATES.json')
    spec = read(FIELD / 'inputs/openwave-tuning/common_dynamics/quadratic_field/FOCK.json')['spectral_certificate']
    u = F(1, 64); R = 95; M = 16; K = 12; gamma2 = F(receipt['coexact_Gamma_squared_upper_exact'])
    assert config['cell_time_exact'] == str(u) and config['physical_time_exact'] == str(u / 8)
    assert config['coexact_cap'] == R and config['harmonic_half_width'] == M and config['Dyson_order'] == K
    assert M > K and config['harmonic_bits'] == 5
    # Independently enclose each displacement using exact numerator data.
    norms = []; coefficient_count = 0
    stored_norms = {(r['cell'], r['axis']): r for r in receipt['native_certification']['native_midpoint_rows']}
    for route in coord['currents']:
        if route['sign'] != 1: continue
        lo2 = hi2 = F(); radius_sum = F()
        for numerator, lam in zip(route['coexact_rational_numerators'], spec['eigenvalues_exact']):
            wlo, whi = sqrt_bounds(F(3, 2) * F(lam), 320)
            low = sqrt_bounds(1 / (4 * whi), 320)[0]
            high = sqrt_bounds(1 / (4 * wlo), 320)[1]
            f = abs(F(numerator)); low *= f; high *= f
            lo2 += low * low; hi2 += high * high
            radius_sum += 2 * high * F(1024, 1 << 256) + F(1024, 1 << 256) ** 2
            coefficient_count += 1
        saved = stored_norms[route['cell'], route['axis']]
        assert (saved['cell'], saved['axis']) == (route['cell'], route['axis'])
        midnorm = F(saved['midpoint_norm_squared_exact'])
        assert lo2 - radius_sum <= midnorm <= hi2 + radius_sum
        assert hi2 < gamma2 and midnorm < gamma2
        norms.append((lo2, hi2))
    # Independent rational upper for exp(-lambda), using a longer positive sum.
    physical = F(receipt['vacuum_error_exact'])
    eta = F(spec['actual_curvature_operator_error_upper'])
    w0 = sqrt_bounds(F(3, 2) * (F(spec['minimum_eigenvalue']) - eta))[0]
    v = 1 / (4 * w0)
    delta = F(spec['frame_residual_upper']) / F(spec['minimum_eigenvalue']) + F(spec['polar_Gram_error_upper'])
    for k, saved in enumerate(receipt['history_rows']):
        weight = (48 * u) ** k / math.factorial(k); lam = k * k * gamma2
        assert F(saved['absolute_history_weight_exact']) == weight
        assert F(saved['lambda_upper_exact']) == lam and F(saved['harmonic_error_exact']) == 0
        geom = 3 * u * eta * v * up((4 * lam + 190)**2 + 16 * lam + 380) / 2
        geom += delta * k * up(5) * up(v * (4 * lam + 1))
        assert geom <= F(saved['geometry_estimate_exact'])
        if k == 0: co = F()
        elif lam <= R + 1:
            term = total = F(1)
            for j in range(1, ceil(3 * lam) + 101):
                term *= lam / j; total += term
            p = lam ** (R + 1) / (math.factorial(R + 1) * total)
            z = p * (R + 2) / (R + 2 - lam)
            co = up(2 * (z + 190 * lam * (R + 1) * p))
            assert co <= F(saved['coexact_estimate_exact'])
        else: co = F(2)
        value = min(F(2), co + geom)
        assert value <= F(saved['complete_history_norm_bound_exact'])
        physical += weight * value
    remainder = 2 * (48 * u) ** 13 / math.factorial(13)
    assert remainder == F(receipt['two_Dyson_remainders_exact'])
    physical += remainder
    assert physical <= F(receipt['error_terms_exact']['physical_history_comparison'])
    # Rebuild the G=2 Coulomb polynomial directly in occupation masks.
    P = [[F(x) for x in sector['particular']] for sector in coord['charge_sectors']]
    gram = [[sum((x*y for x,y in zip(a,b)),F()) for b in P] for a in P]
    zpoly = defaultdict(F); zpoly[0] += sum((x for row in gram for x in row), F())
    for i in range(32):
        c = i // 4; linear = gram[c][c]/4 - sum(gram[c], F())
        zpoly[0] += linear/2; zpoly[1 << i] -= linear/2
        for j in range(i+1, 32):
            coeff = gram[c][j//4]/8
            zpoly[0] += coeff; zpoly[1 << i] -= coeff
            zpoly[1 << j] -= coeff; zpoly[(1 << i) | (1 << j)] += coeff
    zpoly = {k:v for k,v in zpoly.items() if v}
    l1 = sum(map(abs, zpoly.values()), F())
    assert l1 == F(13771,924) == F(config['Coulomb_l1_bound_exact'])
    rng = np.random.default_rng(591721)
    masks = [0, (1 << 32)-1, sum(1 << i for i in range(1,32,2))]
    masks += [int(x) for x in rng.integers(0,1 << 32,size=192)]
    for mask in masks:
        z = sum((val * (-1)**((key & mask).bit_count()) for key,val in zpoly.items()), F())
        q = [sum((mask >> (4*c+s)) & 1 for s in range(4))-2 for c in range(8)]
        target = sum((gram[c][d]*q[c]*q[d] for c in range(8) for d in range(8)), F())/4
        assert z == target
    X0 = 190*11*(R+F(1,2)) + F(3*M*M,48) + l1 + 4
    N = ceil(u*u*(48*X0+1140)*10**6)
    assert N == config['segments'] and F(config['X0_upper_exact']) == X0
    steps = ceil(F(192*N)*gamma2*R*(R+1)*10**6/4)
    assert steps == config['route_subdivisions']
    outer = u*u*(48*X0+1140)/N; pairs = F(192*N)*gamma2*R*(R+1)/(4*steps)
    assert F(config['finite_outer_exact']) == outer and F(config['route_pairs_exact']) == pairs
    coD = 190*steps*R*128*65
    fluxD = sum(2**(j+2)*(8*(j+2)+1) for j in range(5))
    free = 190*128*57 + 3*32*41 + len(zpoly)*257 + 32*2*257
    pulses = N*(free + 96*(2*(coD+fluxD)+514)) + 2000
    assert pulses == config['native_pulses_upper']
    bits = config['precision_bits']; assert bits == 256
    rounding = (192*N*190*(2*R*1024+2*R)+u*190*(R+F(1,2))) / (1 << bits)
    assert rounding == F(receipt['native_certification']['algebraic_rounding_exact'])
    assert F(receipt['native_certification']['angle_upper_exact']) == 4*(R+1)+l1+M*M+100
    assert config['total_register_bits'] == 32+190*7+3*5+3 == 1380
    reset = config['reset_plan']; assert reset['registers'] == 1444
    assert F(reset['complete_timeout_union_bound_exact']) == 1444*F(reset['base_uniform_failure_exact'])**reset['blocks_per_register']
    total = sum(map(F,receipt['error_terms_exact'].values()), F())
    assert total == F(receipt['complete_record_TV_upper_exact']) and total < F(1,100000)
    # Every directed kick changes harmonic flux by a unit; parity survives wrap.
    for p in range(-M,M):
        for sign in (-1,1):
            wrapped = (p+sign+M)%(2*M)-M
            assert (wrapped-p)%2 == 1
            if abs(p) <= K: assert wrapped == p+sign
    return {'coefficient_intervals':coefficient_count,'exact_Coulomb_masks':len(masks),
            'history_orders':13,'segments':N,'native_pulse_upper':pulses,
            'complete_TV_upper_exact':str(total),'complete_TV_display':float(total),
            'full_program_executed':False}


def native_replay():
    raw = read(FIELD/'inputs/reverse-engineering-reality/code/source_selection_model/response.json')
    generators = []
    for blocks in raw['generators']:
        g = np.zeros((6,6),complex)
        for k,block in enumerate(blocks):
            for i,row in enumerate(block):
                for j,x in enumerate(row):
                    g[3*k+i,3*k+j] = float(F(x[0]))+math.sqrt(5)*float(F(x[1]))+1j*(float(F(x[2]))+math.sqrt(5)*float(F(x[3])))
        assert np.linalg.norm(g+g.conj().T) < 1e-13; generators.append(g)
    @lru_cache(None)
    def response(controls):
        return expm(sum((float(F(c))*g for c,g in zip(controls,generators)),np.zeros((6,6),complex)))
    def execute(tape,bits):
        dim=1 << bits; state=np.eye(dim,dtype=complex); failure=np.zeros((dim,dim),complex)
        diagonals=[]; metadata=[]
        for j,event in enumerate(tape):
            targets=event['targets']; local=1 << len(targets); other=[q for q in range(bits) if q not in targets]
            ix=np.array([[sum(((a>>(len(targets)-1-r))&1)<<q for r,q in enumerate(targets))+
                          sum(((b>>r)&1)<<q for r,q in enumerate(other)) for b in range(1<<len(other))] for a in range(local)])
            old=state[ix].copy(); matrix=response(tuple(event['controls_exact']))
            flat=old.reshape(local,-1)
            for c in range(6-local):
                failed=(matrix[[local+c],:local]@flat).reshape(-1,dim)
                failure+=failed.conj().T@failed
                diagonals.append(np.sum(abs(failed)**2,axis=0));metadata.append([j,c])
            state[ix]=(matrix[:local,:local]@flat).reshape(old.shape)
        defect=np.linalg.norm(state.conj().T@state+failure-np.eye(dim),2)
        assert defect<2e-10
        return state,failure,np.array(diagonals),np.array(metadata),float(defect)
    tapes=read(PACK/'LOCAL_PROGRAMS.json'); stored=np.load(PACK/'LOCAL_MAPS.npz')
    controls=read(PACK/'CONTROLS.json'); gamma=float(sum(map(F,controls['gamma_interval_exact']))/2)
    coordinates=read(FIELD/'COORDINATES.json')
    route=next(r for r in coordinates['currents'] if r['sign']==1 and
               all(r[k]==v for k,v in controls['actual_route'].items()))
    mode=controls['actual_coexact_mode']
    spec=read(FIELD/'inputs/openwave-tuning/common_dynamics/quadratic_field/FOCK.json')['spectral_certificate']
    wlo,whi=sqrt_bounds(F(3,2)*F(spec['eigenvalues_exact'][mode]),320)
    slo=sqrt_bounds(1/(4*whi),320)[0];shi=sqrt_bounds(1/(4*wlo),320)[1]
    numerator=F(route['coexact_rational_numerators'][mode])
    lower,upper=sorted((numerator*slo,numerator*shi))
    claimed_lower,claimed_upper=map(F,controls['gamma_interval_exact'])
    assert claimed_lower<=lower<=upper<=claimed_upper
    pointer,pfail,*_=execute(tapes['pointer'],3); select=pointer[:,[0,4]]
    effects={}; states=[(0,0,1,1)+(0,)*8,(0,1,1)+(0,)*9,(1,0,1)+(0,)*9,(0,2)+(0,)*10]
    for law in ('C','E'):
        groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)));meta={}
        for code,state in enumerate(states):
            block=select[[code,code+4]]
            for coin in (0,1):
                diff=state[0]-state[1]; current=(diff-(2*coin-1)*(abs(diff)%2))//2
                out=list(state);out[0]-=current;out[1]+=current
                key=(coin,current if law=='C' else diff)
                groups[key][tuple(out)]+=block/math.sqrt(2)
                meta[key]={'coin':coin,'label':key[1],'current':current,'bit':int(2*current+2*coin-1<0)}
        effects[law]=[(meta[key],sum((K.conj().T@K for K in group.values()),np.zeros((2,2),complex))) for key,group in sorted(groups.items())]
    # Independent endpoint CAR matrices; the middle bit is the full spectator parity.
    def ann(mode):
        out=np.zeros((8,8),complex)
        for mask in range(8):
            if mask>>mode&1:out[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
        return out
    hop=ann(2).conj().T@ann(0); p=np.diag([(x>>2)&1 for x in range(8)])
    ladder=np.zeros((4,4),complex);ladder[0,1]=1;ladder[1,2]=math.sqrt(2)
    W=np.kron(np.roll(np.eye(4),1,axis=0),expm(1j*gamma*(ladder+ladder.conj().T)))
    D=np.kron(np.eye(16),np.eye(8)-p)+np.kron(W,p)
    probabilities=0;first_exports=0;cases=[]
    for saved in controls['cases']:
        label=saved['label']; actual,failed,diag,metadata,defect=execute(tapes[label],7)
        assert np.linalg.norm(actual-stored[label+'_success'],2)<2e-11
        assert np.linalg.norm(failed-stored[label+'_failure_Gram'],2)<2e-11
        assert np.max(abs(diag-stored[label+'_first_failure_diagonal']))<1e-12
        assert np.array_equal(metadata,stored[label+'_failure_metadata'])
        quadrature=label[0];sign=-1 if 'sign-1' in label else 1
        bare=hop+hop.conj().T if quadrature=='X' else -1j*hop+1j*hop.conj().T
        target=expm(-.25j*D@np.kron(np.eye(16),sign*bare/2)@D.conj().T)
        residual=float(np.linalg.norm(actual-target,2))
        assert residual <= float(F(saved['finite_pair_error_upper_exact'])+F(saved['native_error_upper_exact']))+1e-10
        for law,rows in effects.items():
            subtotal=np.zeros((128,128),complex)
            for (rec,effect),arch in zip(rows,saved['fine_pointer_laws'][law]['fine_records']):
                assert rec==arch['record']
                full=actual.conj().T@np.kron(np.eye(16),np.kron(effect,np.eye(4)))@actual
                assert np.max(abs(np.diag(full).real-arch['all128_input_probabilities']))<1e-11
                probabilities+=128;subtotal+=full
            pf=pfail[np.ix_([0,4],[0,4])]
            assert np.linalg.norm(subtotal+actual.conj().T@np.kron(np.eye(16),np.kron(pf,np.eye(4)))@actual+failed-np.eye(128),2)<2e-10
        first_exports+=len(metadata)
        cases.append({'label':label,'pulses':len(tapes[label]),'all_columns':128,'Gram_residual':defect,'exact_W_residual':residual})
    clocks=read(PACK/'CLOCK_PROGRAMS.json'); clock_receipt=read(PACK/'CLOCK.json')
    for saved in clock_receipt['cases']:
        axis=saved['axis'];actual,failed,diag,metadata,defect=execute(clocks[axis],1)
        assert len(metadata)==len(saved['first_source_export_records'])
        for row,arch in zip(diag,saved['first_source_export_records']):
            assert np.max(abs(row-arch['all2_input_probabilities']))<1e-12
        for law,rows in effects.items():
            contrast=0
            for (rec,E),arch in zip(rows,saved['laws'][law]['fine_records']):
                value=float(np.vdot(actual[:,0],E@actual[:,0]).real)
                assert rec==arch['record'] and abs(value-arch['probability'])<1e-12
                contrast+=(-1)**rec['bit']*value;probabilities+=1
            target=math.cos(1/4096) if axis=='X' else math.sin(1/4096)
            assert abs(contrast-target)<1e-12
        first_exports+=len(metadata)
    return {'local_cases':cases,'fine_probabilities_checked':probabilities,
            'first_export_maps_rebuilt':first_exports,'clock_all_column_cases':2,
            'raw_pointer_all_columns':8,'producer_modules_imported':False}


def main():
    pins=custody(); print('CUSTODY',len(pins),flush=True)
    budget=budget_and_source();print('BUDGET_AND_SOURCE_PASS',flush=True)
    native=native_replay();print('RAW_NATIVE_PASS',flush=True)
    assert custody()==pins
    result={'schema':'oph.N16.native.root.audit.v1','checks_pass':True,
            'budget_and_source':budget,'native':native,'source_pins':pins,
            'code_sha256':sha(Path(__file__)),'full_criterion_met':False,
            'scope':'Independent component audit; no native interacting-ground preparation or physical-row award.'}
    (HERE/'autonomous_car_native_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',flush=True)


if __name__=='__main__':main()
