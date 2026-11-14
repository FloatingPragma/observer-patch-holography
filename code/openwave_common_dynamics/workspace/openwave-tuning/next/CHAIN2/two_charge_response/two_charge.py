"""Bounded exterior-square execution of the frozen shared-field native source.

Arrays are blocks of the complete 2048-mode exterior square. Zero blocks
have an implicit exact-zero value; no amplitude threshold is used.
"""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import gzip
import hashlib
import json
import math
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'inputs'
MANIFEST = json.loads((HERE / 'INPUT_SNAPSHOT.json').read_text())


def verify_snapshot():
    pins=dict(MANIFEST['files'])
    annex=HERE/'IMPORT_ANNEX.json'
    if annex.exists():pins.update(json.loads(annex.read_text())['files'])
    for name, expected in pins.items():
        if hashlib.sha256((INPUT / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Immutable input mismatch: ' + name)
    owned_modules={'autonomous','manybody','mode_field','common_walk','compile_walk','current_walk',
                   'cell_complex','mass_coin','massive_coupled','particle_records','composition','native_field'}
    for name,module in tuple(sys.modules.items()):
        path=getattr(module,'__file__',None)
        if path and (name.startswith('model.') or name in owned_modules):
            if not Path(path).resolve().is_relative_to(INPUT):
                raise RuntimeError('Imported source outside immutable snapshot: '+name)


verify_snapshot()
sys.path.insert(0, str(INPUT / 'openwave-tuning/common_dynamics/vector_field'))
import autonomous as af

MODEL = 'native-local12-mass-fermion-plaquette-v3-proposed'
CAP = 120000
MAX_STEPS = 3
RESOURCE_LOG = []


class ResourceRefusal(RuntimeError):
    pass


def check_cap(state, phase):
    if len(state) > CAP:
        raise ResourceRefusal(f'{phase}: {len(state)} blocks exceeds {CAP}; no outcomes pruned')
    RESOURCE_LOG.append({'stage': phase, 'blocks': len(state),
                         'complex_scalar_slots': 256 * len(state),
                         'array_payload_bytes': 4096 * len(state)})


def fresh():
    return np.zeros((16, 16), dtype=np.complex128)


@lru_cache(None)
def coin(mass):
    r = -np.eye(16)
    r[:12, :12] += 1 / 6
    mu = mass / 8
    return math.cos(mu) * r - 1j * math.sin(mu) * np.eye(16)


def initial(x=0, y=1, dark=True):
    """Two orthogonal local orbitals and F=P_x+P_y, G=2 delta_0."""
    if x == y:
        raise ValueError('disjoint starting ports required')
    a, b = af.mb.common.encode_position(x), af.mb.common.encode_position(y)
    u = np.zeros(16, complex); u[:12] = 1 / np.sqrt(12)
    v = np.zeros(16, complex)
    if dark:
        v[0], v[2] = 1 / np.sqrt(2), -1 / np.sqrt(2)
    else:
        v[:12] = 1 / np.sqrt(12)
    field = bytes(((af.current_walk.geometry()['P'][:, x]
                    + af.current_walk.geometry()['P'][:, y]) % 16).astype(np.uint8))
    if a < b:
        return {(a, b, field): np.outer(u, v)}
    return {(b, a, field): -np.outer(v, u)}


def onsite(state, mass=0, car=True):
    # Rank-one formula avoids artificial support from a dense product on
    # the exact held-dark invariant subspace. No threshold is applied.
    diagonal = -np.exp(1j*mass/8)
    uniform = math.cos(mass/8)/6
    def apply(v):
        out=diagonal*v.copy()
        out[:12,:] += uniform*v[:12,:].sum(axis=0)[None,:]
        second=diagonal*out.copy()
        second[:,:12] += uniform*out[:,:12].sum(axis=1)[:,None]
        return second
    result = {}
    for (a, b, f), v in state.items():
        if a == b:
            # Exterior-square action, including two modes at the same port.
            antisymmetric = v - v.T if car else v + v.T
            out = apply(antisymmetric)
            result[a, b, f] = np.triu(out, 1)
        else:
            result[a, b, f] = apply(v)
    return result


@lru_cache(None)
def mode_move(mode):
    permutation = af.mb.common.expected_swap_permutation()
    after = int(permutation[mode])
    address, event = divmod(mode, 16)
    vertex = af.mb.common.decode_position(address)
    if vertex is None or event >= 12 or event % 2 == 0:
        return after, ()
    edge = int(af.mb.common.geometry_tables()[1][vertex, event // 2])
    left = af.mb.common.graph().seams[edge][0]
    return after, ((edge, 1 if vertex == left else -1),)


@lru_cache(None)
def pair_move(a, b, i, j):
    m, c = mode_move(16*a + i); n, d = mode_move(16*b + j)
    changes = defaultdict(int)
    for edge, delta in c + d:
        changes[edge] += delta
    sign = 1 if m < n else -1
    m, n = sorted((m, n))
    return m//16, n//16, m%16, n%16, sign, tuple(sorted(changes.items()))


def shift_field(f, changes):
    if not changes:
        return f
    out = bytearray(f)
    for e, delta in changes:
        out[e] = (out[e] + delta) % 16
    return bytes(out)


def transport(state, car=True):
    result = defaultdict(fresh)
    for (a, b, f), v in state.items():
        for i, j in zip(*np.nonzero(v)):
            x, y, k, l, sign, current = pair_move(a, b, int(i), int(j))
            result[x, y, shift_field(f, current)][k, l] += (sign if car else 1) * v[i, j]
        if len(result) > CAP:
            raise ResourceRefusal('transport cap exceeded; no outcomes pruned')
    result = dict(result)
    check_cap(result, 'transport')
    return result


@lru_cache(None)
def face_data(beta):
    a = af.data()[2][:, 0]
    pivot = int(np.flatnonzero(a)[0]); sign = int(a[pivot])
    q = np.arange(16)
    fourier = np.exp(2j * np.pi * np.outer(q, q) / 16) / 4
    kick = (fourier * np.exp(1j * beta * np.cos(2*np.pi*q/16))[None, :]) @ fourier.conj().T
    return a, pivot, sign, kick


def field_step(state, alpha=.17, beta=.41):
    if beta == 0:
        return {key: value * np.exp(-1j*alpha*af.electric_energy(key[2]))
                for key, value in state.items()}
    a, pivot, sign, kick = face_data(beta)
    groups = defaultdict(lambda: np.zeros((16, 256), complex))
    for (x, y, f), v in state.items():
        coordinate = (sign * f[pivot]) % 16
        canonical = af.add_field(f, -coordinate*a)
        groups[x, y, canonical][coordinate] += v.reshape(-1)
    if 16 * len(groups) > CAP:
        raise ResourceRefusal('plaquette orbit cap exceeded; no outcomes pruned')
    result = {}
    for (x, y, f), v in groups.items():
        out = kick @ v
        for t in range(16):
            ft = af.add_field(f, t*a)
            result[x, y, ft] = out[t].reshape(16, 16) * np.exp(-1j*alpha*af.electric_energy(ft))
    check_cap(result, 'field')
    return result


def step(state, mass=0, alpha=.17, beta=.41, car=True):
    return field_step(transport(onsite(state, mass, car), car), alpha, beta)


def norm(state):
    return math.fsum(float(np.vdot(v, v).real) for v in state.values())


def dark_occupation(state, vertex):
    """Norm of a_d |psi>, with complete contraction signs."""
    address = af.mb.common.encode_position(vertex)
    result = defaultdict(lambda: np.zeros(16, complex))
    d = np.zeros(16); d[0], d[2] = 1/np.sqrt(2), -1/np.sqrt(2)
    for (a, b, f), v in state.items():
        if a == b == address:
            result[a, f] += d @ (v-v.T)
        else:
            if a == address:
                result[b, f] += d @ v
            if b == address:
                result[a, f] -= v @ d
    return math.fsum(float(np.vdot(v, v).real) for v in result.values())


def records(state, dark_vertex=None, last_current=True):
    B = af.current_walk.geometry()['B']
    positions = np.zeros(128); pair_positions = np.zeros((128, 128))
    current = np.zeros(288); mode_loads=np.zeros((128,16)); energy = 0.; gauss = True
    for (a, b, f), v in state.items():
        probs = abs(v)**2; p = float(probs.sum())
        pair_positions[a, b] += p; positions[a] += p; positions[b] += p
        energy += p * af.electric_energy(f)
        x, y = af.mb.common.decode_position(a), af.mb.common.decode_position(b)
        n = np.zeros(96, int)
        if x is not None: n[x] += 1
        if y is not None: n[y] += 1
        # All executed preparations lie in the invariant physical sector.
        if p:
            gauss &= not np.any((n+B@np.frombuffer(f,dtype=np.uint8)-2*np.eye(96,dtype=int)[:,0]) % 16)
        if last_current:
            mode_loads[a]+=probs.sum(axis=1)
            mode_loads[b]+=probs.sum(axis=0)
    if last_current:
        # Linearity permits summing terminal mode occupations before
        # decoding current. This equals the complete pair-law decoder.
        for mode,weight in enumerate(mode_loads.reshape(-1)):
            if weight:
                for edge,delta in mode_move(mode)[1]:
                    current[edge]-=delta*weight
    physical = [float(positions[af.mb.common.encode_position(x)]) for x in range(96)]
    result = {'norm_squared':norm(state),'all128_expected_port_occupations':positions.tolist(),
              'all96_expected_physical_occupations':physical,
              'all128_unordered_pair_position_law':pair_positions.tolist(),
              'last288_native_mean_currents':current.tolist(),
              'electric_energy_expectation':energy,'modular_Gauss_G_equals_2_delta0':bool(gauss),
              'blocks':len(state),'renormalized':False,'discarded_probability':0.0}
    if dark_vertex is not None:
        probe = np.array(physical); probe[dark_vertex] -= norm(state)
        result['probe96_position_law'] = probe.tolist()
        result['dark_orbital_occupation'] = dark_occupation(state,dark_vertex)
    return result


def write_law(name, state, description):
    """Complete sparse law, including complex state for independent replay."""
    path = HERE / (name + '.json.gz')
    with path.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as gz:
            gz.write(json.dumps({'schema':'oph.two_charge.complete_joint_law.v1','description':description,
                'all_omitted_entries':'Exact zero in the declared block representation; no cutoff.',
                'mode_order':'Global canonical order16*address+event; ordered pair m<n.',
                'current_decoder':'Undo expected_swap_permutation on each terminal mode; sum origin-to-destination currents. Field kicks commute with this record.',
                'native_complements':'Retained by CONTRACT; these are unconditional ideal-source laws, compared through a complete native TV bound.'},separators=(',',':'))[:-1].encode())
            gz.write(b',"blocks":[')
            first=True
            for (a,b,f),v in sorted(state.items()):
                if not first:gz.write(b',')
                first=False
                rows=[[int(i),int(j),float(v[i,j].real),float(v[i,j].imag),float(abs(v[i,j])**2)]
                      for i,j in zip(*np.nonzero(v))]
                gz.write(json.dumps([a,b,f.hex(),rows],separators=(',',':')).encode())
            gz.write(b']}')
    return {'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'compressed_bytes':path.stat().st_size}


def embed_dark(probe, vertex, phase=1):
    """Isometry a_d^dagger on the exactly orthogonal probe subspace.

    The probe argument carries the SHARED field P_x+P_y and its subsequent
    autonomous evolution. This is an invariant-sector representation of the
    two-particle source, not a substitution of separate gauge registers.
    """
    result=defaultdict(fresh)
    address=af.mb.common.encode_position(vertex)
    for (x,f),u in probe.items():
        a=af.mb.common.encode_position(x)
        for i,value in enumerate(u):
            for j,weight in ((0,1/np.sqrt(2)),(2,-1/np.sqrt(2))):
                m,n=16*a+i,16*address+j
                if m==n:continue
                if m<n:result[a,address,f][i,j] += phase*value*weight
                else:result[address,a,f][j,i] -= phase*value*weight
    return dict(result)


def reduced_dark(name,y=1,mass=0,alpha=.17,beta=.41,steps=3,save=True,independent=False):
    if type(steps)is not int or not 0<=steps<=MAX_STEPS:
        raise ResourceRefusal('development domain is at most three steps')
    field=next(iter(initial(y=y)))[2]
    if independent:field=bytes(288)
    probe={(0,field):np.ones(12,dtype=complex)/np.sqrt(12)}
    rows=[]; states=[]
    for n in range(steps+1):
        phase=(-np.exp(1j*mass/8))**n
        state=embed_dark(probe,y,phase)
        if independent:
            # Record only the separate probe's field; its two-charge Gauss
            # failure deliberately distinguishes this invalid substitution.
            pass
        row=records(state,y,last_current=n>0);row['step']=n;rows.append(row)
        assert abs(row['norm_squared']-1)<2e-10
        assert abs(row['dark_orbital_occupation']-1)<2e-10
        if not independent:assert row['modular_Gauss_G_equals_2_delta0']
        states.append(state)
        if n<steps:
            r=np.ones((12,12))/6-np.eye(12)
            mu=mass/8
            # The number-preserving kick precedes the frozen R/S source.
            kick=math.cos(mu)*np.eye(12)-1j*math.sin(mu)*r
            probe=af.step({key:kick@v for key,v in probe.items()},alpha,beta)
        print(json.dumps({'experiment':name,'step':n,'blocks':len(state),'norm':row['norm_squared']}),flush=True)
    result={'name':name,'second_vertex':y,'second_orbital':'held dark','mass':mass,'h':'1/8',
            'alpha':alpha,'beta_face0':beta,'other207_betas':0,'rows':rows,
            'execution':'Proved exact invariant two-Fock-particle sector, checked against generic exterior-square dynamics.',
            'shared_field':not independent}
    if save:result['complete_joint_law']=write_law(name,state,name)
    return result,states


def execute(name,y=1,dark=True,mass=0,alpha=.17,beta=.41,steps=3,save=True,car=True):
    if type(steps) is not int or not 0<=steps<=MAX_STEPS:
        raise ResourceRefusal('development domain is at most three steps')
    state=initial(y=y,dark=dark); rows=[];start=time.monotonic()
    for n in range(steps+1):
        row=records(state,y if dark else None,last_current=n>0);row['step']=n;rows.append(row)
        if car:
            assert abs(row['norm_squared']-1)<2e-10,(name,n,row['norm_squared'])
            assert row['modular_Gauss_G_equals_2_delta0']
            if dark:assert abs(row['dark_orbital_occupation']-1)<2e-10
        if n<steps:
            before=np.array(row['all96_expected_physical_occupations'])
            state=step(state,mass,alpha,beta,car)
            after=records(state,None,True)
            residual=np.array(after['all96_expected_physical_occupations'])-before+af.current_walk.geometry()['B']@np.array(after['last288_native_mean_currents'])
            row['continuity_residual']=float(abs(residual).max())
            if car:assert row['continuity_residual']<2e-10
        print(json.dumps({'experiment':name,'step':n,'blocks':len(state),'norm':row['norm_squared']}),flush=True)
    receipt={'name':name,'second_vertex':y,'second_orbital':'held dark' if dark else 'uniform moving',
             'mass':mass,'h':'1/8','alpha':alpha,'beta_face0':beta,'other207_betas':0,
             'rows':rows,'runtime_seconds':time.monotonic()-start}
    if save:receipt['complete_joint_law']=write_law(name,state,receipt['name'])
    return receipt,state


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--steps',type=int,default=2)
    parser.add_argument('--y',type=int,default=1);parser.add_argument('--mobile',action='store_true')
    parser.add_argument('--mass',type=int,default=0);args=parser.parse_args()
    result,_=execute('development',y=args.y,dark=not args.mobile,mass=args.mass,steps=args.steps)
    (HERE/'DEVELOPMENT.json').write_text(json.dumps(result,indent=2)+'\n')
