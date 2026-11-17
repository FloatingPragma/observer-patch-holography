#!/usr/bin/env python3
"""Executed R-work source; record-only characteristic field and native census."""
from __future__ import annotations
from collections import Counter
from dataclasses import asdict
from fractions import Fraction
import argparse
import cmath
import gzip
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RW_DIR = HERE.parent / 'rows08_10_reversible_work'
MEM_DIR = HERE.parent / 'rows08_10_native_packet_memory'


def encode(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()


def sha(value):
    return hashlib.sha256(value).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_sources():
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    for path, expected in pins.items():
        if sha((ROOT / path).read_bytes()) != expected:
            raise ValueError('Source changed: ' + path)
    sys.path.insert(0, str(RW_DIR))
    rw = load_module('characteristic_work', RW_DIR / 'verify.py')
    mem = load_module('characteristic_memory', MEM_DIR / 'verify.py')
    return rw, mem


def symbol(packet):
    if packet is None:
        return 0
    e, c, j = packet
    return 1 + 4 * e + 2 * c + (j + 1) // 2


def packet_from_symbol(s):
    if s == 0:
        return None
    value = s - 1
    return [value // 4, (value % 4) // 2, 2 * (value % 2) - 1]


class Observer(list):
    """This reader receives only the initial public reference and event records."""
    def __init__(self, seams):
        super().__init__()
        self.seams = seams
        self.x, self.slots, self.frames = [], [], []
        self.errors = Counter()
        self.kinds = Counter()
        self.transfers = []

    def field(self):
        u = [0 if p[0] is None else p[0][2] for p in self.slots]
        v = [0 if p[1] is None else p[1][2] for p in self.slots]
        return dict(u=u, v=v, E=[a+b for a,b in zip(u,v)],
                    B=[a-b for a,b in zip(u,v)],
                    epsilon=[a*a+b*b for a,b in zip(u,v)])

    def budget(self):
        return sum(z*z for cell in self.x for z in cell), sum(self.field()['epsilon'])

    def append(self, value):
        rec = json.loads(json.dumps(value))
        super().append(rec)
        kind = rec[0]
        self.kinds[kind] += 1
        if kind == 'reference':
            self.x = [cell[:] for cell in rec[1]]
            self.slots = json.loads(json.dumps(rec[2]))
            return
        if kind == 'observation':
            self.frames.append(dict(t=rec[1], phase=rec[2], field=self.field(),
                load_square=self.budget()[0], charge=sum(map(sum,self.x)),
                packets=json.loads(json.dumps(self.slots)), record_index=len(self)-1))
            return
        before = self.budget()
        before_charge = sum(map(sum, self.x))
        i = rec[1]
        if kind in ('native', 'emit', 'capacity_hold'):
            _, _, e, c, j = rec
            a,b = self.seams[e]
            contrast = self.x[i][a]-self.x[i][b]
            expected = (contrast - (2*c-1)*(abs(contrast)%2))//2
            if kind == 'capacity_hold':
                self.errors['capacity_guard'] += abs(contrast) != 2 or all(p is None for p in self.slots[i]) or j != 0
            else:
                self.errors['native_current'] += j != expected
                self.x[i][a] -= j
                self.x[i][b] += j
                if kind == 'emit':
                    self.errors['emit_guard'] += abs(contrast) != 2 or any(p is not None for p in self.slots[i])
                    self.slots[i] = [[e,c,j], [e,c,j]]
        elif kind == 'absorb':
            _, _, e, c, j = rec
            a,b = self.seams[e]
            self.errors['absorb_guard'] += self.slots[i] != [[e,c,-j],[e,c,-j]] or self.x[i][a] != 0 or self.x[i][b] != 0
            self.x[i][a] -= j
            self.x[i][b] += j
            self.slots[i] = [None,None]
        elif kind == 'swap':
            _, i, lane, k, other = rec
            self.slots[i][lane], self.slots[k][other] = self.slots[k][other], self.slots[i][lane]
        elif kind == 'scramble':
            self.slots[i][rec[2]] = rec[3]
        elif kind in ('detector_idle','detector_single','detector_mismatch','detector_busy'):
            left,right = self.slots[i]
            valid = ((kind == 'detector_idle' and left is None and right is None)
                or (kind == 'detector_single' and ((left is None) != (right is None)))
                or (kind == 'detector_mismatch' and left is not None and right is not None and left != right))
            if kind == 'detector_busy' and left is not None and left == right:
                a,b = self.seams[left[0]]
                valid = bool(self.x[i][a] or self.x[i][b])
            self.errors['detector_outcome_guard'] += not valid
        elif kind not in ('detector_idle','detector_single','detector_mismatch','detector_busy','detector_disabled'):
            self.errors['unknown_kind'] += 1
        after = self.budget()
        self.errors['combined_budget'] += sum(before) != sum(after)
        self.errors['charge'] += sum(map(sum,self.x)) != before_charge
        if kind in ('emit','absorb'):
            change = [after[k]-before[k] for k in range(2)]
            self.transfers.append(dict(kind=kind, site=i, change=change))
            self.errors['pair_field_transfer'] += change != ([-2,2] if kind == 'emit' else [2,-2])


def cycle(world, clock, observer, transport=True, detector=True):
    size = len(world.x)
    for i in range(size):
        world.charge_event(i, (clock+i)%30, (clock+i)%2)
    observer.append(['observation',clock,'pre_transport'])
    if transport:
        for i in range(size):
            world.swap(i,0,i,1)
        for i in range(size):
            world.swap(i,1,(i+1)%size,0)
    observer.append(['observation',clock,'post_transport'])
    for i in range(size):
        world.detect(i,detector)
    observer.append(['observation',clock,'post_detector'])


def flight(rw, size, e, j, coin, mode):
    observer = Observer(rw.SEAMS)
    world = rw.World(size, observer)
    a,b = rw.SEAMS[e]
    if mode != 'no_pair':
        world.x[0][a],world.x[0][b] = j,-j
    world.reference()
    observer.append(['observation',0,'pre_emission'])
    world.charge_event(0,e,coin)
    observer.append(['observation',0,'launch'])
    if mode == 'scrambled':
        old = world.slots[0][1]
        world.slots[0][1] = (old[0],old[1],-old[2])
        world.ds[0][1] = world.slots[0][1]
        observer.append(['scramble',0,1,list(world.slots[0][1])])
    for t in range(1,size//2+1):
        cycle(world,t,observer,mode != 'transport_off',mode != 'detector_off')
    return dict(kind='flight',L=size,e=e,j=j,coin=coin,mode=mode,records=list(observer),
        source_errors=dict(world.errors), source_native_calls=world.native_calls)


def impulse(rw,size,lane,site):
    observer = Observer(rw.SEAMS)
    world = rw.World(size,observer)
    world.slots[site][lane] = (0,0,1)
    world.reference()
    observer.append(['observation',0,'launch'])
    cycle(world,1,observer,True,False)
    return dict(kind='impulse',L=size,lane=lane,site=site,records=list(observer),
        source_errors=dict(world.errors),source_native_calls=world.native_calls)


def planted(rw,u,v):
    reader = Observer(rw.SEAMS)
    slots = [[None if q == 0 else [0,0,q] for q in (u,v)]]
    reader.append(['reference',[[0]*12],slots])
    reader.append(['observation',0,'planted'])
    return dict(kind='planted',u=u,v=v,records=list(reader),source_errors={},source_native_calls=0)


def read_case(case,seams):
    reader = Observer(seams)
    for rec in case['records']:
        reader.append(rec)
    return reader


def transport_errors(before,after):
    n = len(before['u'])
    errors = Counter()
    for i in range(n):
        left,right = (i-1)%n,(i+1)%n
        errors['characteristics'] += after['u'][i] != before['u'][left] or after['v'][i] != before['v'][right]
        flux_left = before['u'][left]**2-before['v'][i]**2
        flux_right = before['u'][i]**2-before['v'][right]**2
        errors['local_continuity'] += after['epsilon'][i]-before['epsilon'][i] != flux_left-flux_right
        errors['electric_update'] += 2*after['E'][i] != before['E'][left]+before['E'][right]+before['B'][left]-before['B'][right]
        errors['magnetic_update'] += 2*after['B'][i] != before['B'][left]+before['B'][right]+before['E'][left]-before['E'][right]
    return errors


def field_report(cases,seams):
    errors = Counter()
    source_errors = Counter()
    kinds = Counter()
    waves = 0
    flights, controls, impulses, planted_rows = [], [], {}, []
    checkpoints = []
    for index,case in enumerate(cases):
        reader = read_case(case,seams)
        errors.update(reader.errors)
        source_errors.update(case['source_errors'])
        kinds.update(reader.kinds)
        frames = {(f['t'],f['phase']):f for f in reader.frames}
        for f in reader.frames:
            ff = f['field']
            errors['field_norm_identity'] += any(E*E+B*B != 2*ep for E,B,ep in zip(ff['E'],ff['B'],ff['epsilon']))
        if case['kind'] == 'planted':
            f = reader.frames[0]['field']
            ok = f['E'] == [case['u']+case['v']] and f['B'] == [case['u']-case['v']]
            errors['planted_reader'] += not ok
            planted_rows.append(dict(u=case['u'],v=case['v'],field=f,correct=ok))
            continue
        size = case['L']
        if case['kind'] == 'impulse':
            f = frames[1,'post_transport']['field']
            impulses.setdefault(size,[]).append((case['lane'],case['site'],f['u']+f['v']))
            errors.update(transport_errors(frames[1,'pre_transport']['field'],f))
            continue
        mode = case['mode']
        actual_checks = Counter()
        for t in range(1,size//2+1):
            e = transport_errors(frames[t,'pre_transport']['field'],frames[t,'post_transport']['field'])
            actual_checks.update(e)
            if mode != 'transport_off':
                errors.update(e)
        sequence = [frames[0,'launch']['field']] + [frames[t,'post_transport']['field'] for t in range(1,size//2+1)]
        if mode == 'native':
            for t in range(1,len(sequence)-1):
                for field in ('E','B'):
                    for i in range(size):
                        bad = sequence[t+1][field][i]+sequence[t-1][field][i] != sequence[t][field][(i-1)%size]+sequence[t][field][(i+1)%size]
                        errors['wave_recurrence'] += bad
                        waves += 1
        first = None
        current_time = 0
        for rec in case['records']:
            if rec[0] == 'observation':
                current_time = rec[1]
            if rec[0] == 'absorb' and rec[1] == size//2 and first is None:
                first = current_time
        a,b = seams[case['e']]
        received = reader.x[size//2][a] == case['j'] and reader.x[size//2][b] == -case['j']
        row = dict(case=index,L=size,e=case['e'],j=case['j'],coin=case['coin'],mode=mode,
            first_response=first,received=received,final_load_square=reader.budget()[0],
            final_field_norm=reader.budget()[1],transfers=reader.transfers,
            unit_transport_defects=dict(actual_checks))
        if mode == 'native':
            flights.append(row)
            errors['arrival'] += first != size//2 or not received
            errors['arrival_budget'] += reader.budget() != (2,0)
            errors['flight_budget'] += any(sum(f['field']['epsilon']) != 2 or f['load_square'] != 0 for f in [frames[0,'launch']] + [frames[t,'post_transport'] for t in range(1,size//2+1)])
            if size == 8:
                for name,t,phase,site in (('source',0,'launch',0),('flight',1,'post_transport',1),
                    ('detector_before',4,'post_transport',4),('detector_after',4,'post_detector',4)):
                    f = frames[t,phase]
                    checkpoints.append(dict(case=index,checkpoint=name,t=t,site=site,
                        symbols=[symbol(p) for p in f['packets'][site]],
                        E=f['field']['E'][site],B=f['field']['B'][site]))
        else:
            controls.append(row)
    mode_rows = []
    for size,rows in impulses.items():
        matrix = np.zeros((2*size,2*size),dtype=np.int64)
        for lane,site,out in rows:
            matrix[:,lane*size+site] = out
        expected = np.zeros_like(matrix)
        S = np.zeros_like(matrix)
        for lane,site in product(range(2),range(size)):
            expected[lane*size+(site+(1 if lane == 0 else -1))%size,lane*size+site] = 1
            for step in (-1,1):
                S[lane*size+(site+step)%size,lane*size+site] += 1
        identity = np.eye(2*size,dtype=np.int64)
        errors['impulse_transfer'] += not np.array_equal(matrix,expected)
        errors['transfer_orthogonality'] += not np.array_equal(matrix.T@matrix,identity)
        errors['wave_matrix'] += int(np.count_nonzero(matrix@matrix-S@matrix+identity))
        for lane,m in product(range(2),range(size)):
            k = 2*np.pi*m/size
            f = np.zeros(2*size,dtype=complex)
            f[lane*size:(lane+1)*size] = [cmath.exp(1j*k*i) for i in range(size)]
            phase = cmath.exp((-1j if lane == 0 else 1j)*k)
            observed = np.vdot(f,matrix@f)/np.vdot(f,f)
            residual = float(np.max(np.abs(matrix@f-phase*f)))
            errors['mode_phase'] += residual >= 1e-10
            mode_rows.append(dict(L=size,lane=lane,m=m,k=float(k),
                observed_phase=[float(observed.real),float(observed.imag)],
                source_phase=[phase.real,phase.imag],residual=residual))
    control_checks = {
        'all_remote_controls_silent': all(r['first_response'] is None and not r['received'] for r in controls),
        'no_pair_zero_field': all(r['final_field_norm'] == 0 for r in controls if r['mode'] == 'no_pair'),
        'transport_off_rejects_unit_transfer': all(r['unit_transport_defects']['characteristics'] > 0 for r in controls if r['mode'] == 'transport_off'),
        'detector_off_retains_field': all(r['final_field_norm'] == 2 for r in controls if r['mode'] == 'detector_off'),
        'scramble_retains_field': all(r['final_field_norm'] == 2 for r in controls if r['mode'] == 'scrambled'),
    }
    return dict(errors=dict(errors),source_errors=dict(source_errors),event_kinds=dict(kinds),
        cases=len(cases),native_calls=sum(c['source_native_calls'] for c in cases),
        wave_equalities=waves,flights=flights,controls=controls,control_checks=control_checks,
        modes=mode_rows,planted=planted_rows,checkpoints=checkpoints)


def bits(h,s,r,stage):
    prefix = f'characteristic-field-v1:{h}:{s}:{r}:{stage}'.encode()
    counter = 0
    while True:
        for byte in hashlib.sha256(prefix+b':'+str(counter).encode()).digest():
            for bit in range(8):
                yield (byte >> bit)&1
        counter += 1


def memory_service(mem,h,s,request,disabled=False):
    C = mem.C
    binary,all_states = C.native.states(binary=True),C.native.states()
    blank = binary[0]
    events,reads = [],[]
    for r in range(2):
        initial = all_states[(s+13*h+17*r)%78]
        context = dict(service=h,symbol=s,carrier=r,write_disabled=disabled)
        result = C.execute(initial,bits(h,s,r,'reset'),full=True,mode='reset',target=blank)
        mem.record_execution('reset',initial,result,events,context)
        if not result.outcome.success:
            reads.append(None)
            continue
        state = result.output
        if not disabled:
            target = binary[s//66 if r == 0 else s%66]
            written = mem.write(state,target,bits(h,s,r,'write'))
            events.append(dict(written,stage='write',**context))
            if not written['success']:
                reads.append(None)
                continue
            state = tuple(written['output'])
        result = C.execute(state,bits(h,s,r,'read'),full=False,mode='measure')
        mem.record_execution('read',state,result,events,context)
        reads.append(binary.index(result.outcome.initial) if result.outcome.success else None)
    conclusive = all(r is not None for r in reads)
    return dict(service=h,symbol=s,request=request,write_disabled=disabled,events=events,
        read_digits=reads,decoded=reads[0]*66+reads[1] if conclusive else None,
        conclusive=conclusive)


def memory_controls(mem):
    C = mem.C
    binary = C.native.states(binary=True)
    blank = binary[0]
    target = next(t for t in binary if C.route(blank,t))
    seam = C.route(blank,target)[0]
    hold = int(blank[seam[0]] > blank[seam[1]])
    written = mem.write(blank,target,iter([hold]*32))
    write_event = dict(written,stage='forced_write_timeout')
    coins = [int(blank[a] > blank[b]) for _ in range(4) for a,b in C.quantum.seams()]
    result = C.execute(blank,iter(coins),mode='measure')
    read_event = dict(stage='read',initial=list(blank),output=list(result.output),
        outcome=asdict(result.outcome),records=result.records,forced_holding=True)
    return [write_event,read_event]


def memory_report(mem,services,controls,checkpoints):
    events = [e for s in services for e in s['events']] + controls
    replay = mem.replay(events)
    binary = mem.C.native.states(binary=True)
    errors = Counter()
    stage_outcomes = Counter()
    normal = [s for s in services if not s['write_disabled']]
    disabled = [s for s in services if s['write_disabled']]
    for service in services:
        reads = []
        for carrier in (0,1):
            rows = [e for e in service['events'] if e['carrier'] == carrier]
            for row in rows:
                outcome = row.get('outcome')
                stage_outcomes[row['stage']+':'+(outcome['reason'] if outcome else row['reason'])] += 1
            reading = next((e for e in rows if e['stage'] == 'read'),None)
            reads.append(binary.index(tuple(reading['outcome']['initial'])) if reading and reading['outcome']['success'] else None)
        conclusive = all(r is not None for r in reads)
        decoded = reads[0]*66+reads[1] if conclusive else None
        errors['stored_outcome_matches_replay'] += reads != service['read_digits'] or conclusive != service['conclusive'] or decoded != service['decoded']
        expected = 0 if service['write_disabled'] else service['symbol']
        errors['conclusive_symbol'] += conclusive and decoded != expected
    by_checkpoint = {(s['request']['checkpoint_index'],s['request']['lane']):s for s in normal}
    field_reads = []
    for i,checkpoint in enumerate(checkpoints):
        pair = [by_checkpoint[i,lane] for lane in (0,1)]
        conclusive = all(s['conclusive'] for s in pair)
        EB = None
        if conclusive:
            packets = [packet_from_symbol(s['decoded']) for s in pair]
            u,v = [0 if p is None else p[2] for p in packets]
            EB = [u+v,u-v]
            errors['conclusive_field'] += EB != [checkpoint['E'],checkpoint['B']]
        field_reads.append(dict(checkpoint=i,conclusive=conclusive,EB=EB))
    reset_budget = mem.C.budget(full=True,mode='reset',target=binary[0])
    read_budget = mem.C.budget(mode='measure')
    bound = 2*reset_budget['failure_upper'] + Fraction(8,2**32) + 2*read_budget['failure_upper']
    max_calls = 2*reset_budget['native_attempts_max']+8*32+2*read_budget['native_attempts_max']
    return dict(errors=dict(errors),record_replay=replay,stage_outcomes=dict(stage_outcomes),
        normal_services=len(normal),normal_conclusive=sum(s['conclusive'] for s in normal),
        normal_inconclusive=sum(not s['conclusive'] for s in normal),
        disabled_services=len(disabled),disabled_conclusive=sum(s['conclusive'] for s in disabled),
        disabled_inconclusive=sum(not s['conclusive'] for s in disabled),
        field_reads=field_reads,forced_controls=[{k:v for k,v in e.items() if k != 'records'} for e in controls],
        bounds=dict(per_service_failure_upper=str(bound),normal_union_failure_upper=str(len(normal)*bound),
            all_services_union_failure_upper=str(len(services)*bound),native_attempts_max_per_service=max_calls,
            native_event_tape_bits_per_service=8*max_calls),
        checks=dict(disabled_write_detected=any(s['conclusive'] and s['decoded'] != s['symbol'] for s in disabled),
            forced_write_timeout=not controls[0]['success'] and controls[0]['attempts'] == 32,
            forced_census_timeout=not controls[1]['outcome']['success'] and controls[1]['outcome']['initial'] is None
                and controls[1]['outcome']['reason'] == 'seed_timeout' and len(controls[1]['records']) == 120))


def validate(data,rw,mem):
    field = field_report(data['cases'],rw.SEAMS)
    memory = memory_report(mem,data['memory_services'],data['memory_controls'],field['checkpoints'])
    checks = dict(field_exact=not any(field['errors'].values()),source_exact=not any(field['source_errors'].values()),
        controls=all(field['control_checks'].values()),contact_denominator=len(field['flights']) == 360,
        full_mode_denominator=len(field['modes']) == 112,planted_denominator=len(field['planted']) == 9,
        native_memory_exact=not any(memory['errors'].values()) and not any(memory['record_replay']['errors'].values()),
        native_memory_controls=all(memory['checks'].values()),
        memory_denominator=memory['normal_services'] == 960 and memory['disabled_services'] == 240)
    return dict(schema='oph.characteristic_field.v1',candidate='R-work-characteristic-field',class_='N+',
        field=field,memory=memory,checks=checks,all_checks_pass=all(checks.values()),
        full_physical_rows_solved=[],physical_energy_identified=False,
        execution_scope='Unchanged symbolic R-work plus separately executed native checkpoint storage/readout; no full native swap/inverse implementation.')


def run(rw,mem):
    cases = [flight(rw,L,e,j,c,'native') for L,e,j,c in product((8,16,32),range(30),(-1,1),(0,1))]
    cases += [flight(rw,L,0,1,0,mode) for L,mode in product((8,16,32),('no_pair','transport_off','detector_off','scrambled'))]
    cases += [impulse(rw,L,lane,site) for L in (8,16,32) for lane in (0,1) for site in range(L)]
    cases += [planted(rw,u,v) for u,v in product((-1,0,1),repeat=2)]
    field = field_report(cases,rw.SEAMS)
    services = []
    for index,checkpoint in enumerate(field['checkpoints']):
        for lane,s in enumerate(checkpoint['symbols']):
            services.append(memory_service(mem,len(services),s,dict(checkpoint_index=index,lane=lane)))
    for index,checkpoint in enumerate(field['checkpoints']):
        if checkpoint['checkpoint'] == 'source':
            for lane,s in enumerate(checkpoint['symbols']):
                services.append(memory_service(mem,len(services),s,dict(checkpoint_index=index,lane=lane),True))
    data = dict(cases=cases,memory_services=services,memory_controls=memory_controls(mem))
    return data,validate(data,rw,mem)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    rw,mem = load_sources()
    data,receipt = run(rw,mem)
    raw = encode(data)
    receipt['records'] = dict(file='RECORDS.json.gz',uncompressed_sha256=sha(raw))
    if args.write:
        (HERE/'RECORDS.json.gz').write_bytes(gzip.compress(raw,mtime=0))
        (HERE/'RECEIPT.json').write_bytes(encode(receipt))
    print(json.dumps(dict(checks=receipt['checks'],field_errors=receipt['field']['errors'],
        source_native_calls=receipt['field']['native_calls'],memory_errors=receipt['memory']['errors'],
        memory_native_calls=receipt['memory']['record_replay']['attempts'],
        normal_conclusive=receipt['memory']['normal_conclusive'],normal_inconclusive=receipt['memory']['normal_inconclusive'],
        receipt_sha256=sha(encode(receipt))),indent=2,sort_keys=True))
    return 0 if receipt['all_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
