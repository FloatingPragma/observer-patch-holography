"""Site-local decoding and regeneration of authenticated dual-rail words."""
import hashlib
import hmac
import math

import numpy as np

from .chart import addressed
from ..kernel import reflect, require


def legs(route):
    return [(axis,d) for axis,d in enumerate(route) if d]


def front(route,phase=1):
    """Earliest-arrival probability, summing coherent paths before measuring the slot."""
    directions = legs(route)
    require(bool(directions),'empty signal route')
    slots = [(2*a+int(d<0),abs(d)) for a,d in directions]
    limit = tuple(n for _,n in slots)
    start = (0,)*len(slots)
    prepared = reflect(np.eye(8,dtype=complex)[slots[0][0]])
    layer = {start:prepared}
    for _ in range(sum(limit)):
        following = {}
        for at,amplitudes in layer.items():
            scattered = reflect(amplitudes,phase)
            for j,(slot,n) in enumerate(slots):
                if at[j]==n:
                    continue
                target = at[:j]+(at[j]+1,)+at[j+1:]
                if target not in following:
                    following[target] = np.zeros(8,complex)
                following[target][slot^1] += scattered[slot]
        layer = following
    p = float(np.vdot(layer[limit],layer[limit]).real)
    require(0<=p<=1+1e-12,'unphysical boundary probability')
    return min(1.,p)


def itinerary(route,member,segment):
    if member in ('addressed','passive'):
        return [tuple(route)]
    result = []
    for axis,distance in legs(route):
        n,rem = divmod(abs(distance),segment)
        require(rem==0,'relay outside the enrolled menu')
        r = [0]*len(route)
        r[axis] = int(np.sign(distance))*segment
        result += [tuple(r)]*n
    return result


def front_formula(route,phase):
    counts = [abs(d) for d in route if d]
    hops = sum(counts)
    if hops==1:
        return abs(1+(phase-1)/8)**2
    total = 0.
    for last in range(len(counts)):
        amplitude = 0j
        for first in range(len(counts)):
            remaining = counts.copy()
            remaining[first] -= 1
            remaining[last] -= 1
            if min(remaining)<0:
                continue
            paths = math.factorial(hops-2)//math.prod(math.factorial(c) for c in remaining)
            amplitude += paths*(int(first==0)+(phase-1)/8)
        total += abs(amplitude*((1+phase)/8)**(hops-1))**2
    return total


def tag(key,session,role,payload,length):
    return hmac.digest(key,session+role+payload,'sha256')[:length]


def frame(key,session,role,nonce,spec):
    return nonce+tag(key,session,role,nonce,spec['tag_bytes'])


def valid(data,key,session,role,spec):
    if data is None or len(data)!=spec['nonce_bytes']+spec['tag_bytes']:
        return False
    return hmac.compare_digest(data,frame(key,session,role,data[:spec['nonce_bytes']],spec))


def decode(hits):
    require(hits.ndim==2 and hits.shape[1]==2 and hits.shape[0]>0 and hits.shape[0]%8==0 and
            hits.dtype==np.int16 and np.all(hits>=0),
            'malformed receiver histograms')
    occupied = hits>0
    if not np.all(occupied.sum(axis=1)==1):
        return None
    return np.packbits(occupied[:,1]).tobytes()


def forecast(spec,signal):
    require(math.isfinite(signal) and signal>0,'a positive measured signal speed is required')
    result = {}
    for i,route in enumerate(spec['routes']):
        measured = []
        for axis,d in legs(route):
            at,forward,_ = addressed((0,0,0),[(axis,d)],256,'range')
            end,back,_ = addressed(at,[(axis,-d)],256,'echo')
            require(end==(0,0,0),'range echo missing')
            measured.append((forward+back)/2)
        distance = math.hypot(*measured)
        for member in spec['access']:
            segments = itinerary(route,member,spec['segment'])
            hops = sum(abs(d) for d in route)
            budget = 1 if member=='addressed' else spec['copies']
            t = 2*hops+(2*len(segments)-1)*spec['processing']
            both = segments+[tuple(-d for d in s) for s in reversed(segments)]
            probs = [1. if member=='addressed' else front_formula(s,1j if member=='phase-i' else 1) for s in both]
            bits = 8*(spec['nonce_bytes']+spec['tag_bytes'])
            def detected(p):
                return 1. if p==1 else -math.expm1(budget*math.log1p(-p))
            success = math.prod(detected(p)**bits for p in probs)
            result[f'{i}_{member}'] = dict(segments=len(segments),front_probabilities=probs,
                copies=budget,elapsed=t,range=distance,scalar_roundtrip=2*distance/signal,
                beams_per_completed_trip=len(both)*bits*budget,
                ideal_acceptance=success if member not in ('tamper','stale','silent','blocked') else 0.)
    return result


def acquire(spec,sealed,seed):
    rng = np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(b'counts'+seed).digest(),'big')))
    key = hashlib.sha256(b'paired endpoints'+seed).digest()
    raw,summary = {},{}
    for i,route in enumerate(spec['routes']):
        for member in spec['access']:
            label = f'{i}_{member}'
            f = sealed[label]
            copies = 1 if member=='addressed' else spec['copies']
            segments = itinerary(route,member,spec['segment'])
            all_legs = segments+[tuple(-d for d in s) for s in reversed(segments)]
            gates = {leg:(1. if member=='addressed' else front(leg,1j if member=='phase-i' else 1))
                     for leg in set(all_legs)}
            middle = len(segments)
            bits = 8*(spec['nonce_bytes']+spec['tag_bytes'])
            tape = np.zeros((spec['trials'],len(all_legs),bits,2),dtype=np.int16)
            launches = 0
            stamps = []
            for trial in range(spec['trials']):
                session = f'{i}/{trial}'.encode()+b':'
                # Matched access controls get the same independently chosen message.
                nonce = hashlib.sha256(b'sender only'+seed+session).digest()[:spec['nonce_bytes']]
                word = frame(key,session if member!='stale' else b'old:'+session,b'?',nonce,spec)
                if member=='silent':
                    word = None
                elapsed = 0
                site = (0,0,0)
                for j,leg in enumerate(all_legs):
                    at,travel,delivered = addressed(site,legs(leg),256,word if member=='addressed' else None)
                    # On the wave channels only detected counts carry the word to the next station.
                    # The address traversal above determines arrival time and station, not its payload.
                    site = at
                    elapsed += travel
                    if member=='addressed':
                        word = delivered
                    if word is not None:
                        rails = np.unpackbits(np.frombuffer(word,dtype=np.uint8))
                        launches += bits*copies
                        p = 0 if member=='blocked' and j==middle-1 else gates[leg]
                        hits = rng.binomial(copies,p,size=bits)
                        tape[trial,j,np.arange(bits),rails] = hits
                    word = decode(tape[trial,j])
                    if member=='tamper' and j==middle-1 and word is not None:
                        word = word[:-1]+bytes([word[-1]^1])
                    if j==middle-1:
                        # The remote endpoint has just its local receipt, setup key and session.
                        word = (frame(key,session,b'!',word[:spec['nonce_bytes']],spec)
                                if valid(word,key,session,b'?',spec) else None)
                    if j<len(all_legs)-1:
                        elapsed += spec['processing']
                require(site==(0,0,0) and elapsed==f['elapsed'],'invalid reunion schedule')
                accepted = valid(word,key,session,b'!',spec) and word[:spec['nonce_bytes']]==nonce
                stamps.append(elapsed if accepted else -1)
            raw[label] = tape
            summary[label] = dict(accepted=sum(t>=0 for t in stamps),ahead=sum(0<=t<f['scalar_roundtrip'] for t in stamps),
                                  stamps=stamps,launched=launches)
    return raw,summary
