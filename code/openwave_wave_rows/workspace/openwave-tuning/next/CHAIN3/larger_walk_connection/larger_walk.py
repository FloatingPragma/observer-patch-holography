"""Padded arbitrary-size LOCAL12 massive walk and declared phase controls.

The array evaluator implements the certified target; native_word.py emits
actual rational native controls with an additive full-response error bound.
No read occurs inside step. This is the scalar mass-coin source, not Dirac.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'openwave-tuning/common_dynamics/particles'),
                str(ROOT/'openwave-tuning/next/CHAIN2/local12_limit_audit')]
import mass_coin as mass_source
import refinement_adapter as adapter

SOURCE_ID = 'chain3-larger-local12-mass-phase-v1'


def geometry_spec():
    """Immutable translation data shared with static-field consumers."""
    return (tuple(tuple(edge) for row in adapter.common.INTERNAL_MATCHINGS for edge in row),
            tuple(tuple(tuple(pair) for pair in axis) for axis in adapter.spatial.links()))


def verify_sources():
    pins = json.loads((HERE/'SOURCES.json').read_text())['files_sha256']
    for path, expected in pins.items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != expected:
            raise ValueError('frozen source changed: '+path)
    return pins


@dataclass(frozen=True)
class Geometry:
    size: int
    graph: object
    position_dimension: int
    addresses: np.ndarray
    invalid_addresses: np.ndarray
    endpoints: np.ndarray
    slots: np.ndarray
    left: np.ndarray
    right: np.ndarray
    neighbours: np.ndarray
    edge_ids: np.ndarray

    @property
    def shape(self): return (self.position_dimension, 16)


@lru_cache(None)
def geometry(size):
    info = adapter.layout(size)
    g = adapter.spatial.graph(size)
    pdim = 1 << (info['data_bits']-4)
    addresses = 16*(np.arange(g.ports)//12)+np.arange(g.ports)%12
    valid = np.zeros(pdim, bool); valid[addresses] = True
    endpoints = np.asarray(g.seams, dtype=np.int64)
    colour = {tuple(edge): c for c,row in enumerate(adapter.common.INTERNAL_MATCHINGS) for edge in row}
    directions = np.array([colour[(int(a%12),int(b%12))] if a//12==b//12 else 5
                           for a,b in endpoints], dtype=np.int64)
    slots = 2*directions+1
    left = 16*addresses[endpoints[:,0]]+slots
    right = 16*addresses[endpoints[:,1]]+slots
    neighbours = np.empty((g.ports,6),dtype=np.int64)
    edge_ids = np.empty_like(neighbours)
    for e, ((a,b),d) in enumerate(zip(endpoints,directions)):
        neighbours[a,d]=b; neighbours[b,d]=a
        edge_ids[a,d]=e; edge_ids[b,d]=e
    arrays = (addresses,np.flatnonzero(~valid),endpoints,slots,left,right,neighbours,edge_ids)
    for value in arrays: value.flags.writeable=False
    return Geometry(size,g,pdim,*arrays)


def _real_vector(value, length, name):
    if value is None:
        out=np.zeros(length);out.flags.writeable=False;return out
    out = np.asarray(value)
    if out.shape != (length,) or np.iscomplexobj(out):
        raise ValueError(name+' requires one real value per graph entry')
    out = np.array(out,dtype=float,copy=True)
    if not np.isfinite(out).all(): raise ValueError(name+' must be finite')
    out.flags.writeable=False
    return out


class Walk:
    """W=S_theta C_mu P_V; phase theta is positive from seam a to b (a<b).

    P_V=exp(-i*h*V) acts on all16 event slots at each valid port. It is
    identity on invalid position addresses. S_theta has conjugate phases
    on reverse moves, and is identity on hold, event-padding and source-
    padding slots. C_mu acts on every position address, including padding.
    """
    def __init__(self,size,h=F(1),mass=F(0),seam_phase=None,port_potential=None):
        self.h,self.mass=mass_source.parameters(h,mass)
        self.geometry=geometry(size)
        self.mu=float(self.h*self.mass)
        self.seam_phase=_real_vector(seam_phase,len(self.geometry.graph.seams),'seam phase')
        self.port_potential=_real_vector(port_potential,self.geometry.graph.ports,'port potential')
        self.link=np.exp(1j*self.seam_phase)
        self.potential=np.ones(self.geometry.position_dimension,complex)
        self.potential[self.geometry.addresses]=np.exp(-1j*float(self.h)*self.port_potential)
        self.link.flags.writeable=False;self.potential.flags.writeable=False

    def _state(self,state):
        out=np.asarray(state,dtype=complex)
        if out.shape[-2:]!=self.geometry.shape or not np.isfinite(out).all():
            raise ValueError('complete finite padded state required; final axes must be '+str(self.geometry.shape))
        return out

    def coin(self,state,adjoint=False):
        value=self._state(state)
        mu=-self.mu if adjoint else self.mu
        result=(-np.cos(mu)-1j*np.sin(mu))*value
        result[...,:12]+=np.cos(mu)*np.sum(value[...,:12],axis=-1,keepdims=True)/6
        return result

    def route(self,state):
        value=self._state(state)
        out=value.copy().reshape(value.shape[:-2]+(-1,))
        src=value.reshape(value.shape[:-2]+(-1,))
        g=self.geometry
        out[...,g.right]=self.link*src[...,g.left]
        out[...,g.left]=self.link.conjugate()*src[...,g.right]
        return out.reshape(value.shape)

    def step(self,state):
        value=self._state(state)
        return self.route(self.coin(value*self.potential[:,None]))

    def adjoint(self,state):
        return self.coin(self.route(state),adjoint=True)*self.potential.conjugate()[:,None]

    def current(self,state):
        """Mean oriented move current, evaluated after the local coin."""
        value=self.coin(self._state(state)*self.potential[:,None])
        flat=value.reshape(value.shape[:-2]+(-1,)); g=self.geometry
        return abs(flat[...,g.left])**2-abs(flat[...,g.right])**2

    def divergence(self,current):
        value=np.asarray(current); g=self.geometry
        if value.shape!=(len(g.graph.seams),): raise ValueError('one current per seam required')
        out=np.zeros(g.graph.ports,dtype=value.dtype)
        np.add.at(out,g.endpoints[:,0],value); np.add.at(out,g.endpoints[:,1],-value)
        return out

    def field_step(self,state,modulus=None):
        """Coherent current-register lift; field changes E_e -> E_e-j_e.

        Keys are (padded position, complete integer seam field), values are
        length16 event vectors. With B[a,e]=+1, G=B E-rho is conserved.
        A finite modulus preserves G modulo q. This register has no field
        Hamiltonian, force law, or source-independent Gaussian closure.
        """
        if modulus is not None and (type(modulus)is not int or modulus<2):raise ValueError('integer field modulus at least2 required')
        g=self.geometry;out={}
        for (position,field),events in state.items():
            if type(position)is not int or not 0<=position<g.position_dimension or len(field)!=len(g.graph.seams) or any(type(x)is not int for x in field):raise ValueError('complete padded position and integer seam field required')
            events=np.asarray(events,complex)
            if events.shape!=(16,) or not np.isfinite(events).all():raise ValueError('complete16 event amplitudes required')
            v=(-np.cos(self.mu)-1j*np.sin(self.mu))*events
            v[:12]+=np.cos(self.mu)*np.sum(events[:12])/6
            v*=self.potential[position]
            for event,amplitude in enumerate(v):
                record=self.raw_record(position,event)
                destination=position;target_field=field;phase=1
                if record['admitted'] and event%2:
                    # raw_record at the input reverses the same involutive
                    # route, so its current sign is opposite the outgoing one.
                    edge=record['edge_index'];current=-record['current']
                    destination=int(g.addresses[record['source_before']])
                    changed=list(field);changed[edge]-=current
                    if modulus is not None:changed[edge]%=modulus
                    target_field=tuple(changed);phase=np.exp(1j*current*self.seam_phase[edge])
                key=(destination,target_field)
                if key not in out:out[key]=np.zeros(16,complex)
                out[key][event]+=phase*amplitude
        return out

    def records(self,state):
        p=abs(self._state(state))**2; g=self.geometry
        if p.ndim!=2: raise ValueError('one state per record law required')
        return {'joint_source_event_probabilities':p,
                'valid_position_probabilities':p[g.addresses,:12].sum(axis=1),
                'source_padding_probability':float(p[g.invalid_addresses].sum()),
                'event_padding_probability_on_valid_source':float(p[g.addresses,12:].sum()),
                'original_norm_squared':float(p.sum()),'renormalized':False}

    def raw_record(self,final_position,event,law='C'):
        g=self.geometry
        if type(final_position)is not int or not 0<=final_position<g.position_dimension or type(event)is not int or not 0<=event<16 or law not in ('C','E'):
            raise ValueError('complete padded final address, event and C/E law required')
        cell,port=divmod(final_position,16)
        if cell>=g.graph.cells or port>=12 or event>=12:
            return {'admitted':False,'reason':'source_padding' if cell>=g.graph.cells or port>=12 else 'event_padding'}
        after=12*cell+port; direction,sigma=divmod(event,2)
        before=int(g.neighbours[after,direction]) if sigma else after
        edge=int(g.edge_ids[after,direction]); a,b=map(int,g.endpoints[edge])
        return {'admitted':True,'edge_index':edge,'seam':[a,b],
                'source_before':before,'source_after':after,'direction':direction,
                'recoded_sigma':sigma,'law':law,'coin':sigma^int(before==a),
                'current':int(before==a)-int(after==a),
                'original_E_contrast':(1 if before==a else -1) if law=='E' else None,
                'fresh_uniform_attempt_claimed':False}

    def initial_cell(self,cell=0):
        if type(cell)is not int or not 0<=cell<self.geometry.graph.cells: raise ValueError('valid cell required')
        out=np.zeros(self.geometry.shape,complex);out[16*cell:16*cell+12,:12]=1/12
        return out

    def gaussian_packet(self,center,width,momentum=(0,0,0)):
        if width<=0:raise ValueError('positive width required')
        g=self.geometry; center=np.asarray(center,float); momentum=np.asarray(momentum,float)
        if center.shape!=(3,) or momentum.shape!=(3,) or not np.isfinite(center).all() or not np.isfinite(momentum).all(): raise ValueError('finite three-vectors required')
        delta=(np.asarray(g.graph.cell_coordinates)-center+g.size/2)%g.size-g.size/2
        amplitudes=np.exp(-np.sum(delta**2,axis=1)/(4*width**2)+1j*delta@momentum)
        amplitudes/=np.linalg.norm(amplitudes)
        out=np.zeros(g.shape,complex)
        out[g.addresses,:12]=np.repeat(amplitudes,12)[:,None]/12
        return out

    def trajectory(self,state,steps,read_steps=(),channel='terminal_only',seed=None):
        """Selective trajectory, with every requested pointer retained.

        event: project event, retain it; position: project padded position;
        joint: project both. No implicit event reset. A fresh-register
        channel is a different instrument and is deliberately explicit below.
        """
        if type(steps)is not int or steps<0:raise ValueError('nonnegative steps required')
        schedule=tuple(read_steps)
        if len(set(schedule))!=len(schedule) or any(type(n)is not int or not 1<=n<=steps for n in schedule): raise ValueError('distinct integer read steps in1..steps required')
        if channel not in ('terminal_only','event','position','joint','event_reset'):raise ValueError('declared measurement channel required')
        if channel=='terminal_only' and schedule:raise ValueError('terminal_only has no intermediate read steps')
        state=self._state(state).copy()
        if state.ndim!=2 or abs(np.linalg.norm(state)-1)>1e-10:raise ValueError('normalized pure trajectory input required')
        rng=np.random.default_rng(seed); records=[]
        for n in range(1,steps+1):
            state=self.step(state)
            if n not in schedule:continue
            p=abs(state)**2
            weights=p.sum(axis=0) if channel in ('event','event_reset') else p.sum(axis=1) if channel=='position' else p.reshape(-1)
            choice=int(rng.choice(len(weights),p=weights/weights.sum())); probability=float(weights[choice])
            selected=np.zeros_like(state)
            if channel in ('event','event_reset'):selected[:,choice]=state[:,choice]
            elif channel=='position':selected[choice]=state[choice]
            else:selected.reshape(-1)[choice]=state.reshape(-1)[choice]
            state=selected/np.sqrt(probability)
            if channel=='event_reset':
                position=state[:,choice].copy();state[:]=0;state[:,:12]=position[:,None]/np.sqrt(12)
            record={'step':n,'channel':channel,'outcome':choice,'probability':probability}
            if channel=='joint':record['native_C_record']=self.raw_record(choice//16,choice%16,'C')
            records.append(record)
        return {'state':state,'intermediate_records':records,'terminal_law':self.records(state),
                'read_steps':list(schedule),'read_rate_per_macro':len(schedule)/steps if steps else 0,
                'channel':channel,'seed':seed}


def resource_estimate(size):
    info=adapter.layout(size);dim=1<<info['data_bits']
    return {'size':size,'state_complex_amplitudes':dim,'complex128_state_bytes':16*dim,
            'physical_ports':12*size**3,'physical_seams':36*size**3,
            'matrix_free_arithmetic_per_step':'O(256*2^ceil(log2(L^3))); no dense walk matrix',
            'native_free_word':adapter.resource_upper_bound(size),
            'graph_construction':'O(L^3) graph storage and one breadth-first connectivity check',
            'native_compiler_materializes_reversible_gates':True,
            'native_word_execution_claimed_for_large_sizes':False}
