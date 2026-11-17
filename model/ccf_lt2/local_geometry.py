#!/usr/bin/env python3
"""Bounded local lapse coupled to complete gauge-invariant source energy shards."""
from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
import sympy as sy

try:
    from .geometry import diagnostic_source
except ImportError:
    from .geometry import diagnostic_source

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TOL=1e-10


@dataclass(frozen=True)
class LocalGeometryConfig:
    phi_values: tuple[float,...]=(-1.0,0.0,1.0)
    hopping: float=.5
    stiffness: float=.5
    lapse_coupling: float=.25
    regulator: float=1.0
    boundary_weights: tuple[float,...]=()
    source_energy_references: tuple[float,...]=()


@dataclass
class LocalGeometryExtension:
    hamiltonian: sparse.csr_matrix
    geometry_hamiltonian: sparse.csr_matrix
    source_shards: tuple[sparse.csr_matrix,...]
    lapses: tuple[sparse.csr_matrix,...]
    phi: tuple[sparse.csr_matrix,...]
    labels: tuple[tuple[float,...],...]
    graph_laplacian: np.ndarray
    config: LocalGeometryConfig

    @property
    def auxiliary_dimension(self):
        return len(self.labels)


def _norm(a):
    return float(sparse.linalg.norm(a)) if sparse.issparse(a) else float(np.linalg.norm(a))


def _exact(value):
    z=complex(value)
    return sy.Rational(z.real)+sy.I*sy.Rational(z.imag)


def _graph(graph):
    n=int(graph['nvertices'] if isinstance(graph,dict) else graph.nvertices)
    edges=graph['edges'] if isinstance(graph,dict) else graph.edges
    normalized=[]
    for edge in edges:
        if len(edge) not in (2,3):
            raise ValueError('graph edges must give two vertices and optional nonnegative weight')
        u,v=int(edge[0]),int(edge[1]);weight=float(edge[2]) if len(edge)==3 else 1.0
        if not 0<=u<n or not 0<=v<n or u==v or not np.isfinite(weight) or weight<0:
            raise ValueError('invalid geometry graph edge')
        normalized.append((min(u,v),max(u,v),weight))
    if n<1 or len({(u,v) for u,v,_ in normalized})!=len(normalized):
        raise ValueError('nonempty graph with no duplicate undirected edges required')
    return n,tuple(sorted(normalized))


def _references(n,cfg):
    result=cfg.source_energy_references or (0.0,)*n
    if len(result)!=n or not np.isfinite(result).all():
        raise ValueError('one finite energy reference per source vertex required')
    return tuple(result)


@lru_cache(maxsize=16)
def _auxiliary(n,edges,cfg):
    values=np.asarray(cfg.phi_values,float)
    params=(cfg.hopping,cfg.stiffness,cfg.lapse_coupling,cfg.regulator)
    if (len(values)<2 or not np.isfinite(values).all() or not np.all(np.diff(values)>0)
        or not all(np.isfinite(params)) or min(params)<0 or cfg.stiffness<=0
        or min(1+cfg.lapse_coupling*values)<=0):
        raise ValueError('finite ordered coordinates, positive lapse and nonnegative geometry couplings required')
    weights=cfg.boundary_weights or (0.0,)*n
    if len(weights)!=n or not np.isfinite(weights).all() or min(weights)<0:
        raise ValueError('one nonnegative boundary weight per geometry vertex required')
    _references(n,cfg)
    lap=np.diag(np.array(weights)+cfg.regulator)
    for u,v,w in edges:
        lap[u,u]+=w;lap[v,v]+=w;lap[u,v]-=w;lap[v,u]-=w
    labels=tuple(product(tuple(float(x) for x in values),repeat=n))
    dimension=len(labels);table=np.asarray(labels)
    potential=(cfg.stiffness/2)*np.einsum('ai,ij,aj->a',table,lap,table)
    hg=sparse.diags(potential,format='csr')
    count=len(values);local=np.zeros((count,count))
    for i in range(count-1):
        local[i,i]+=1;local[i+1,i+1]+=1;local[i,i+1]=local[i+1,i]=-1
    for vertex in range(n):
        operator=sparse.csr_matrix([[1.]])
        for k in range(n):
            operator=sparse.kron(operator,local if k==vertex else sparse.eye(count),format='csr')
        hg=hg+cfg.hopping*operator
    phi=tuple(sparse.diags(table[:,v],format='csr') for v in range(n))
    identity=sparse.eye(dimension,format='csr')
    lapse=tuple(identity+cfg.lapse_coupling*q for q in phi)
    def columns(matrix):
        matrix=matrix.tocsc();result=[]
        for a in range(dimension):
            lo,hi=matrix.indptr[a:a+2]
            result.append(tuple((int(b),_exact(value)) for b,value in
                zip(matrix.indices[lo:hi],matrix.data[lo:hi]) if value!=0))
        return tuple(result)
    return hg.tocsr(),phi,lapse,labels,lap,columns(hg),tuple(columns(p) for p in lapse)


def geometry_factors(graph,config=LocalGeometryConfig()):
    """Finite geometry factors only; no source state space is enumerated."""
    n,edges=_graph(graph)
    hg,phi,lapse,labels,lap,_,_=_auxiliary(n,edges,config)
    return {'H_geometry':hg,'phi':phi,'lapses':lapse,'labels':labels,
            'graph_laplacian':lap,'dimension':len(labels),'nvertices':n}


def _ordered(shards,n):
    if hasattr(shards,'keys'):
        if set(shards)!=set(range(n)):
            raise ValueError('energy shards must include every graph vertex exactly once')
        return tuple(shards[v] for v in range(n))
    if len(shards)!=n:
        raise ValueError('one complete energy shard per graph vertex required')
    return tuple(shards)


def apply_local_source_column(h_shard_actions,source_key,auxiliary_index,graph,
                              config=LocalGeometryConfig()):
    """Apply all full native energy-shard columns without source projection."""
    n,edges=_graph(graph);actions=_ordered(h_shard_actions,n)
    _,_,_,labels,_,acols,ncols=_auxiliary(n,edges,config)
    if not isinstance(auxiliary_index,(int,np.integer)) or not 0<=auxiliary_index<len(labels):
        raise ValueError('auxiliary index outside the complete local geometry carrier')
    hash(source_key);answer={};refs=_references(n,config)
    def add(key,value):
        hash(key);answer[key]=answer.get(key,sy.Integer(0))+value
    for row,value in acols[auxiliary_index]:
        add((row,source_key),value)
    for vertex,action in enumerate(actions):
        source=action(source_key)
        if not hasattr(source,'items'):
            raise TypeError('each source shard must return a complete coefficient mapping')
        for row,value in ncols[vertex][auxiliary_index]:
            for target,coefficient in source.items():
                add((row,target),value*sy.sympify(coefficient))
            add((row,source_key),-_exact(refs[vertex])*value)
    add((int(auxiliary_index),source_key),sum((_exact(x) for x in refs),sy.Integer(0)))
    answer={key:sy.expand(value) for key,value in answer.items()}
    return {key:value for key,value in answer.items() if value!=0}


def _hermitian(matrix,name):
    matrix=sparse.csr_matrix(matrix,dtype=complex)
    if matrix.shape[0]!=matrix.shape[1] or not matrix.shape[0] or not np.isfinite(matrix.data).all() or _norm(matrix-matrix.getH())>TOL:
        raise ValueError(name+' must be a finite Hermitian matrix')
    return matrix


def build_local_geometry_extension(source_shards,graph,gauss_operators=(),
                                    config=LocalGeometryConfig(),source_hamiltonian=None):
    n,edges=_graph(graph)
    shards=tuple(_hermitian(x,'source energy shard') for x in _ordered(source_shards,n))
    dimension=shards[0].shape[0]
    if any(h.shape!=(dimension,dimension) for h in shards):
        raise ValueError('source shards must act on the same complete carrier')
    total=sum(shards,sparse.csr_matrix((dimension,dimension)))
    if source_hamiltonian is not None and _norm(total-_hermitian(source_hamiltonian,'full source'))>TOL:
        raise ValueError('energy shards omit or alter the supplied full source')
    gauss=tuple(_hermitian(g,'Gauss generator') for g in gauss_operators)
    if any(g.shape!=shards[0].shape or _norm(g@h-h@g)>TOL for g in gauss for h in shards):
        raise ValueError('each local energy shard must preserve every Gauss generator')
    hg,phi,lapses,labels,lap,_,_=_auxiliary(n,edges,config)
    identity=sparse.eye(dimension,format='csr');refs=_references(n,config)
    h=sparse.kron(hg,identity,format='csr')
    for lapse,shard,ref in zip(lapses,shards,refs):
        h=h+sparse.kron(lapse,shard-ref*identity,format='csr')
    h=h+sum(refs)*sparse.eye(h.shape[0],format='csr')
    return LocalGeometryExtension(h.tocsr(),hg,shards,lapses,phi,labels,lap,config)


def complete_records(extension,state,source_projectors):
    d=extension.source_shards[0].shape[0]
    psi=np.asarray(state,complex).reshape(-1)
    if len(psi)!=extension.hamiltonian.shape[0] or abs(np.vdot(psi,psi)-1)>TOL:
        raise ValueError('normalized full joint state required')
    projectors=tuple((str(label),_hermitian(p,'source projector')) for label,p in source_projectors)
    if _norm(sum((p for _,p in projectors),sparse.csr_matrix((d,d)))-sparse.eye(d))>TOL:
        raise ValueError('source instrument is incomplete')
    for i,(_,p) in enumerate(projectors):
        if _norm(p@p-p)>TOL or any(_norm(p@q)>TOL for _,q in projectors[i+1:]):
            raise ValueError('orthogonal source projectors required')
    records=[];source_marginal={label:0. for label,_ in projectors}
    geometry_marginal=np.zeros(extension.auxiliary_dimension);after=0.
    for a,label in enumerate(extension.labels):
        select=np.zeros(extension.auxiliary_dimension);select[a]=1
        for outcome,p in projectors:
            branch=sparse.kron(sparse.diags(select),p,format='csr')@psi
            probability=float(np.vdot(branch,branch).real)
            after+=float(np.vdot(branch,extension.hamiltonian@branch).real)
            geometry_marginal[a]+=probability;source_marginal[outcome]+=probability
            records.append({'geometry':list(label),'source_outcome':outcome,'probability':probability})
    energy=float(np.vdot(psi,extension.hamiltonian@psi).real)
    return {'records':records,'source_marginal':source_marginal,
        'geometry_marginal':geometry_marginal.tolist(),
        'phi_means':(geometry_marginal@np.asarray(extension.labels)).tolist(),
        'charged_plus_probability':sum(p for label,p in source_marginal.items() if label.endswith('_plus')),
        'probability_sum':sum(source_marginal.values()),'total_energy':energy,
        'terminal_measurement_work':after-energy}


def static_comparison():
    lap=sy.Matrix([[2,-1],[-1,2]]);source=sy.Matrix([1,2])
    coupling=sy.Rational(1,4);stiffness=sy.Rational(1,2)
    phi=-(coupling/stiffness)*lap.inv()*source
    energy=(stiffness*(phi.T*lap*phi)[0]/2+coupling*(phi.T*source)[0])
    predicted=-coupling**2*(source.T*lap.inv()*source)[0]/(2*stiffness)
    cross=-coupling**2*lap.inv()[0,1]/stiffness
    return {'scope':'Continuous unbounded commuting-source static comparison; not exact finite quantum elimination.',
        'laplacian_inverse':str(lap.inv()),'stationary_phi':[str(x) for x in phi],
        'energy':str(energy),'predicted_energy':str(predicted),'cross_coefficient':str(cross),
        'stationarity_exact':stiffness*lap*phi+coupling*source==sy.zeros(2,1),
        'energy_identity_exact':energy==predicted,'attractive_cross_term':cross<0}


def witness():
    graph={'nvertices':2,'edges':((0,1),)}
    detached={'nvertices':2,'edges':()}
    h,pieces,q,projectors=diagnostic_source()
    shards={0:pieces['neutral_rest'],1:pieces['charged_rest']}
    configs=(('coupled',graph,LocalGeometryConfig()),
             ('coupling_deleted',graph,replace(LocalGeometryConfig(),lapse_coupling=0.)),
             ('edge_deleted_diagonal_matched',detached,replace(LocalGeometryConfig(),regulator=2.)))
    runs=[];matrices={};max_energy=0.;max_gauss=0.;local_steps=True
    for name,carrier,cfg in configs:
        extension=build_local_geometry_extension(shards,carrier,[q],cfg,h)
        matrices[name]=extension
        # Full-operator locality: geometry kinetic changes one adjacent label.
        coo=extension.geometry_hamiltonian.tocoo()
        for a,b,value in zip(coo.row,coo.col,coo.data):
            if a!=b and value!=0:
                delta=np.abs(np.array(extension.labels[a])-np.array(extension.labels[b]))
                local_steps &= np.count_nonzero(delta)==1 and max(delta)==1
        g=sparse.kron(sparse.eye(extension.auxiliary_dimension),q)
        max_gauss=max(max_gauss,_norm(g@extension.hamiltonian-extension.hamiltonian@g))
        for neutral in (0,1):
            src=np.zeros(8,complex);src[4*neutral]=src[4*neutral+3]=1/np.sqrt(2)
            aux=np.zeros(extension.auxiliary_dimension);aux[extension.labels.index((0.,0.))]=1
            initial=np.kron(aux,src)
            energy=float(np.vdot(initial,extension.hamiltonian@initial).real)
            for time in (0.,.5,1.,2.,4.):
                state=sparse.linalg.expm_multiply(-1j*time*extension.hamiltonian,initial)
                result=complete_records(extension,state,projectors)
                result.update(control=name,neutral_preparation=neutral,time=time,
                              energy_residual=result['total_energy']-energy)
                max_energy=max(max_energy,abs(result['energy_residual']));runs.append(result)
    lookup={(r['control'],r['neutral_preparation'],r['time']):r for r in runs}
    a,b=[lookup['coupled',n,4.] for n in (0,1)]
    clock_response=abs(a['charged_plus_probability']-b['charged_plus_probability'])
    phi_response=[abs(x-y) for x,y in zip(a['phi_means'],b['phi_means'])]
    deleted_clock=max(abs(lookup['coupling_deleted',0,t]['charged_plus_probability']-lookup['coupling_deleted',1,t]['charged_plus_probability']) for t in (0.,.5,1.,2.,4.))
    deleted_geometry=max(abs(x-y) for t in (0.,.5,1.,2.,4.) for x,y in zip(lookup['coupling_deleted',0,t]['phi_means'],lookup['coupling_deleted',1,t]['phi_means']))
    edge_clock=max(abs(lookup['edge_deleted_diagonal_matched',0,t]['charged_plus_probability']-lookup['edge_deleted_diagonal_matched',1,t]['charged_plus_probability']) for t in (0.,.5,1.,2.,4.))
    edge_receiver=max(abs(lookup['edge_deleted_diagonal_matched',0,t]['phi_means'][1]-lookup['edge_deleted_diagonal_matched',1,t]['phi_means'][1]) for t in (0.,.5,1.,2.,4.))
    extension=matrices['coupled'];dim=extension.hamiltonian.shape[0]
    column_error=0.;calls=[0,0]
    def action(vertex):
        def apply(key):
            calls[vertex]+=1
            return {row:_exact(shards[vertex][row,key]) for row in range(8) if shards[vertex][row,key]!=0}
        return apply
    callbacks={v:action(v) for v in range(2)}
    full=extension.hamiltonian.toarray();column_records=[]
    for aux in range(extension.auxiliary_dimension):
        for src in range(8):
            out=apply_local_source_column(callbacks,src,aux,graph)
            column=np.zeros(dim,complex)
            for (b,t),value in out.items():column[b*8+t]=complex(value)
            column_error=max(column_error,float(max(abs(column-full[:,aux*8+src]))))
            column_records.append({'auxiliary':aux,'source':src,'outputs':[
                {'auxiliary':b,'source':t,'coefficient':str(value)} for (b,t),value in out.items()]})
    c=.5;shifted={0:shards[0]+c*np.eye(8),1:shards[1]-c*np.eye(8)}
    fixed=build_local_geometry_extension(shifted,graph,[q],source_hamiltonian=h)
    common=build_local_geometry_extension(shifted,graph,[q],replace(LocalGeometryConfig(),source_energy_references=(c,-c)),h)
    shift_operator_error=_norm(common.hamiltonian-extension.hamiltonian)
    aux=np.zeros(extension.auxiliary_dimension);aux[extension.labels.index((0.,0.))]=1
    src=np.zeros(8,complex);src[0]=src[3]=1/np.sqrt(2);initial=np.kron(aux,src)
    fixed_record=complete_records(fixed,sparse.linalg.expm_multiply(-4j*fixed.hamiltonian,initial),projectors)
    redistribution_change=max(abs(x['probability']-y['probability']) for x,y in zip(a['records'],fixed_record['records']))
    try:
        build_local_geometry_extension({0:shards[0],1:np.zeros((8,8))},graph,[q],source_hamiltonian=h)
        omission=False
    except ValueError:
        omission=True
    static=static_comparison()
    gates={name:bool(value) for name,value in {
        'complete_source_shards_sum':np.array_equal(shards[0]+shards[1],h),
        'missing_source_piece_rejected':omission,
        'positive_lapses':all(min(l.diagonal())>0 for l in extension.lapses),
        'full_hermiticity':max(_norm(e.hamiltonian-e.hamiltonian.getH()) for e in matrices.values())<=TOL,
        'full_gauss_preserved':max_gauss<=TOL,
        'local_geometry_transitions':local_steps,
        'all_symbolic_columns_match':column_error<=TOL,
        'each_shard_called_once_per_column':calls==[dim,dim],
        'complete_record_normalization':max(abs(r['probability_sum']-1) for r in runs)<=TOL,
        'total_energy_conserved':max_energy<=TOL,
        'neutral_changes_distant_clock':clock_response>1e-7,
        'neutral_changes_source_geometry':phi_response[0]>1e-6,
        'neutral_changes_receiving_geometry':phi_response[1]>1e-7,
        'coupling_deletion_removes_response':max(deleted_clock,deleted_geometry)<=TOL,
        'edge_deletion_removes_distant_response':max(edge_clock,edge_receiver)<=TOL,
        'zero_sum_reference_shift_preserves_joint':shift_operator_error<=TOL,
        'zero_sum_fixed_reference_redistribution_changes_records':redistribution_change>1e-7,
        'conditional_static_square_completion':static['stationarity_exact'] and static['energy_identity_exact'] and static['attractive_cross_term'],
        'all_outcomes_retained':all(len(r['records'])==72 for r in runs)}.items()}
    return {'scope':'Bounded local energy/geometry construction with complete tiny-source witness; no physical gravity or master-source execution.',
        'gates':gates,'all_passed':all(gates.values()),'joint_dimension':dim,
        'neutral_distant_clock_response':clock_response,'neutral_phi_responses':phi_response,
        'coupling_deleted_clock_response':deleted_clock,'edge_deleted_clock_response':edge_clock,
        'edge_deleted_receiving_phi_response':edge_receiver,'max_energy_residual':max_energy,
        'max_symbolic_column_error':column_error,'max_gauss_residual':max_gauss,
        'zero_sum_common_reference_shift_residual':shift_operator_error,
        'zero_sum_fixed_reference_record_change':redistribution_change,
        'conditional_static_comparison':{k:bool(v) if isinstance(v,(bool,sy.logic.boolalg.BooleanAtom)) else v for k,v in static.items()},
        'complete_columns':column_records,'runs':runs,'redistributed_complete_record':fixed_record}

