"""Exact certification of a proposed fair-source spectral dictionary."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

import numpy as np

from tools import strict_json
from model.native_energy_identification import energy
from model.native_effect_synthesis import synthesis

DIM=66
GRID=2**44
FLOOR=F(2,3)
GRAM_GATE=F(1,10**9)
RESIDUAL_GATE=F(1,10**9)
PROPOSAL=Path(__file__).with_name('proposal.json')


@lru_cache(maxsize=1)
def source_counts():
    """Native fair population columns in the existing coherent source chart."""
    source=energy.exchange.sector.source
    family=source.local_family(source.local.Parameters(port_weights=(1,)*12))
    states=source.local.basis()
    if tuple(states)!=tuple(synthesis.tomo.compiled.quantum.STATES):
        raise ArithmeticError('fair occupation and coherent source coordinates differ')
    masks=tuple(sum(int(x)<<i for i,x in enumerate(state)) for state in states)
    where={mask:i for i,mask in enumerate(masks)}
    columns=[]
    for mask in masks:
        column=source.population_column(family,mask)
        row=[]
        for target in masks:
            count=60*column.get(target,F())
            if count.denominator!=1:
                raise ArithmeticError('fair native endpoint probabilities left the60-count chart')
            row.append(int(count))
        columns.append(row)
    rows=tuple(tuple(columns[j][i] for j in range(DIM)) for i in range(DIM))
    if any(sum(row)!=60 or min(row)<0 for row in rows):
        raise ArithmeticError('native fair transfer ceased to be stochastic')
    if any(rows[i][j]!=rows[j][i] for i in range(DIM) for j in range(DIM)):
        raise ArithmeticError('native fair transfer ceased to be symmetric')
    if min(2*F(rows[i][i],60)-1 for i in range(DIM))<FLOOR:
        raise ArithmeticError('native diagonal lower bound lost its positive transfer floor')
    if len(where)!=DIM or len(set(synthesis.tomo.codes()))!=DIM:
        raise ArithmeticError('the exact source chart lost injectivity')
    return rows


def load(path=PROPOSAL):
    value=strict_json.load(path)
    if (set(value)!={'schema','denominator','columns','eigenvalues'}
            or value['schema']!='oph.native_energy_measurement.proposal.v1'
            or type(value['denominator']) is not int or value['denominator']!=GRID):
        raise ValueError('the exact fair spectral proposal schema is required')
    columns=value['columns']; values=value['eigenvalues']
    if (type(columns) is not list or len(columns)!=DIM or type(values) is not list
            or len(values)!=DIM or any(type(c) is not list or len(c)!=DIM
            or any(type(x) is not int or abs(x)>GRID for x in c) for c in columns)
            or any(type(x) is not int or not FLOOR*GRID<=x<=GRID for x in values)):
        raise ValueError('all66 bounded exact source rays and positive transfer labels required')
    return tuple(tuple(F(x,GRID) for x in c) for c in columns),tuple(F(x,GRID) for x in values)


def certify(columns=None,values=None):
    if columns is None and values is None:columns,values=load()
    if (type(columns) is not tuple or type(values) is not tuple
            or len(columns)!=DIM or len(values)!=DIM
            or any(type(c) is not tuple or len(c)!=DIM or any(type(x) is not F for x in c) for c in columns)
            or any(type(x) is not F or not FLOOR<=x<=1 for x in values)):
        raise ValueError('the full exact source spectral components are required')
    counts=source_counts()
    gram=[[sum((a*b for a,b in zip(columns[i],columns[j],strict=True)),F())-int(i==j)
           for j in range(DIM)] for i in range(DIM)]
    gram_bound=max(sum(map(abs,row),F()) for row in gram)
    residual=[[sum((F(counts[i][k],60)*c[k] for k in range(DIM) if counts[i][k]),F())-lam*c[i]
               for c,lam in zip(columns,values,strict=True)] for i in range(DIM)]
    residual_bound=max(max(sum(map(abs,row),F()) for row in residual),
                       max(sum(map(abs,col),F()) for col in zip(*residual,strict=True)))
    if gram_bound>GRAM_GATE or residual_bound>RESIDUAL_GATE:
        raise ArithmeticError('rational proposal fails its fair-source residual certificate')
    return {'dimension':DIM,'gram_operator_upper':str(GRAM_GATE),
            'source_residual_operator_upper':str(RESIDUAL_GATE),
            'transfer_floor':str(FLOOR),
            'nearby_transfer_operator_error_upper':str(RESIDUAL_GATE+2*GRAM_GATE),
            'counted_energy_operator_error_upper':str(60*(RESIDUAL_GATE+2*GRAM_GATE)/FLOOR),
            'ideal_complete_measurement_halfdiamond_upper':str(6*DIM*GRAM_GATE),
            'polar_basis_is_executed':False,'eigenvalue_gap_required':False}


def energy_labels(values=None):
    if values is None:_,values=load()
    return tuple(tuple(60*x for x in energy.negative_log_interval(value,96)) for value in values)


def proposed_arrays():
    columns,values=load()
    return np.array(columns,float).T,np.array(values,float)
