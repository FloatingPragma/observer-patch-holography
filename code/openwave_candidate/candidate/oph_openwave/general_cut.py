"""Exact general-cut record identity on the actual canonical ring geometry.

Pure incidence arithmetic, without quantum evolution or acceptance execution.
A coherent modular state has record residues, not a measured integer trajectory.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp

from .geometry import flat_ring


def objects():
    c=flat_ring(1,1)
    # Every half-edge is oriented from port to occupied-transit vertex.
    columns=[]
    for edge,(a,b) in enumerate(c.seams):
        columns.extend(((a,c.ports+edge),(b,c.ports+edge)))
    incidence=sp.zeros(c.ports+len(c.seams),len(columns))
    for col,(source,target) in enumerate(columns):incidence[source,col]=1;incidence[target,col]=-1
    return c,columns,incidence


def cut_coefficients(incidence,region):
    indicator=sp.Matrix([[int(i in region) for i in range(incidence.rows)]])
    return tuple(int(x) for x in indicator*incidence)


def raw_gauss(incidence,position,records):
    occupation=sp.zeros(incidence.rows,1);occupation[position]=1
    return occupation+incidence*sp.Matrix(records)


def run():
    c,columns,a=objects();n=a.rows;N=6
    # Symbolic coefficients for arbitrary region S and arbitrary F: summing
    # each column yields s_source-s_target, so internal edges cancel exactly.
    symbols=sp.symbols('s:'+str(n))
    s=sp.Matrix([symbols]);symbolic=s*a
    symbolic_exact=all(sp.expand(symbolic[e]-(symbols[source]-symbols[target]))==0
                       for e,(source,target) in enumerate(columns))
    if not symbolic_exact:raise ValueError('symbolic cut cancellation failed')
    regions=[set(range(12)),set(range(c.ports)),{0},set(range(0,n,2))]
    checks=[]
    for region in regions:
        coeff=cut_coefficients(a,region)
        if any(value not in (-1,0,1) for value in coeff):raise ValueError('cut orientation failed')
        # Exact coefficient identity proves Q_S=G_S-(s A)F for all one-load
        # basis states, not just a sampled register vector or special bridge.
        checks.append({'region_vertices':sorted(region),'boundary_half_edges':sum(bool(x) for x in coeff),
                       'internal_cancellation_exact':True})
    cycle=[0]*(2*len(c.seams))
    for source,target in ((0,1),(1,5),(5,0)):
        edge=next(i for i,pair in enumerate(c.seams) if set(pair)=={source,target})
        first,second=c.seams[edge];sign=1 if first==source else -1
        cycle[2*edge]+=sign;cycle[2*edge+1]-=sign
    zero=[0]*len(cycle);wound=[N*x for x in cycle]
    g0=raw_gauss(a,0,zero);gw=raw_gauss(a,0,wound)
    if g0!=gw or any(x%N for x in wound):raise ValueError('winding alias control failed')
    cut_loop_sums=[sum(x*y for x,y in zip(cut_coefficients(a,region),wound)) for region in regions]
    if any(cut_loop_sums):raise ValueError('closed-boundary loop cancellation failed')
    # Aliasing of distinct integer Gauss sectors is a separate issue.
    shifted=zero.copy();shifted[0]=N;g_shift=raw_gauss(a,0,shifted)
    if g_shift==g0 or tuple(int(x)%N for x in g_shift)!=tuple(int(x)%N for x in g0):
        raise ValueError('sector alias control failed')

    return {'schema':'oph.general_cut.exact.v1','geometry':'canonical flat_ring(1,1)',
      'ports':c.ports,'seams':len(c.seams),'half_incidence_shape':list(a.shape),
      'symbolic_region_identity':symbolic_exact,'identity':'Q_S = G_S - sum_h (s_source-s_target) F_h; modulo N when only residues are specified',
      'arbitrary_quantum_state_scope':'A diagonal operator identity extends to every density matrix supported in the declared sector; probabilities of modular charge outcomes are pushforwards, not inferred histories.',
      'region_controls':checks,
      'cycle_alias':{'modulus':N,'same_modular_record':True,'same_integer_gauss_sector':g0==gw,
                    'integer_half_record_difference':wound,'closed_boundary_current_differences':cut_loop_sums},
      'sector_alias':{'same_modular_gauss':True,'different_integer_gauss':g_shift!=g0,
                      'scope':'Excluded only if a fixed integer-lift sector is independently specified.'},
      'integer_boundary_current':'In a fixed one-load integer Gauss sector initialized with zero records, net complete-boundary current is initial regional occupancy minus final regional occupancy. It is determined by endpoints; seam circulation winding is not.',
      'scope':'Exact general-cut operator identity on the canonical ring; modular alias controls, no physical charge identification.',
      'checks':{'symbolic_identity':bool(symbolic_exact), 'cycle_boundary_cancels':not any(cut_loop_sums)},
      'mutations':{'integer_sector_alias':bool(g_shift!=g0), 'integer_winding_alias':any(wound)},
      'full_criterion_met':False,
      'quantum_evolutions_executed':0,'reserved_configurations_executed':False,'official_rows_promoted':0}

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
