#!/usr/bin/env python3
"""Small field/headline checks and optional full CCF-2 common-source join."""
from __future__ import annotations
import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import tempfile
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
import numpy as np
from scipy.sparse.linalg import expm_multiply
import sympy as sy
from . import core
from . import spatial

HERE = Path(__file__).resolve().parent


def field_check():
    graph = spatial.triangle()
    runs = []
    for magnetic in (1, 0):
        d = spatial.field_operators(graph, cap=3, electric=1, magnetic=magnetic)
        initial = np.zeros(d['dimension'], dtype=complex)
        zero = (d['dimension']-1)//2
        initial[zero] = 1
        state = expm_multiply(-.5j*d['H'], initial)
        probabilities = np.abs(state)**2
        assert abs(float(probabilities.sum())-1) < 1e-10
        energy0 = np.vdot(initial, d['H']@initial).real
        energy = np.vdot(state, d['H']@state).real
        assert abs(energy-energy0) < 1e-10
        for G in d['gauss']:
            assert np.linalg.norm(G@state) < 1e-10
        if magnetic:
            assert 1-probabilities[zero] > 1e-5
        else:
            assert abs(1-probabilities[zero]) < 1e-12
        post_energy = float(np.dot(probabilities, d['H'].diagonal().real))
        runs.append({'magnetic_coefficient': magnetic, 'all_outcome_probabilities': probabilities.tolist(),
                     'nonzero_loop_probability': float(1-probabilities[zero]),
                     'initial_energy': float(energy0), 'energy_residual': float(abs(energy-energy0)),
                     'complete_meter_work': post_energy-float(energy),
                     'all_gauss_complements_retained': True})
    shell = spatial.radial_shell()
    B, C = shell.incidence()
    assert not np.any(B@C.T)
    assert shell.nvertices == 8 and len(shell.edges) == 16 and len(shell.faces) == 14
    # A same-letter corruption must fail incidence, not be projected away.
    bad = C.copy()
    bad[0, 0] += 1
    assert np.any(B@bad.T)
    reversed_graph = replace(graph, edges=((1, 0), (0, 2), (1, 2)))
    original_B, original_C = graph.incidence()
    reversed_B, reversed_C = reversed_graph.incidence()
    assert np.array_equal(reversed_B[:, 0], -original_B[:, 0])
    assert np.array_equal(reversed_C[:, 0], -original_C[:, 0])
    assert not np.any(reversed_B @ reversed_C.T)
    return {'all_passed': True, 'full_flux_states': 343, 'runs': runs,
            'radial_shell_incidence': {'vertices': 8, 'edges': 16, 'faces': 14,
                                      'boundary_of_boundary_zero': True},
            'reversed_directed_edge_incidence_exact': True,
            'scope': 'Compact field-loop evolution, not a Maxwell continuum or emitted photon identification.'}


def actual_local_capacity_columns(source):
    from .local_capacity import geometry_shard_column, phi_labels
    from .joint_geometry import auxiliary_columns
    from .local_geometry import LocalGeometryConfig, _exact
    from . import source as master
    graph = source['graph'];n=graph.nvertices;local_index=1
    config=LocalGeometryConfig();labels=phi_labels(local_index,n,config)
    acols,mcols=auxiliary_columns();results=[]
    for label,mask in [('native_empty',0),('native_neutral',1<<core.index(0,0,10,0))]:
        native_key=(mask,(0,)*len(graph.edges))
        key=((0,)*n,(local_index,native_key))
        full=master.apply_total_column(source,key)
        assert full and all(c!=0 for c in full.values())
        expected_diagonal=sy.Integer(0)
        for v in range(n):
            lapse=1+_exact(config.lapse_coupling)*_exact(config.phi_values[labels[v]])
            shard=master.apply_source_column(source,mask,native_key[1],vertex=v)
            kv=geometry_shard_column(local_index,v,graph,config).get(local_index,0)+lapse*shard.get(native_key,0)
            expected_diagonal+=dict(acols[0]).get(0,0)+dict(mcols[0]).get(0,0)*kv
        assert sy.expand(full.get(key,0)-expected_diagonal)==0
        charges,_=core.mode_metadata(n);B=source['field']['B']
        for _,(_,(outmask,outflux)) in full:
            local_charge=[sum(charges[i] for i in range(96*v,96*(v+1)) if outmask>>i&1) for v in range(n)]
            assert local_charge==[sum(int(B[v,e])*outflux[e] for e in range(len(outflux))) for v in range(n)]
        for edge in range(len(graph.edges)):
            phase=master.apply_total_phase_column(source,key,edge,matter_only=False)
            expected=core.clean({out:sy.I*(out[1][1][1][edge]-native_key[1][edge])*c for out,c in full.items()})
            assert phase==expected
            for other in range(len(graph.edges)):
                contact=master.apply_total_phase_column(source,key,edge,other,matter_only=False)
                expected=core.clean({out:-(out[1][1][1][edge]-native_key[1][edge])*(out[1][1][1][other]-native_key[1][other])*c for out,c in full.items()})
                assert contact==expected
        # B times the matter link derivative equals i times the native charge change.
        currents=[master.apply_total_phase_column(source,key,e) for e in range(len(graph.edges))]
        for v in range(n):
            lhs=core.add(*(core.scale(current,int(B[v,e])) for e,current in enumerate(currents)))
            initial_q=sum(charges[i] for i in range(96*v,96*(v+1)) if mask>>i&1)
            rhs=core.clean({out:sy.I*(sum(charges[i] for i in range(96*v,96*(v+1)) if out[1][1][0]>>i&1)-initial_q)*c for out,c in full.items()})
            assert lhs==rhs
        results.append({'preparation':label,'complete_nonzero_outputs':len(full),
            'diagonal_exact':str(full.get(key,0)),'independent_diagonal_matches':True,
            'all_outputs_gauss_consistent':True,'local_phase_current_and_contact_exact':True,
            'exact_total_charge_continuity':True,'no_native_occupation_projection':True})
    assert results[0]['diagonal_exact']!=results[1]['diagonal_exact']
    return {'all_passed':True,'native_modes':96*n,'matter_digest':core.polynomial_digest(source['matter']),
            'capacity_auxiliary_dimension_per_vertex':len(acols),'capacity_auxiliary_count':n,
            'local_response_dimension':len(config.phi_values)**n,'columns':results,
            'all_source_components_retained':list(source['components'])+['electric','magnetic'],
            'scope':'Literal complete local-capacity native columns; no arbitrary full-state propagation.'}


def run(full=False):
    started = time.perf_counter()
    from . import dirac
    from . import geometry
    from . import local_capacity_checks
    result = {'field': field_check(), 'dirac': dirac.verify(),
              'per_vertex_capacity': geometry.witness(), 'local_capacity': local_capacity_checks.run()}
    if full:
        from . import local_geometry
        from . import source
        from . import interactions, audit_core, gauging_compiler
        result['interactions'] = interactions.checks()
        result['inherited_charge_and_preparation'] = {
            'charge': audit_core.charge_audit(),
            'topology_and_repair': audit_core.topology_repair_audit(),
            'quantum': audit_core.quantum_audit(),
            'full_fock': audit_core.fock_audit(),
            'compiler': gauging_compiler.tests(),
        }
        result['local_geometry'] = local_geometry.witness()
        result['source_edge'] = source.audit(spatial.edge())
        result['source_triangle'] = source.audit(spatial.triangle())
        result['actual_local_capacity_join'] = actual_local_capacity_columns(source.assemble(spatial.triangle()))
    for name, value in result.items():
        if isinstance(value, dict) and 'all_passed' in value:
            assert value['all_passed'], 'Failed component: '+name
    result.update({'all_passed': True, 'full_source_included': full,
                   'elapsed_seconds': time.perf_counter()-started})
    return result

