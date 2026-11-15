"""Exact implicit uniform-fiber compact transit mode; finite modular controls.

Large record fibers are represented by their full linear-constraint definition,
not enumerated or approximated. No physical particle or acceptance is claimed.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

PACKAGE=Path(__file__).resolve().parent
from .geometry import flat_ring
from .construction import Construction
from .modular_native import LocalCoarseModular
from .uniform_native import UniformFiberReduction


def require(ok,message):
    if not ok:raise ValueError(message)


def canonical_control(lx,ly):
    c=flat_ring(lx,ly,.23)
    integer_cycle=np.zeros(len(c.seams),dtype=int)
    for edge,sign in [((0,11),1),((5,11),-1),((0,5),-1)]:
        integer_cycle[c.seams.index(edge)]=sign
    divergence=c.incidence.astype(int)@integer_cycle
    require(np.array_equal(divergence,np.zeros(c.ports,dtype=int)),'cycle not closed')
    require(int(integer_cycle@integer_cycle)==3,'wrong support normalization')
    psi=np.r_[np.zeros(c.ports),integer_cycle/np.sqrt(3)]
    h=c.hamiltonian
    residual=float(np.max(abs(h@psi+.23*psi)))
    evolved=expm(-.37j*h)@psi
    phase_error=float(np.linalg.norm(evolved-np.exp(.37j*.23)*psi))
    require(residual<1e-12 and phase_error<1e-12,'bare eigenstate identity failed')
    broken=integer_cycle.copy();broken[c.seams.index((0,11))]=0
    broken_psi=np.r_[np.zeros(c.ports),broken/np.linalg.norm(broken)]
    broken_evolved=expm(-.37j*h)@broken_psi
    broken_leak=float(np.sum(abs(broken_evolved[:c.ports])**2))
    require(broken_leak>.01,'broken-cycle control failed to leak')
    beta=len(c.seams)-c.ports+1
    fiber_size=6**beta
    reduction=UniformFiberReduction(c,6)
    require(reduction.fiber_cardinality==fiber_size,'reduction fiber normalization mismatch')
    local_records=[];max_kraus_error=0.
    for e,operations in enumerate(reduction.maps):
        for record,src,dst in operations:
            image=np.zeros_like(psi);image[dst]=.5*psi[src]
            expected=.5*psi if record['record_value']==0 else np.zeros_like(psi)
            max_kraus_error=max(max_kraus_error,float(np.max(abs(image-expected))))
            local_records.append({**record,
                'probability_conditional_on_commit':str(F(1,4*len(c.seams)) if record['record_value']==0 else F(0)),
                'post_state':'same normalized uniform-fiber compact transit state' if record['record_value']==0 else 'zero-probability outcome'})
    require(max_kraus_error<1e-12,'selected reduced local Kraus action differs from proof')
    require(sum(F(r['probability_conditional_on_commit']) for r in local_records)==1,
            'incomplete local records')
    return {'ring_shape':[lx,ly],'ports':c.ports,'seams':len(c.seams),
            'mass_exact':'23/100','energy_exact':'-23/100',
            'integer_cycle_divergence':divergence.tolist(),
            'spatial_support':[{'seam':list(c.seams[e]),'amplitude':f'{int(integer_cycle[e])}/sqrt(3)',
                                'probability_exact':'1/3'} for e in np.flatnonzero(integer_cycle)],
            'record_fiber':{'modulus':6,'cycle_rank':beta,'fiber_cardinality_exact':str(fiber_size),
                            'full_connected_sector_dimension_exact':str((c.ports+len(c.seams))*fiber_size),
                            'prepared_nonzero_joint_basis_count_exact':str(3*fiber_size),
                            'constraint':'For each supported transit position z: subdivided-incidence F = delta_initial_port0 - delta_z modulo6',
                            'amplitude_on_each_joint_basis':'signed 1/sqrt(3*fiber_cardinality); all other amplitudes zero',
                            'enumerated':False},
            'complete_position_readout':{'port_probabilities':['0']*c.ports,
                                         'transit_probabilities':['1/3' if x else '0' for x in integer_cycle]},
            'complete_native_outcomes_conditional_on_commit':local_records,
            'every_reduced_kraus_selected_state_max_error':max_kraus_error,
            'no_commit_outcome_at_p_029':{'probability':'71/100','post_state':'same state'},
            'conditional_commit_probability':'29/100',
            'numerical_phase_error':phase_error,'eigen_residual':residual,
            'broken_cycle_port_leakage_t037':broken_leak,
            'proof_domain':'Implicit complete modular Gauss fibers; exact algebra, not an executable preparation of the enormous state'}


def explicit_triangle_control():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    model=LocalCoarseModular(c,6)
    v=np.zeros((model.n,6))
    fibers=[]
    for pos in range(6):
        fiber=[i for i,state in enumerate(model.states) if state[0]==pos]
        fibers.append(fiber);v[fiber,pos]=1/np.sqrt(len(fiber))
    bare=np.r_[np.zeros(3),np.ones(3)/np.sqrt(3)]
    lifted=v@bare
    # Cumulative orientation of the triangle is (0,1),(1,2),(2,0).
    intertwining=float(np.max(abs(model.hamiltonian@v-v@c.hamiltonian)))
    require(intertwining<1e-12,'uniform-fiber generator control failed')
    residual=float(np.linalg.norm(model.hamiltonian@lifted+.23*lifted))
    rho=np.outer(lifted,lifted)
    out,ledger=model.step(rho,.37,.29)
    stability=float(np.linalg.norm(out-rho,ord='nuc')/2)
    require(residual<1e-12 and stability<1e-12,'actual modular compact-state stability failed')
    definite=np.zeros(model.n)
    for pos in range(3,6):
        chosen=min(fibers[pos],key=lambda i:model.states[i][1])
        definite[chosen]=1/np.sqrt(3)
    definite_evolved=expm(-.37j*model.hamiltonian)@definite
    leakage=float(sum(abs(definite_evolved[i])**2 for i,s in enumerate(model.states) if s[0]<3))
    require(leakage>.01,'definite record lift did not reject bare cancellation transfer')
    fully,_=model.step(rho,0.,1.,fully_monitored_control=True)
    fully_evolved=expm(-.37j*model.hamiltonian)@fully@expm(.37j*model.hamiltonian)
    monitor_leak=float(sum(fully_evolved[i,i].real for i,s in enumerate(model.states) if s[0]<3))
    require(monitor_leak>.01,'global monitor control failed to destroy dark cancellation')
    return {'dimension':model.n,'fiber_sizes':[len(x) for x in fibers],
            'intertwiner_max_error':intertwining,'lifted_eigenstate_residual':residual,
            'local_coarse_step_trace_distance':stability,
            'complete_joint_state':[{'position':int(s[0]),'record_residues':list(s[1]),
                                     'amplitude':float(lifted[i])} for i,s in enumerate(model.states)],
            'complete_local_step_outcomes':ledger,
            'definite_per_position_record_lift_port_leak_t037':leakage,
            'full_monitor_then_coherent_port_leak_t037':monitor_leak}


def run():
    controls=[canonical_control(1,1),canonical_control(2,2)]
    triangle=explicit_triangle_control()
    files=['construction.py','geometry.py','modular_native.py','uniform_native.py','compact_native.py']
    return {'schema':'oph.openwave.c6.uniform_fiber.compact_transit.v1',
            'scope':'Spatially compact transit mode with nonlocal uniform modular record preparation; no electron identity, physical clock, propagation, magnetic moment or official criterion',
            'canonical_flat_controls':controls,
            'explicit_actual_modular_triangle':triangle,
            'checks':{'closed_canonical_cycles':all(not any(c['integer_cycle_divergence']) for c in controls),
                      'all_selected_kraus_actions':all(c['every_reduced_kraus_selected_state_max_error']<1e-12 for c in controls),
                      'complete_exact_native_probabilities':all(sum(F(r['probability_conditional_on_commit']) for r in c['complete_native_outcomes_conditional_on_commit'])==1 for c in controls),
                      'actual_modular_eigenstate':triangle['lifted_eigenstate_residual']<1e-12,
                      'same_instrument_persistence':triangle['local_coarse_step_trace_distance']<1e-12},
            'mutations':{'broken_cycle_rejected':all(c['broken_cycle_port_leakage_t037']>.01 for c in controls),
                         'definite_record_transfer_rejected':triangle['definite_per_position_record_lift_port_leak_t037']>.01,
                         'global_monitor_stability_rejected':triangle['full_monitor_then_coherent_port_leak_t037']>.01},
            'supplied':{'modulus':6,'mass':.23,'preparation':'uniform complete modular Gauss fiber on a signed three-seam transit cycle',
                        'transit_policy':'local C/E records treat all transit positions as elsewhere, without a transit flag',
                        'canonical_graphs':[[1,1],[2,2]],'development_triangle':True},
            'physical_rows':{'4':'Exact persistence for the selected uniform-fiber transit state in the same local-coarse instrument; no identified physical particle.',
                             '3':'State density and all local readouts are stationary; global phase is not an internal clock.',
                             '5':'No field B or intrinsic magnetic moment follows from a flat transit mode.',
                             '6':'No physical angular-momentum dictionary follows from a compact graph cycle.'},
            'full_criterion_met':False,'official_rows_changed':[],
            'source_sha256':{f'oph_openwave/{p}':hashlib.sha256((PACKAGE/p).read_bytes()).hexdigest() for p in files}}


if __name__=='__main__':print(json.dumps(run(),indent=2))
