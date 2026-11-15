"""Native spectral/interference readouts with an explicit optional GeV anchor."""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np
from scipy.linalg import expm

BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.geometry import flat_ring
from oph_openwave.uniform_native import UniformFiberReduction


def gev_scale(reference_native, reference_gev):
    if not all(isinstance(x,(int,float,np.integer,np.floating)) and not isinstance(x,(bool,np.bool_)) and np.isfinite(x) and x>0 for x in (reference_native,reference_gev)):
        raise ValueError('positive finite native and GeV reference energies required')
    scale=float(reference_gev)/float(reference_native)
    if not math.isfinite(scale) or scale<=0:
        raise ValueError('GeV scale is outside finite positive floating-point range')
    return scale


def run(mass=.23,anchor_gev=.000511):
    if not isinstance(mass,(int,float,np.integer,np.floating)) or isinstance(mass,(bool,np.bool_)) or not np.isfinite(mass) or mass<=0:
        raise ValueError('strictly positive finite mass parameter required')
    gev_scale(mass,anchor_gev)  # Validate supplied unit inputs before numerical work.
    model=UniformFiberReduction(flat_ring(1,1,mass),6)
    h=model.hamiltonian;energy,vectors=np.linalg.eigh(h)
    reference_gap=float(energy[np.argmin(abs(energy-mass))]-energy[np.argmin(abs(energy+mass))])
    reference_half_gap=reference_gap/2
    scale=gev_scale(reference_half_gap,anchor_gev)
    lap=np.linalg.eigvalsh(model.construction.incidence@model.construction.incidence.T)
    distinct=[]
    for value in lap:
        if not distinct or abs(value-distinct[-1])>1e-8:distinct.append(max(0,float(value)))
    pairs=[];max_error=0.;trace_error=0.
    for lam in distinct:
        target=math.sqrt(mass*mass+lam)
        ip=int(np.argmin(abs(energy-target)));im=int(np.argmin(abs(energy+target)))
        plus,minus=vectors[:,ip],vectors[:,im]
        gap=float(energy[ip]-energy[im])
        # Explicit two-level interferometric apparatus with a retained complement.
        initial=(plus+minus)/math.sqrt(2);opposite=(plus-minus)/math.sqrt(2)
        pplus=np.outer(initial,initial);pminus=np.outer(opposite,opposite)
        complement=np.eye(len(energy))-pplus-pminus
        samples=[]
        for t in np.linspace(0,4*np.pi/gap,9):
            state=expm(-1j*h*t)@initial
            probabilities=[float(abs(np.vdot(initial,state))**2),float(abs(np.vdot(opposite,state))**2),float(np.vdot(state,complement@state).real)]
            predicted=(1+math.cos(gap*t))/2
            max_error=max(max_error,abs(probabilities[0]-predicted))
            trace_error=max(trace_error,abs(sum(probabilities)-1))
            samples.append({'native_time':float(t),'plus_minus_complement_probabilities':probabilities})
        rho=np.outer(initial,initial)
        measured,records=model.step(rho,.37,.29)
        native_energy_probs=np.einsum('ij,ji->i',vectors.T@measured,vectors).real
        pairs.append({'laplacian_eigenvalue':lam,'positive_energy_native':float(energy[ip]),
                      'negative_energy_native':float(energy[im]),'gap_native':gap,
                      'gap_gev':gap*scale,'half_gap_gev':gap*scale/2,
                      'half_gap_to_supplied_m_ratio':gap/(2*mass),
                      'positive_position_ipr':float(np.sum(abs(plus)**4)),
                      'positive_state_amplitudes':plus.tolist(),'negative_state_amplitudes':minus.tolist(),
                      'interference_samples':samples,'native_step_energy_probabilities':native_energy_probs.tolist(),
                      'native_step_complete_record_probabilities':records})
    shifted=np.linalg.eigvalsh(h+7*np.eye(len(h)))
    gaps=np.array([p['gap_native'] for p in pairs])
    shifted_gaps=np.array([shifted[np.argmin(abs(energy-p['positive_energy_native']))]-shifted[np.argmin(abs(energy-p['negative_energy_native']))] for p in pairs])
    paths=[Path(__file__).resolve(),BASE/'candidate/oph_openwave/geometry.py',BASE/'candidate/oph_openwave/construction.py',BASE/'candidate/oph_openwave/uniform_native.py',BASE/'candidate/oph_openwave/modular_native.py']
    return {'schema':'oph.native.energy_gev.v1','scope':'Finite canonical native spectra and complete declared interferometric readouts. GeV values are conditional calibration outputs, not identified particle masses.',
            'parameters':{'mass_native':mass,'ring_shape':[1,1],'modulus':6,'preparation':'uniform Gauss fibers plus explicit eigenvector superpositions','spectral_measurement':'declared complete H eigenbasis; all 87 real eigenvectors archived as columns; resolving a degenerate basis is an added apparatus choice','interference_measurement':'two rank-one superposition projectors and their complement; supplied cross-branch phase reference and nonlocal spectral POVM, not implemented by native C/E alone','native_readout_scope':'energy PVM and event-record probabilities are separate marginals; no joint record-by-energy law or sequential measurement is claimed','native_comparison_step':{'time':.37,'probability':.29,'law':'equal coarse C/E'},'vacuum_state_identified':False},
            'calibration':{'reference_native':reference_half_gap,'supplied_rest_parameter_native':mass,'reference_gev':anchor_gev,'reference_interpretation':'half of the measured +m/-m spectral gap is mapped to the row-2 511 keV benchmark; this shift-invariant calibration does not identify a physical rest mass','gev_per_native_unit':scale,'electron_mass_predicted':False,'physical_particle_identifications':[]},
            'energy_eigenvectors_columns':vectors.tolist(),'all_energy_levels_native':energy.tolist(),'all_energy_levels_gev':(scale*energy).tolist(),'pairs':pairs,
            'checks':{'all_branches_retained':len(energy)==87,'complete_energy_basis_orthonormal':float(np.max(abs(vectors.T@vectors-np.eye(len(energy)))))<1e-11,'interference_agrees_with_gap':max_error<1e-11,'complete_apparatus_normalized':trace_error<1e-11,'scalar_shift_preserves_gaps':float(np.max(abs(gaps-shifted_gaps)))<1e-11,'gev_anchor_scaling':gev_scale(reference_half_gap,2*anchor_gev)==2*scale},
            'errors':{'interference_max':max_error,'apparatus_normalization_max':trace_error,'scalar_shift_gap_max':float(np.max(abs(gaps-shifted_gaps)))},
            'not_established':['physical particle species','vacuum-to-particle mass gap','infinite-volume rest pole or Lorentz mass','electron/muon/proton mass prediction','self-starting native clock'],
            'full_physical_criterion_met':False,'official_rows_promoted':0,
            'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}

if __name__=='__main__':
    result=run()
    if not all(result['checks'].values()):raise SystemExit('energy measurement check failed')
    Path(__file__).with_name('ENERGY_GEV.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':result['checks'],'errors':result['errors'],'gap_gev':[p['gap_gev'] for p in result['pairs']]},indent=2))
