"""Same-C10 global-reference charged excitation, complete energy/charge readout."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
C10=ROOT/'openwave-tuning/next/C10';sys.path.insert(0,str(C10))
from joint_source import JointNativeSource
from spectral_readout import GroupedEnergyInstrument

def run():
 paths=[Path(__file__),C10/'joint_source.py',C10/'charged_reservoir.py',C10/'modular_field_adapter.py',C10/'spectral_readout.py',C10/'joint_source_states.npz',HERE.parent/'gauge_loop/GAUGE_LOOP.json',HERE.parent/'gauge_loop/gauge_loop.py']
 pins={str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 model=JointNativeSource();archive=np.load(C10/'joint_source_states.npz');assert np.array_equal(model.hamiltonian,archive['hamiltonian'])
 instrument=GroupedEnergyInstrument(model);reference=archive['global_allN_reference'];ground=model.ground
 reports=[];candidates=[];total=0.;angle=np.linspace(0,2*np.pi,33);q=model.numbers-3
 for meta,k in model.source.reservoir_instrument():
  label=f"global_allN_reservoir_{meta['reservoir_initial']}_{meta['reservoir_final']}";state=archive[label]
  assert np.max(abs(k@reference@k.T-state))<1e-12
  rows=instrument.measure(state,label);p=float(np.trace(state));total+=p
  for row in rows:
   row['relative_charge']=row['source_number']-3;row['excitation_energy']=row['energy']-ground
   row['gauge_orbit_degree']=row['relative_charge']
   if row['relative_charge']==1 and row['probability']>1e-10 and row['excitation_energy']>1e-10:candidates.append(row)
  reports.append({**meta,'preparation_joint_probability':p,'complete_joint_number_energy_charge_outcomes':rows})
 selected=max(candidates,key=lambda r:r['probability']);branch=instrument.reconstruct(selected)
 # Integral gauge degree derives from exact boundary label, not a Z6 residue.
 degree=selected['relative_charge'];ids=np.flatnonzero(np.diag(branch).real>1e-14)
 assert all(q[i]==degree for i in ids)
 loop=np.exp(1j*angle*degree);measured=float(np.angle(loop[1:]/loop[:-1]).sum()/(2*np.pi))
 continued,ledger=instrument.continue_branch(branch);prob=selected['probability']
 native_work=float(np.trace(model.hamiltonian@continued).real-selected['energy']*prob)
 assert abs(total-1)<1e-12 and abs(sum(x['joint_probability'] for x in ledger)-prob)<1e-12
 assert abs(measured-degree)<1e-12 and np.linalg.norm(model.hamiltonian@np.diag(q)-np.diag(q)@model.hamiltonian)<1e-12
 np.savez_compressed(HERE/'joint_charge_energy_states.npz',selected_poststate=branch,continued_poststate=continued,**instrument.arrays)
 result={'schema':'oph.same_source.charge_energy_bridge.v1','checks':{'bit_exact_frozen_interacting_H':True,'actual_global_reference_reservoir_poststates_recovered':True,'all_joint_charge_energy_outcomes_retained':True,'positive_charged_excitation_selected_from_complete_law':True,'exact_integer_boundary_charge_not_mod6_lift':True,'complete_native_continuation_on_same_H':True},'source':'C10 384-state all-number CAR/shared-Z6-record periodic-field source; unchanged parameters and clock','global_reference_energy':ground,'complete_preparation_and_readout':reports,'selected_excitation':selected,'selected_complete_native_continuation':ledger,'continued_unnormalized_energy_work':native_work,'gauge_orbit':{'degree':measured,'definition':'exp(i theta (N-3)); Q=-B exact supplied boundary register','readout_service':'Existing gauge_loop source-compiled character/pointer receipt; this bridge reuses its proven Q character, does not claim a newly executed phase grid.','topology_scope':'Gauge-orbit degree only; literal-native spatial dual-cell winding has not been intertwined with this quantum source.'},'selection':'Largest-probability positive-energy Q=+1 outcome after complete charged preparation; selected poststate is conditional, all rejected/zero outcomes remain above.','scope':'Joint number/energy Luders PVM is declared additional apparatus, not synthesized by native C/E. Excitation energy is relative to actual interacting all-number global reference. CAR statistics, boundary-energy-zero choice, Z6 field regulator and addressed clock remain supplied. No spatial-defect/electromagnetic/species identity or physical mass prediction.','state_archive':'joint_charge_energy_states.npz','state_archive_sha256':hashlib.sha256((HERE/'joint_charge_energy_states.npz').read_bytes()).hexdigest(),'source_sha256':pins,'official_rows_changed':[]}
 assert pins=={str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 return result
if __name__=='__main__':
 r=run();(HERE/'JOINT_CHARGE_ENERGY.json').write_text(json.dumps(r,indent=2)+'\n');print(r['checks']);print(r['selected_excitation'])
