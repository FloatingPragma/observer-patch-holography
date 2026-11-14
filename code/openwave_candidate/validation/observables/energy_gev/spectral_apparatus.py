"""Same-native-H complete two-energy/local-record apparatus and work ledger.

Reuses the energy-shell isometry method of mass_instrument.complete_experiment;
no CAR reservoir, filled vacuum or mass identity is transferred.
"""
from pathlib import Path
import sys
import hashlib
import numpy as np
BASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.uniform_native import UniformFiberReduction
try:from .finite_spectral_adapter import spectral_coupling
except ImportError:from finite_spectral_adapter import spectral_coupling

def run():
    model=UniformFiberReduction(Construction(3,((0,1),(1,2),(2,0)),.23),6)
    h=model.hamiltonian;e,v=np.linalg.eigh(h);d=len(e);rho=np.eye(d)/d
    # This is a separate sequential experiment: initial energy is measured.
    # Its canonical conjugate reference remains explicitly correlated to i.
    outcomes=[];max_energy_error=0.;normalization=0.;source_delta=0.;battery_delta=0.
    for meta,k in model.kraus_for(0):
        matrix=v.conj().T@k@v
        for i in range(d):
            for j in range(d):
                probability=float(abs(matrix[j,i])**2/d)
                change=float(e[j]-e[i]);battery=-change
                max_energy_error=max(max_energy_error,abs(change+battery))
                normalization+=probability;source_delta+=probability*change;battery_delta+=probability*battery
                outcomes.append({'initial_energy_index':i,'final_energy_index':j,'reference_energy_label':i,
                                 'native_record':meta,'probability':probability,'source_energy_change':change,
                                 'battery_energy_change':battery})
    p=np.zeros((d,d));p[0,0]=1
    spectral=spectral_coupling(h,rho,p)
    # Actual supported shell columns: Gram identity follows from sum K* K=I,
    # with different initial energies in orthogonal total-energy shells.
    gram=sum((k.conj().T@k for _,k in model.kraus_for(0)),np.zeros_like(h))
    err=float(np.max(abs(gram-np.eye(d))))
    return {'schema':'oph.native.spectral_apparatus.v1','dimension':d,'checks':{'complete_joint_records':abs(normalization-1)<1e-12,
            'supported_shell_isometry':err<1e-12,'energy_work_balance':max_energy_error<1e-12 and abs(source_delta+battery_delta)<1e-12},
            'complete_records':outcomes,'total_probability':normalization,'mean_source_energy_change':source_delta,
            'mean_battery_energy_change':battery_delta,'shell_gram_error':err,
            'unmeasured_initial_energy_alternative_spectral_coupling':spectral,
            'resources':{'source':'actual modular uniform-fiber triangle H_m and every localcoarse Kraus at seam0','initial_density':'I/d, stationary; retained conjugate purification vec(I)/sqrt(d)','reference':'initial spectral label i retained after initial-energy PVM; actual measurement changes its correlations','battery':'declared energy labels B0+Ei−Ej, all initial-supported shell columns; no autonomous finite reusable battery implementation claimed','apparatus':'complete initial H eigenbasis, native C/E outcome, final H eigenbasis; degenerate basis resolution supplied','preparation_work':'stationary maximally mixed source/reference preparation supplied, not claimed free','mass_instrument_reuse':'same energy-shell isometry method; original CAR reservoir law and36outcome probabilities are not transplanted'},
            'scope':'Complete sequential experiment distinct from original unmeasured-reference spectral response and from energy_gev separate marginals. No species, vacuum or GeV calibration inferred.',
            'full_criterion_met':False,'official_rows_changed':[],
            'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('finite_spectral_adapter.py'),BASE/'candidate/oph_openwave/construction.py',BASE/'candidate/oph_openwave/uniform_native.py']}}
if __name__=='__main__':
    import json
    Path(__file__).with_name('SPECTRAL_APPARATUS.json').write_text(json.dumps(run(),indent=2)+'\n')
