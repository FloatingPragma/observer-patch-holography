"""Independent supported-shell construction and same-source continuity audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.linalg import expm
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'candidate'))
sys.path.insert(0,str(BASE/'validation/observables/energy_gev'))
from oph_openwave.uniform_native import UniformFiberReduction
from oph_openwave.construction import Construction
from oph_openwave import native_currents
import spectral_apparatus

def run():
    m=UniformFiberReduction(Construction(3,((0,1),(1,2),(2,0)),.23),6)
    h=m.hamiltonian;e,v=np.linalg.eigh(h);ops=m.kraus_for(0);d=len(e)
    # Explicit finite supported-shell output map, independent of helper ledger.
    # Equal battery energies are identified; arbitrary splitting would hide overlaps.
    labels=[];rows={};triples=[]
    for letter,(_,k) in enumerate(ops):
        a=v.conj().T@k@v
        for i in range(d):
            for j in range(d):
                delta=float(e[i]-e[j]);key=(letter,j,round(delta,10))
                if key not in rows:rows[key]=len(rows);labels.append((j,delta))
                triples.append((rows[key],i,a[j,i]))
    w=np.zeros((len(rows),d),complex)
    for row,col,value in triples:w[row,col]+=value
    gram=float(np.max(abs(w.conj().T@w-np.eye(d))))
    energy=float(max(abs((e[j]+delta)*w[row,col]-e[col]*w[row,col]) for row,(j,delta) in enumerate(labels) for col in range(d)))
    receipt=spectral_apparatus.run()
    law={}
    for letter,(meta,k) in enumerate(ops):
        a=v.conj().T@k@v
        for i in range(d):
            for j in range(d):law[(letter,i,j)]=abs(a[j,i])**2/d
    assert np.isclose(sum(law.values()),1)
    records=receipt['complete_records'];assert len(records)==len(law)
    err=max(abs(row['probability']-law[(q//(d*d),row['initial_energy_index'],row['final_energy_index'])]) for q,row in enumerate(records))
    p,j=native_currents.regional_current(h,[0]);psi=np.array([1,2j,0,1-1j,0,0]);psi=psi/np.linalg.norm(psi)
    dt=1e-5
    charge=lambda t:float(np.vdot(expm(-1j*t*h)@psi,p@(expm(-1j*t*h)@psi)).real)
    continuity=abs((charge(dt)-charge(-dt))/(2*dt)-np.vdot(psi,j@psi).real)
    assert max(gram,energy,err)<1e-10 and continuity<1e-8
    sources=[Path(__file__),BASE/'candidate/oph_openwave/native_currents.py',BASE/'validation/observables/energy_gev/spectral_apparatus.py',BASE/'validation/observables/energy_gev/finite_spectral_adapter.py']
    return {'outcome':'PASS_WITH_SCOPE_NOTES','errors':{'actual_supported_shell_gram':gram,'actual_total_energy_intertwiner':energy,'complete_sequential_outcome_probabilities':float(err),'finite_difference_same_H_continuity':float(continuity)},'scope':['supported initial battery shell only, not autonomous reusable finite battery','initial and final energy PVM supplied, degenerate basis resolution supplied','same-H whole-position charge is not yet physical electric current','arbitrary idle reference preserved by explicit isometry/Kraus maps'],'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'full_physical_criterion_met':False}
if __name__=='__main__':
    r=run();Path(__file__).with_name('CODEX_REVIEW_NATIVE_APPARATUS.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['errors']))
