"""Numerical coordinates only; the consumer rechecks full native residuals."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh
from .cube import recoupled_matrices,matrices_on_states,trial_state_set
from model.openwave_coupled_haar_spectrum.certificate import source_digest
import sympy as sy


def main():
    records={};base_states,E,faces=recoupled_matrices(32)
    trial_states=trial_state_set();Et,Wt=matrices_on_states(trial_states)
    for g in (4,8,16):
        alpha=sy.Rational(g,2);beta=sy.Rational(2,g)
        H=alpha*E+beta*(6*sy.eye(E.rows)-sum(faces,sy.zeros(E.rows)))
        Ht=alpha*Et+beta*(6*sy.eye(Et.rows)-sum(Wt,sy.zeros(Et.rows)))
        values,frame=eigh(np.asarray(H.evalf(18),float))
        energies,trial=eigh(np.asarray(Ht.evalf(18),float),subset_by_index=(0,0))
        D=2**40
        records[str(g)]={'g_squared_exact':str(g),'complete_base_quarter_energy':32,
            'exact_complete_base_matrix_sha256':source_digest(H),
            'ordered_complete_base_spins':[list(x)for x in base_states],
            'complete_finite_frame':{'dyadic_denominator':D,
                'dyadic_eigenvalues':[str(round(float(x)*D))for x in values],
                'dyadic_eigenframe':[[str(round(float(x)*D))for x in row]for row in frame]},
            'trial_quarter_energy':56,'ordered_trial_spins':[list(x)for x in trial_states],
            'ground_denominator':2**48,
            'ground_integer_coordinates':[str(round(float(x)*2**48))for x in trial[:,0]]}
        print(g,len(base_states),len(trial_states),energies[0],flush=True)
    Path(__file__).with_name('numerical_witnesses.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':main()
