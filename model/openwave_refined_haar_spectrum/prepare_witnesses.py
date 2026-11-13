"""Propose fixed integer frames; the independent producer rechecks every bound."""
from pathlib import Path
import json
import numpy as np
from .certificate import data,D
from model.openwave_coupled_haar_spectrum.certificate import source_digest


def main():
    result={}
    for which in('coarse','refined'):
        H=data(which)['H'];matrix=np.array(H,float);e,v=np.linalg.eigh(matrix);d=2**40
        result[which]={'dimension':H.rows,'exact_matrix_sha256':source_digest(H),'ground_denominator':D,
            'ground_integer_coordinates':[str(round(float(x)*D))for x in v[:,0]],
            'complete_finite_frame':{'dyadic_denominator':d,'dyadic_eigenframe':[[str(round(float(x)*d))for x in row]for row in v],
                'dyadic_eigenvalues':[str(round(float(x)*d))for x in e]}}
    Path(__file__).with_name('numerical_witnesses.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
if __name__=='__main__':main()
