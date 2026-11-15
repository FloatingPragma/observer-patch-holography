"""Explicit native subtraction of the known cosine constant phase."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
import quantum_state as qs

HERE=Path(__file__).resolve().parent


def word(beta,steps,bit=0):
    if type(steps) is not int or steps<0 or type(bit) is not int or bit<0:
        raise ValueError('nonnegative count and allocated data wire required')
    return [qs.nf.compiler.pulse('I',-float(F(beta)*steps),[bit])]


def run():
    paths=[Path(__file__),Path(qs.nf.compiler.__file__),Path(qs.nf.cert.__file__)]
    pins={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    rows=[]
    for eps2,N in ((F(1,3),8),(F(1,12),8),(F(1,192),64)):
        beta=qs.parameters(eps2,N)['beta'];tape=word(beta,N)
        bound=qs.nf.root_certificate.exact_template_bound(tape,[('I',(0,),-beta*N)])
        actual=qs.nf.compiler.execute(tape,1)
        residual=float(np.linalg.norm(actual-np.eye(2)*np.exp(-1j*float(beta*N)),2))
        assert residual<float(bound)+2e-12
        rows.append({'epsilon_squared':str(eps2),'substeps':N,'source_scalar_phase':str(beta*N),
                     'native_counterphase_word':tape,'native_error_exact':str(bound),
                     'numeric_full_column_residual':residual})
    assert pins=={str(p.relative_to(qs.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.local12.explicit_field_reference_phase.v1',
       'checks':{'literal_global_phase_compiled':True,'full_native_columns_checked':True,
                 'operator_error_certified':True,'source_bytes_unchanged':True},
       'cases':rows,'scope':'Optional explicit reference intervention. Original cosine-source records and controlled phase retain exp(i sum beta). Appending this native word implements the stated reference phase; its error and pulse/service duration are additional. It is not silently dropped from any preceding receipt.',
       'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'REFERENCE_PHASE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'cases':len(rows)}),flush=True)


if __name__=='__main__':run()
