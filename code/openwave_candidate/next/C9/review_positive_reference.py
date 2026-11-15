"""Independent scalar-cost and exact exterior-instrument audit; no producer import."""
from pathlib import Path
import sys,json,hashlib,itertools
import numpy as np
import sympy as sp
BASE=Path(__file__).resolve().parents[2];sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.uniform_native import UniformFiberReduction
from oph_openwave.geometry import flat_ring


def run():
    source=BASE/'validation/bridges/positive_reference_attachment.py';receipt=source.with_name('POSITIVE_REFERENCE_ATTACHMENT.json');r=json.loads(receipt.read_text())
    model=UniformFiberReduction(flat_ring(1,1,.23),6)
    eigen=np.linalg.eigvalsh(model.hamiltonian);costs=np.abs(eigen)
    saved=np.array(r['attachment']['reference_covariance_real'])+1j*np.array(r['attachment']['reference_covariance_imag'])
    projector_error=float(np.linalg.norm(saved@saved-saved));commutator=float(np.linalg.norm(model.hamiltonian@saved-saved@model.hamiltonian))
    # Exhaust all 64 Fock occupations for triangle, directly summing occupied levels.
    tiny=UniformFiberReduction(Construction(3,((0,1),(1,2),(2,0)),.23),6)
    values,vectors=np.linalg.eigh(tiny.hamiltonian);sea=sum(x for x in values if x<0)
    actual=[sum(values[j] for j in range(6) if mask>>j&1)-sea for mask in range(64)]
    identity=[sum(abs(values[j]) for j in range(6) if bool(mask>>j&1)!=(values[j]<0)) for mask in range(64)]
    n=sum(values<0);indices=list(itertools.combinations(range(6),n));gram=sp.zeros(len(indices))
    for _,k in tiny.kraus_for(0):
        exact=sp.Matrix([[sp.Rational(int(round(2*x)),2) for x in row] for row in k])
        lift=sp.Matrix([[exact.extract(rows,cols).det() for cols in indices] for rows in indices])
        gram+=lift.T*lift
    saved_diagonal=list(map(sp.Rational,r['exact_rejecting_control']['sum_exterior_KdagK_diagonal']))
    exact_diagonal_match=list(gram.diagonal())==saved_diagonal
    checks={'finite_reference_projector':projector_error<1e-12,'reference_commutes_with_same_H':commutator<1e-12,
      'charged_gap_matches_same_H':abs(min(costs)-r['attachment']['charged_excitation_gap'])<1e-12,
      'full64_fock_cost_identity':max(abs(np.array(actual)-identity))<1e-12,
      'exact_naive_lift_diagonal_reproduced':exact_diagonal_match,
      'both_raw_and_coin_corrected_lifts_not_TP':gram!=sp.eye(len(indices)) and 16*gram!=sp.eye(len(indices))}
    checks={k:bool(value) for k,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    return {'schema':'oph.c9.positive_reference_independent_review.v1','checks':checks,
      'reference_projector_error':projector_error,'same_H_commutator_error':commutator,
      'charged_gap':float(min(costs)),'neutral_fixed_number_gap':float(min(eigen[eigen>0])-max(eigen[eigen<0])),
      'exact_naive_gram_diagonal':list(map(str,gram.diagonal())),
      'findings':['Positive reference and excitation identity are genuine finite CAR results on unchangedH, not merely absent services.',
         'Complete native CAR extension remains missing: neither raw wedgeK nor coin-corrected wedgeK forms a complete instrument.',
         'One-body covariance mapping is explicitly diagnostic; it need not preserve fermionic0<=C<=I and is not a many-body channel.',
         'No physical vacuum/species or calibrated physical electron mass follows. Select existing compatible CAR native apparatus or specify a new complete architecture before native stability claims.'],
      'full_criterion_met':False,'reviewer':'Codex independent worker; not Fable opposite-team',
      'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,receipt,Path(__file__),BASE/'candidate/oph_openwave/uniform_native.py',BASE/'candidate/oph_openwave/geometry.py',BASE/'candidate/oph_openwave/construction.py')}}

if __name__=='__main__':
    r=run();Path(__file__).with_name('POSITIVE_REFERENCE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
