"""Independent consumption of the frozen common normalization."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'units'))
import units


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    manifest_path=HERE.parent/'units/MANIFEST.json'
    assert sha(manifest_path)=='fd3abf81be6147577b31e36a00f339664ad8f287e44052ccbc3adafc07e3876f'
    manifest=json.loads(manifest_path.read_text())
    pins=dict(manifest['source_pins'])
    for name,digest in manifest['files'].items():pins[str((manifest_path.parent/name).relative_to(ROOT))]=digest
    pins[str(manifest_path.relative_to(ROOT))]=sha(manifest_path)
    pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    assert all(sha(ROOT/path)==digest for path,digest in pins.items())
    receipt=json.loads((manifest_path.parent/'UNITS.json').read_text())
    W=s.Matrix([[s.Rational(x) for x in row] for row in receipt['geometry']['unit_flux_basis']])
    A0,B0,D,terms,labels=units.field.bloch_data()
    crossings=s.zeros(3,36)
    for axis in range(3):
        crossings[axis,30+2*axis]=crossings[axis,31+2*axis]=1
    assert B0*W==s.zeros(12,3) and A0.T*W==s.zeros(26,3)
    assert crossings*W==s.eye(3) and W.T*W==s.Rational(2,3)*s.eye(3)
    theta=s.Rational(3,2)*W
    assert theta.T*W==s.eye(3)
    canonical_pairs=[]
    for alpha,beta,N in ((s.Rational(1,8),s.Rational(1,4),1),
                         (s.Rational(1,16),s.Rational(3,8),8),
                         (s.Rational(1,32),s.Rational(3,4),8)):
        electric=N*alpha*s.Rational(2,3);magnetic=N*beta
        independent={'electric':electric,'magnetic':magnetic,
                     'speed_squared':electric*magnetic,
                     'charge_coupling_squared':s.sqrt(electric/magnetic)}
        got=units.field_normalization(F(alpha),F(beta),N,homogeneous=True)
        assert all(got[k]==v for k,v in independent.items())
        canonical_pairs.append({k:str(v) for k,v in independent.items()})
    assert canonical_pairs[1]['speed_squared']==canonical_pairs[2]['speed_squared']=='1'
    assert canonical_pairs[1]['charge_coupling_squared']=='1/3'
    assert canonical_pairs[2]['charge_coupling_squared']=='1/6'
    unit=units.Units(F(7,3),F(5,2),F(11,7),F(13,11))
    assert unit.clock_phase(3,F(1,8))==F(3,64)
    assert unit.elapsed(3,F(1,8))==3*unit.time/8
    assert unit.eigenphase_energy(F(-1,8),F(1,8))==unit.rest_energy(1)
    u2=units.Units(unit.length,unit.time,unit.action,2*unit.charge)
    assert u2.rest_energy(1)==unit.rest_energy(1)
    assert all(sha(ROOT/path)==digest for path,digest in pins.items())
    output={'schema':'oph.common.units.opposing_review.v1',
            'checks':{'source_bytes_match':True,'unit_flux_columns_on_actual_cell_symbol':True,
                      'canonical_symplectic_pair':True,'shared_clock_and_energy':True,
                      'speed_does_not_select_charge_coupling':True},
            'canonical_settings':canonical_pairs,'joint_clock_phase':'3/64',
            'UV_stability_boundary':'The scalar speed-matched setting has conservative curvature upper product79/16>4, so that upper bound alone does not certify its full ultraviolet stability. It does not prove instability. The Dirac8substep setting has237/64<4.',
            'preparation_boundary':'Definite F0, finite cycle-character states and a canonical Maxwell vacuum are different preparations. Coefficient/geometry normalization is conditional on a controlled common quantum field state and limit.',
            'physical_scope':'No choice of external length/time/action/charge units changes the dimensionless gauge coupling. The conditional Coulomb prefactor is not an executed two-charge force or a renormalized electromagnetic alpha.',
            'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'UNITS_REVIEW.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'checks':output['checks'],'settings':canonical_pairs}))


if __name__=='__main__':run()
