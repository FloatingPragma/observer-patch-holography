"""Publish scoped source claims without changing physical-premise statuses."""
from pathlib import Path
from hashlib import sha256
import argparse
import importlib
import json
from tools.native_bridge_results import NAMES,project

ROOT=Path(__file__).resolve().parents[1]
SPECS=(
 ('Narrow Lambda spectral weight','unchanged strong source; enlarged full-operator trial frame','Spectral projector residual inequality and exact translated-packet Gram forms.'),
 ('Repeated charged-source clock records','declared apparatus and finite reference/reset resources','Finite quantum-instrument composition and the energy-covariance obstruction.'),
 ('Autonomous fractional-charge string records','declared gauge-covariant hopping extension','Exact CAR hole identity and a reducing Jacobi chain inside the full Gauss carrier.'),
 ('Nuclear field attraction and complete-input work','declared neutral oscillator extension','Complete-square lower bound, coherent-state upper bound and full-singlet quadratic Casimir trace.'),
 ('Full matter with an integer U(1) cycle','declared spatial and field extension','Exact Gauss-kernel unitary, bounded perturbation and full-tail Temple inequality.'),
 ('Local gauge memory on three refinements','declared full-Haar sources and common native times','Local interaction-picture propagation and complete temporal-instrument anticommutator identity.'),
 ('Finite quantum stress and tensor acquisition','declared native coordinate source and metric readout','Ten-mean factorization, complete quantum instruments, sharp null inverse and staggered transport.'),
 ('Geometry covariance and physical unit normalization','declared matter/photon/geometry source and resources','Complete paired covariance records, finite-sample bounds and source-dependent unit conversion.'),
 ('Certified finite-wavelength cosmological response','declared classical Einstein/scalar/fluid model and late preparation','Exact comoving reduction, interval Taylor transport and a complete Volterra remainder.'),
 ('OpenWave M10 status and arithmetic certification','authenticated external status and exact source replay','Complete seven-column matrix/tally check and three original exact implementations.'),
)

def generate():
    data={};rows=[]
    for name,(title,relation,shortcut)in zip(NAMES,SPECS,strict=True):
        prefix='model/openwave_'+name;path=ROOT/prefix/'results.json'
        value=json.loads(path.read_text());owner=importlib.import_module('model.openwave_'+name+'.__main__')
        if value.get('source_sha256')!=owner.pins()or not value.get('verdict','').startswith('PASS_'):
            raise ValueError('source result lacks current original-owner custody: '+name)
        data[name]=value
        rows.append({'id':'source_record_'+name,'label':title,'owner_medium':'protocol_record',
            'status':'verified_in_stated_source_domain'if name!='m10_review_gate'else'authenticated_external_status_and_exact_replay',
            'construction_relation':relation,'finite_theorem_shortcut':shortcut,
            'proof':prefix+'/PROOF.md','receipt':prefix+'/results.json','receipt_sha256':sha256(path.read_bytes()).hexdigest(),
            'native_twelve_port_implementation_proved':False,'independent_Nature_match_added':False,
            'physical_particle_mass_identification_added':False,'official_criterion_awarded_by_this_work':False})
    project(data)
    canonical=('tracking/premise_register.json','claims/physical_identification_registry.json','claims/selection_ledger.json')
    return {'schema':'oph.native_process_bridge_claims.v1',
        'mission':'One specified self-referential executable process reproduces compatible observer records for all physical observables.',
        'architecture_uniqueness_required':False,'changes_existing_physical_statuses':False,
        'canonical_source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest()for p in canonical},
        'claims':rows,'physical_premises_closed':[],
        'independent_physical_masses_GeV_added':[],
        'official_status':data['m10_review_gate']['official_status'],
        'native_process_gate':data['m10_review_gate']['native_process_acceptance_gate'],
        'next_acceptance':'Prove an apparatus/state/clock/readout implementation on the actual twelve-port process for one criterion, retaining complete outcomes and a quantitative source-to-record error. Physical identification and held-out comparison follow on that same implementation.',
        'scope':'Finite and operator-theoretic source results retain their stated hypotheses. Installed fields, calibrated masses, arithmetic certification and compatible labels do not establish the full source-to-Universe mission.'}

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args()
    encoded=json.dumps(generate(),sort_keys=True,indent=2,ensure_ascii=False)+'\n';path=ROOT/'claims/native_process_bridge.json'
    if a.write:path.write_text(encoded)
    else:assert path.read_text()==encoded,'native process claim ledger differs'
    print('PASS_TEN_SOURCE_RECORD_CLAIMS_WITHOUT_PHYSICAL_PROMOTION')
if __name__=='__main__':main()
