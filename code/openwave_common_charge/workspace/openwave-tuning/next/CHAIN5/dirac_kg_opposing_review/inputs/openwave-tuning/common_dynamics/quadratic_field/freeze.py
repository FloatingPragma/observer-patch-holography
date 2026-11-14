"""Freeze this bounded proposal without selecting it or promoting any row."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RECEIPTS=('SOURCE','PHOTONS','FOCK','NATIVE_PHOTON')
PROVIDERS={
 'openwave-tuning/common_dynamics/vector_field/MANIFEST.json':'aa4a18d89c1a8893ee86eb5e9151bc75e8a0b3a527640cee66aa25382ca55e76',
 'openwave-tuning/common_dynamics/quantum_field/MANIFEST.json':'38cd7609e1be37ab5bb3751fb0cbfc6b35d6c393d52fbe9e2e5e40817c12ec83',
 'openwave-tuning/common_dynamics/spin/MANIFEST.json':'65075e0bc14f883c5cdb877171650b08ee57c8d1bdffce281e68937bb07eb401',
 'openwave-tuning/common_dynamics/units/MANIFEST.json':'fd3abf81be6147577b31e36a00f339664ad8f287e44052ccbc3adafc07e3876f',
 'openwave-tuning/auditor/harmonic_stability/MANIFEST.json':'9b5723b1b9dfc48e3b38e3480ca308c5f1f33af61409142c9242065fb1b4da2e'}


def run():
    pins={};counts={}
    def add(relative,digest):
        assert len(digest)==64 and hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest,relative
        assert relative not in pins or pins[relative]==digest,relative
        pins[relative]=digest
    for name in RECEIPTS:
        receipt=json.loads((HERE/(name+'.json')).read_text())
        assert receipt['checks'] and all(receipt['checks'].values()),name
        counts[name]=len(receipt['checks'])
        for relative,digest in receipt['source_sha256'].items():add(relative,digest)
    for relative,digest in PROVIDERS.items():
        add(relative,digest)
        provider=json.loads((ROOT/relative).read_text())
        for key in ('files','files_sha256','source_pins','source_sha256'):
            for source,expected in provider.get(key,{}).items():
                if isinstance(expected,str) and len(expected)==64:
                    target=ROOT/source
                    if not target.is_file():target=(ROOT/relative).parent/source
                    add(str(target.relative_to(ROOT)),expected)
    core='openwave-tuning/auditor/native_walk_bridge/COMMON_SOURCE_CUSTODY.json'
    add(core,hashlib.sha256((ROOT/core).read_bytes()).hexdigest())
    for source,digest in json.loads((ROOT/core).read_text())['files_sha256'].items():add(source,digest)
    files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
           for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
    manifest={'schema':'oph.common.quadratic_field.freeze.v1','worker':'DYN16',
      'candidate_id':'native-local12-dirac-noncompact-quadratic-v4-proposed',
      'files':files,'source_pins':pins,'receipt_checks':counts,
      'checks':{'four_receipts_pass':True,'all_source_and_provider_pins_match':True,
        'actual190coexact_modes_and3harmonic_packets_retained':True,
        'fixed_coupling_photon_theorem_and_native_spectral_readout':True,
        'vacuum_phase_and_cap_boundary_explicit':True,'no_source_selection_or_row_promotion':True},
      'interfaces':{'source_geometry':'source.geometry','Gaussian_reference':'source.reference',
        'free_mean_evolution':'source.propagate_mean','actual_routed_current':'source.transverse_current',
        'ordered_Gaussian_overlap':'source.gaussian_translation_overlap',
        'same_field_interference':'source.two_route_records','photon_limit':'photons.frequency_bound',
        'integer_geometry_certificate':'fock.spectral_certificate',
        'native_free_phase':'native_photon.photon_phase','native_free_photon_clock':'native_photon.word'},
      'source_amendment':'Noncompact local quadratic field on the actual face complex, retaining native code-service, Gauss and integer matter-current ancestry; not equivalent to the frozen compact cosine field.',
      'scope':'Full continuous Gaussian free field and exact charge-sector energy; controlled actual-face fixed-coupling two-polarization excitation frequencies; bounded cap-one native spectral-code/clock/readout. Harmonic native discretization, local field encoder, interacting Weyl-current accuracy, moving massive particles, calibrated forces and official validation remain open.',
      'selected_by_this_worker':False,'independent_review_pending':True,
      'reserved_acceptance_runs':False,'official_rows_changed':[]}
    path=HERE/'MANIFEST.json';path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'manifest_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'checks_total':sum(counts.values()),'receipt_checks':counts,'files':len(files),'source_pins':len(pins)}),flush=True)


if __name__=='__main__':run()
