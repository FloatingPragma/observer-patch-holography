"""Reproducible conditional pixel/mass closure and an excluded scale design."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from .closure import ROOT, DESIGNS, certificate, Point
from . import census, exclusion, memory, replay
from .sources import DEPENDENCIES

PACKAGE=Path(__file__).resolve().parent
OUTPUT=PACKAGE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths=set(DEPENDENCIES)
    paths.update(str(p.relative_to(ROOT)) for p in PACKAGE.iterdir() if p.suffix in ('.py','.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def report():
    pins=source_bindings()
    closures=[]
    for design in DESIGNS:
        row=certificate(design)
        row['memory_admissibility']=memory.certify(row)
        row['complete_record_controls']=memory.point_controls(row)
        # Same-coordinate exclusion of the top term is a retained control.
        point=Point(design,40,18,include_top=False)
        a=point.b.num(row['display']['g'])
        row['same_coordinate_top_omission_control']={
            'omitted_charged_fermion':'top',
            'inverse_response_without_top':point.b.ctx.nstr(point.inner(a)['alpha_inv'],28),
            'is_a_self_consistent_root':False}
        closures.append(row)
    result={'schema':'oph.pixel_mass_closure.v1','all_pass':True,'issues':[1143],
            'census':census.coefficients(),'excluded_scale_design':exclusion.certificate(),
            'closures':closures,'generic_memory_identities':memory.old.exact_identities(),
            'units':{'E_star':'declared energy unit, no numerical laboratory value supplied',
                     'v_mz_MU':'ratios to E_star','mass_coefficients':'m_i/v',
                     'P_alpha_alpha_U':'dimensionless'},
            'response_contract':{
                'masses':'same MCPR charged leptons and Stage5 quark coefficients as the chosen field design',
                'fermions':'three charged leptons and six quarks, including the finite top contribution',
                'kernel':'2*N*Q^2/pi times integral_0^1 x(1-x) log(1+(mZ/m)^2*x(1-x)) dx',
                'quark_factor':'declared factor 1-3*alpha3/pi',
                'omissions':['charged-vector and charged-scalar terms','supersymmetric thresholds and matching',
                             'confining hadronic spectral density','renormalized physical mass and clock maps'],
                'sm_color_beta7':'same minimal census, chosen Higgs-scale transmutation coefficient |b3|=7',
                'mssm_beta4_proxy':'MSSM running ansatz with beta4, extra MSSM fields not realized by native carrier'},
            'scope':{'conditional_roots_and_inverse_memory_certified':True,
                     'same_mass_law_reenters_declared_feedback':True,
                     'runtime_measured_masses_or_alpha_read':False,
                     'MCPR_architecture_target_informed':True,
                     'complete_physical_SM_Thomson_response':False,
                     'native_current_spectral_response_derived':False,
                     'physical_particle_map_complete':False,'integrated_gate_accepted':False,
                     'source_law_uniquely_selected':False},
            'canonical_physical_premises_discharged':[],
            'source_bindings':pins}
    result['independent_simultaneous_root_and_quadrature_replay']=replay.verify_report(result)
    if pins != source_bindings():
        raise ValueError('source changed during report')
    return result


def validate(value):
    if value.get('source_bindings') != source_bindings():
        raise ValueError('source bindings are incomplete or stale')
    if canonical(value) != canonical(report()):
        raise ValueError('receipt differs from complete recomputation')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',const=str(OUTPUT))
    group.add_argument('--check',nargs='?',const=str(OUTPUT))
    args=parser.parse_args()
    if args.write:
        Path(args.write).write_text(canonical(report()))
        print('pixel/mass closure written')
    else:
        validate(strict_json.load(Path(args.check)))
        print('pixel/mass closure: interval proof and independent numerical replay pass')


if __name__=='__main__':main()
