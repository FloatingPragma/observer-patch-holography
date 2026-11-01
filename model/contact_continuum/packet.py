"""Source-bound exact continuum-exclusion certificate and diagnostics."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from . import continuum as c

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCES=('model/contact_continuum/continuum.py',
         'model/contact_continuum/packet.py',
         'model/contact_continuum/test_continuum.py',
         'model/contact_continuum/README.md',
         'model/contact_scattering/onshell.py',
         'model/contact_scattering/ONSHELL.md',
         'model/interacting_field/lattice.py',
         'model/bound_pair/fiber.py')


def exact(value):
    if isinstance(value,F): return str(value)
    if isinstance(value,dict):return {k:exact(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [exact(v) for v in value]
    return value


def run():
    velocity,energy=F(1,8),F(1,64)
    rows=[]
    for spacing in (F(1),F(1,2),F(1,4),F(1,16),F(1,64)):
        rows.append({'spacing':spacing,'speed':velocity/spacing,
                     'fixed_U_one':c.scattering_bound(energy,speed=velocity/spacing,interaction=1),
                     'diverging_U_a_to_minus_ten':c.scattering_bound(
                         energy,speed=velocity/spacing,interaction=spacing**-10)})
    diagnostic=[]
    for x in (.125,.25,.5):
        value=c.broadened_green(x,.5,order=96)
        diagnostic.append({'energy':x,'eta':.5,'speed':1.,
                           'green_real':value.real,'green_imaginary':value.imag})
    return exact({'schema':'oph.repulsive_contact_continuum.v1',
        'issues':[1142,1143], 'source_bindings':{
            name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        'sign_certificate':c.principal_value_certificate(),
        'scaling_rows':rows,
        'uniform_window':c.spacing_certificate(F(1,16),F(1,1000)),
        'kinetic_example':{'momentum':[F(1,2),F(-2,3),F(3,4)],'spacing':F(1,32),
            'all_eight_nodes_error':c.kinetic_error([F(1,2),F(-2,3),F(3,4)],F(1,32))},
        'broadened_diagnostic_not_boundary_value':diagnostic,
        'verdict':'POSITIVE_ENERGY_CONTACT_SCATTERING_IS_UNIFORMLY_TRANSPARENT_AS_SPACING_VANISHES',
        'quantified_domain':{'total_momentum':[0,0,0], 'energy':'0<E<=Emax, same model zero',
            'kinetic_scaling':'c=v/a, fixed v>0', 'contact_scaling':'every finite U(a)>=0',
            'energy_ratio':'E/c<=1/2', 'all_eight_Weyl_nodes_retained':True,
            'S_domain':'absolutely continuous two-body scattering subspace'},
        'above_band_pair':'if present, E_pair>2 sqrt(3) v/a, so it escapes every bounded positive-energy window',
        'exclusions_not_claimed':['negative-energy resonances','filled or shifted physical vacuum',
            'nonzero total-momentum scaling','added interactions or counterterms','many-body continuum',
            'whole OPH realization class'],
        'physical_particle_identifications_added':0, 'physical_gate_accepted':False})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=run()
    if args.check is not None:
        compare_receipt(strict_json.loads(args.check.read_text()),result)
        print('Repulsive-contact continuum certificate reproduced')
    elif args.write:
        OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(OUTPUT.relative_to(ROOT))
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
