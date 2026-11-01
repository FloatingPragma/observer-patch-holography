"""Retained exact certificate for the interacting two-body band-edge limit."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from . import limit as c
ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCES=('model/band_edge_contact/limit.py','model/band_edge_contact/packet.py',
         'model/band_edge_contact/test_limit.py','model/band_edge_contact/README.md',
         'model/contact_wavepacket/continuum.py','model/contact_continuum/continuum.py',
         'model/contact_scattering/onshell.py','model/interacting_field/lattice.py',
         'model/bound_pair/fiber.py','paper/tex_fragments/BAND_EDGE_CONTACT.tex')


def exact(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:exact(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [exact(v) for v in x]
    return x


def run():
    return exact({'schema':'oph.band_edge_contact.v1','issues':[1142,1143],
        'source_bindings':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
        'all_pass':True,'threshold_certificate':c.green_certificate(),
        'scattering_schedule':c.spacing_certificate(1,1,F(1,1000)),
        'scattering_rows':[c.scattering_certificate(F(1,n),1,1) for n in [256,1024,4096]],
        'pole_rows':[c.pole_certificate(F(1,n),1) for n in [256,1024,4096]],
        'dispersion_example':c.dispersion_certificate(F(1,32),(F(1,3),F(-2,5),F(3,7))),
        'quantified_domain':{'total_momentum':[0,0,0],'spatial_dimensions':3,
          'free_speed':'sqrt(3)/a^2','energy_edge':'B=6/a^2',
          'inverse_contact':'a^2 W0/2-a^3 gamma/pi, W0=<1/sum cos(k_j)^2>',
          'gamma':'declared fixed positive renormalized parameter',
          'comparison':'J=B-H; explicit energy shift and reversed Hamiltonian time',
          'scattering_representation':'cyclic contact ac energy representation, compact shifted energy intervals',
          'all_edge_valleys_retained':True,'orthogonal_free_channels_retained':True},
        'verdict':'CONTROLLED_NONIDENTITY_BAND_EDGE_SCATTERING_AND_BOUND_POLE',
        'physical_particle_identifications_added':0,'physical_clock_claimed':False,
        'norm_resolvent_limit_claimed':False,'many_body_QFT_limit_claimed':False,
        'remaining':['physical vacuum and time/energy interpretation','nonzero total momentum',
                     'many-body stability and field limit','massive/massless relativistic sectors',
                     'physical detector comparison and units']})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--write',action='store_true')
    args=p.parse_args();result=run()
    if args.check is not None:
        compare_receipt(strict_json.loads(args.check.read_text()),result)
        print('Band-edge contact limit certificate reproduced')
    elif args.write:
        OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(OUTPUT.relative_to(ROOT))
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
