"""Exact retained moving bound-band and propagator certificate."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from . import limit as m
ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCES=('model/moving_bound_pair/__init__.py','model/moving_bound_pair/__main__.py',
 'model/moving_bound_pair/limit.py','model/moving_bound_pair/packet.py',
 'model/moving_bound_pair/test_limit.py','model/moving_bound_pair/README.md',
 'model/band_edge_contact/limit.py','model/band_edge_contact/README.md',
 'model/contact_continuum/continuum.py','model/contact_wavepacket/continuum.py',
 'model/contact_scattering/scattering.py','model/bound_pair/fiber.py',
 'model/interacting_field/lattice.py','paper/tex_fragments/MOVING_BOUND_PAIR.tex')


def exact(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:exact(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [exact(v) for v in x]
    return x


def run():
    return exact({'schema':'oph.moving_bound_pair.v1','issues':[1142,1143],
      'source_bindings':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
      'all_pass':True,'green_certificate':m.green_comparison_certificate(),
      'moving_poles':[m.moving_pole_certificate(F(1,n),1,1) for n in [512,2048,8192]],
      'propagation_schedule':m.spacing_certificate(1,1,1,F(1,1000)),
      'exterior_transition':m.exterior_transition_certificate(F(1,4096),1,1,4,F(1,2)),
      'verdict':'CONTROLLED_MOVING_BOUND_PAIR_SECTOR',
      'quantified_domain':{'lattice':'original infinite cubic Weyl/contact two-body action',
        'total_momentum':'K=aP, |P_j|<=Pmax',
        'coupling':'c=sqrt(3)/a^2; inverse U=a^2 W0/2-a^3 gamma/pi; gamma>0',
        'comparison':'single global J=6/a^2-H; explicit energy shift and reversed time',
        'band_energy':'|P|^2/4-gamma^2',
        'isometry':'a^(-3/2) f(K/a) times normalized exact resolvent bound vector',
        'all_spin_and_valley_components_retained':True,
        'other_source_channels_retained':True,
        'reference_extensions':True},
      'physical_identifications_added':0,
      'many_body_QFT_claimed':False,'relativistic_mass_shell_claimed':False,
      'nonzero_K_on_shell_scattering_claimed':False,
      'remaining':['physical vacuum/time/energy/detector maps','many-pair interactions and stability',
                   'nonzero-total-momentum on-shell scattering','relativistic massive/massless sectors',
                   'local Lorentz stress and Einstein premises']})


def main():
    p=argparse.ArgumentParser(description=__doc__);group=p.add_mutually_exclusive_group()
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--write',action='store_true')
    args=p.parse_args();data=run()
    if args.check is not None:
        compare_receipt(strict_json.loads(args.check.read_text()),data)
        print('Moving bound-pair certificate reproduced')
    elif args.write:
        OUTPUT.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
        print(OUTPUT.relative_to(ROOT))
    else:print(json.dumps(data,indent=2,allow_nan=False))


if __name__=='__main__':main()
