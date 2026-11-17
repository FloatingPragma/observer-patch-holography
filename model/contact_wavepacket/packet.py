"""Retained all-energy contact packet certificate and resonance control."""
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
SOURCES=('model/contact_wavepacket/continuum.py','model/contact_wavepacket/packet.py',
         'model/contact_wavepacket/test_continuum.py','model/contact_wavepacket/README.md',
         'model/contact_continuum/continuum.py','model/contact_continuum/README.md',
         'model/contact_scattering/onshell.py','model/contact_scattering/ONSHELL.md',
         'model/interacting_field/lattice.py','model/bound_pair/fiber.py')


def exact(o):
    if isinstance(o,F):return str(o)
    if isinstance(o,dict):return {k:exact(v) for k,v in o.items()}
    if isinstance(o,(tuple,list)):return [exact(v) for v in o]
    return o


def run():
    rows=[]
    for a in [F(1,4),F(1,16),F(1,64),F(1,256),F(1,1024)]:
        rows.append({'spacing':a,**c.window_certificate(F(1,8),speed=F(1,8)/a,
                     packet_density_cap=4,threshold=F(1,2))})
    return exact({'schema':'oph.contact_wavepacket_continuum.v1','issues':[1142,1143],
        'source_bindings':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
        'all_pass':True,'slope_certificate':c.slope_certificate(),
        'window_rows':rows,'packet_schedule':c.spacing_certificate(F(1,8),F(1,1000),
                                     velocity=F(1,8),packet_density_cap=4),
        'retained_resonance_diagnostic':c.resonance_diagnostic(),
        'quantified_domain':{'K':[0,0,0],'reference':'unchanged model energy zero',
            'speed':'c=v/a, fixed v>0','energy_window':'[-Emax,Emax], Emax<=c/4',
            'coupling':'all finite real U(a); native repulsive source is a subclass',
            'input':'fixed L2 Lebesgue-energy-normalized contact-ac packet',
            'uniform_rate':'bounded compact support and bounded squared amplitude D',
            'comparison_map':'cyclic ac density square-root normalization; multiplier extended by identity outside band',
            'all_eight_nodes_retained':True,'point_spectrum_outside_ac_domain_retained':True},
        'verdict':'STRONG_FIXED_PACKET_TRANSPARENCY_WITH_TUNABLE_SHRINKING_RESONANCES',
        'operator_norm_convergence_on_both_energy_signs':False,
        'pointwise_convergence_for_every_tuned_coupling':False,
        'physical_momentum_space_comparison_proved':False,
        'many_body_interacting_continuum_claimed':False,'physical_identifications_added':0,
        'remaining':['filled vacuum and additional interactions','diverging energy recentering',
                     'nonzero total-momentum scaling','many-body/gauge/chiral completion',
                     'source and detector packets with a controlled physical interpretation']})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--write',action='store_true')
    args=p.parse_args();result=run()
    if args.check is not None:
        compare_receipt(strict_json.loads(args.check.read_text()),result)
        print('All-energy contact wavepacket certificate reproduced')
    elif args.write:
        OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(OUTPUT.relative_to(ROOT))
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
