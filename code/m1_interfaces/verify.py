"""Fail-closed, producer-free replay of all committed interface evidence."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

from . import check
import pin_scope

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def pins():
    paths=[p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')]
    paths += ['Lean/Geometry/M1Interfaces.lean','Lean/Geometry/M1InterfacesAxiomAudit.lean',
              'Lean/Geometry/FlatDiamondNormalization.lean',
              'Lean/Geometry/OrderingFractionFourDimensional.lean',
              'paper/tex_fragments/M1_INTERFACES.tex','.github/workflows/m1-interfaces.yml']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(pin_scope.bound(paths))}


def parent_claims():
    claims=check.strict_load(ROOT/'claims/claim_registry.yaml')['claims']
    selected={
        'OPH-BH-CROSSING-READ-AREA-LAW':'variational_cut_extension_without_physical_entropy_identification',
        'OPH-GR-D4B-SOURCE-CAUSAL-CONTINUUM':'named_kinematics_retained_without_native_source_identification'}
    matches=[row for row in claims if row['claim_id'] in selected]
    check.require(len(matches)==len(selected),'parent claims must occur exactly once')
    found={row['claim_id']:row for row in matches}
    check.require(set(found)==set(selected),'parent claim catalog')
    return {name:{'original_claim_sha256':check.digest(found[name]),'disposition':meaning}
            for name,meaning in selected.items()}


@lru_cache(maxsize=None)
def graph(config, family, topology):
    masks,cuts,reads,action=check.graph_cuts(config,family,topology)
    q=config[0]
    vectors=check.stencil(config,family)
    rows={}
    for name,mask in masks.items():
        rows[name]=dict(volume_sites=sum(mask),cut_pairs=cuts[name],
                        ordered_crossing_reads=2*cuts[name],
                        pi_times_normalized_entropy=str(F(4*cuts[name],q**4)),
                        mask_sha256=hashlib.sha256(bytes(mask)).hexdigest())
    changed=sum(x!=y for x,y in zip(masks['localized'],masks['connected']))
    check.require(changed>0 and changed<=q,'axial bridge size')
    check.require(not check.connected(config,masks['localized']),'localized control must be disconnected')
    check.require(check.connected(config,masks['connected']),'axial bridge fails to connect')
    check.require(abs(cuts['localized']-cuts['connected'])<=len(vectors)*changed,'cut change bound')
    if topology=='periodic':
        m=config[1]
        if family=='T': expected=2*q**3*(m-1)//m
        else:
            fine=check.ball_columns(config[3],q)
            extra=sum(1<<b for b in range(m.bit_length()-1) if 1<<b>config[3])
            expected=2*q**3*(fine['flux']+extra)//m
        check.require(cuts['residue']==expected,'exact residue cut identity')
        check.require(rows['residue']['volume_sites']==q**3//2,'residue volume')
    return dict(config=list(config),family=family,topology=topology,
                stencil_sha256=check.digest(vectors),degree=len(vectors),
                ordered_reads=reads,masks=rows,uniform_action=action,
                bridge=dict(changed_sites=changed,cut_change=cuts['connected']-cuts['localized'],
                            degree_bound=len(vectors)*changed,nearest_neighbor_connected=True))


def verify(packet):
    check.require(type(packet) is dict and set(packet)=={'schema','sources','parent_claims','scales','moments','graphs'},
                  'exact interface receipt schema')
    check.require(packet['schema']=='oph-m1-interfaces-v2','interface schema version')
    check.require(check.canonical(pin_scope.bound_pins(packet['sources']))==check.canonical(pins()),'source custody mismatch')
    check.require(check.canonical(packet['parent_claims'])==check.canonical(parent_claims()),'parent claim scope mismatch')
    check.require(type(packet['scales']) is dict and set(packet['scales'])=={str(t) for t in range(1,13)},
                  'scale catalog')
    for t in range(1,13):
        check.require(check.canonical(packet['scales'][str(t)])==check.canonical(check.scale(t)),'false scale certificate')
    check.require(type(packet['moments']) is dict and set(packet['moments'])=={'raw','critical','repaired'},
                  'moment regime catalog')
    for kind in ('raw','critical','repaired'):
        check.require(check.canonical(packet['moments'][kind])==check.canonical(check.moment_case(kind)),
                      'false complete moment or alias certificate')
    keys={f'{c[0]}:{f}:{b}' for c in check.CONFIGS for f in ('T','U') for b in ('periodic','clipped')}
    check.require(type(packet['graphs']) is dict and set(packet['graphs'])==keys,'finite graph catalog')
    total_reads=total_sites=total_cases=0
    for config in check.CONFIGS:
        for family in ('T','U'):
            for topology in ('periodic','clipped'):
                key=f'{config[0]}:{family}:{topology}'
                expected=graph(config,family,topology)
                check.require(check.canonical(packet['graphs'][key])==check.canonical(expected),'false finite cut execution')
                total_reads+=expected['ordered_reads']
                total_sites+=config[0]**3
                total_cases+=len(expected['masks'])
    return dict(graphs=len(keys),cut_cases=total_cases,site_instances=total_sites,
                ordered_read_incidences=total_reads,scale_levels=12,moment_regimes=3)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,default=HERE/'receipt.json')
    args=parser.parse_args()
    print(check.canonical(verify(check.strict_load(args.receipt))).decode())


if __name__=='__main__': main()
