"""Independent boundary, K2-support and neutral variational review."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from dataclasses import replace
import hashlib
import json
import random
import sympy as s

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=ROOT/'model/openwave_neutral_candidate'


def run():
    saved=json.loads((PACKAGE/'results.json').read_text())
    for path,digest in saved['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms,cells,cell_index,shift_cell
    from model.coherent_defect_wave.chiral_joint.peter_weyl_links import model as native
    from model.openwave_neutral_candidate.source import translate_key
    from model.defect_mass_join.source import car,one_cell_polynomial
    cfg=Configuration();terms=kinetic_terms(cfg)
    positions=cells(2);L=(-0,0,0,0,0,0,0,0,0,-1,1,1,0,0,0)
    B=(-1,-1,1,1,-1,1,1,1,1,0,0,0,-1,-1,-1)
    assert all((a+b)%2==1 for a,b in zip(L,B))
    for weights in (L,B):
        charge=car.one_body(s.diag(*(weights[i%15] for i in range(30))))
        assert not car.commutator(charge,one_cell_polynomial())
    for spin in (0,1):
        probe=car.neutral_composite(spin)
        for generator in range(12):assert not car.commutator(car.gauge_current(generator),probe)
        assert all(len(word)%2 for word,_,_ in probe)

    # Independently implement the fermionic exterior permutation, including
    # the seam phase, and compare every mode and generic dense occupied keys.
    def independent_translate(key,axis,periodic):
        mask,links,scalars=key;destinations=[];seam=0
        for mode in range(240):
            if mask>>mode&1:
                site,local=divmod(mode,30);point=positions[site]
                target=cell_index(2,shift_cell(2,point,axis))
                destinations.append(30*target+local);seam+=point[axis]
        inversions=sum(a>b for j,a in enumerate(destinations)for b in destinations[j+1:])
        phase=(-1)**(inversions+(0 if periodic else seam))
        outlinks=[0]*24;outscalars=[0]*8
        for site,point in enumerate(positions):
            target=cell_index(2,shift_cell(2,point,axis));outscalars[target]=scalars[site]
            outlinks[3*target:3*target+3]=links[3*site:3*site+3]
        return(sum(1<<mode for mode in destinations),tuple(outlinks),tuple(outscalars)),phase
    rng=random.Random(7221);keys=[(1<<mode,(0,)*24,(0,)*8)for mode in range(240)]
    keys +=[(rng.getrandbits(240),tuple(rng.randrange(119)for _ in range(24)),tuple(rng.randrange(5)for _ in range(8)))for _ in range(32)]
    checks=0
    for key in keys:
        for periodic in (False,True):
            for axis in range(3):
                first,phase=independent_translate(key,axis,periodic)
                assert (first,phase)==translate_key(key,axis,periodic=periodic)
                second,again=independent_translate(first,axis,periodic)
                assert second==key and phase*again==(1 if periodic else (-1)**key[0].bit_count())
                checks+=1
    for key in keys[-32:]:
        for axis,other in ((0,1),(0,2),(1,2)):
            a,p=independent_translate(key,axis,True);ab,q=independent_translate(a,other,True)
            b,r=independent_translate(key,other,True);ba,t=independent_translate(b,axis,True)
            assert ab==ba and p*q==r*t

    # Actual blank-link PW coordinates land strictly outside the blank row.
    # H0 is local matter/scalar plus diagonal link number, so every edge
    # support is preserved also by its resolvent. Disjoint single-edge
    # supports force every cross-edge second-order contraction to zero.
    for adjoint in (False,True):
        blank=native.single_link_columns(adjoint)[0]
        assert blank and all(row!=0 for _,_,row,_ in blank)
    signs=[]
    for term in terms:
        sign=(-1)**positions[term.left_cell][term.axis]
        assert sign*term.coefficient==s.I/4
        assert len(term.path)==1
        assert sum(positions[term.left_cell])%2!=sum(positions[term.right_cell])%2
        signs.append(sign)
    assert signs.count(-1)==12 and signs.count(1)==12
    # Test actual first-hop support with a populated native coordinate at
    # the right endpoint; no abstract hopping matrix replaces the source.
    for term in terms:
        key=(1<<(30*term.right_cell+9),(0,)*24,(0,)*8)
        image=native.hopping_column(key,(term,));assert image
        assert all([j for j,entry in enumerate(out[1])if entry]!=[] and
                   set(j for j,entry in enumerate(out[1])if entry)==set(term.path)for out in image)
        changed=replace(term,coefficient=signs[terms.index(term)]*term.coefficient)
        variant=native.hopping_column(key,(changed,))
        assert variant=={out:signs[terms.index(term)]*value for out,value in image.items()}

    # A native odd-L defect changes one local charge label. Its eight
    # translated vectors therefore have an orthogonal regular group orbit.
    # Verify all characters on actual signed fermionic translation orbits.
    masks=[]
    base=sum(1<<(30*site+10)for site in range(8))
    for change in ('remove','add'):
        mask=base^(1<<10) if change=='remove' else base|(1<<11)
        seed=(mask,(0,)*24,(0,)*8);orbit={}
        for address in product((0,1),repeat=3):
            value=seed;phase=1
            for axis,power in enumerate(address):
                if power:value,factor=independent_translate(value,axis,True);phase*=factor
            assert value not in orbit;orbit[value]=phase
        assert len(orbit)==8
        for axis in range(3):
            translated={}
            for key,phase in orbit.items():
                out,factor=independent_translate(key,axis,True);translated[out]=phase*factor
            assert translated==orbit
        patterns=[]
        for key in orbit:
            patterns.append(tuple(sum(L[mode%15]for mode in range(30*site,30*(site+1))if key[0]>>mode&1)for site in range(8)))
        assert len(set(patterns))==8
        assert sum(s.Rational(abs(value)**2,64)for value in orbit.values())==s.Rational(1,8)

    result=saved['threshold'];gap=F(result['complete_atomic_gap_lower']);v=F(result['perturbation_norm_upper'])
    assert v==360*F(1,2**34) and gap>v
    assert F(result['variant_quartic_effective_H_remainder_upper'])==25*v**4/gap**3
    assert F(result['inherited_ground_isolation_positive_check'])>0
    scale=F(result['calibration']['GeV_per_native_energy'])
    for row in result['sectors']:
        low,high=map(F,row['periodic_variant_zero_momentum_sector_threshold_interval'])
        assert low==gap-v and high==F(row['local_ground_relative_trial_gap_interval'][1])+v
        assert F(row['local_lowest_line_probability_for_fixed_setting_interval'][0])>0
        assert list(map(F,row['conditional_GeV_energy_threshold_interval']))==[low*scale,high*scale]
        assert F(row['remaining_threshold_under_that_perturbation_lower'])==low/2
    assert not saved['boundary']['original_odd_sector_zero_momentum_exists']
    assert not result['scope']['relativistic_mass_shell_or_asymptotic_particle_proved']
    assert not result['scope']['active_or_sterile_neutrino_identified']
    return {'schema':'oph.independent_neutral_candidate_review.v1','all_pass':True,
            'independent_actual_fermion_translation_comparisons':checks,
            'actual_single_edge_PW_hop_support_and_seam_sign_checks':24,
            'K2_transfer_reason':'Complete H0 preserves every nonblank edge support; cross-edge resolvent contractions vanish and equal-edge seam signs square to one.',
            'nonzero_zero_momentum_trial_reason':'Native charge patterns distinguish the eight translated localized vectors; their signed CAR orbit gives exactly1/8 projection success.',
            'scope':'Original odd sector lacks zero source momentum. A separately declared periodic boundary has a protected positive zero-momentum threshold, not a known relativistic particle mass.',
            'reviewed_sha256':{str((PACKAGE/name).relative_to(ROOT)):hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest()for name in ('source.py','threshold.py','PROOF.md','results.json')}}


if __name__=='__main__':
    value=run();target=Path(__file__).resolve().parents[1]/'neutral_review.json'
    target.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    print('PASS: native neutral charge, fermionic boundary, K2 support and momentum threshold')
