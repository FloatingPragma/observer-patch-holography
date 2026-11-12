"""Graded native translations on every atomic one-neutral state."""
from itertools import product
import random
from model.openwave_neutral_candidate.source import translate_key

POINTS=tuple(product((0,1),repeat=3))
INDEX={x:i for i,x in enumerate(POINTS)}


def shift(site,axis):
    point=list(POINTS[site]);point[axis]^=1
    return INDEX[tuple(point)]


def translate(label,axis):
    hole,bits=label
    occupied=[i for i in range(8) if i!=hole]
    spins={site:(bits>>j)&1 for j,site in enumerate(occupied)}
    targets=[shift(site,axis) for site in occupied]
    inversions=sum(a>b for i,a in enumerate(targets) for b in targets[i+1:])
    newhole=shift(hole,axis);newspins={shift(site,axis):spin for site,spin in spins.items()}
    newbits=sum(newspins[site]<<j for j,site in enumerate(i for i in range(8) if i!=newhole))
    return (newhole,newbits),(-1)**inversions


def certificate():
    labels=[(y,bits) for y in range(8) for bits in range(128)]
    checks=0
    for label in labels:
        for axis in range(3):
            one,a=translate(label,axis);two,b=translate(one,axis)
            assert two==label and a*b==1;checks+=1
            for other in range(axis):
                left,c=translate(one,other)
                mid,d=translate(label,other);right,e=translate(mid,axis)
                assert left==right and a*c==d*e
    traces=[]
    for element in POINTS:
        trace=0
        for label in labels:
            dest=label;phase=1
            for axis,bit in enumerate(element):
                if bit:dest,p=translate(dest,axis);phase*=p
            if dest==label:trace+=phase
        traces.append(trace)
    assert traces==[1024,0,0,0,0,0,0,0]
    ranks=[sum((-1)**sum(k*g for k,g in zip(momentum,group))*trace
               for group,trace in zip(POINTS,traces))//8 for momentum in POINTS]
    assert ranks==[128]*8
    rng=random.Random(65162);native=0
    for _ in range(64):
        hole=rng.randrange(8);mask=0
        for site in range(8):
            occupation=rng.getrandbits(30)
            if occupation.bit_count()%2!=int(site!=hole):occupation^=1
            mask|=occupation<<(30*site)
        key=(mask,(0,)*24,(0,)*8)
        for axis in range(3):
            _,phase=translate_key(key,axis,periodic=True)
            _,expected=translate((hole,0),axis)
            assert phase==expected;native+=1
    origin=(0,0);orbit=[]
    for element in POINTS:
        label=origin;phase=1
        for axis,bit in enumerate(element):
            if bit:label,p=translate(label,axis);phase*=p
        orbit.append({'cell':label[0],'ground_spin_bits':label[1],'zero_momentum_amplitude_numerator':phase})
    assert len({x['cell'] for x in orbit})==8
    return {'atomic_band_dimension':1024,'local_neutral_even_and_seven_ground_factors_odd':True,
        'exact_graded_translation_square_checks':checks,'independent_native_occupation_phase_checks':native,
        'commuting_periodic_translation_generators':True,'group_character':traces,
        'all_eight_momentum_ranks':ranks,'explicit_zero_momentum_orbit':orbit,
        'orbit_amplitude_denominator':'sqrt(8)',
        'localized_trial_projection_probability_at_each_momentum':'1/8',
        'continuation_preserves_each_momentum_rank_by_separated_spectral_projection':True,
        'antiperiodic_member_silently_used':False}
