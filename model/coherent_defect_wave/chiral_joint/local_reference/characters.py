"""Complete exact Sp4 occupation characters and the excluded-sector bound."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb

from .algebra import quark_clifford


@lru_cache(maxsize=1)
def inventory():
    fundamental=quark_clifford()[-1]
    character={(0,0,0):1}
    for _ in range(6):
        for u,v in fundamental:
            for (n,a,b),multiplicity in list(character.items()):
                key=n+1,a+u,b+v
                character[key]=character.get(key,0)+multiplicity
    denominator={(0,0):1}
    for u,v in ((2,-2),(0,2),(2,0),(2,2)):
        for (a,b),coefficient in list(denominator.items()):
            key=a-u,b-v;denominator[key]=denominator.get(key,0)-coefficient
        denominator={key:value for key,value in denominator.items() if value}
    alternating=defaultdict(int)
    for (n,a,b),multiplicity in character.items():
        for (u,v),coefficient in denominator.items():
            alternating[n,a+u,b+v]+=multiplicity*coefficient
    alternating={key:value for key,value in alternating.items() if value}
    modules={key:value for key,value in alternating.items() if key[1]>=key[2]>=0}
    assert all(value>0 for value in modules.values())
    reconstructed=defaultdict(int)
    for (n,a,b),multiplicity in modules.items():
        for swap in (0,1):
            coordinates=(a+3,b+1) if not swap else (b+1,a+3)
            for signs in product((-1,1),repeat=2):
                weight=(coordinates[0]*signs[0]-3,coordinates[1]*signs[1]-1)
                reconstructed[n,*weight]+=multiplicity*((-1)**swap)*signs[0]*signs[1]
    reconstructed={key:value for key,value in reconstructed.items() if value}
    assert reconstructed==alternating
    rows=[]
    for n in range(25):
        dimension=0;entries=[]
        for (occupation,a,b),multiplicity in sorted(modules.items()):
            if occupation!=n:continue
            size=Fraction((a+3)*(b+1)*(a-b+2)*(a+b+4),24)
            casimir=Fraction(a*(a+6)+b*(b+2),4)
            assert size.denominator==1 and a>=b>=0 and (a-b)%2==0
            assert (a+b)//2<=6
            dimension+=multiplicity*int(size)
            entries.append({'twice_highest_weight':[a,b],'dimension':int(size),
                            'multiplicity':multiplicity,'Casimir':str(casimir)})
        assert dimension==comb(24,n)
        rows.append({'transformed_quark_number':n,'raw_B3':n-12,
                     'full_dimension':dimension,'modules':entries})
    assert sum(row['full_dimension'] for row in rows)==2**24
    return rows


def certificate():
    rows=inventory();modules=[(row['transformed_quark_number'],entry) for row in rows for entry in row['modules']]
    maximum=max(Fraction(entry['Casimir']) for _,entry in modules)
    largest=[(n,entry) for n,entry in modules if Fraction(entry['Casimir'])==maximum]
    assert maximum==54 and len(largest)==1
    number,entry=largest[0]
    assert number==12 and entry['twice_highest_weight']==[12,0]
    assert entry['dimension']==140 and entry['multiplicity']==1
    admissible=[Fraction(entry['Casimir']) for n,entry in modules if n%3==0 and Fraction(entry['Casimir'])!=54]
    assert max(admissible)==42
    # (sqrt(84)+sqrt(6))^2=90+12 sqrt(14)<135, checked without a square root.
    assert Fraction(135-90,12)**2>14
    lower=Fraction(-1113,100)
    assert lower*(lower-1)>135
    return {'full_quark_Fock_dimension':2**24,'complete_occupation_characters':rows,
        'distinct_collective_irreps':len({tuple(entry['twice_highest_weight']) for _,entry in modules}),
        'Weyl_alternant_reconstructed_at_every_weight':True,
        'quark_Casimir_operator_upper_bound':'54',
        'maximal_module':{'transformed_quark_number':12,**entry},
        'all_other_color_admissible_quark_Casimir_upper_bound':'42',
        'excluded_physical_sectors_squared_transfer_upper_bound':'90+12 sqrt(14)',
        'strict_rational_squared_transfer_upper_bound':'135',
        'excluded_physical_sectors_energy_lower_bound':str(lower),
        'exclusion_includes_every_lepton_occupation_and_spin':True,
        'global_ground_proved_by_Casimir_bound_alone':False}
