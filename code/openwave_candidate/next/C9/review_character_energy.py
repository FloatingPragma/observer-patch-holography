"""Independent analytic cycle spectrum and all fixed-number CAR energies.

No producer, CAR implementation or Fable Hamiltonian is imported. The three
port Laplacian eigenvalues at flux theta are 2-2 cos((2 pi r+theta)/3).
The six port/transit energies are both square-root signs; complete three-mode
occupation energies are their twenty subset sums.
"""
from pathlib import Path
from itertools import combinations
import math,json,hashlib
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]


def levels(theta,mass=.23):
    positive=[math.sqrt(mass*mass+2-2*math.cos((2*math.pi*r+theta)/3)) for r in range(3)]
    single=sorted([-x for x in positive]+positive)
    many=sorted(sum(x) for x in combinations(single,3))
    return single,many


def run():
    reports=[];paths=[Path(__file__)]
    for modulus,name in ((3,'CHARACTER_REFERENCE.json'),(6,'CHARACTER_MINIMA.json')):
        path=HERE/name;paths.append(path);value=json.loads(path.read_text())
        actual=value['sector_spectra'] if modulus==3 else value['characters']
        errors=[];all_levels=[];minima=[];signed=[];positive=[]
        for chi,row in enumerate(actual):
            one,three=levels(2*math.pi*chi/modulus)
            saved=row['spectrum'] if modulus==3 else row['CAR_spectrum']
            errors.append(max(abs(a-b) for a,b in zip(three,saved)))
            all_levels.extend(three);minima.append(three[0]);signed.append(one[0]);positive.append(min(x for x in one if x>0))
        all_levels.sort();ground=all_levels[0];multiplicity=sum(abs(x-ground)<1e-10 for x in all_levels)
        gap=next(x-ground for x in all_levels if x>ground+1e-10)
        recorded=value['global_ground_energy'] if modulus==3 else value['global_ground']
        recorded_gap=value['fixed_number_gap_above_ground_multiplet'] if modulus==3 else value['global_fixed_number_gap']
        if max(errors)>1e-12 or abs(ground-recorded)>1e-12 or abs(gap-recorded_gap)>1e-12:raise ValueError('analytic all-character spectrum differs')
        argmin=lambda xs:[i for i,x in enumerate(xs) if abs(x-min(xs))<1e-10]
        reports.append({'modulus':modulus,'all_twenty_levels_per_character_max_error':max(errors),
          'global_ground':ground,'ground_multiplicity':multiplicity,'gap_above_ground_multiplet':gap,
          'filled_CAR_minimizing_characters':argmin(minima),'lowest_signed_one_load_minimizing_characters':argmin(signed),
          'lowest_positive_one_load_minimizing_characters':argmin(positive),
          'sector_ground_energies':minima,'sector_positive_edges':positive})
    return {'schema':'oph.c9.analytic_character_energy_review.v1','outcome':'PASS','reviewer':'Codex coordinator; independent analytic circulant and subset-sum calculation',
      'results':reports,'scope':'Full finite spectral comparison only. Mod6 class3 minimizes the filled CAR and lowest signed one-load energies; class0 minimizes the positive one-load edge. No dynamic cooling, character transfer, spin derivation, physical mass or independent state preparation follows. Native-accessible response gaps require actual outcome probabilities and are distinct from global spectral gaps.',
      'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}

if __name__=='__main__':
    result=run();HERE.joinpath('ANALYTIC_CHARACTER_ENERGY_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['results'],indent=2))
