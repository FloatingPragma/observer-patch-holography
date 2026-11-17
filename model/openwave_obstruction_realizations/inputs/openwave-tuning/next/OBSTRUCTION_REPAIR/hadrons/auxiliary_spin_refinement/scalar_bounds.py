"""Spin-resolved full625-source bounds using complete representation sectors."""
from pathlib import Path
from fractions import Fraction as F
import sys,json
HERE=Path(__file__).resolve().parent;OLD=HERE.parent;sys.path.insert(0,str(OLD/'shared_auxiliary'))
from feshbach_bounds import sqrt_bounds,lower_pair

def main():
    branching=json.loads((HERE/'BRANCHING.json').read_text())['results'];source=json.loads((OLD/'SPIN_FLAVOUR_SPECTRUM.json').read_text());full=json.loads((OLD/'shared_auxiliary/FESHBACH_BOUNDS.json').read_text());local=json.loads((OLD/'shared_auxiliary/LOCAL_ALGEBRA.json').read_text());inp=json.loads((OLD/'RAW_COEFFICIENT_INPUTS.json').read_text());m=F(inp['432']['declared_light_midpoint']);omega=sqrt_bounds(F(1,64)+m*m);norms={r['N']:r for r in local['norms']};results=[]
    for old,bound in zip(source['results'],full['results']):
        kappa=F(old['kappa']);records=[]
        for group in branching:
            N=group['N'];branches={tuple(r['U4_highest_weight']):r['spin_isospin_branches'] for r in group['source_SU4_branching']}
            for sector in group['spin_isospin_sectors']:
                s,i=sector['twice_spin'],sector['twice_isospin'];eligible=[]
                for r in old['all_competing_irrep_records']:
                    if r['N']!=N:continue
                    if any(z['twice_spin']==s and z['twice_isospin']==i for z in branches[tuple(r['U4_highest_weight'])]):eligible.append(r)
                assert eligible
                plo=min(F(r['lowest_energy_exact_interval'][0]) for r in eligible);phi=min(F(r['lowest_energy_exact_interval'][1]) for r in eligible)
                diagonals=[]
                for auxnum in [1,2]:
                    q=sector['Q0_free_band_bounds'][str(auxnum)]
                    if not q['sector_present']:continue
                    coefficient=q['exact_free_energy_coefficient'];free=coefficient*(omega[1] if coefficient<0 else omega[0]);diagonals.append(auxnum*3*kappa+free)
                a2=F(sector['exact_P0_Q0_coupling_norm_squared'])
                if not sector['Q0_free_band_bounds']['1']['sector_present']:assert a2==0
                if not diagonals:
                    assert a2==0;qlo=None;lower=plo
                elif len(diagonals)==1:
                    qlo=diagonals[0];lower=lower_pair(plo,qlo,a2)
                else:
                    qlo=lower_pair(diagonals[0],diagonals[1],F(norms[N]['Q_one_to_two_aux_norm_squared_bound']));lower=lower_pair(plo,qlo,a2)
                p0hi=phi
                if N==12 and s==i==0:phi=min(phi,F(bound['vacuum_trial']['exact_full_source_variational_upper']))
                records.append(dict(N=N,twice_spin=s,twice_isospin=i,P0_highest_weight_dimension=sector['complete_P0_highest_weight_dimension'],
                    P0_ground_interval_exact=list(map(str,[plo,p0hi])),P0_max_Gsum_squared=sector['exact_P0_max_Gsum_squared'],
                    Q0_lower_exact=None if qlo is None else str(qlo),P0_Q0_norm_squared=str(a2),
                    full_source_spin_sector_interval_exact=list(map(str,[lower,phi])),
                    source_energy_interval_diagnostic=[float(lower),float(phi)]))
            globalQ=F(next(r['Q0_lower_exact'] for r in bound['all_N_sector_bounds'] if r['N']==N));assert globalQ>0
        results.append(dict(kappa=str(kappa),all_P0_present_spin_isospin_sectors=records,
            Q0_only_sectors_excluded_below_zero_by_frozen_global_Q_lower=True))
        for N in [9,12]:
            lowest=sorted([r for r in records if r['N']==N],key=lambda r:F(r['full_source_spin_sector_interval_exact'][0]))
            print('SPIN BOUNDS',kappa,N,[{k:r[k] for k in ['twice_spin','twice_isospin','source_energy_interval_diagnostic']} for r in lowest[:5]],flush=True)
    (HERE/'SCALAR_SPIN_BOUNDS.json').write_text(json.dumps(dict(schema='oph.auxiliary_spin_refinement.scalar_bounds.v1',
        actual_full625_source=True,all_spin_isospin_sectors_included=True,exact_character_and_radical_bounds=True,
        complete_full_source_baryon_spin_identification=False,results=results),indent=2)+'\n')
if __name__=='__main__':main()
