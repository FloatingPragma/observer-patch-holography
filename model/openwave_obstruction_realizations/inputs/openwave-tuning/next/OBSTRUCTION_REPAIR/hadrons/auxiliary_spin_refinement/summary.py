"""Same-source vacuum, baryon spin and mass consequences of exact bounds."""
from pathlib import Path
from fractions import Fraction as F
import json
HERE=Path(__file__).resolve().parent;OLD=HERE.parent

def main():
    qb=json.loads((HERE/'Q_COLOUR_BOUNDS.json').read_text())['results'];mb=json.loads((HERE/'MATRIX_FESHBACH.json').read_text())['results'];vb=json.loads((OLD/'shared_auxiliary/FESHBACH_BOUNDS.json').read_text())['results'];ph=json.loads((HERE/'PARTICLE_HOLE.json').read_text());assert ph['actual_full625_source_intertwiner']
    inp=json.loads((OLD/'RAW_COEFFICIENT_INPUTS.json').read_text());unit=F(inp['electron_anchor_GeV'])/F(inp['432']['declared_light_midpoint']);out=[]
    for q,m,v in zip(qb,mb,vb):
        vlo=min(F(r['refined_full_source_lower_exact']) for r in q['records'] if r['N']==12);vhi=F(v['vacuum_trial']['exact_full_source_variational_upper']);
        other_vac=min(F(r['refined_full_source_lower_exact']) for r in q['records'] if r['N']==12 and (r['twice_spin'],r['twice_isospin'])!=(0,0))
        blo=min([F(r['full_source_lower_exact']) for r in m['records']]+[F(r['refined_full_source_lower_exact']) for r in q['records'] if r['N']==9 and (r['twice_spin'],r['twice_isospin']) not in [(1,1),(3,3)]])
        bhi=min(F(r['trial']['full_source_upper_exact']) for r in m['records']);gap=[blo-vhi,bhi-vlo]
        # Q-only physical spin sectors have positive energy by the uniform
        # frozen Q0 bound; they cannot undercut any of these negative trials.
        for N in [9,12]:assert F(next(r['Q0_lower_exact'] for r in v['all_N_sector_bounds'] if r['N']==N))>0
        assert v['N12_vacuum_sector_certified'] and vlo<=vhi and blo<=bhi
        out.append(dict(kappa=q['kappa'],N12_vacuum_sector_certified=True,
            vacuum_energy_interval_exact=list(map(str,[vlo,vhi])),vacuum_energy_interval_diagnostic=[float(vlo),float(vhi)],
            all_other_vacuum_spin_isospin_lower_exact=str(other_vac),vacuum_spin0_isospin0_certified=other_vac>vhi,
            baryon_energy_interval_exact=list(map(str,[blo,bhi])),baryon_energy_interval_diagnostic=[float(blo),float(bhi)],
            lowest_Bplusminus1_spin_half_isospin_half_certified=m['full_source_baryon_spin_half_isospin_half_certified'],
            Bplus_and_Bminus_equal_by_exact_full_source_intertwiner=True,
            same_source_baryon_gap_exact_interval=list(map(str,gap)),conditional_baryon_gap_GeV_interval=[float(x*unit) for x in gap],
            conditional_proton_comparison=('excluded_above_1_GeV' if gap[0]*unit>1 else 'excluded_below_half_GeV' if gap[1]*unit<F(1,2) else 'broad_interval_only_no_mass_prediction')))
    result=dict(schema='oph.auxiliary_spin_refinement.summary.v1',unchanged_source_parameters=dict(t='1/8',g='1/8',Delta='3*kappa',kappa=['1/2','1','2']),
        conditional_energy_unit_GeV_exact=str(unit),conditional_energy_unit_GeV=float(unit),common_spatial_Dirac_transfer=False,
        no_empirical_parameter_selection=True,actual_full625_source_not_effective_eigenvectors=True,results=out,
        unresolved=['weakest-coupling baryon spin','weakest-coupling vacuum spin','precise absolute mass at kappa1','isospin-symmetric neutron/proton degeneracy','full-source B2 binding and mobility','native preparation and readout','common Dirac clock and interacting pole map'],
        physical_rows_solved=False)
    (HERE/'SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
