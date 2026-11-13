"""Exact complete low nuclear bands from local characters and Schur bounds."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
from .source import PARAMETERS,UNIT
from .chain_bound import full_complement
from .complement import local_rows
from .exact import certified_levels
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.positivity import joint_stability

INPUTS=json.loads((Path(__file__).parent/'inputs.json').read_text())


def upper_P(cells,B,Iz):
    p=PARAMETERS
    scalar=3*cells*p['mbar']+cells*p['A']+F(3*cells,4)*(p['Js']+p['Jw'])+3*cells*p['t']-p['delta']*Iz
    scalar+=3*(cells-2)*degree_coefficients()[1]
    return scalar+(cells-1)*(3*p['spin_gradient']+2*p['eta'])-F(3*(cells-2),4)*(p['spin_gradient']+p['eta'])+2*B*(cells-2)


def singlet_dimension(cells):
    # Exterior Cauchy: SU3 singlet at N=3A is GL(4L) Young diagram(3^A).
    A=cells;n=4*cells;value=F(1)
    for i in range(A):
        for j in range(3):value*=F(n+j-i,(A-i)+(3-j)-1)
    assert value.denominator==1
    return int(value)


def graph(cells,B):
    Iz=F(-1,2)if cells%2 else F(0);p=PARAMETERS
    Q=full_complement(cells,B,True)
    floor=Q['floor_relative_to_minus_delta_Iz']-p['delta']*Iz
    top=upper_P(cells,B,Iz);gap=floor-top;assert gap>0
    # Full N=3L hopping norm: sqrt(2)*9t for path3 and <=2*12t for path4.
    b2=F(81,32)if cells==3 else F(9)
    error=b2/gap;tail=b2/(gap*gap)
    rows=certified_levels(cells,B)
    for row in rows:
        row['full_source_GeV_intervals']=[((lo-error)*UNIT,hi*UNIT)for lo,hi in row['source_energy_intervals']]
        row['all_other_magnetic_components']='Same energies for every Sz; other Iz shift exactly by-delta*(Iz-Iz_reference).'
    count=sum(r['rank']*(2*F(r['spin'])+1)*(2*F(r['isospin'])+1)for r in rows)
    assert count==4**cells
    dimension=singlet_dimension(cells);success=F(4**cells,dimension)
    return {'cells':cells,'complete_local_quark_Fock_dimension':4096,
        'local_irrep_rows':len(local_rows()),'complete_physical_N3A_colour_singlet_dimension':dimension,
        'complete_low_magnetic_dimension':4**cells,'low_multiplets':rows,
        'exact_complement':Q,'complement_floor_at_reference_Iz':floor,'P_upper':top,
        'full_hopping_off_block_norm_squared_upper':b2,
        'uniform_low_energy_error_source':error,'uniform_low_energy_error_GeV':error*UNIT,
        'complete_original_preparation_and_spectral_records':{
            'input':'I/d on the complete physical fixed-N Gauss singlet carrier, with an opposite-charge reference on every charge block.',
            'preparation_success':success,'preparation_failure':1-success,
            'retained_failed_state':'The full complement of the local N=3,colour-singlet,spin/isospin1/2 product filter.',
            'success_and_low_spectrum_probability_interval':[success*(1-tail),success],
            'success_and_high_complement_probability_interval':[0,success*tail],
            'arbitrary_external_quantum_references':True,
            'mean_work':'Zero on I/d for the complete unital preparation filter; zero for the complete commuting spectral readout.',
            'outside_fixed_N_inputs':'Identity failure; prepared weight zero.'}}


def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:serial(v)for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v)for v in value]
    return value


@lru_cache(None)
def calculate():
    B=F(INPUTS['three_body_coefficient']);assert B>0
    stability=joint_stability(B);rows=[graph(cells,B)for cells in(3,4)]
    triton=next(r for r in rows[0]['low_multiplets']if r['spin']==r['isospin']=='1/2')['full_source_GeV_intervals'][0]
    target=F(INPUTS['nuclear_rest_energy_GeV']['triton']);assert triton[0]<target<triton[1]
    for row in rows[0]['low_multiplets']:
        if row['spin']!='1/2' or row['isospin']!='1/2':assert row['full_source_GeV_intervals'][0][0]>triton[1]
    helium3=[x-PARAMETERS['delta']*UNIT for x in triton]
    alpha=next(r for r in rows[1]['low_multiplets']if r['spin']==r['isospin']=='0')['full_source_GeV_intervals'][0]
    assert helium3[1]<F(INPUTS['nuclear_rest_energy_GeV']['helion'])
    assert alpha[0]>F(INPUTS['nuclear_rest_energy_GeV']['alpha'])
    # The N<=6 operator has no triple-occupation window, not merely a small expectation.
    return serial({'schema':'oph.light_nuclear_curve.v1',
        'verdict':'PASS_COMPLETE_A3_A4_LOW_BANDS_TRITON_CALIBRATION_AND_HELIUM_FAILURES',
        'three_body_coefficient':B,'same_degree_coefficients':degree_coefficients(),
        'full_Fock_joint_stability':stability,'complete_graphs':rows,
        'mass_comparisons':{
            'triton':{'source_GeV_interval':triton,'target_GeV':str(target),'role':'calibration','spin':'1/2','isospin':'1/2'},
            'helion':{'source_GeV_interval':helium3,'target_GeV':INPUTS['nuclear_rest_energy_GeV']['helion'],'role':'held_out','matches':False},
            'alpha':{'source_GeV_interval':alpha,'target_GeV':INPUTS['nuclear_rest_energy_GeV']['alpha'],'role':'held_out','matches':False,'spin':'0','isospin':'0'}},
        'inherited_calibration_domain_unchanged':'The occupancy-gated range-three interaction is exactly zero for total N<=6. The degree terms are exactly zero on two cells. All seven two-cell light calibration domains are retained.',
        'additional_calibration':'One triton mass sets the range-three coefficient. Its mathematical spectral enclosure, rather than a Ritz value alone, contains the target.',
        'scope':'Complete finite-box native low nuclear spectra and original observer records on 3 and 4 cells, using one shared local action. The helium comparisons fail. No asymptotic breakup threshold, observed halo, complete binding curve or infinite-volume nuclear mass is inferred from these boxes.'})
