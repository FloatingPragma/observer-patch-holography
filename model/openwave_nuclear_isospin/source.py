"""Exact isospin obstruction and a source-bound electrostatic comparison."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
import sympy as sy
from model.openwave_light_nuclear_curve.source import projected_source, PARAMETERS, UNIT, lift
from model.openwave_light_nuclear_curve.certificate import INPUTS, graph
from model.openwave_light_nuclear_curve.exact import rational_matrix
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals
from model.openwave_capacity_stress.source import enclosure
import model.openwave_light_nuclear_curve.source as parent
ROOT=Path(parent.__file__).resolve().parents[2]


def coefficient():
    alpha_input=json.loads((ROOT/'model/openwave_lepton_widths/inputs.json').read_text())
    length=json.loads((ROOT/'model/openwave_string_positive_cone/results.json').read_text())['common_length_attachment']['a_fm']
    hc=json.loads((ROOT/'model/openwave_transverse_gravity/inputs.json').read_text())['hbar_c_GeV_fm']
    return F(hc)/(F(alpha_input['alpha_inverse'])*F(length)*UNIT),{'alpha_inverse':alpha_input['alpha_inverse'],'length_fm':length,'hbar_c_GeV_fm':hc}


@lru_cache(None)
def electrostatic(cells,charge):
    B=F(INPUTS['three_body_coefficient']);d=projected_source(cells,B,True);g,inputs=coefficient();Iz=F(charge)-F(cells,2)
    reference=graph(cells,B);refiz=F(-1,2)if cells%2 else F(0)
    shift=-PARAMETERS['delta']*(Iz-refiz)
    Qfloor=reference['complement_floor_at_reference_Iz']+shift-2*g*(3*cells)**2/9
    top=reference['P_upper']+shift+g*charge*(charge-1)/2
    assert Qfloor>top
    error=reference['full_hopping_off_block_norm_squared_upper']/(Qfloor-top)
    Q=[np.diag(lift(np.diag([1,0,1,0]),x,cells)).real for x in range(cells)]
    electro=[sum(F(int(Q[x][j]*Q[y][j]),y-x)for x in range(cells)for y in range(x+1,cells))for j in range(d['dimension'])]
    rows=[]
    for twice_s in range(cells%2,cells+1,2):
        spin=F(twice_s,2)
        ids=np.flatnonzero((np.diag(d['Sz'])==float(spin))&(np.diag(d['Iz'])==float(Iz)))
        S=rational_matrix(d['S2'][ids][:,ids]);cols=(S-sy.Rational(spin*(spin+1))*sy.eye(len(ids))).nullspace()
        if not cols:continue
        P=sy.Matrix.hstack(*cols);G=P.T*P;H=sy.zeros(len(ids))
        for c,a in d['rational_terms']:H+=sy.Rational(c)*rational_matrix(a[ids][:,ids])
        H+=sy.diag(*[sy.Rational(g*electro[j])for j in ids])
        all_intervals=eigen_intervals(P.T*H*P,G)
        rows.append({'spin':str(spin),'rank':len(cols),'full_GeV_intervals':[[str((lo-error)*UNIT),str(hi*UNIT)]for lo,hi in all_intervals]})
    assert sum(r['rank']*(2*F(r['spin'])+1)for r in rows)==__import__('math').comb(cells,charge)*2**cells
    ground=(min(F(r['full_GeV_intervals'][0][0])for r in rows),min(F(r['full_GeV_intervals'][0][1])for r in rows))
    rank=__import__('math').comb(cells,charge)*2**cells
    dimension=reference['complete_physical_N3A_colour_singlet_dimension']
    success=F(rank,dimension);tail=reference['full_hopping_off_block_norm_squared_upper']/(Qfloor-top)**2
    return {'cells':cells,'charge':charge,
       'complete_preparation_success_exact':str(success),
       'complete_preparation_failure_exact':str(1-success),
       'original_success_high_complement_probability_upper':str(success*tail),'ground_GeV_interval':[str(x)for x in ground],'ground_GeV_display':enclosure(ground),
       'all_spin_multiplets':rows,'full_complement_floor_source':str(Qfloor),'P_upper_source':str(top),
       'uniform_full_source_error_GeV':str(error*UNIT),'entire_high_complement_retained':True,
       'complete_low_magnetic_dimension':sum(int(r['rank']*(2*F(r['spin'])+1))for r in rows),
       'input_and_energy_work':'The parent complete uniform Gauss-carrier filter with every failure; full commuting spectral readout has zero mean work.'}


def run():
    p=PARAMETERS;target=INPUTS['nuclear_rest_energy_GeV'];gap=F(target['triton'])-F(target['helion']);native=p['delta']*UNIT;needed=native-gap
    assert needed>0
    rows=[electrostatic(3,1),electrostatic(3,2),electrostatic(4,2)]
    for name,row in zip(('triton','helion','alpha'),rows):
        row['name']=name;row['target_GeV']=target[name]
        lo,hi=map(F,row['ground_GeV_interval']);row['target_in_interval']=lo<=F(target[name])<=hi
        row['target_role']='inherited triton calibration'if name=='triton'else'held-out comparison'
    g,inputs=coefficient()
    return {'schema':'oph.native_nuclear_isospin_electrostatics.v1',
      'verdict':'PASS_EXACT_ISOSPIN_OBSTRUCTION_AND_COMPLETE_ELECTROSTATIC_LOW_BANDS',
      'exact_isospin_obstruction':{'strong_triton_minus_helion_GeV':str(native),
        'observed_gap_GeV':str(gap),'additional_charge_breaking_gap_needed_GeV':str(needed),
        'minimum_charge_breaking_operator_norm_GeV':str(needed/2),
        'arbitrary_isospin_scalar_nuclear_counterterms_can_fix_both':False,
        'proof':'H=H_isospin_scalar-delta*I_z; the native raising map relates all corresponding multiplet eigenvalues exactly. Weyl gives at most2||V|| change in a level difference.'},
      'declared_electrostatic_member':{'operator':'gamma sum_(x<y) Q_x Q_y/(y-x), Q_x=N_x/6+I_z,x',
        'gamma_source_exact':str(g),'inputs':inputs,
        'geometric_domain':'three/four collinear cells at the independently strong-calibrated length',
        'pair_kernel_role':'declared Coulomb comparison; no transfer from the full dynamical field claimed',
        'additional_nuclear_mass_fits':0,'full_positive_Gauss_field_equivalence_proved':False,
        'seven_parent_mass_calibrations_transferred':False,
        'full_fixed_N_negative_perturbation_bound':'V>=-2 gamma N^2/9',
        'P_to_Q_electrostatic_coupling_exact':'0'},
      'complete_finite_spectra':rows,
      'infinite_volume_breakup_or_halo_proved':False,
      'physical_nuclear_mass_spectrum_closed':False}
