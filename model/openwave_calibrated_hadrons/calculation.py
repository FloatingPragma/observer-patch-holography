"""Source-bound spectra, explicit fits, withheld comparisons and complete records."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from model.openwave_six_quark.certificate import spectral_enclosure
from .source import INPUTS,parameters,sector,exact_matrix,components


def source_charge(kind,label):
    return sum(x[2]==0 for x in label)-1 if kind=='baryon'else label[1][2]-label[0][2]


@lru_cache(None)
def spectra():
    out={}
    for kind in ('baryon','meson'):
        d=sector(kind);H=exact_matrix(kind);rows=[]
        for charge in sorted({source_charge(kind,x)for x in d['labels']}):
            ids=[i for i,x in enumerate(d['labels'])if source_charge(kind,x)==charge]
            outside=[i for i in range(H.rows)if i not in ids]
            assert all(H[i,j]==0 for i in ids for j in outside)
            witness=spectral_enclosure(H.extract(ids,ids))
            rows.append({'charge':charge,'source_basis_indices':ids,'certificate':witness})
        out[kind]=rows
    return out


def selected(kind,charge,spin):
    d=sector(kind);n=len(d['h']);q=d['q'];s=d['spin'];z=d['z'];I=np.eye(n)
    f=(q-charge*I)@(q-charge*I)+(s-spin*(spin+1)*I)@(s-spin*(spin+1)*I)+(z-spin*I)@(z-spin*I)
    e,v=eigh(f);P=v[:,e<1e-8];energy,U=eigh(P.conj().T@d['h']@P)
    # The entire nearly degenerate parity band is retained. A declared exact
    # basis-seed filter prepares its numerical spectral projection at t=0.
    line=P@U[:,0];full=d['B']@line;prob=abs(full)**2
    assert abs(prob.sum()-1)<2e-12
    return {'kind':kind,'charge':charge,'internal_spin':spin,'selected_source_dimension':P.shape[1],
       'all_selected_eigenvalues':energy.tolist(),'lowest_energy':float(energy[0]),
       'conditional_GeV':float(F(INPUTS['electron_unit_GeV']))*float(energy[0]),
       'complete_uniform_input_dimension':len(d['basis']),
       'rank_one_filter_success_probability':1/len(d['basis']),
       'complete_success_native_occupation_probabilities':(prob/len(d['basis'])).tolist(),
       'complete_failure_native_occupation_probabilities':((1-prob)/len(d['basis'])).tolist(),
       'orthogonal_full_Fock_complement_probability':0.,
       'nonselective_filter_and_full_occupation_mean_work':0.,
       'initial_uniform_sector_energy':float(d['H'].diagonal().real.mean()),
       'spectral_residual_diagnostic':float(np.linalg.norm(d['h']@line-energy[0]*line)),
       'record_scope':'Complete t=0 uniform-input/rank-one-filter occupation instruments. Floating spectral vectors are diagnostic; the separately stored full-spectrum enclosures certify energies. No unique-parity-line or arbitrary-time rounded-vector stationarity is asserted.'}


def certificate():
    p=parameters();unit=F(INPUTS['electron_unit_GeV']);rows=[]
    for kind,charges,spin in [('baryon',(0,1),.5),('baryon',(-1,0,1,2),1.5),('meson',(-1,0,1),0.),('meson',(-1,0,1),1.)]:
        rows.extend(selected(kind,q,spin)for q in charges)
    targets=INPUTS['calibration_targets_GeV'];fits=[];spec=spectra()
    difference=(F(targets['proton'])+F(targets['neutron']))/(2*unit)-F(3,2)*F(targets['charged_pion'])/unit
    bar_trial=F(3,4)*(p['Js']+p['Jw']+p['eta'])+3*p['t']
    pion_trial=F(3,4)*(p['Jw']+p['eta'])+2*p['t']
    left=bar_trial-difference;right=F(5,9)*4096-F(3,2)*pion_trial-difference
    mass_lower=(F(targets['charged_pion'])/unit-pion_trial)/2
    assert left<0<right and mass_lower>p['delta']/2
    for name,kind,q,s in [('proton','baryon',1,.5),('neutron','baryon',0,.5),('charged_pion','meson',1,0.)]:
        row=next(x for x in rows if(x['kind'],x['charge'],x['internal_spin'])==(kind,q,s));target=float(targets[name]);error=row['conditional_GeV']-target
        assert abs(error)<3e-11
        enclosure=next(x for x in spec[kind]if x['charge']==q)['certificate']['all_eigenvalue_intervals'][0]
        interval=[unit*F(x)for x in enclosure]
        assert interval[0]<F(targets[name])<interval[1]
        fits.append({'comparison':name,'target_GeV':targets[name],'computed_GeV':row['conditional_GeV'],
                     'computed_energy_interval_GeV':list(map(str,interval)),
                     'residual_GeV':error,'status':'CALIBRATION_INPUT_REPRODUCED_NOT_INDEPENDENT_PREDICTION'})
    comparisons=[]
    for name,kind,q,s in [('neutral_pion','meson',0,0.),('rho_vector','meson',0,1.),('Delta_baryon','baryon',2,1.5)]:
        row=next(x for x in rows if(x['kind'],x['charge'],x['internal_spin'])==(kind,q,s));value=float(INPUTS['held_out_comparisons_GeV'][name])
        comparisons.append({'comparison':name,'reference_GeV':value,'computed_candidate_GeV':row['conditional_GeV'],
          'fractional_disagreement':row['conditional_GeV']/value-1,'physical_assignment_supported':False})
    assert p['delta']>1
    return {'schema':'oph.positive_native_calibrated_hadron_source.v1','verdict':'PASS_POSITIVE_NATIVE_HADRON_CALIBRATION_WITH_EXPLICIT_HELD_OUT_FAILURES',
      'full_source':'60 native CAR modes on two cells; Q/u^c/d^c native colour currents, positive quark/conjugate number masses, colour and spin/weak costs, baryon-density square and finite-edge weak-current square. Charged lepton energy1, neutral energy0; their full complementary Fock sectors are retained.',
      'exact_action':'sum m_f N_f+t sum_edges(c_x-c_y)^dagger(c_x-c_y)+U sum C_x^2+Js sum S_x^2+Jw sum I_x^2+A sum B_x^2+eta sum_edges(I_x+I_y)^2',
      'native_baryon_density':'B_x=(N_Q,x-N_(u^c,d^c),x)/3; representation-defined one-body current, not a composite projector',
      'parameters_rational':{k:str(v)for k,v in p.items()},'positive_quark_mass_coefficients':{'up':str(p['mbar']-p['delta']/2),'down':str(p['mbar']+p['delta']/2)},
      'calibration_inputs':INPUTS,'three_calibrated_coefficients':['mbar','delta','A'],
      'declared_unfitted_coefficients':['U=2048','Js=Jw=1/4','t=1/8','eta=1'],
      'unique_calibration_shortcut':'Ebar(A)-3Epi(A)/2 has slope at least2/9: baryon B-density square in[5/9,1], meson in[0,2/9]. mbar then follows linearly and delta equals the fitted neutron-proton difference.',
      'exact_calibration_existence_bounds':{'objective_upper_at_A_zero':str(left),'objective_lower_at_A_4096':str(right),'positive_mbar_lower':str(mass_lower),'monotonic_increment_lower':'2h/9'},
      'fits':fits,'complete_colour_spectra':spec,'complete_particle_records':rows,'held_out_comparisons':comparisons,
      'exact_beta_available_energy_GeV':str(unit*(p['delta']-1)),
      'beta_channel_kinematically_allowed':True,'beta_matrix_element_rate_or_parity_proved_here':False,
      'number_mass_only_control':'A free positive number-mass fit gives baryon/meson ratio3/2, incompatible with the fitted nucleon/pion ratio. The added B-density square supplies an independent interacting baryon scale.',
      'positive_full_Fock_floor':'H >= min(m_up,m_down) N_colourful+N_charged_lepton; massless neutral spectators prevent a claimed positive global vacuum gap.',
      'physical_boundary':'Three empirical hadron masses fix three action coefficients. No independent prediction of those masses, complete observed spectrum, relativistic pole dictionary, physical clock, local gauge completion or continuum nuclear limit follows. Near-degenerate parity states, vector mesons and spin-three-halves partners are retained.'}
