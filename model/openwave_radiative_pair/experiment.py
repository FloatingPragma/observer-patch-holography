"""Complete quantum radiation experiment, finite cutoffs and countercontrols."""
from .radiation import run_probe,source_controls,construction
from .exact import certificate as exact_certificate
from model.lepton_mass_compiler.compiler import electron_calibration
import numpy as np

def run():
    controls=source_controls();exact=exact_certificate()
    rows={}
    for cap in(1,2):
        rows[str(cap)]={'pair':run_probe(cap,times=(0,5,20,50)),
            'vacuum_control':run_probe(cap,times=(20,50),initial_kind='vacuum')}
    zero=run_probe(1,eta=0,times=(20,));detuned=run_probe(1,times=(20,50),omega_ratio=1.2)
    assert zero[0]['two_photon_no_fermion_probability']==0
    for cap,values in rows.items():
        assert values['pair'][-1]['two_photon_no_fermion_probability']>.02
        assert abs(values['pair'][-1]['full_energy_before_readout']-2)<2e-8
        assert values['pair'][-1]['forbidden_total_momentum_probability']<1e-24
    return {'schema':'oph.native_current_quantum_radiation.v1','native_source_controls':controls,'exact_short_time_certificate':exact,
        'declared_source':{'parent':'native charged Dirac frame with analyzer mass; one family on a periodic three-site ring',
            'all_other_native_modes':'retained as exact free spectator factors; the coupling acts only on the declared family/ring factor',
            'ring_spacing':'sqrt(3)/m; an exposed resonance choice giving cochain photon frequency m',
            'transverse_photon_modes':'q=+/-2pi/3, both x/y polarizations; finite occupation caps1 and2',
            'interaction':'eta*m sum_q,pol (a_q,pol^dagger dGamma(J_-q,pol)+h.c.), eta=1/20',
            'current':'J is the actual derivative -tau3 tensor sigma_pol in the native charged frame, with literal spatial Fourier momentum transfer',
            'not_inserted_annihilation_projector':True,'not_a_full_PW119_or_local_Gauss_QED_transfer':True,
            'no_bosonic_infinite_CCR_claim':True,'preparation':'one scalar rest particle-hole pair over the filled free negative bands, photon vacuum',
            'same_action_actual_interacting_ground_prepared':False},
        'complete_cap_experiments':rows,'zero_coupling_control':zero,'detuned_photon_control':detuned,
        'conditional_calibration':electron_calibration(),
        'readout':'joint full photon occupation with source free-ground/complement; all other fermion outcomes summed explicitly, retaining quantum branch states',
        'reference_domain':'the finite source-field Hamiltonian is Hermitian on the entire Fock direct sum, so its unitary and complete instruments extend to arbitrary idle references',
        'energy_boundary':'The initial pair/free-vacuum expectation difference is2m. The full H expectation is conserved. Number readout can remove interaction energy; its work is printed. Free ground is not called the interacting ground.',
        'interpretation':'Executed native-current conversion to two transverse cochain radiation quanta in a declared finite quantum extension. No observed positron identification, continuum annihilation rate, asymptotic outgoing scattering or fine-structure constant is derived.',
        'numerical_scope':'Longer-time sparse exponential outputs and cap comparisons are numerical diagnostics, not outward interval bounds. The positive short-time certificate is analytic.',
        'verdict':'PASS_NATIVE_CURRENT_PAIR_TO_TRANSVERSE_RADIATION'}
