from pathlib import Path
from functools import lru_cache
import numpy as np
from model.openwave_intrinsic_magnetic.solver import dirac,source_current,moment,SPIN,GAMMA,BETA,SIGMA,native_reduction
from model.openwave_field_observables.lattice import curl,curl_adjoint,divergence_adjoint
from .certified_probe import run as certified_probe
HERE=Path(__file__).resolve().parent


def check_witness(path):
    data=np.load(path,allow_pickle=False)
    p=dict(zip(data['parameter_names'].tolist(),data['parameters'].tolist(),strict=True))
    for key in ('n','field_n','spin_sign'):
        assert p[key]==int(p[key]);p[key]=int(p[key])
    psi=data['state'];A=data['field'];a=p['spacing'];k=p['coupling']
    assert abs(np.vdot(psi,psi)-1)<1e-12
    h,ph=dirac(p['n'],a,p['field_n'],A,p['mass'],p['trap'],p['wilson'],p['external_B'])
    energy=float(np.vdot(psi,h@psi).real);residual=float(np.linalg.norm(h@psi-energy*psi))
    assert residual<1e-10 and energy>p['mass']
    j=source_current(psi,ph,p['field_n'],p['wilson']);b=curl(A,a)
    maxwell=float(np.linalg.norm(curl_adjoint(b,a)-k*j/a**3))
    continuity=float(np.linalg.norm(divergence_adjoint(j,a)))
    assert maxwell<1e-9 and continuity<1e-10
    gauge=float(np.linalg.norm(divergence_adjoint(A,a)))
    harmonic=float(np.linalg.norm(A.sum(axis=(0,1,2))))
    assert gauge<1e-10 and harmonic<1e-10
    spin=float(np.sum(SPIN.diagonal()[:,None]*abs(psi.reshape(4,-1))**2).real)
    mu=moment(j,a);field_energy=a**3*np.sum(b*b)/(2*k);rest=energy+field_energy
    return {'witness':path.name,'parameters':p,'matter_boundary_extent':(p['n']+1)*a,
        'periodic_field_extent':p['field_n']*a,'positive_D_eigenvalue':energy,
        'actual_rest_energy':float(rest),'field_energy':float(field_energy),'spin':spin,'moment':mu.tolist(),
        'g_using_supplied_mass':float(-2*p['mass']*mu[2]/spin),
        'g_using_actual_rest_energy':float(-2*rest*mu[2]/spin),
        'g_using_doublet_half_generator_and_actual_energy':float(-4*rest*mu[2]),
        'source_equation_residuals':{'eigenstate':residual,'Maxwell':maxwell,'current_continuity':continuity,
            'transverse_gauge':gauge,'zero_harmonic_field':harmonic}}


@lru_cache(maxsize=1)
def run():
    algebra=native_reduction()
    chirality=GAMMA[0]@GAMMA[1]@GAMMA[2]@BETA
    theta=np.block([[1j*SIGMA[1],np.zeros((2,2))],[np.zeros((2,2)),1j*SIGMA[1]]])
    assert np.array_equal(chirality@chirality,np.eye(4))
    assert np.array_equal(theta@theta.conj(),-np.eye(4))
    for gamma in (*GAMMA,BETA):assert not np.any(chirality@gamma+gamma@chirality)
    for gamma in GAMMA:assert np.array_equal(theta@gamma.conj()@theta.conj().T,-gamma)
    assert np.array_equal(theta@BETA.conj()@theta.conj().T,BETA)
    values={p.stem:check_witness(p) for p in sorted((HERE/'witnesses').glob('*.npz')) if not p.stem.startswith('probe_B') and p.stem!='probe_zero'}
    refinements=[values['stationary_n'+str(n)] for n in (5,9,13,17)]
    assert all(abs(r['matter_boundary_extent']-6)<1e-12 and abs(r['periodic_field_extent']-21)<1e-12 for r in refinements)
    a=np.array([r['parameters']['spacing'] for r in refinements]);g=np.array([r['g_using_actual_rest_energy'] for r in refinements])
    assert np.all(np.diff(g)>0)
    fits=[]
    for label,degree,index in (('linear_finest_two',1,slice(-2,None)),('quadratic_first_three',2,slice(None,3)),('quadratic_finest_three',2,slice(1,None)),('cubic_all_four',3,slice(None,None))):
        coeff=np.polyfit(a[index],g[index],degree)
        fits.append({'fit':label,'coefficients':coeff.tolist(),'zero_spacing_intercept':float(coeff[-1]),
            'finest_grid_residual':float(np.polyval(coeff,a[-1])-g[-1])})
    base=values['stationary_n9'];probes=[]
    for B,label in ((.02,'0p02'),(.01,'0p01'),(.005,'0p005')):
        plus=values['probe_n9_'+label];minus=values['probe_n9_minus'+label]
        derivative=(plus['actual_rest_energy']-minus['actual_rest_energy'])/(2*B)
        probes.append({'absolute_B':B,'two_sided_energy_derivative':derivative,
            'minus_zero_field_moment':-base['moment'][2],
            'difference':derivative+base['moment'][2],
            'g_from_symmetric_probe':2*base['actual_rest_energy']*derivative/base['spin']})
    assert abs(probes[-1]['difference'])<abs(probes[0]['difference'])
    return {'schema':'oph.finite_magnetic_refinement.v1',
        'complete_native_algebra':algebra,'zero_probe_Kramers_and_all_probe_chiral_identities_checked':True,
        'witness_method':'Numerical states and fields are immutable witnesses. The producer rebuilds the declared operator, current and Maxwell equations and verifies their residuals. Regeneration uses the same full sparse eigenproblem and no target g.',
        'fixed_both_volume_spacing_sequence':refinements,
        'independent_fixed_spacing_matter_volume_sequence':[values['volume_field_N31'],values['volume_matter_n7'],values['volume_matter_n9']],
        'independent_field_volume_sequence':[values['volume_field_N13'],values['stationary_n5'],values['volume_field_N31']],
        'external_field_amplitude_controls':probes,'all_probe_witnesses':[v for k,v in values.items() if k.startswith('probe_n9')],
        'spacing_extrapolation_diagnostics':fits,
        'continuum_intercept_is_certified':False,
        'fit_order_and_volume_limitations':'The finite-g increase and held-out prediction are measured; intercepts depend on fit order. Fixed-spacing matter-volume changes are material. No fit is promoted to a physical or continuum g prediction.',
        'same_refined_source_precession':certified_probe(),
        'limits':{'finite_bound_state_g_and_actual_precession_on_one_declared_member':True,
            'common_physical_electron_dictionary':False,'mass_or_alpha_derived':False,
            'full_propagating_quantum_electromagnetism':False,
            'reference_uniform_linear_CP_claim_for_mean_field_dynamics':False},
        'verdict':'PASS_MAGNETIC_SPACING_VOLUME_PROBE_AND_ALL_MODE_LARMOR_CONTROLS'}
