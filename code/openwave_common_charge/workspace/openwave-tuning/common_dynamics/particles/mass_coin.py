"""One explicit local mass-coin amendment of the frozen carrier walk.

The parameter m is supplied architecture. No observed particle mass is input
or inferred. The same event reflection supplies the mass generator.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WALK = ROOT / 'openwave-tuning/auditor/native_walk_bridge'
REVIEW = ROOT / 'openwave-tuning/auditor/common_walk_review'
sys.path[:0] = [str(WALK), str(REVIEW), str(HERE.parent/'clock'),
               str(ROOT/'reverse-engineering-reality')]
import common_walk as common
import compile_walk as compiled
import source_error
from clock_contract import effective_time
from model.observer.source_scalar import compiler
from model.certified_reference import certificate as cert, interval as iv
from model.native_spatial_refinement import source as spatial

MASS_ID = 'native-local12-mass-coin-v1-proposed'


def rational(value, name, positive=False):
    if type(value) not in (int, F) or (positive and value <= 0):
        raise ValueError(name+' must be an exact '+('positive ' if positive else '')+'rational')
    return F(value)


def parameters(h, mass):
    h, mass = rational(h, 'spacing', True), rational(mass, 'mass parameter')
    if abs(h*mass) > F(1,2):
        raise ValueError('the declared mass-coin window is |h*m|<=1/2')
    return h, mass


def kick_word(h, mass):
    h,mass = parameters(h,mass)
    mu = h*mass
    prepare = compiled.controls()['preparation_word']
    diagonal = compiler.compile_diagonal([-float(mu)]+[float(mu)]*15)
    return compiler.inverse(prepare)+diagonal+prepare


def coin_word(h, mass):
    """Fuse the mass kick and R; keep the complete scalar/padding phase."""
    h,mass=parameters(h,mass)
    mu=h*mass
    prepare=compiled.controls()['preparation_word']
    diagonal=compiler.compile_diagonal([-float(mu)]+[math.pi+float(mu)]*15)
    return compiler.inverse(prepare)+diagonal+prepare


@lru_cache(None)
def coin_matrix(h,mass):
    return compiler.execute(coin_word(h,mass),4)


@lru_cache(None)
def kick_matrix(h, mass):
    return compiler.execute(kick_word(h,mass),4)


def step(state, h=F(1), mass=F(0)):
    """Actual source-word kick, then the same frozen W; no spectral input."""
    value = np.asarray(state,complex)
    if value.shape != (128,16) or not np.isfinite(value).all():
        raise ValueError('complete finite128-by16 state required')
    reflected=value@coin_matrix(h,mass).T
    result=np.empty(2048,complex)
    result[compiled.controls()['swap_permutation']]=reflected.reshape(-1)
    return result.reshape(128,16)


def native_error(h, mass):
    h,mass = parameters(h,mass)
    mu=h*mass
    prepare=compiled.controls()['preparation_word']
    Uerror,_=source_error.uniform_bound(prepare)
    word=kick_word(h,mass)
    middle=word[len(prepare):-len(prepare)]
    template=[]
    for mask in range(16):
        coefficient=7*mu/8 if mask==0 else -mu/8
        label=''.join('Z' if mask>>q&1 else 'I' for q in range(4))
        template+=cert.pauli(label,coefficient)
    diagonal=source_error.exact_template_bound(middle,template)
    extra=2*Uerror+diagonal
    base=json.loads((REVIEW/'SOURCE_ERROR.json').read_text())
    fused=coin_word(h,mass)
    middle=fused[len(prepare):-len(prepare)]
    offset,fixed=0,F()
    for mask in range(16):
        label=''.join('Z' if mask>>q&1 else 'I' for q in range(4))
        offset,error,variables=cert.consume(middle,offset,cert.pauli(label))
        assert len(variables)==1
        pcoef=F(15,16) if mask==0 else F(-1,16)
        mcoef=7*mu/8 if mask==0 else -mu/8
        fixed+=error+(iv.Interval.point(variables[0])-pcoef*iv.pi()-mcoef).absolute_upper()
    assert offset==len(middle)
    fused_error=2*Uerror+fixed+cert.controls_bound(middle,source_error.CONTROL_TABLE)
    return {'mass_parameter_exact':str(mass),'spacing_exact':str(h),'mu_exact':str(mu),
            'kick_native_pulses':len(word),'extra_operator_error_exact':str(extra),
            'whole_step_operator_error_exact':str(fused_error+F(base['swap_operator_error_exact'])),
            'fused_coin_native_pulses':len(fused),
            'fused_whole_step_native_pulses':compiled.controls()['source_swap_native_pulses']+len(fused),
            'word':word,'fused_coin_word':fused,
            'scope':'Complete source/export bound under the same ideal coherent services; initial and detector resets are additional.'}


def continuum_bound(h, mass, bandwidth, steps, tau0=F(1)):
    h,mass=parameters(h,mass)
    K=rational(bandwidth,'bandwidth')
    tau0=rational(tau0,'time scale',True)
    if K<0 or h*K>F(1,2) or type(steps) is not int or steps<0:
        raise ValueError('nonnegative bandwidth/steps with h*K<=1/2 required')
    L=K+abs(mass); b=K/4+abs(mass); x=5*K/4
    coefficient=x*(L+b)+F(13,16)*L*L+b*b/2
    error=5*h*K/2+steps*h*h*coefficient
    dispersion=(F(3,4)*(mass*mass+5*K*K/64)**2+mass**4/12
                +5*mass*mass*K*K/128+K**4/16)
    return {'time_exact':str(effective_time(steps,h,tau0)),
            'complete_terminal_TV_bound_exact':str(min(F(1),error)),
            'uncapped_vector_bound_exact':str(error),'one_step_coefficient_exact':str(coefficient),
            'squared_frequency_error_upper_exact':str(h*h*dispersion/(tau0*tau0)),
            'domain':'Bandlimited P0 preparation; arbitrary reference and terminal instrument. Retains all26 low modes and actual optical outputs. Native implementation errors and reset laws add separately.'}


def run():
    paths=[Path(__file__),Path(common.__file__),Path(compiled.__file__),Path(source_error.__file__),
           REVIEW/'SOURCE_ERROR.json',HERE.parent/'continuum/MANIFEST.json',
           HERE.parent/'clock/MANIFEST.json',Path(spatial.__file__)]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    u=np.r_[np.ones(12)/np.sqrt(12),np.zeros(4)]
    R=2*np.outer(u,u)-np.eye(16)
    T=common.native_transfer();values,vectors=np.linalg.eigh(T)
    swap=common.expected_swap_permutation()
    def S(value):
        out=np.empty(2048,complex);out[swap]=value.reshape(-1)
        return out.reshape(128,16)
    source=common.initial_cell_state();rows=[]
    errors=[]
    for h,mass in ((F(1),F(0)),(F(1),F(1,4)),(F(1,4),F(1)),(F(1,8),F(-1))):
        mu=float(h*mass)
        actual=kick_matrix(h,mass)
        target=np.cos(mu)*np.eye(16)-1j*np.sin(mu)*R
        error=float(np.linalg.norm(actual-target,2));errors.append(error)
        assert error<1e-11
        fused=coin_matrix(h,mass)
        assert np.linalg.norm(fused-(np.cos(mu)*R-1j*np.sin(mu)*np.eye(16)))<1e-11
        assert np.linalg.norm(actual[:,12:]-np.exp(1j*mu)*np.eye(16)[:,12:])<1e-11
        residual=0.
        for j,t in enumerate(values):
            A=np.zeros((128,16),complex)
            for x in range(96):A[common.encode_position(x)]=vectors[x,j]*u
            SA=S(A)
            if t<1-1e-11:
                q=(SA-t*A)/math.sqrt(1-t*t)
                basis=np.column_stack((A.reshape(-1),q.reshape(-1)))
                matrix=basis.conj().T@np.column_stack((step(A,h,mass).reshape(-1),step(q,h,mass).reshape(-1)))
                residual=max(residual,float(abs(np.linalg.det(matrix)-1)),float(abs(np.trace(matrix)-2*t*np.cos(mu))))
                phases=np.linalg.eigvals(matrix)
                omega=np.arccos(t*np.cos(mu))
                residual=max(residual,float(np.max(np.min(abs(phases[:,None]-np.exp(1j*np.array([omega,-omega]))[None,:]),axis=1))))
            else:
                residual=max(residual,float(np.linalg.norm(step(A,h,mass)-np.exp(-1j*mu)*A)))
        assert residual<2e-11
        dark=np.zeros((128,16),complex);dark[0,0]=1/np.sqrt(2);dark[0,2]=-1/np.sqrt(2)
        assert np.linalg.norm(step(dark,h,mass)+np.exp(1j*mu)*dark)<1e-11
        invalid=[x for x in range(128)if common.decode_position(x)is None]
        assert len(invalid)==32
        for address in invalid:
            padding=np.zeros((128,16),complex);padding[address]=u
            assert np.linalg.norm(step(padding,h,mass)-np.exp(-1j*mu)*padding)<1e-11
        for address in range(128):
            padding=np.zeros((128,16),complex);padding[address,12]=1
            assert np.linalg.norm(step(padding,h,mass)+np.exp(1j*mu)*padding)<1e-11
        assert abs(np.vdot(step(source,h,mass),step(source,h,mass)).real-1)<1e-11
        if not mass:assert np.linalg.norm(step(source,h,mass)-common.apply_walk(source))<1e-11
        rows.append({'h':str(h),'mass':str(mass),'full16_column_kick_error':error,
                     'all96_source_spectral_identity_residual':residual,'native_error':native_error(h,mass)})
    bounds=[]
    for h in (F(1,8),F(1,16),F(1,32),F(1,64)):
        m=F(1);K=F(1);n=int(1/h)
        law=continuum_bound(h,m,K,n)
        for p in ((1.,0.,0.),(.5,.5,0.),(.2,.3,.4)):
            lam=float(np.linalg.eigvalsh(spatial.symbol([float(h)*x for x in p]))[0])
            omega=math.acos((1-lam/12)*math.cos(float(h*m)))
            difference=abs((omega/float(h))**2-(float(m*m)+sum(x*x for x in p)/48))
            assert difference<=float(F(law['squared_frequency_error_upper_exact']))+1e-12
        bounds.append({'h':str(h),**law})
    assert continuum_bound(F(1,8),F(0),F(1),8)['one_step_coefficient_exact']=='77/32'
    # Consume the frozen exact continuum matrices, then check the mass
    # anticommutation with integer numerators, not a fitted dispersion.
    sys.path.insert(0,str(HERE.parent/'continuum'))
    import continuum
    from scipy.linalg import expm
    bloch=continuum.rational_bloch_data()
    Pd=bloch['exact_projector_denominator']
    Pn=np.rint(Pd*bloch['P']).astype(np.int64)
    Rn=np.rint(6*bloch['R']).astype(np.int64)
    assert np.array_equal(Rn@Pn,Pn@Rn)
    assert np.array_equal(12*Rn@Pn,Pd*np.ones((144,144),dtype=np.int64)-72*Pn)
    for H,denominator in zip(bloch['H'],bloch['exact_gradient_matrix_denominators']):
        Hn=np.rint(denominator*H).astype(np.int64)
        assert not np.any(Rn@Hn+Hn@Rn)
    complete_controls=[]
    P=bloch['P'];R0=bloch['R']
    for h in (F(1,16),F(1,32)):
        m=F(1);n=int(1/h)
        for p in ((0.,0.,0.),(.5,.25,.125)):
            W=continuum.bloch_walk([float(h)*x for x in p])
            V=W@(math.cos(float(h*m))*np.eye(144)-1j*math.sin(float(h*m))*R0)
            K=continuum.effective_generator(p)-1j*float(m)*(R0@P)
            error=float(np.linalg.norm(np.linalg.matrix_power(V,n)@P-P@expm(float(n*h)*K),2))
            upper=float(F(continuum_bound(h,m,F(1),n)['uncapped_vector_bound_exact']))
            assert error<=upper
            complete_controls.append({'h':str(h),'mass':str(m),'p':list(p),'steps':n,
                                      'full144_component_operator_error':error,'certified_bound':upper})
    try:parameters(F(1),F(1))
    except ValueError:pass
    else:raise AssertionError('outside mass window accepted')
    result={'schema':'oph.local12.mass_coin.v1','candidate_id':MASS_ID,
            'selection':'explicit common amendment proposal, not silently substituted for the frozen base',
            'checks':{'full_native_mass_kick_columns':True,'same_source_zero_mass_restriction':True,
                      'exact_cyclic_dispersion_on_all96_modes':True,'flat_and_padding_phases_retained':True,
                      'massive_continuum_error_bounds':True,'massless_bound_recovered':True,
                      'exact_R_anticommutation_on_all26_low_modes':True,
                      'complete144_component_limit_controls':True,
                      'invalid_mass_window_rejected':True,'source_bytes_unchanged':True},
            'definition':'W_mu=S R exp(-i mu R), mu=h*m; one source-defined reflection, no independent Hamiltonian',
            'dispersion':'cos(omega)=t*cos(mu), t=1-lambda/12',
            'physical_sector_spectrum':{'constant_branch':'exp(-i mu), multiplicity1',
                'cycle_branches':'exp(+i mu), multiplicity193',
                'held_and_other_minus_branches':'-exp(+i mu), multiplicity768',
                'dispersive_branches':'exp(±i arccos(t cos(mu))),95 pairs'},
            'padding_sector_spectrum':{'invalid_address_uniform':'exp(-i mu), multiplicity32',
                'other_invalid_address_and_event_padding':'-exp(+i mu), multiplicity864',
                'scope':'All896 padded columns retained; these are decoder flags, not physical species.'},
            'continuum':'Two branches E=±sqrt(m²+|p|²/48)/tau0 and24 flat E=-m/tau0 branches onP0. All modes retained.',
            'Gauss':'The kick acts on event labels at fixed port. Its Fock lift preserves each total port occupation and commutes with every finite/integer Gauss operator.',
            'parameter_ancestry':'m is a supplied architecture parameter. It is not derived from the carrier, fitted here, or identified with any observed species mass.',
            'native_cases':rows,'refinement_bounds':bounds,'complete_state_controls':complete_controls,
            'scope':'Controlled free massive amendment with exact local source kick. Physical vacuum, interacting continuum, spin and species identification require their same-source constructions. Larger native member compilation and full field/Fock apparatus are separate obligations.',
            'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    (HERE/'MASS_COIN.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'maximum_kick_error':max(errors),
                      'refinement_bounds':[r['complete_terminal_TV_bound_exact']for r in bounds]}))


if __name__=='__main__':run()
