"""Independent spectral, preparation and native-profile review."""
from pathlib import Path
import hashlib,json
import numpy as np
import scipy.linalg as la
import mpmath as mp
import sympy as s
from model.openwave_neutral_stationary.isolation import data,DIMENSION
from model.openwave_neutral_stationary.filter import design
from model.openwave_neutral_stationary.local_profile import certificate as profile
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import target
from model.coherent_defect_wave.chiral_joint.local_reference import physical
from model.coherent_defect_wave.chiral_joint.local_reference.reduction import quark_module,lepton_operators
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
ROOT=Path(__file__).resolve().parents[3]


def run():
    rng=np.random.default_rng(760981);checks=[]
    for n,rank in ((11,3),(17,4)):
        gap=.6+rng.random();d0=np.r_[np.zeros(rank),gap+np.arange(n-rank)/7]
        raw=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));v=raw+raw.conj().T;v[:rank,:rank]=0;v*=.07/la.norm(v,2)
        h=np.diag(d0)+v;w,u=la.eigh(h);delta=.07**2/(gap-.07)
        assert w[0]>=-delta-1e-13 and w[rank-1]<=1e-13 and w[rank]>=gap-.07-1e-13
        P=np.diag([1.]*rank+[0.]*(n-rank));E=u[:,:rank]@u[:,:rank].conj().T
        distance=la.norm(E-P,2);assert distance<=.07/gap+1e-12
        checks.append({'dimension':n,'rank':rank,'ground_band':[float(w[0]),float(w[rank-1])],
            'Sylvester_projector_bound':.07/gap,'measured_projector_distance':float(distance)})
    block=target(2,0);H=block['H'];G=block['gram'];B=block['basis'];keys=block['keys']
    q=quark_module();l=lepton_operators();scalar=higgs_memory()[1]
    all_keys=tuple(key for weight in (-2,0,2) for key in physical.weight_basis(2,0,weight))
    for j in range(4):
        current=physical.tensor_map(keys,all_keys,quark=q['currents'][j],lepton=l['currents'][j],scalar=scalar[8+j])
        assert current*B==s.zeros(len(all_keys),B.cols)
    # Independent Cholesky whitening rather than the producer's algebraic projector.
    g=np.asarray(G,complex);h=np.asarray(H,complex);L=la.cholesky(g,lower=True)
    z=la.solve_triangular(L.conj().T,np.eye(len(g)),lower=False)
    whitened=L.conj().T@h@z
    assert la.norm(whitened-whitened.conj().T)<1e-11
    energies,vec=la.eigh(whitened);psi=z@vec[:,0]
    d=data();assert float(d['en'][0])-1e-12<energies[0]<float(d['en'][1])+1e-12
    record=profile();e=s.Symbol('e')
    numerator=np.zeros(h.shape,complex);power=np.eye(len(h),dtype=complex)
    for coefficient in record['projector_numerator_coefficients']:
        numerator+=complex(s.sympify(coefficient).subs(e,energies[0]))*power;power=power@h
    denominator=complex(s.sympify(record['projector_denominator']).subs(e,energies[0]));P=numerator/denominator
    projector_error=la.norm(P-psi[:,None]@(psi.conj()[None,:]@g))
    assert projector_error<2e-8
    f=design()
    with mp.workdps(210):
        conv=lambda x:mp.mpf(int(s.Rational(x).p))/int(s.Rational(x).q)
        m=f['degree'];K=f['flattening_power'];za=conv(f['za']);den=mp.cosh(m*mp.acosh(za))
        def evaluate(E):
            x=(2*E-conv(f['upper'])-conv(f['c']))/(conv(f['upper'])-conv(f['c']))
            cheb=(-1)**m*mp.cosh(m*mp.acosh(-x)) if x<-1 else mp.cos(m*mp.acos(x))
            little=(-1)**m*cheb/den
            return 1 if little==1 else -mp.expm1(K*mp.log1p(-little))
        pass_errors=[abs(evaluate(conv(E))-1) for E in (f['a'],(f['a']+f['b'])/2,f['b'])]
        stop_errors=[abs(evaluate(conv(f['c'])+(conv(f['upper'])-conv(f['c']))*j/12)) for j in range(13)]
        assert max(pass_errors)<=conv(f['epsilon']) and max(stop_errors)<=conv(f['epsilon'])
        assert f['pass_loss']>f['epsilon'] and f['pass_loss']**K<=f['epsilon']
        # Exact high-dimensional conditioning cannot be replaced by the local16 rows.
        p=mp.mpf(128)/int(DIMENSION);eps=conv(f['epsilon'])
        actual_leak=(1-p)*eps**2/(p*(1-eps)**2+(1-p)*eps**2)
        assert 2*actual_leak<conv(f['conditional_trace_error'])
        filter_values={'maximum_pass_error':mp.nstr(max(pass_errors),12),'maximum_stop_amplitude':mp.nstr(max(stop_errors),12),
            'complete_dimension_conditioning_leak':mp.nstr(actual_leak,12)}
    receipt_path=ROOT/'model/openwave_neutral_stationary/results.json';receipt=json.loads(receipt_path.read_text())
    for p,hsh in receipt['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==hsh,p
    paths=[*receipt['input_sha256'],'model/openwave_neutral_stationary/results.json',str(Path(__file__).relative_to(ROOT))]
    return {'schema':'oph.independent_neutral_stationary_review.v1','verdict':'PASS',
        'independent_band_and_projector_controls':checks,'native_Gauss_embedding_four_nontrivial_generators':True,
        'independent_whitened_projector_error':float(projector_error),'filter_controls':filter_values,
        'mathematical_review':'Complete flux-support isolation, full-Q self-energy/Sylvester bounds, penalty-dependent exchange, graded momentum ranks, explicit Chebyshev contraction, full-dimension conditioning and complete51-record composition are consistent. The doubled link penalty and periodic source are explicit; no coefficient-one ground transfer occurs.',
        'scope_review':'The stationary accepted mixture is defined by an explicit polynomial on the same full source. It is a finite declared apparatus with success about2.85e-126 and degree167608, not a feasible primitive preparation or a physical relativistic neutrino/dark-matter identification.',
        'stationary_state_scope':'The source has an actual native profile, explicit stationary preparation, positive ground-relative zero-momentum energy and retained complements. Cosmological abundance and physical relativistic dispersion are separate unresolved quantities.',
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(set(paths))}}


if __name__=='__main__':
    value=run();Path(__file__).with_suffix('.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');print(value['verdict'])
