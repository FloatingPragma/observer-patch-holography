"""Same-kernel Majorana contraction, original failures and phase envelopes."""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sympy as s
from model.common_source import calibrated_matter as matter
from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:raise ValueError(message)


def clean(a):return a.applyfunc(s.expand)
def zero(a):return all(x==0 or s.expand(x)==0 or s.simplify(x)==0 for x in a)
def sparse(a):return dict(shape=list(a.shape),entries=[[i,j,str(v)]for(i,j),v in sorted(s.ImmutableSparseMatrix(a).todok().items())])
def display(x):return dict(exact=str(x),numeric=float(s.re(s.N(x,20))))


@lru_cache(maxsize=1)
def source():
    masses=matter.neutrino_mass_coordinates();u=matter.pmns_matrix();m=matter.neutrino_family_kernel()
    require(masses[0]==0 and 0<masses[1]<masses[2],'selected normal rank-two domain required')
    require(zero(u.H*u-s.eye(3)),'actual PMNS is not unitary')
    neutral=[3*row['index']+family for row in action.ordering()['all_left_components']
             if row['internal_type']=='L'and row['electric_Q']=='0'for family in range(3)]
    require(len(neutral)==3,'three actual neutral fields required')
    k=matter.mass_kernel(include_neutrinos=True)
    require(zero(k.extract(neutral,neutral)-m),'full K45 neutral extraction differs')
    outside=[i for i in range(45)if i not in neutral]
    require(zero(k.extract(outside,neutral))and zero(k.extract(neutral,outside)),'neutral source block leaks')
    delta=matter.pairing_kernel(include_neutrinos=True)
    indices=neutral+[90+45+i for i in neutral]
    # The actual full p=0 BdG generator is [[0,Delta90],[Delta90*,0]].
    # H means conjugate transpose in the second block.
    h=s.zeros(6);h[:3,3:]=delta.extract(neutral,[45+i for i in neutral]);h[3:,:3]=h[:3,3:].H
    require(zero(h[:3,3:]-m.conjugate())and zero(h[3:,:3]-m),'native pairing orientation differs')
    # Weak charged-current contraction uses the actual mass-frame current.
    w=matter.mass_basis_change(include_neutrinos=True);weak=w.H*matter.weak_raising_current()*w
    leptons=[3*row['index']+family for row in action.ordering()['all_left_components']
             if row['internal_type']=='L'and row['electric_Q']=='-1'for family in range(3)]
    current=weak.extract(neutral,leptons)
    require(zero(current-u.H),'actual mass-frame weak raising block differs')
    return dict(masses=masses,u=u,m=m,h=h,indices=indices,current=current)


@lru_cache(maxsize=1)
def instrument():
    data=source();masses,u,m,h=(data[k]for k in ('masses','u','m','h'));scale=masses[2]
    # Takagi frame proves all identities on arbitrary states and idle references.
    frame=s.diag(u,u.conjugate());c=clean(h/scale)
    q=s.diag(*(s.sqrt(1-(x/scale)**2)for x in masses))
    defect=clean(frame*s.diag(q,q)*frame.H)
    require(zero(c.H-c),'mass contraction is not Hermitian')
    # Certify in the exact Takagi frame rather than numerically testing a square root.
    cframe=s.zeros(6);cframe[:3,3:]=s.diag(*(x/scale for x in masses));cframe[3:,:3]=cframe[:3,3:]
    require(zero(h*frame-frame*cframe*scale),'full native neutral intertwiner failed')
    require(zero(cframe*cframe+s.diag(q,q)**2-s.eye(6)),'Julia defect identity failed')
    require(zero(cframe*s.diag(q,q)-s.diag(q,q)*cframe),'Julia defect commutator failed')
    initial=s.eye(6)[:,0];success=clean(c*initial)
    probabilities=[s.simplify(s.conjugate(x)*x)for x in success]
    fail=s.simplify(1-sum(probabilities))
    mbeta2=s.simplify(sum(abs(u[0,i])**2*masses[i]**2 for i in range(3)))
    require(s.simplify(sum(probabilities)-mbeta2/scale**2)==0,'weak-current weighted mass norm differs')
    require(s.simplify(probabilities[3]-abs(m[0,0])**2/scale**2)==0,'electron LNV record differs')
    return dict(full_Nambu_dimension=180,neutral_invariant_dimension=6,spectator_dimension=174,
        embedded_indices=data['indices'],contraction_norm_exact='1',normalization_GeV=str(scale),
        neutral_contraction=sparse(c),neutral_failure_operator=sparse(defect),Takagi_frame=sparse(frame),
        complete_record_count=8,
        record_definition=['six success records K_z=|z><z| C on the neutral block',
                           'one failure record D=sqrt(I-C^2) on that block',
                           'one identity record on the full 174-dimensional spectator complement'],
        Julia_unitary='[[C,D],[D,-C]] on neutral block tensor flag, identity on spectators',
        full_effect_identity_exact=True,arbitrary_input_reference_preserved=True,
        electron_input_original_probabilities=[dict(record=str(i),probability=display(p))for i,p in enumerate(probabilities)]+
            [dict(record='failure',probability=display(fail)),dict(record='spectator',probability=display(s.Integer(0)))],
        electron_success_sum=display(mbeta2/scale**2),
        inverse=dict(m_beta_squared_GeV2='m3^2 * sum_all_six_success_probabilities',
                     m_ee_squared_GeV2='m3^2 * P(success, conjugate-electron-family)'),
        literal_native_pulse_compilation_claimed=False,nuclear_or_beta_decay_rate_claimed=False)


def phase_controls():
    data=source();u,masses,m=(data[k]for k in ('u','masses','m'))
    inputs=matter.neutrino_parameters()['oscillation_inputs']
    delta=s.Rational(inputs['delta_degrees'])*s.pi/180
    a=s.simplify(masses[1]*abs(u[0,1])**2);b=s.simplify(masses[2]*abs(u[0,2])**2)
    rows=[]
    for name,z,target in [('maximum',s.exp(s.I*delta),(a+b)**2),
                          ('minimum',s.I*s.exp(s.I*delta),(a-b)**2)]:
        rotated=u*s.diag(1,1,z)
        modified=clean(rotated.conjugate()*s.diag(*masses)*rotated.H)
        electron=s.simplify(abs(modified[0,0])**2)
        require(s.simplify(electron-target)==0,'attained Majorana envelope failed')
        # The unit-modulus column phase leaves U diag(m^2) U* exactly invariant.
        require(s.simplify(z*s.conjugate(z))==1,'Majorana control is not a phase')
        rows.append(dict(name=name,third_column_phase=str(z),Mnu_ee_GeV=str(modified[0,0]),
                         electron_original_success_probability=display(electron/masses[2]**2),
                         beta_mass_and_oscillation_kernel_unchanged=True))
    return rows


def mass_observables():
    data=source();u,masses,m=(data[k]for k in ('u','masses','m'))
    a=s.simplify(masses[1]*abs(u[0,1])**2);b=s.simplify(masses[2]*abs(u[0,2])**2)
    require(a>b>0,'phase-envelope ordering changed')
    mbeta2=s.simplify(sum(abs(u[0,i])**2*masses[i]**2 for i in range(3)))
    mee2=s.simplify(abs(m[0,0])**2)
    return dict(neutrino_masses_GeV=[display(x)for x in masses],mass_sum_GeV=display(sum(masses)),
        m_beta_squared_GeV2=display(mbeta2),m_beta_GeV=display(s.sqrt(mbeta2)),
        m_ee_at_packet_phases_GeV=display(s.sqrt(mee2)),
        m_ee_phase_envelope_GeV=[display(a-b),display(a+b)],
        envelope_terms_GeV=dict(a=display(a),b=display(b)),
        phase_law='m_ee^2=a^2+b^2+2ab cos(alpha31-alpha21-2 delta)',
        phase_envelope_attained=dict(maximum='alpha31-alpha21=2 delta modulo 2pi',
                                    minimum='alpha31-alpha21=2 delta+pi modulo 2pi'),
        oscillation_and_beta_mass_invariant_under_Majorana_phase=True,
        one_remaining_relative_Majorana_phase=True,
        all_masses_are_conditional_calibration=True)


def bindings():
    names={p.relative_to(ROOT).as_posix()for p in HERE.iterdir()if p.suffix in('.py','.md')}
    names.update(matter.source_bindings())
    names.update(('model/source_flavour_identifiability/PROOF.md','model/conditional_particle_spectrum/PROOF.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(import_closure(names))}


def certificate():
    pins=bindings();out=dict(schema='oph.neutrino_mass_observables.v1',source_bindings=pins,
        all_checks_pass=True,mass_observables=mass_observables(),complete_instrument=instrument(),phase_sensitive_controls=phase_controls(),
        calibration_ancestry=matter.neutrino_parameters(),physical_promotions=[])
    require(bindings()==pins,'inputs changed during calculation');return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
    r=certificate();text=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if a.write:(HERE/'results.json').write_text(text)
    if a.check:require((HERE/'results.json').read_text()==text,'neutrino receipt differs')
    print(json.dumps(r['mass_observables']))


if __name__=='__main__':main()
