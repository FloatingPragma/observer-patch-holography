"""Independent finite-carrier and actual-record audit of current protection.

This audit uses the source Hamiltonians and representation as data, independently
constructs the exterior weights, enumerates every one-family occupation and
builds all actual polar record maps. No full Fock matrix is needed.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy import linalg as la
from model.compact_defect_transport import source, transport
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.particle_common_contract.provenance import import_closure

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=Path(__file__).resolve().parents[1]


def audit():
    subsets,_,currents,_=representation()
    # Independent fundamental weights, then the actual exterior basis.
    qfund=[Fraction(-1,3)]*3+[Fraction(1),Fraction(0)]
    cfund=[Fraction(1,3),Fraction(1,3),Fraction(-2,3),Fraction(0),Fraction(0)]
    q=[sum((qfund[j] for j in row),Fraction()) for row in subsets]
    c=[sum((cfund[j] for j in row),Fraction()) for row in subsets]
    q_native=list((currents[10]+currents[11]).diagonal())
    c_native=list((Fraction(2,3)*currents[6]+Fraction(4,3)*currents[7]).diagonal())
    assert all(str(a)==str(b) for a,b in zip(q,q_native,strict=True))
    assert all(str(a)==str(b) for a,b in zip(c,c_native,strict=True))
    q3=[int(3*x) for x in q]; c3=[int(3*x) for x in c]
    assert all((a+b)%3==0 for a,b in zip(q3,c3,strict=True))
    states=[(0,0)]
    for a,b in zip(q3,c3,strict=True):
        states += [(x+a,y+b) for x,y in states]
    assert len(states)==2**15
    assert all((a+b)%3==0 and 2*a!=3 for a,b in states)
    neutral=[(a,b) for a,b in states if b%3==0]
    assert all(a%3==0 for a,b in neutral)
    one=Counter((a,b%3) for a,b in states)
    total=Counter({(0,0):1})
    for _ in range(3):
        next_total=Counter()
        for (a,b),count in total.items():
            for (x,y),multiplicity in one.items():
                next_total[a+x,(b+y)%3]+=count*multiplicity
        total=next_total
    assert sum(total.values())==2**45
    assert all((a+b)%3==0 and (b!=0 or a%3==0) for a,b in total)
    # Halving the charge convention creates an actual half-unit eigenvalue
    # on a centre-neutral occupied mode: this is not normalization-free.
    halved_half_states=sum(a==3 and b%3==0 for a,b in states)
    assert halved_half_states>0
    assert any((Fraction(a,2)+b)%3!=0 for a,b in zip(q3,c3,strict=True))

    texture_states,records,_=source.continuing_protocol()
    hs=[transport.hamiltonian(state)[0] for state in texture_states]
    frames=[]
    for h in hs:
        values,vectors=la.eigh(h,driver='evr')
        frame=vectors[:,abs(values)<2]
        assert frame.shape[1]==4
        frames.append(frame)
    n=len(hs[0]); identity=np.eye(n)
    maps={}; errors=[]
    for a,b in ((0,1),(0,2),(1,2)):
        p=frames[a]@frames[a].conj().T
        qproj=frames[b]@frames[b].conj().T
        s=qproj@p+(identity-qproj)@(identity-p)
        # Full SVD polar factor is independent of the producer's active-span
        # QR implementation; all 686 spatial-spin coordinates are retained.
        left,_,right=la.svd(s,full_matrices=False)
        u=left@right
        producer=transport.polar_transport(frames[a],frames[b])
        errors.extend([la.norm(u-producer),la.norm(u.conj().T@u-identity),
                       la.norm(u@p@u.conj().T-qproj)])
        maps[a,b]=u; maps[b,a]=u.conj().T
    charge45=np.tile(np.asarray([float(x) for x in q]),3)
    rng=np.random.default_rng(20261112)
    probe=rng.normal(size=(n,45))+1j*rng.normal(size=(n,45))
    probe/=la.norm(probe)
    branch_rows=[]
    complete_errors=[]
    for before in range(3):
        complete=np.zeros((n,n),complex)
        exact_probability=Fraction()
        for record in records:
            if record['source_state']!=before:continue
            after=record['target_state']
            probability=Fraction(record['probability'])
            exact_probability+=probability
            u=identity if before==after else maps[before,after]
            if not record['accepted']:
                assert before==after
            # For the actual U, test the tensor action without building its
            # 30870-dimensional Kronecker matrix. The exact identity follows
            # for every entry from U acting only on the spatial factor.
            err=la.norm(u@(probe*charge45)-(u@probe)*charge45)
            errors.append(err)
            complete+=float(probability)*u.conj().T@u
            branch_rows.append({'source':before,'target':after,
                'proposal':record['proposal_label'],'accepted':record['accepted'],
                'probability':str(probability),'charge_commutator_below_1e_minus_11':bool(err<1e-11)})
        assert exact_probability==1
        complete_errors.append(la.norm(complete-identity))
    assert len(branch_rows)==12 and max(errors+complete_errors)<1e-10
    # Mixing a neutral and charged internal coordinate destroys the current
    # identity, despite retaining a perfectly valid spatial unitary.
    neutral_index=q.index(Fraction(0)); charged_index=q.index(Fraction(1))
    bad=np.eye(15); bad[[neutral_index,charged_index]]=bad[[charged_index,neutral_index]]
    bad_commutator=la.norm(bad@np.diag([float(x) for x in q])-np.diag([float(x) for x in q])@bad)
    assert bad_commutator>1

    inputs=['model/openwave_charge_protection/protection.py',
        'model/openwave_charge_protection/PROOF.md',
        'model/openwave_charge_protection/results.json',
        'model/compact_defect_transport/source.py',
        'model/compact_defect_transport/transport.py',
        'model/compact_defect_transport/PROOF.md',
        'model/coherent_defect_wave/physical_attachment/chiral.py',
        'model/common_source/chiral_candidate/action.py',
        'model/openwave_charge_protection/audits/charge_review.py']
    inputs=sorted(import_closure(set(inputs)|set(SOURCE_FILES)))
    return {'schema':'oph.current_protection_independent_review.v1','verdict':'PASS',
        'reviewed_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
        'one_family_all_occupation_states':len(states),
        'one_family_centre_neutral_states':len(neutral),
        'three_family_exact_occupation_multiplicity':sum(total.values()),
        'three_family_distinct_charge_centre_pairs':len(total),
        'charge_times_three_range':[min(a for a,b in total),max(a for a,b in total)],
        'colour_centre_identity_all_occupations':True,
        'centre_neutral_implies_integer_charge':True,
        'centre_neutral_alone_proves_colour_singlet':False,
        'centre_neutral_is_necessary_for_colour_singlet':True,
        'colour_singlet_charge_integer_follows':True,
        'actual_source_spatial_spin_dimension':n,'native_internal_family_dimension':45,
        'full_spatial_polar_oracle_compared_with_producer_pairs':3,
        'all_actual_records':branch_rows,
        'every_actual_source_texture_has_complete_instrument':True,
        'maximum_polar_transport_and_current_residual':float(max(errors)),
        'maximum_complete_instrument_residual':float(max(complete_errors)),
        'negative_controls':{'halved_charge_has_centre_neutral_half_unit_occupation_states':halved_half_states,
            'halved_charge_breaks_phase_identity':True,
            'internal_neutral_charge_swap_commutator_norm':float(bad_commutator)},
        'exact_extension':'Exterior products preserve the one-particle phase identity; Gamma(U tensor I45) commutes with dGamma(I tensor Q) on every occupation sector. No finite occupation cutoff or selected-record postselection is used.',
        'scope':{'fixed_declared_Q_normalization':True,'physical_electric_charge_identified':False,
            'local_Gauss_transport_proved_by_uniform_commutation':False,
            'winding_current_identity_required':False,
            'full_interacting_repair_dynamics_covered':False},
        'local_Gauss_shortcut':'Use the same native Q as the charge representation in an explicitly gauge-covariant link transport and retain its link electric-field response. The native current and exact colour-centre relation remove the need for Q=f(winding), but uniform global commutation does not prove preservation of spatial Gauss constraints.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=PACKAGE/'charge_review.json')
    args=parser.parse_args()
    value=audit();args.output.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:value[k] for k in ('verdict','one_family_all_occupation_states','three_family_exact_occupation_multiplicity','maximum_polar_transport_and_current_residual')}))
