"""Independent CAR, native-current and complete-record controls."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
import sympy as sy
from scipy.linalg import eigh, expm
from model.openwave_current_averaging import averaging as avg
from model.openwave_meson_weak_current import decay as meson
from model.openwave_larmor_dynamics.certificate import equivalent

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
TARGETS=('openwave_current_averaging','openwave_meson_weak_current')


def pins():
    files={str(Path(__file__).resolve().relative_to(ROOT))}
    for package in TARGETS:
        result=ROOT/'model'/package/'results.json'
        files.add(str(result.relative_to(ROOT)))
        files.update(json.loads(result.read_text())['source_sha256'])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(files)}


def independent_word(mask,word):
    """Sorted-list Fock action, independent of producer's bitwise CAR action."""
    occupied=[i for i in range(mask.bit_length())if mask&(1<<i)]
    coefficient=1
    for mode,create in word[::-1]:
        position=sum(x<mode for x in occupied)
        if create:
            if mode in occupied:return None,0
            occupied.insert(position,mode)
        else:
            if mode not in occupied:return None,0
            occupied.remove(mode)
        coefficient*=(-1)**position
    return sum(1<<x for x in occupied),coefficient


def exact_four_fock_controls():
    a=[]
    for mode in range(2):
        matrix=sy.zeros(4)
        for state in range(4):
            target,sign=independent_word(state,((mode,False),))
            if sign:matrix[target,state]=sign
        a.append(matrix)
    trials=0
    for mass in (sy.Rational(3,5),sy.Integer(1),sy.Rational(17,3)):
        for momentum in (-sy.Rational(2,3),sy.Integer(0),sy.Rational(2,3)):
            for spin in (-1,1):
                epsilon=-spin*momentum;delta=spin*mass
                energy=sy.sqrt(epsilon**2+mass**2)
                u=sy.sqrt((energy+epsilon)/(2*energy))
                v=spin*sy.sqrt((energy-epsilon)/(2*energy))
                left=u*a[0]-v*a[1].T;right=u*a[1]+v*a[0].T
                source=epsilon*(left.T*left+right.T*right)+delta*(left.T*right.T+right*left)
                target=energy*(a[0].T*a[0]+a[1].T*a[1])+(epsilon-energy)*sy.eye(4)
                assert (source-target).applyfunc(sy.simplify)==sy.zeros(4)
                assert (left*right+right*left).applyfunc(sy.simplify)==sy.zeros(4)
                trials+=1
    return trials


def exact_resonant_columns():
    samples=((F(3,5),F(4,5),F(5,13),F(12,13),F(8,17),F(15,17)),
             (F(5,13),F(12,13),F(8,17),F(15,17),F(3,5),F(4,5)))
    basis=(avg.parent(),avg.receiver());expected=((F(173,2),F(1)),(F(1),F(85)))
    checked=0
    for sample in samples:
        for column,mask in enumerate(basis):
            output=defaultdict(lambda:[F(0),F(0)])
            for operations,(powers,root,coefficient),frequency in avg.words()[1]:
                assert frequency==(0,)*6
                target,sign=independent_word(mask,operations)
                if not sign:continue
                value=coefficient*sign
                for x,p in zip(sample,powers):value*=x**p
                output[target][root]+=value
            output={k:v for k,v in output.items()if any(v)}
            assert set(output)==set(basis)
            for row,target in enumerate(basis):assert output[target]==[expected[row][column],F(0)]
            checked+=1
    frequency=avg.frequency_certificate();L=F(frequency['native_current_norm_upper'])
    primitive=8*L/F(1,50000);assert primitive==16780800000
    error=F(1,2**72)*primitive*(1+2*L*F(44,35))
    assert error<F(1,2500000)
    return {'independent_exact_rational_resonant_column_checks':checked,
            'exact_resonant_matrix':[['173/2','1'],['1','85']],
            'primitive_bound':str(primitive),'first_peak_full_source_unitary_error_upper':str(error),
            'analytic_review':'The full nonzero-frequency primitive and bounded integration by parts give G*A*(1+2*G*L*t); all CAR complements and external references are retained.'}


def meson_controls():
    d=meson.construction();B,h,residual=meson.spin_parity_reduction(d)
    # A direct matrix exponential supplies a propagation algorithm independent
    # of the producer's diagonal spectral expansion.
    reduced=expm(-720j*h)@(B.conj().T@d['initial']);state=B@reduced
    row=meson.complete_records(d,720.)
    occupied=np.array([v%49!=0 for v in d['keep']])
    pair=float(np.vdot(state[occupied],state[occupied]).real)
    target=float(abs(np.vdot(d['target'],state))**2)
    assert abs(pair-row['conditional_total_lepton_pair_probability'])<2e-9
    assert abs(target-row['conditional_light_bound_spin_zero_receiver_probability'])<2e-9
    # Full sector intertwining is stronger than comparing only the prepared line.
    fullH=d['H'].toarray()
    assert np.linalg.norm(fullH@B-B@h)<5e-11
    energy=float(np.vdot(state,fullH@state).real)
    initial_energy=float(np.vdot(d['initial'],fullH@d['initial']).real)
    assert abs(energy-initial_energy)<3e-7
    for key in ('complete_position_colour_triplet_probability_each','complete_momentum_colour_triplet_probability_each'):
        assert abs(3*sum(row[key])+row['preparation_failure_probability']-1)<3e-12
    momentum=meson.momentum_readout(d).toarray()
    assert np.linalg.norm(momentum@momentum.conj().T-np.eye(888))<1e-12
    form=np.vdot(d['target'],d['V']@d['initial'])
    assert abs(form-d['C'][0,1]/np.sqrt(2))<2e-12
    assert abs(np.vdot(d['initial'],d['V']@d['initial'])+3)<1e-12
    assert abs(np.vdot(d['target'],d['V']@d['target'])-2.5)<1e-12
    # Closure follows from native quark-up number minus charged lepton number.
    for value in d['keep']:
        hi,li=divmod(int(value),49);w=d['hadron_labels'][hi][3]
        assert int(li!=0)==1-w
    full_ref=(1<<12)-1
    for mask in d['masks']:
        for family in range(3):
            count=sum(bool(mask&(1<<(weak*12+cell*6+spin*3+family)))
                      for weak in range(2)for cell in range(2)for spin in range(2))
            assert count==4
        assert (mask&full_ref).bit_count()>=11
    return {'full_sector_dimension':888,'spin_zero_even_dimension':150,
            'independent_matrix_exponential_time':720,
            'independent_pair_probability':pair,'independent_bound_receiver_probability':target,
            'independent_energy_drift_diagnostic':abs(energy-initial_energy),
            'full_sector_intertwiner_residual':residual,
            'exact_current_parent_receiver_contact_expectations':['-3','5/2'],
            'source_resonance_q_is_architecture_input':True,
            'scope':'Exact native invariant-sector reasoning and complete instruments; floating current/evolution controls. No physical kaon identity, lifetime, outward probability certificate, or forced resonance/port parameters.'}


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true');args=parser.parse_args()
    before=pins();result={'schema':'oph.independent_native_flavour_current_review.v1',
       'exact_four_fock_Bogoliubov_controls':exact_four_fock_controls(),
       'full_native_current_averaging':exact_resonant_columns(),'native_bound_meson':meson_controls(),
       'verdict':'PASS_INDEPENDENT_NATIVE_FLAVOUR_CURRENT_REVIEW'}
    assert pins()==before;result['source_sha256']=before;path=HERE/'flavour_current_review.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'])

if __name__=='__main__':main()
