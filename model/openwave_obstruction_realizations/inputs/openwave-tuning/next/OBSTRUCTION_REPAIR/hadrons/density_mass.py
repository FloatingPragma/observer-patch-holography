"""Test the existing baryon-density mechanism without calibrating its scale.

U=1 is the original bridge; U=2048 is archived architecture from the
calibrated-hadron package, whose independence from target knowledge is NOT
established here. A=U adds no independent fitted coefficient. No mass target
is read until all source coefficients and operator inequalities are fixed.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json
import audit

HERE=Path(__file__).resolve().parent


def diagonal_addition(d,U):
    values=[]
    for ns in d['patterns']:
        local_colour=F(2,3)*sum(n*(3-n) for n in ns)
        B2=sum(F(n*n,9) for n in ns)
        values.append((U-1)*local_colour+U*B2)
    return values+[F()]*(d['code_dimension']-d['dimension'])


def run():
    # Every three-quark occupation pattern, not only the low band.
    patterns=[]
    for L in [2,3,4,5]:
        for labels in itertools.combinations_with_replacement(range(L),3):
            ns=[labels.count(i) for i in range(L)]
            C2=F(2,3)*sum(n*(3-n) for n in ns)
            B2=sum(F(n*n,9) for n in ns)
            electric=F(4,3)*(labels[-1]-labels[0])
            assert C2+B2>=1
            assert electric-C2/2>=0
            patterns.append(dict(L=L,occupations=ns,colour_plus_density=str(C2+B2),nonnegative_string_remainder=str(electric-C2/2)))
    rows=[]
    for U in [F(1),F(2048)]:
        # Spin half entails sum_x S_x^2>=3/4. The local endpoint trial attains
        # this, density=1, colour=0, string=0, and kinetic expectation=3/8.
        lo=U+F(51,16)
        hi=U+F(57,16)
        d=audit.source.data(2,F(1),1)
        A,G=audit.rational_matrix(d,diagonal_addition(d,U))
        for x,expected in [(lo,0),(hi,1)]:
            B=[row[:] for row in A]
            for i,g in enumerate(G): B[i][i]-=x*g
            n=audit.inertia(B)
            assert n>=expected if expected else n==0
        # Translation-independent bounds also retain all spin/weak states in
        # this three-quark domain, not merely this checked projection.
        rows.append(dict(U=str(U),A_equals_U=True,source_interval=[str(lo),str(hi)],
            archive_ancestry='bridge fixed U=1' if U==1 else 'pre-existing declared U=2048; independent target ancestry unproved'))
    unit=F('0.00051099895069');target=F('0.93827208943')
    for row in rows:
        lo,hi=map(F,row['source_interval'])
        row['mass_interval_GeV']=list(map(str,[lo*unit,hi*unit]))
        row['mass_interval_GeV_decimal']=list(map(float,[lo*unit,hi*unit]))
        row['proton_target_in_interval']=lo*unit<=target<=hi*unit
        row['distance_to_proton_GeV']=str(max(lo*unit-target,target-hi*unit))
        assert not row['proton_target_in_interval']
    out=dict(schema='oph.obstruction_repair.density_mass.v1',
        mechanism='H -> H + U sum_x B_c,x^2, and same U on the existing local colour cost; B_c signed native baryon current',
        archive='reverse-engineering-reality/model/openwave_calibrated_hadrons/PROOF.md',
        new_fitted_coefficient=False,native_particle_projector_used=False,
        vacuum_and_electron_sector_unchanged=True,full_Fock_positivity_preserved=True,
        zero_span_trial_obstruction_removed=True,
        scale_selected_by_vacuum_or_Gauss_law=False,
        ground_reference_shift_changes_excitation_gap=False,
        exact_pattern_count=len(patterns),patterns=patterns,rows=rows,
        full_physical_mass_solution=False,
        unresolved='The positive density mechanism generates a nonzero colour-singlet rest cost, but its coefficient is still a dimensionless architecture input. The vacuum and Gauss constraints admit every positive U. Both unretuned inherited values miss the proton comparison.')
    (HERE/'DENSITY_MASS_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['exact_pattern_count','zero_span_trial_obstruction_removed','scale_selected_by_vacuum_or_Gauss_law','rows']},indent=2))


if __name__=='__main__':run()
