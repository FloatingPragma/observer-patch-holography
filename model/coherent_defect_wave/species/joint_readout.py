"""Complete translation and phase records on the actual two-cell CAR clock."""
import sympy as s

from .fock import MODES,mode
from .spatial_fock import PairFockState,exchange,gamma,onsite,positive_basis


def combine(left,right,a=1,b=1):
    if left.norm_denominator!=right.norm_denominator:
        raise ValueError('one common normalization required')
    out={key:a*value for key,value in left.amplitudes}
    for key,value in right.amplitudes:
        out[key]=s.expand(out.get(key,0)+b*value)
    return PairFockState(out,left.norm_denominator)


def translate(state):
    """Swap complete cells, with the CAR sign for swapping mode blocks."""
    out={}
    for (links,occupied),value in state.amplitudes:
        left=sum(i<MODES for i in occupied);right=len(occupied)-left
        moved=tuple(sorted((i+MODES if i<MODES else i-MODES) for i in occupied))
        out[tuple(reversed(links)),moved]=(-1)**(left*right)*value
    return PairFockState(out,state.norm_denominator)


def tick(state,inverse=False):
    """U_pair=Hopping U_right U_left, with the same shared-link onsite U_F."""
    if inverse:
        return onsite(onsite(exchange(state,inverse=True),1,inverse=True),0,inverse=True)
    return exchange(onsite(onsite(state,0),1))


def record_branches(state,evolved=None):
    """P_eta=(I+eta*T)/2 then K_zeta=(I+zeta*U_pair)/2.

    The proved commutation permits one actual U evaluation. Every detector
    alternative is retained; zero-probability branches are included.
    """
    evolved=tick(state) if evolved is None else evolved
    shifted=translate(state);shifted_evolved=translate(evolved)
    output=[]
    for eta in (1,-1):
        prepared=combine(state,shifted,s.Rational(1,2),s.Rational(eta,2))
        advanced=combine(evolved,shifted_evolved,s.Rational(1,2),s.Rational(eta,2))
        for zeta in (1,-1):
            branch=combine(prepared,advanced,s.Rational(1,2),s.Rational(zeta,2))
            output.append(({'translation_eigenvalue':eta,
                            'dimensionless_momentum':'0' if eta==1 else 'pi',
                            'phase_bin':'0' if zeta==1 else 'pi'},branch))
    return tuple(output)


def constraint_preparation(state):
    """All four local Gauss-filter outcomes, without renormalizing branches."""
    out=[]
    for left in (1,-1):
        first=combine(state,gamma(state,0),s.Rational(1,2),s.Rational(left,2))
        for right in (1,-1):
            branch=combine(first,gamma(first,1),s.Rational(1,2),s.Rational(right,2))
            out.append(((left,right),branch))
    return tuple(out)


def certificate():
    index=mode(7,3,0,0)
    vacuum=positive_basis((1,1),())
    source=positive_basis((1,1),(index,))
    assert vacuum.norm_squared()==source.norm_squared()==1
    assert tick(vacuum)==vacuum and translate(vacuum)==vacuum
    assert translate(translate(source))==source
    evolved=tick(source)
    assert evolved.norm_squared()==1
    moved_flux_weight=s.expand(sum((s.conjugate(v)*v for (links,_),v in evolved.amplitudes
                                    if any(g not in (1,6553635) for g in links)),s.Integer(0))/evolved.norm_denominator)
    assert moved_flux_weight>0
    assert tick(translate(source))==translate(evolved)
    assert all(gamma(source,c)==source and gamma(evolved,c)==evolved for c in (0,1))
    preparation=constraint_preparation(PairFockState({((1,1),(index,)):1}))
    assert all(branch.norm_squared()==s.Rational(1,4) for _,branch in preparation)
    assert sum(branch.norm_squared() for _,branch in preparation)==1
    plus_branch=dict(preparation)[1,1]
    assert plus_branch.as_dict()=={key:value/4 for key,value in source.amplitudes}
    assert plus_branch.norm_denominator==1 and source.norm_denominator==4
    rows=[]
    for record,branch in record_branches(source,evolved):
        eigen=record['translation_eigenvalue']
        assert translate(branch)==PairFockState({k:eigen*v for k,v in branch.amplitudes},branch.norm_denominator)
        assert all(gamma(branch,c)==branch for c in (0,1))
        assert all(len(occupation)==1 for (_,occupation),_ in branch.amplitudes)
        rows.append({**record,'probability':str(branch.norm_squared()),'support':len(branch.amplitudes)})
    assert sum(s.Rational(row['probability']) for row in rows)==1
    assert s.Rational(3,4)+sum(s.Rational(row['probability'])/4 for row in rows)==1
    assert all(sum(s.Rational(row['probability']) for row in rows if row['translation_eigenvalue']==eta)==s.Rational(1,2)
               for eta in (1,-1))
    assert len({row['probability'] for row in rows})>1
    # Explicit fermionic block-swap sign, beyond the one-wave execution domain.
    pair=PairFockState({((1,3),(index,MODES+index)):1})
    expected=PairFockState({((3,1),(index,MODES+index)):-1})
    assert translate(pair)==expected and translate(translate(pair))==pair
    return {'clock':'U_pair=Hopping U_right U_left',
            'translation':'whole-cell swap of link words and CAR mode blocks',
            'source_mode':[7,3,0,0],'source_occupation':1,
            'source':'positive-Gauss localized one-wave state over the identical local link pairs',
            'reference':'positive-Gauss empty-wave state over the identical local link pairs',
            'reference_fixed_by_clock_and_translation':True,
            'source_orthogonal_to_reference_by_wave_number':True,
            'clock_translation_commutation_executed':True,
            'fermionic_translation_sign_tested':True,
            'complete_joint_records':rows,
            'joint_records_conditional_on_both_Gauss_preparation_success':True,
            'Gauss_preparation_outcome_probabilities':[
                {'local_eigenvalues':list(labels),'probability':str(branch.norm_squared())}
                for labels,branch in preparation],
            'Gauss_preparation_success_probability':'1/4',
            'joint_record_probabilities_per_original_ambient_trial':[
                str(s.Rational(row['probability'])/4) for row in rows],
            'retained_preparation_failure_probability':'3/4',
            'probability_at_least_one_flux_leaves_its_initial_face':str(moved_flux_weight),
            'both_constraints_preserved_in_every_branch':True,
            'phase_kernel':'(1 + zeta*cos(phi))/2 on a U_pair eigenphase phi',
            'records_are_exact_energy_eigenprojectors':False,
            'physical_momentum_units_attached':False,
            'three_dimensional_dispersion_established':False,
            'rest_mass_identified':False}
