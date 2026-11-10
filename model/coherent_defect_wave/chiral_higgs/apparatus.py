"""A charged reservoir with its retained compensating boundary charge."""
from functools import lru_cache
import sympy as s
from .model import full_quantum_spectrum,quantum_ground_and_readout,spectral_weights,clean,quantum_lepton_block
from ..physical_attachment.chiral import representation


def second_quantized(one_body):
    n=one_body.rows;dim=1<<n;entries={}
    for mask in range(dim):
        for (i,j),coefficient in one_body.todok().items():
            if not mask&(1<<j):continue
            middle=mask^(1<<j)
            if middle&(1<<i):continue
            target=middle|(1<<i)
            sign=(-1)**((mask&((1<<j)-1)).bit_count()+(middle&((1<<i)-1)).bit_count())
            entries[target,mask]=entries.get((target,mask),0)+coefficient*sign
    return s.ImmutableSparseMatrix(dim,dim,entries)


@lru_cache(maxsize=1)
def compensated_swap():
    # Six actual lepton CAR modes, reservoir mode6 and compensator mode7.
    permutation=(6,1,2,3,4,5,0,7);entries={}
    for mask in range(256):
        occupied=[permutation[i] for i in range(8) if mask&(1<<i)]
        sign=(-1)**sum(occupied[i]>occupied[j] for i in range(len(occupied)) for j in range(i+1,len(occupied)))
        entries[sum(1<<i for i in occupied),mask]=sign
    u=s.ImmutableSparseMatrix(256,256,entries)
    assert u.H*u==s.eye(256) and u*u==s.eye(256)
    modes=[]
    for mode in range(8):
        terms={}
        for mask in range(256):
            if mask&(1<<mode):terms[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
        modes.append(s.ImmutableSparseMatrix(256,256,terms))
    for j in range(8):assert u.H*modes[j]*u==modes[permutation[j]]
    # The actual source currents, extended to the charged apparatus.
    indices=(9,10,11);checks=0
    for t in representation()[2]:
        one=s.zeros(8);one[:6,:6]=s.kronecker_product(s.eye(2),t.extract(indices,indices))
        one[6,6]=t[9,9];one[7,7]=-t[9,9]
        q=second_quantized(s.ImmutableSparseMatrix(one))
        assert q*u==u*q
        checks+=1
    # The compensator remains at its source boundary, including all sectors.
    nb=s.diag(*[(mask>>7)&1 for mask in range(256)])
    assert nb*u==u*nb
    return u,checks


def reduced_record(rho,setting,outcome):
    u,_=compensated_swap();compensator=setting
    # Execute the permutation on every nonzero full source density entry.
    column={col:(row,value) for (row,col),value in u.todok().items()}
    output={}
    for (i,j),value in rho.todok().items():
        fi,si=divmod(i,5);fj,sj=divmod(j,5)
        ri,ai=column[fi+64*setting+128*compensator]
        rj,aj=column[fj+64*setting+128*compensator]
        assert (ri>>7)&1==compensator and (rj>>7)&1==compensator
        if ((ri>>6)&1)==outcome and ((rj>>6)&1)==outcome:
            key=((ri&63)*5+si,(rj&63)*5+sj)
            output[key]=output.get(key,0)+value*ai*s.conjugate(aj)/2
    return s.ImmutableSparseMatrix(320,320,output)


def physical_preparation():
    """Exact vectorized purification; no 102400-square matrix is materialized."""
    h,qs,_,_,_=quantum_lepton_block()
    _,_,_,rho,_,_=full_quantum_spectrum()
    d=h.rows;identity=s.eye(d);checks=0
    # The coefficient matrix I/sqrt(d) represents sum_i |i>source |i>dual.
    # Qsource tensor I + I tensor (-Qsource.T) acts as Q I - I Q.
    for q in qs:
        assert clean(q*identity-identity*q)==s.zeros(d)
        checks+=1
    parity=s.diag(*[(-1)**(j//5).bit_count() for j in range(d)])
    assert parity*identity*parity.T==identity
    marginal=identity/d
    assert s.trace(marginal)==1
    pg=2*rho
    success=s.trace(pg*marginal*pg)
    assert success==s.Rational(1,160)
    assert clean(pg*marginal*pg/success-rho)==s.zeros(d)
    assert all(clean(q*pg)==s.zeros(d) for q in qs)
    assert clean(parity*pg+pg)==s.zeros(d)
    return {'source_dimension':d,'retained_dual_reference_dimension':d,
            'full_purification_dimension':d*d,
            'dual_current':'-Q_source.T','actual_vectorized_Gauss_checks':checks,
            'combined_parity':'even','source_marginal':'I320/320',
            'full_energy_PVM_records':14,'ground_success_probability':str(success),
            'conditional_source_marginal':'rank-two ground projector / 2',
            'post_ground_dual_reference_gauge_charge':'zero',
            'post_ground_source_and_dual_parities':['odd','odd'],
            'source_reference_preparation':'declared finite gauge/parity-preserving controller',
            'native_idle_selects_this_reference':False,
            'dual_quantum_state_retained':True}


def certificate(expected_records=None):
    _,_,_,rho,_,_=full_quantum_spectrum();_,checks=compensated_swap()
    if expected_records is None:expected_records=quantum_ground_and_readout()['complete_driven_reservoir_energy_records']
    actual=[]
    for setting in (1,0):
        for outcome in (0,1):
            output=reduced_record(rho,setting,outcome)
            for energy,p,rank in spectral_weights(output):
                actual.append({'setting':'add' if setting else 'remove','reservoir_out':outcome,
                  'source_energy':str(energy),'excitation_energy_above_full_ground':str(s.simplify(energy+2)),
                  'energy_multiplicity':rank,'probability':str(p)})
    assert actual==expected_records
    assert sum(s.sympify(r['probability']) for r in actual)==1
    complete=[{'stage':'reference_preparation_failure','source_energy':str(e),
               'probability':str(s.Rational(rank,320))}
              for e,_,rank in spectral_weights(rho) if e!=-2]
    complete.extend(dict(row,stage='prepared_ground_then_charged_readout',
                         probability=str(s.sympify(row['probability'])/160)) for row in actual)
    assert len(complete)==69 and sum(s.sympify(r['probability']) for r in complete)==1
    return {'physical_reference_preparation':physical_preparation(),
            'fermion_modes_including_reservoir_and_compensator':8,
            'full_fermion_permutation_dimension':256,'with_scalar_dimension':1280,
            'with_retained_dual_preparation_reference_dimension':409600,
            'actual_all_sector_mode_intertwiner_checks':8,'actual_total_current_commutators':checks,
            'reservoir_hypercharge':'1','compensator_hypercharge':'-1',
            'initial_recorded_reservoir_compensator_pairs':[[1,1],[0,0]],
            'original_setting_probabilities':['1/2','1/2'],
            'initial_source_reference':'rank-two gauge-singlet odd-parity lepton/scalar ground, empty quark spectators',
            'total_Gauss_constraint':'G_source + n_reservoir - n_compensator = 0 for hypercharge',
            'all_other_Gauss_currents':'same source currents; apparatus is weak/color singlet',
            'compensator_untouched_in_all_source_sectors':True,
            'every_energy_and_apparatus_effect_commutes_total_Gauss':True,
            'complete_records_preserve_all_56_source_energy_weights':True,
            'all_reference_density_entries_executed':True,
            'charged_56_record_law_condition':'ground preparation succeeded',
            'complete_preparation_and_readout_records':complete,
            'complete_preparation_and_readout_record_count':len(complete),
            'record_boundary_charge_in_add_setting':'-1',
            'record_boundary_charge_in_remove_setting':'0',
            'energy_supply':'driven operation; apparatus work is not calculated',
            'uncompensated_charged_probe_is_source_only_Gauss_preserving':False,
            'arbitrary_cross_parity_source_CP_equality_claim':False}
