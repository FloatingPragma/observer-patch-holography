"""Recorded probe settings with the full native failure and parity reference."""
import sympy as s
from ..neutral_composite.model import instrument_certificate,local_instrument
from ..higgs_composites import model as car
from ...chiral_higgs.model import PAULI


def certificate():
    inherited=instrument_certificate()
    assert len(inherited['branches'])==2
    assert all(row['two_sided_unitarity_checked'] and
        row['total_reservoir_source_parity_and_lepton_number_conserved']
        for row in inherited['branches'])
    for spin in (0,1):
        for axis in range(3):
            expected=car.add(*(car.scale(car.neutral_composite(beta),-PAULI[axis][spin,beta]/2)
                               for beta in (0,1)))
            assert car.commutator(car.spin_current(axis),car.neutral_composite(spin))==expected
    rows=[]
    for kind in ('remove','add'):
        for spin in (0,1):
            j,success,failure,projectors=local_instrument(kind,spin)
            assert car.add(car.multiply(car.adjoint(success),success),
                           car.multiply(car.adjoint(failure),failure))==car.monomial()
            assert all(not car.commutator(car.gauge_current(a),op)
                       for a in range(12) for op in (j,failure))
            # The selected spin component exchanges the corresponding S_z
            # with the retained fermionic reference.
            charge=1 if kind=='remove' else -1
            twice_spin=charge*(1 if spin==0 else -1)
            assert car.commutator(car.spin_current(2),j)==car.scale(j,-s.Rational(twice_spin,2))
            rows.append({'kind':kind,'probe_spin':spin,'setting_probability':'1/4',
                'success_operator':'J_kind,spin/2',
                'failure_operator':'P0+sqrt(3)/2 P1+sqrt(2)/2 P2',
                'full_effect_eigenvalues':['0','1','2'],
                'full_native_CAR_Kraus_identity':True,'all_twelve_Gauss_currents_preserved':True,
                'reference_success_raw_L':charge,'reference_success_twice_Sz':twice_spin,
                'reference_success_parity':'odd','reference_failure_parity':'even',
                'reference_gauge_charge':'0','retained_failure_has_no_energy_postselection':True})
    return {'source_preparation':'selected normalized global physical ground projector P_ground/2',
        'ground_selection_is_declared_preparation_law':True,
        'arbitrary_input_ground_filter_success_probability_supplied':False,
        'settings':rows,'local_parity_preserving_unitary_dilation_dimension':40,
        'two_sided_unitarity_and_number_parity_conservation_replayed':True,
        'spin_partner_control_obtained_by_exact_native_spin_rotation':True,
        'reference_and_setting_registers_retained':True,
        'success_energy_measurement':'all nine spectral projectors in the full L0 or L2 maximal physical module',
        'outside_target_projector_retained_with_zero_probability':True,
        'complete_record_count':44,'controller_is_declared_finite_apparatus':True,
        'controller_is_idle_Hamiltonian':False,'conditional_success_renormalization_applied':False}
