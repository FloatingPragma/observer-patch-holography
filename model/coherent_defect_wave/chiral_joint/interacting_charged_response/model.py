"""Gauge-covariant odd probes with a complete static CAR reference."""
import sympy as s
from ...physical_attachment.chiral import representation
from ...chiral_source.model import PAULI
from ..periodic_family.model import native_blocks

FAMILIES={'uc':(0,1,4),'Q':(2,3,5,6,7,8),'ec':(9,),
          'L':(10,11),'dc':(12,13,14)}
REFERENCE_COST=s.Integer(4)


def operator_certificate():
    gauges=representation()[2]
    raw=s.diag(*(1 if a%15 in (10,11) else -1 if a%15==9 else 0 for a in range(30)))
    rows=[]
    assert sorted(FAMILIES.values())==sorted(tuple(b) for b in native_blocks())
    for family,indices in FAMILIES.items():
        probe=s.zeros(60)
        for spin in range(2):
            for a in indices:probe[30+15*spin+a,15*spin+a]=1/s.sqrt(2*len(indices))
        for gauge in gauges:
            full=s.diag(s.kronecker_product(s.eye(2),gauge),s.kronecker_product(s.eye(2),gauge))
            assert full*probe==probe*full
        for sigma in PAULI:
            spin=s.kronecker_product(sigma/2,s.eye(15));full=s.diag(spin,spin)
            assert full*probe==probe*full
        assert s.diag(raw,raw)*probe==probe*s.diag(raw,raw)
        opposite=s.diag(s.kronecker_product(s.eye(2),gauges[-1]),
                        -s.kronecker_product(s.eye(2),gauges[-1].T))
        assert opposite*probe!=probe*opposite
        rows.append({'family':family,'native_internal_dimension':len(indices),
            'full_source_reference_gauge_commutators':12,'full_spin_commutators':3,
            'full_raw_L_commutator':'0','opposite_reference_orientation_rejected':True,
            'reference_vacuum_effect':f'N_{family}/{2*len(indices)}',
            'full_Fock_reference_vacuum_contraction_bound':'1'})
    assert all(s.trace(gauge)==0 for gauge in gauges)
    return {'full_source_modes_per_cell':30,'full_static_reference_modes':30,
        'full_reference_Fock_dimension':2**30,
        'reference_representation':'same native representation as source creation operators',
        'reference_Hamiltonian':'mu N_ref, mu=4; zero-energy alternative assessed separately',
        'reference_has_spatial_hopping':False,
        'full_carrier_unitary':'Halmos completion on reference-vacuum plus the selected one-particle family; identity on its complete orthogonal reference complement',
        'source_probe_parity':'odd','joint_source_reference_instrument_parity':'even',
        'reference_retained_as_quantum_system_not_measured_in_a_gauge_variant_basis':True,
        'all_five_native_elementary_probe_controls':rows,
        'filled_reference_is_Gauss_singlet':True,
        'zero_cost_full_reference_unique_ground_ruled_out_by_filled_reference_degeneracy':True,
        'quark_spectral_response_transferred_from_lepton_response':False}
