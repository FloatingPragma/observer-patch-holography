"""Independent full-source controls of complete collision coefficient readout."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from model.source_family_compatibility import source_action_identification as c
from model.source_family_compatibility.defect_budget import exact_source
from model.source_strong_coupling.readout import Interval


def test_actual_full_source_block_keeps_original_dark_probability():
    row = c.dimer_algebra()
    assert row['full_Fock_dimension'] == 16
    assert row['original_dark_weight'] == row['original_coupled_weight'] == '1/2'
    assert len(row['complete_records']) == 64
    assert sum(r['probability'] == '0' for r in row['complete_records']) == 56
    assert not row['postselection_used']


def test_original_weyl_axes_have_exact_same_instrument_on_number_two_only():
    row = c.original_weyl_attachment()
    assert row['original_initial_masks'] == [5, 5, 9]
    assert row['independent_hopping_and_contact_coefficients_checked']
    assert not row['full_domain_bare_repair_covariance_asserted']


@pytest.mark.parametrize('coupling,contact,time', [(F(1,8),F(1),F(1,64)),
                                                (F(3,2),F(4),F(1,16)),
                                                (F(2),F(1,3),F(1,32))])
def test_independent_full_matrix_exponential_all_fine_probabilities(coupling, contact, time):
    h, _, _ = exact_source(2, ((0,1,1),), coupling, contact)
    psi = expm(-1j*float(time)*np.array(h, complex))[:,9]
    p = sum(abs(psi[i])**2 for i in (3,12))
    enclosed = c.collision_probability(coupling, contact, time)
    assert abs(p-float((enclosed.low+enclosed.high)/2)) < 2e-15
    total = 0
    for law in ('C','E'):
        for (_,raw), matrix in c.raw_matrices(2,(0,1),law).items():
            out = np.array(matrix,complex) @ psi/2
            value = float(np.vdot(out,out).real)
            expected = (1-p)/4 if raw == 0 else p/4 if raw == 8 else 0
            assert abs(value-expected) < 2e-15
            total += value
    assert abs(total-1) < 2e-15


def test_exact_source_time_pi_retains_nine_fiftieths_collision():
    h, _, _ = exact_source(2,((0,1,1),),F(3,2),F(4))
    values = (-s.Rational(1,2),s.Integer(0),s.Integer(4),s.Rational(9,2))
    # Restrict only to compute the exact source evolution of a known N=2 input.
    indices = [i for i in range(16) if i.bit_count()==2]
    block=h.extract(indices,indices); identity=s.eye(6); evolution=s.zeros(6)
    for energy,phase in zip(values,(s.I,1,1,-s.I),strict=True):
        projector=identity
        for other in values:
            if other != energy: projector=projector*(block-other*identity)/(energy-other)
        evolution += phase*projector
    psi=evolution[:,indices.index(9)]
    assert sum(s.expand(s.conjugate(psi[indices.index(i)])*psi[indices.index(i)])
               for i in (3,12)) == s.Rational(9,50)


@pytest.mark.parametrize('coupling,contact',[(F(1,8),F(1)),(F(3,8),F(2)),(F(1),F(1,2))])
def test_rational_readout_contracts_with_nonzero_probability_errors(coupling,contact):
    widths=[]
    for n in (7,8,9):
        t=F(1,2**n); error=t**6
        ps=[c.collision_probability(coupling,contact,m*t).expanded(error) for m in (1,2)]
        result=c.identify_action(t,*ps,2*coupling,2*contact)
        for key,expected in [('coupling_squared',coupling**2),('contact_squared',contact**2),
                             ('contact_over_coupling_squared',(contact/coupling)**2)]:
            assert result[key].contains(expected)
        widths.append(result['contact_over_coupling_squared'].high-result['contact_over_coupling_squared'].low)
    assert widths[2] < widths[1] < widths[0]


def test_common_clock_rescaling_is_an_exact_record_identity():
    a,b,t,scale=F(1,8),F(1),F(1,64),F(7,3)
    first=c.collision_probability(a,b,t)
    second=c.collision_probability(scale*a,scale*b,t/scale)
    assert first == second


def test_precision_and_wrong_dark_normalization_are_not_silently_accepted():
    t=F(1,64)
    with pytest.raises(ValueError):
        c.identify_action(t,Interval(F(0),F(1,10)),Interval(F(0),F(1,10)),F(1),F(2))
    with pytest.raises(ValueError):
        c.identify_action(t,Interval(F(3,4),F(4,5)),Interval(F(0),F(1,2)),F(1),F(2))
    p1,p2=[c.collision_probability(F(1,8),F(1),k*t) for k in (1,2)]
    wrong=c.identify_action(t,Interval(2*p1.low,2*p1.high),Interval(2*p2.low,2*p2.high),F(1,4),F(2))
    assert not wrong['coupling_squared'].contains(F(1,64))


def test_complete_record_reference_marginal_cannot_read_doublon_phase():
    assert len(c.branch_invariants()) == 4
    phase=s.diag(*(s.I**sum(((mask>>(2*x))&3)==3 for x in (0,1)) for mask in range(16)))
    initial=s.Matrix(16,2,lambda i,j:s.Integer(i+1)+s.I*j)
    changed_source=False
    for law in ('C','E'):
        for k in c.raw_matrices(2,(0,1),law).values():
            a=k*initial; b=k*phase*initial
            assert (a.H*a-b.H*b).applyfunc(s.expand) == s.zeros(2)
            changed_source |= (a*a.H-b*b.H).applyfunc(s.expand) != s.zeros(16)
    assert changed_source  # record/reference equality is not source-channel equality


def test_fifo_relocation_preserves_all_actual_history_blocks():
    result=c.fifo_archive_equivalence()
    assert result['allocated_history_blocks']==625
    assert result['complete_source_reference_record_equality']


def test_sampling_budget_uses_original_counts_and_two_time_confidence():
    row=c.sampling_budget(F(1,100),F(1,20))
    # ceil(log2(80))=7 >= ln(80).
    assert row['trials_per_time']==35000
    assert row['total_trials']==70000
    assert row['independent_identically_prepared_trials_required']


@pytest.mark.parametrize('arguments',[(0,F(1,10)),(F(1,10),1),(0.1,F(1,10)),(True,F(1,10))])
def test_sampling_domain_is_exact_and_nonvacuous(arguments):
    with pytest.raises(ValueError): c.sampling_budget(*arguments)


def test_configured_source_is_read_without_capacity_or_physical_promotion():
    row=c.certificate()
    assert (row['configured_c'],row['configured_U'])==('1/8','1')
    for key in ('capacity_selects_action_or_physical_clock','physical_particle_or_mass_identified',
                'physical_strong_coupling_or_QCD_alpha_selected',
                'thermodynamic_ground_selected_by_this_readout',
                'coherent_clock_pulse_tape_executed_here'):
        assert row[key] is False
