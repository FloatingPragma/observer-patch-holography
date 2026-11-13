from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from .certificate import measure,pointer
from .source import construction,exact_parity_matrix,parity_basis

def test_full_contact_lines_and_source_domain():
    n=measure(0);p=measure(1)
    assert [n['complete_charge_sector_dimension'],p['complete_charge_sector_dimension']]==[1400,360]
    assert [n['exact_even_parity_dimension'],p['exact_even_parity_dimension']]==[199,51]
    assert len([x for x in n['complete_spectral_clusters']if x['conditional_probability_diagnostic']>.1])==3
    assert max(x['conditional_probability_diagnostic']for x in p['complete_spectral_clusters'])>.9999
    assert F(n['reference_expectation_shift_GeV'])<0<F(p['reference_expectation_shift_GeV'])

def test_complete_register_against_direct_controlled_unitaries():
    d=construction(1);exact,_=exact_parity_matrix(d);H=np.array(exact.evalf(18),float)
    v=np.array(measure(1)['complete_records']['actual_preparation_even_basis_integer_coefficients'],float);v/=np.linalg.norm(v)
    M=8;a=1792.;W=256.;U=expm(-2j*np.pi*(H-a*np.eye(len(H)))/W);powers=[np.eye(len(H),dtype=complex)]
    for n in range(1,M):powers.append(U@powers[-1])
    Kraus=[sum(np.exp(2j*np.pi*k*n/M)*powers[n]for n in range(M))/M for k in range(M)]
    assert np.linalg.norm(sum(K.conj().T@K for K in Kraus)-np.eye(len(H)))<1e-11
    assert max(np.linalg.norm(K@H-H@K)for K in Kraus)<1e-9
    actual=np.array([np.linalg.norm(K@v)**2 for K in Kraus]);pred=np.array(pointer(measure(1),M)['all_conditional_pointer_probabilities'])
    assert np.sum(abs(actual-pred))/2<float(F(pointer(measure(1),M)['exact_instrument_total_variation_error_upper']))

def test_resolution_refinement_retains_all_outcomes():
    a=pointer(measure(0),1024);b=pointer(measure(0),16384)
    assert len(a['all_conditional_pointer_probabilities'])==1024 and len(b['all_conditional_pointer_probabilities'])==16384
    assert F(a['energy_bin_width_GeV'])==16*F(b['energy_bin_width_GeV'])
    for x in(a,b):
        assert abs(sum(x['all_conditional_pointer_probabilities'])-1)<1e-9
        assert abs(sum(x['all_original_input_pointer_probabilities'])+float(F(x['preparation_failure_probability']))-1)<1e-9
