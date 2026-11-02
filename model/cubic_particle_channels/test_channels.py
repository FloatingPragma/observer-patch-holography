from copy import deepcopy
from fractions import Fraction as F
import json
import numpy as np
import pytest
from . import group as g, channels as c, packet
from model.cube_eta_protection import cube
from model.field_symmetry import symmetry


def complex_spinor(q):
    real,imag=g.spinor_parts(q)
    value=lambda x:float(x.a)+float(x.b)*np.sqrt(2.)
    return np.array([[value(real[i][j])+1j*value(imag[i][j]) for j in range(2)] for i in range(2)])


def test_actual_source_rotation_and_negative_control():
    h=cube.GEOMETRY.single_h().toarray()
    bad=[]
    for q in g.elements():
        p=g.permutation(q);P=np.zeros((8,8));P[p,np.arange(8)]=1
        D=np.kron(P,complex_spinor(q))
        assert np.max(np.abs(D@h-h@D))<1e-15
        bad.append(np.max(np.abs(np.kron(P,np.eye(2))@h-h@np.kron(P,np.eye(2)))))
        S,W=g.hidden(q)
        frame=np.zeros((16,16),complex)
        for i,site in enumerate(cube.GEOMETRY.sites):frame[2*i:2*i+2,2*i:2*i+2]=symmetry.onsite_frame(site)
        assert np.max(np.abs(frame@D@frame.conj().T-np.kron(S,complex_spinor(W))))<1e-15
    assert max(bad)>.1


def test_full_fock_central_parity():
    q=g.negative(g.ONE)
    assert g.permutation(q)==tuple(range(8))
    assert np.array_equal(complex_spinor(q),-np.eye(2))
    # Exterior powers of -I on every occupation basis vector, not just N8.
    for mask in range(65536):
        sign=1
        for mode in range(16):
            if mask>>mode&1:sign=-sign
        assert sign==(-1)**mask.bit_count()


def test_complete_group_and_characters():
    result=g.certificate()
    assert result['binary_group_order']==48
    assert result['proper_spatial_group_order']==24
    assert result['hidden_spin_vector']=={'A2':1,'E':1}
    assert result['charged_one_particle']=={'G1':2,'G2':2,'H':2}
    assert result['SU2_restriction_by_twice_spin']['2']=={'T1':1}
    assert result['SU2_restriction_by_twice_spin']['7']=={'G1':1,'G2':1,'H':1}


def test_independent_charged_filter_eigenvalues():
    for name,(real,imag) in c.exact_filters().items():
        Q=(real+1j*imag)/48
        eig=np.linalg.eigvalsh(Q)
        assert np.max(np.minimum(abs(eig),abs(eig-1)))<1e-14
        expected={'G1':4,'G2':4,'H':8}[name]
        assert np.count_nonzero(eig>.5)==expected
        for q in g.elements():
            S,W=g.hidden(q);D=np.kron(S,complex_spinor(W))
            assert np.max(np.abs(Q@D-D@Q))<1e-14


def test_ground_and_commutants():
    s=c.structure()
    assert s['ground_invariance']['all_ground_block_columns_fixed_exactly']
    assert s['ground_invariance']['ground_block_dimension']==32
    assert s['signed_commutant_dimension']==2
    assert s['ordinary_commutant_dimension']==4


def test_every_retained_line_and_trace_weight():
    out=c.spectrum()
    assert len(out['charged'])==228
    for rows in [out['charged']]+[v['rows'] for v in out['neutral']]:
        for row in rows:
            assert sum(map(F,row['complete_probe_trace_fractions'].values()))==1
    assert out['neutral'][0]['rows'][0]['cubic_irreps']=={'A1':1,'E':1}
    assert out['neutral'][1]['rows'][0]['cubic_irreps']=={'T1':1}
    assert all(len(x['rows'])==135 for x in out['neutral'])


@pytest.mark.parametrize('mutation',['rank','weight','sign','residue','dark'])
def test_malformed_spectral_premise_rejected(mutation):
    charged=json.loads((c.ROOT/packet.PRIORS[0]).read_text())
    neutral=json.loads((c.ROOT/packet.PRIORS[1]).read_text())
    if mutation=='rank':charged['spectrum']['rows'][0]['visible_spatial_multiplicity']=3
    if mutation=='weight':charged['spectrum']['rows'][0]['local_branch_weight_interval']=['0','0']
    if mutation=='sign':charged['spectrum']['rows'][0]['normalized_spatial_overlap_intervals'][1]=[-1,1]
    if mutation=='residue':neutral['spectrum']['responses'][0]['rows'][0]['spatial_residue_intervals'][0][1]=[1,2]
    if mutation=='dark':neutral['spectrum']['responses'][0]['rows'][-1]['exactly_dark']=True
    with pytest.raises(ValueError):c.classify(charged,neutral)


@pytest.mark.parametrize('value',[.5,True,'1/2'])
def test_nonexact_coefficient_rejected(value):
    with pytest.raises(ValueError):g.Q(value)


def test_receipt_cannot_promote_physical_scope():
    value=packet._produce()
    altered=deepcopy(value);altered['theorem']['physical_particle_species_identified']=True
    with pytest.raises(ValueError):packet.validate(altered)
    altered=deepcopy(value);altered['spectrum']['charged'][0]['energy_interval']=['1','2']
    with pytest.raises(ValueError):packet.validate(altered)
    altered=deepcopy(value);altered['dependency_custody_pass']=not value['dependency_custody_pass']
    with pytest.raises(ValueError):packet.validate(altered)
