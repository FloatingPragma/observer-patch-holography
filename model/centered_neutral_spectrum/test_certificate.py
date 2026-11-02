import copy,csv,io
from fractions import Fraction as F
import numpy as np
import pytest
from scipy.sparse.linalg import eigsh
from . import certificate,packet,source
from ..centered_particle_spectrum import certificate as old_certificate,packet as old_packet


def test_full_gap_and_exact_dark_counts(result):
    r=result['spectrum'];gap=list(map(F,r['global']['neutral_gap_interval']))
    assert F(15,1000)<gap[0]<gap[1]<F(16,1000)
    assert r['global']['gap_certifies_full_65536_state_source'] is True
    for p in r['responses']:
        assert p['positive_lines']==135 and p['exact_dark_lines']==1
        assert sum(x['strictly_positive'] for x in p['rows'])==135
        assert p['fixed_component_visible_dimension']==357
        for x in p['rows']:
            if x['strictly_positive']:assert F(x['weight_interval'][0])>0
            else:assert all(F(v)==0 for v in x['weight_interval'])
    assert r['all_spin_components_and_density_reachable_N8_dimension']==1428
    assert r['excited_N8_dimension_not_reached_by_these_probes']==11441


def test_independent_unreduced_sparse_gap_agrees(result):
    h=source.parent.centered(4,4).astype(float)/16
    values=np.sort(eigsh(h,k=2,which='SA',tol=1e-12,v0=np.arange(1,h.shape[0]+1,dtype=float),return_eigenvectors=False))
    lo,hi=map(float,result['spectrum']['global']['neutral_gap_interval'])
    assert lo-2e-13<=values[1]-values[0]<=hi+2e-13


def test_spatial_residue_and_first_moment_against_direct_ground_vectors(result):
    w=old_packet.load_witness();g=old_certificate.eigensystem('ground',w['ground_columns'],w['ground_energies'],w['ground_vectors'])
    vector=np.asarray(g['block']['full_basis']@g['vector_numerators'][:,0],dtype=float)/(1<<60)
    vector/=np.linalg.norm(vector);masks=source.parent.masks(4,4)
    energy=float(g['energies'][0]);h=source.parent.centered(4,4).astype(float)/16
    for p in result['spectrum']['responses']:
        qs=[]
        for x in range(8):
            up=(masks>>(2*x))&1;down=(masks>>(2*x+1))&1
            qs.append((up-down)/2 if p['spin'] else up+down-1)
        images=vector[:,None]*np.asarray(qs).T
        expected=images.T@images
        for x in range(8):
            for y in range(8):
                lo=sum(float(r['spatial_residue_intervals'][x][y][0]) for r in p['rows'])
                hi=sum(float(r['spatial_residue_intervals'][x][y][1]) for r in p['rows'])
                assert lo-2e-13<=expected[x,y]<=hi+2e-13
        moment=float(images[:,0]@(h@images[:,0]-energy*images[:,0]))
        lo,hi=map(float,p['first_moment_interval'])
        assert lo-2e-13<=moment<=hi+2e-13


def test_all_outcome_table_and_continuity_weights(result):
    rows=list(csv.DictReader(io.StringIO(packet.table(result))))
    assert len(rows)==272
    assert sum(r['exactly_dark']=='True' for r in rows)==2
    for p in result['spectrum']['responses']:
        for r in p['rows']:
            e0,e1=map(F,r['energy_interval']);w0,w1=map(F,r['weight_interval'])
            c0,c1=map(F,r['current_divergence_weight_interval'])
            # Outward rounded input boxes can exceed the tighter stored current box.
            assert c0<=e1*e1*w1 and c1>=e0*e0*w0


def test_missing_full_block_rejected():
    w=packet.load_witness();del w['global_4_0_0_energies']
    with pytest.raises(ValueError):certificate.verify(w)


def test_tampered_basis_fails_residual():
    w=packet.load_witness();w['global_0_0_0_vectors']=w['global_0_0_0_vectors'].copy()
    w['global_0_0_0_vectors'][:,0]=0
    with pytest.raises(ArithmeticError):certificate.eigensystem(0,0,0,'global',w)


def test_dark_pullback_cannot_be_faked():
    w=packet.load_witness();w['probe_1_0_0_pullback']=w['probe_1_0_0_pullback'].copy();w['probe_1_0_0_pullback'][0]+=1
    old=old_packet.load_witness();g=old_certificate.eigensystem('ground',old['ground_columns'],old['ground_energies'],old['ground_vectors'])
    with pytest.raises(ArithmeticError,match='pullback'):certificate.neutral_response(1,0,w,g)


@pytest.mark.parametrize('field',['finite_gap_is_relativistic_mass','physical_species_identified','native_detector_executed','continuum_mass_gap_or_gravity_proved'])
def test_scope_rejects_boolean_numeric_laundering(result,field):
    bad=copy.deepcopy(result);bad['spectrum'][field]=0
    with pytest.raises(ValueError):packet.validate(bad)


def test_count_cannot_be_boolean(result):
    bad=copy.deepcopy(result);bad['spectrum']['responses'][0]['positive_lines']=True
    with pytest.raises(ValueError):packet.validate(bad)


def test_eta_transport_character_and_nested_scope(result):
    for row in result['spectrum']['global']['blocks']:
        assert row['N8_translation_character']==row['lowest_eta_translation_representative']^(7 if row['eta']%2 else 0)
    bad=copy.deepcopy(result);bad['spectrum']['global']['continuum_gap_claimed']=0
    with pytest.raises(ValueError):packet.validate(bad)
    bad=copy.deepcopy(result);bad['spectrum']['responses'][0]['exact_dark_lines']=True
    with pytest.raises(ValueError):packet.validate(bad)
