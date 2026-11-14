"""Independent Bloch and eigenstate tests; demo-specific tests added when ready."""
import numpy as np
from independent_stacked_reference import incidence,positive_acoustic,hessian_squared_energy,A


def test_vertical_laplacian_and_positive_Hm_spinor():
    q=np.array([.07,.11,.13]);energy,state,b=positive_acoustic(q)
    horizontal=incidence([q[0],q[1],0])
    assert np.allclose(b@b.conj().T,horizontal@horizontal.conj().T+2*(1-np.cos(q[2]))*np.eye(24))
    h=np.block([[.23*np.eye(24),b],[b.conj().T,-.23*np.eye(87)]])
    assert abs(np.vdot(state,state)-1)<1e-12
    assert np.max(abs(h@state-energy*state))<1e-12


def test_physical_acoustic_metric_rank3_and_missing_z_links_control():
    expected=np.diag([15/352,15/352,1.])
    q_hessian=hessian_squared_energy()
    physical=A@q_hessian@A.T/2
    assert np.allclose(physical,expected,atol=2e-6)
    missing=hessian_squared_energy(0.)
    physical_missing=A@missing@A.T/2
    assert np.linalg.matrix_rank(physical,tol=1e-5)==3
    assert np.linalg.matrix_rank(physical_missing,tol=1e-5)==2
    assert abs(physical_missing[2,2])<1e-8


def test_demo_source_matches_independent_bloch_and_eigenstate():
    import stacked_wavepacket as demo
    for q in ((.07,.11,.13),(.75,.45,.10),(.2,-.13,.31)):
        reference=incidence(q)
        assert np.allclose(demo.bloch_incidence(q),reference)
        e,spinor=demo.positive_mode(q)
        e_reference,state,_=positive_acoustic(q)
        assert abs(e-e_reference)<1e-12
        assert np.allclose(np.outer(spinor,spinor.conj()),np.outer(state,state.conj()),atol=1e-12)
        assert demo.internal_offsets().shape==(111,3)
    assert np.allclose(demo.internal_offsets()[84:87],np.array([(demo.D-shift)/2 for shift in (np.zeros(3),demo.A[:,0],demo.A[:,1])]))
    assert np.allclose(demo.internal_offsets()[87:99],np.array([0,0,.5]))
    assert np.allclose(demo.internal_offsets()[99:111],demo.D+np.array([0,0,.5]))


def test_demo_Hellmann_Feynman_velocity_and_FFT_normalization():
    import stacked_wavepacket as demo
    grid=demo.band_grid(8)
    idx=(1,2,1);q=np.array([grid['q'][a] for a in idx]);step=1e-5
    independent=[]
    for axis in range(3):
        shift=np.eye(3)[axis]*step
        independent.append((positive_acoustic(q+shift)[0]-positive_acoustic(q-shift)[0])/(2*step))
    assert np.allclose(grid['velocity_cell'][idx],independent,atol=1e-8)
    assert np.allclose(demo.A@grid['velocity_cell'][idx],A@np.array(independent),atol=1e-8)
    rng=np.random.default_rng(65)
    amplitude=rng.normal(size=(8,8,8))+1j*rng.normal(size=(8,8,8));amplitude/=np.linalg.norm(amplitude)
    momentum=amplitude[...,None]*grid['spinor']
    position=np.fft.ifftn(momentum*np.exp(-.37j*grid['energy'])[...,None],axes=(0,1,2),norm='ortho')
    assert abs(np.sum(abs(position)**2)-1)<1e-12
    assert grid['minimum_phase_pivot']>1e-3


def test_absent_vertical_links_keep_complete_z_marginal_fixed():
    import stacked_wavepacket as demo
    grid=demo.band_grid(8,vertical=False);q=grid['q']
    assert np.allclose(grid['energy'][:,:,0],grid['energy'][:,:,-1])
    assert np.all(grid['velocity_cell'][...,2]==0)
    qq=np.stack(np.meshgrid(q,q,q,indexing='ij'),axis=-1)
    amplitude=np.exp(-np.sum((qq-np.array([.75,.45,.1]))**2,axis=-1)/(4*.5**2)).astype(complex)
    amplitude/=np.linalg.norm(amplitude)
    initial=np.fft.ifftn(amplitude[...,None]*grid['spinor'],axes=(0,1,2),norm='ortho')
    final=np.fft.ifftn(amplitude[...,None]*grid['spinor']*np.exp(-.37j*grid['energy'])[...,None],axes=(0,1,2),norm='ortho')
    assert np.allclose(np.sum(abs(initial)**2,axis=(0,1,3)),np.sum(abs(final)**2,axis=(0,1,3)),atol=1e-12)


def test_saved_packet_displacements_and_wrap_diagnostic():
    from pathlib import Path
    import json
    root=Path(__file__).resolve().parent
    if not (root/'position_distributions.npz').exists():
        return # Source-only checks are runnable before the producer emits data.
    result=json.loads((root/'STACKED_WAVEPACKET.json').read_text())
    with np.load(root/'position_distributions.npz') as data:
        probability=data['probability'];coordinates=data['coordinates'];times=data['times']
    assert np.allclose(np.sum(probability,axis=(1,2,3)),1,atol=1e-12)
    cell_centers=np.einsum('tijk,ijkd->td',probability,coordinates)
    observed=np.array([s['centroid'] for s in result['observables']['samples']])
    # Global component populations are constant in a single eigenband; thus
    # the fixed internal-offset contribution cancels from displacements.
    assert np.allclose(cell_centers-cell_centers[0],observed-observed[0],atol=1e-11)
    cell_indices=coordinates@np.linalg.inv(A).T
    n=result['supplied']['cells_per_axis'];edge=np.any(abs(cell_indices)>=n//2-2-1e-10,axis=-1)
    shell=probability[:,edge].sum(axis=1)
    assert np.allclose(shell,[s['edge_shell_probability'] for s in result['observables']['samples']],atol=1e-12)
    expected=observed[0]+(times-times[0])[:,None]*np.array(result['observables']['weighted_group_velocity'])
    assert np.isclose(np.max(np.linalg.norm(observed-expected,axis=1)),result['observables']['max_centroid_ray_distance'])


def test_saved_complete_offset_bins_reconstruct_position_apparatus():
    from pathlib import Path
    import json
    import stacked_wavepacket as demo
    root=Path(__file__).resolve().parent
    result=json.loads((root/'STACKED_WAVEPACKET.json').read_text())
    with np.load(root/'position_distributions.npz') as data:
        spatial=data['spatial_probability'];cell=data['probability'];offsets=data['internal_offsets']
        component_bins=data['component_to_offset_bin'];coordinates=data['coordinates']
    assert spatial.shape==(*cell.shape,7)
    assert np.all(spatial>=0) and np.allclose(spatial.sum(axis=-1),cell,atol=1e-14)
    assert np.allclose(offsets[component_bins],demo.internal_offsets())
    groups=spatial.sum(axis=(1,2,3))
    centers=np.einsum('tijk,ijkd->td',cell,coordinates)+groups@offsets
    reported=np.array([s['centroid'] for s in result['observables']['samples']])
    assert np.allclose(centers,reported,atol=1e-11)
    for axis in range(3):
        second=np.einsum('tijk,ijk->t',cell,coordinates[...,axis]**2)
        second+=2*np.einsum('tijkb,ijk,b->t',spatial,coordinates[...,axis],offsets[:,axis])
        second+=groups@(offsets[:,axis]**2)
        std=np.sqrt(np.maximum(0,second-centers[:,axis]**2))
        assert np.allclose(std,[s['position_std'][axis] for s in result['observables']['samples']],atol=1e-10)


def test_demo_rejects_invalid_positive_band_and_grid_inputs():
    import pytest
    import stacked_wavepacket as demo
    for mass in (0.,-.1,float('nan'),float('inf')):
        with pytest.raises(ValueError):demo.positive_mode([.1,.2,.3],mass)
    for q in ([.1,.2],[.1,.2,.3,.4],[.1,float('nan'),.3],[.1,float('inf'),.3]):
        with pytest.raises(ValueError):demo.positive_mode(q)
    for n in (0,3,4.5,float('nan'),True):
        with pytest.raises(ValueError):demo.band_grid(n)


def test_demo_rejects_invalid_packet_width_and_time_samples():
    import pytest
    import stacked_wavepacket as demo
    for sigma in (0.,-.1,float('nan'),float('inf')):
        with pytest.raises(ValueError):demo.simulate(n=4,sigma=sigma,times=[0.])
    for times in ([],[float('nan')],[float('inf')],[0.,float('nan')]):
        with pytest.raises(ValueError):demo.simulate(n=4,times=times)


def test_reference_ray_anchors_at_first_sample_time():
    import stacked_wavepacket as demo
    result,_,_=demo.simulate(n=4,times=[1.,1.1])
    samples=result['observables']['samples'];first=np.array(samples[0]['centroid'])
    velocity=np.array(result['observables']['weighted_group_velocity'])
    assert np.allclose(samples[0]['weighted_group_velocity_ray'],first,atol=1e-12)
    assert np.allclose(samples[1]['weighted_group_velocity_ray'],first+.1*velocity,atol=1e-12)
