"""Joint photon, native-clock, geometry and absorbing-detector records."""
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from .field import field,block,packet,packet_certificate,detector
from model.openwave_common_gravity_records.records import SETTINGS,LEVELS,certified_phases,clock_ratios,weak_clock_bound
from model.openwave_universal_energy_geometry.source import rest_block,profile,level_values
from .units import comparison as unit_comparison


def initial(source,polarization=2):
    _,V=eigh(rest_block(source));g=V[:,0]
    if g.sum()<0:g=-g
    p,_=packet(polarization)
    return np.stack([g/np.sqrt(2)]*4),np.stack([np.kron(g,p)/np.sqrt(2)]*4)


def evolve(b,v,w,time):
    return (b['Uv']@(certified_phases(tuple(b['ev']),time)*(b['Uv'].T@v)),
            b['V']@(certified_phases(tuple(b['e']),time)*(b['V'].T@w)))


def energy(vac,waves,blocks):
    return float(sum(np.vdot(v,b['vacuum_H']@v).real+np.vdot(w,b['H']@w).real
                     for v,w,b in zip(vac,waves,blocks,strict=True))/4)


def absorbed_energy(edges,blocks):
    return float(sum(np.einsum('ia,ij,ja->',w.conj(),b['vacuum_H'],w).real
                     for w,b in zip(edges,blocks,strict=True))/4)


def fine(vac,edges,blocks,before):
    # The photon is absorbed into the complete edge/direction pointer. Its
    # vacuum and the retained no-photon reference flag have separate labels.
    branches=np.stack([np.column_stack([v,w]).ravel() for v,w in zip(vac,edges,strict=True)])
    p=np.stack([abs(A@branches)**2/16 for A in SETTINGS])
    ecount=edges.shape[-1]
    diagonal=np.stack([np.repeat(np.diag(b['vacuum_H']),ecount+1) for b in blocks])
    after=sum(float(np.sum(p[k]*(abs(A)**2@diagonal))) for k,A in enumerate(SETTINGS))
    qstate=branches.reshape(2,2,-1)
    coherence=[sum(np.vdot(qstate[0,j],qstate[1,j])for j in range(2))/2,
               sum(np.vdot(qstate[j,0],qstate[j,1])for j in range(2))/2]
    return {'all4_setting4_native2_geometry_vacuum_and_every_edge_probabilities':p.reshape(4,4,2,ecount+1).tolist(),
            'total_probability':float(p.sum()),
            'native_clock_phases':[float(-np.angle(z))for z in coherence],
            'native_clock_visibilities':[float(abs(z))for z in coherence],
            'premeasurement_source_energy':before,'postmeasurement_source_energy':after,
            'complete_absorption_and_clock_readout_work':after-before,
            'all_other_native_clock_spectator_outcomes_probability_exact':'0'}


def edge_amplitudes(waves):
    d=field();m=d['frame'].shape[1]
    return np.einsum('ea,bga->bge',d['frame'],waves.reshape(4,2,m))


def fourier(edges):
    d=field();cells=d['cells']; probabilities=[];momenta=[]
    start=0
    for axis in range(3):
        shape=[cells+1]*3;shape[axis]=cells;count=int(np.prod(shape))
        a=np.fft.fftn(edges[:,:,start:start+count].reshape(4,2,*shape),axes=(2,3,4),norm='ortho')
        probabilities.extend(np.mean(abs(a)**2,axis=0).sum(axis=0).ravel().tolist())
        grid=np.array(np.meshgrid(*(2*np.pi*np.fft.fftfreq(n)for n in shape),indexing='ij')).reshape(3,-1).T
        momenta.extend(grid.tolist());start+=count
    p=np.array(probabilities);mean=np.array(momenta).T@p/p.sum()
    return {'all_directed_edge_Fourier_probabilities':probabilities,'mean_momentum':mean.tolist(),
            'direction_angle_xy':float(np.arctan2(mean[1],mean[0])),
            'separate_from_edge_position_instrument':True}


@lru_cache(None)
def execute(source,gamma=1,polarization=2):
    blocks=[block(source,j,gamma)for j in LEVELS];v0,w0=initial(source,polarization)
    d=field();positions=d['edge_positions'];m=d['frame'].shape[1]
    terminal=[]
    for t in (2.,4.,6.):
        states=[evolve(b,v,w,t)for b,v,w in zip(blocks,v0,w0,strict=True)]
        vac=np.stack([s[0]for s in states]);wave=np.stack([s[1]for s in states]);edges=edge_amplitudes(wave)
        r=fine(vac,edges,blocks,energy(vac,wave,blocks));p=np.mean(abs(edges)**2,axis=0).sum(axis=0)
        assert abs(r['total_probability']-1)<1e-10
        r.update({'source_parameter':t,'conditional_photon_position_mean':(positions.T@p/p.sum()).tolist(),
                  'photon_probability':float(p.sum()),'separate_complete_Fourier_instrument':fourier(edges)})
        terminal.append(r)
    mask,K0,dc=detector();vac,wave=v0.copy(),w0.copy();histories=[];coarse=readout=0.
    for step in range(1,13):
        states=[evolve(b,v,w,.5)for b,v,w in zip(blocks,vac,wave,strict=True)]
        vac=np.stack([s[0]for s in states]);wave=np.stack([s[1]for s in states]);before=energy(vac,wave,blocks)
        arrived=edge_amplitudes(wave)/np.sqrt(2);arrived[:,:,~mask]=0
        wave=np.einsum('ab,igb->iga',K0,wave.reshape(4,2,m)).reshape(4,2*m)
        absorbed=absorbed_energy(arrived,blocks)
        work=energy(vac,wave,blocks)+absorbed-before
        r=fine(np.zeros_like(vac),arrived,blocks,absorbed)
        r.update({'step':step,'source_parameter':step/2,'coarse_detector_work':work})
        histories.append(r);coarse+=work;readout+=r['complete_absorption_and_clock_readout_work']
    failure=fine(vac,edge_amplitudes(wave),blocks,energy(vac,wave,blocks))
    readout+=failure['complete_absorption_and_clock_readout_work']
    success=sum(r['total_probability']for r in histories)
    assert abs(success+failure['total_probability']-1)<1e-9
    initial_energy=energy(v0,w0,blocks)
    final_energy=sum(r['postmeasurement_source_energy']for r in histories)+failure['postmeasurement_source_energy']
    residual=final_energy-initial_energy-coarse-readout
    assert abs(residual)<1e-9
    # Twelve unitary/receiver steps, complete instrument contraction and a
    # conservative sum of all frame/root/eigensystem errors.
    certs=[b['certificate']for b in blocks];pc=packet_certificate(polarization)
    allocation={'normalized_packet':pc['normalized_packet_error_upper'],
                'terminal_Gauss_frame_and_Fourier_instruments':1e-8,
                'all_propagator_application_and_probability_roundoff':3e-8}
    # At most1080 terms per propagator row and at most882 per terminal
    # isometry row. Twelve steps and the finite3-axis FFT stay below the
    # displayed arithmetic allocation with binary64 gamma_n bounds.
    unit_roundoff=np.finfo(float).eps/2
    arithmetic_bound=12*4*(1080*unit_roundoff/(1-1080*unit_roundoff))*np.sqrt(1080)+1e-10
    assert arithmetic_bound<allocation['all_propagator_application_and_probability_roundoff']
    defect=max(6*(c['generator_error_upper']+c['source_operator_error_upper'])+36*c['gram_error_upper']for c in certs)+12*5e-9+5e-8
    TV=2*defect+defect*defect
    return {'source_native_level':source,'gamma_declared':gamma,'prepared_polarization_axis':polarization,
            'terminal_joint_records':terminal,'complete_first_arrival_histories':histories,
            'retained_no_arrival_and_photon_vacuum':failure,
            'arrival_success_probability':success,
            'conditional_arrival_mean_parameter':sum(r['source_parameter']*r['total_probability']for r in histories)/success,
            'complete_history_TV_error_upper':TV,'unitary_certificates':certs,'detector_certificate':dc,
            'packet_certificate':pc,'extra_vector_error_allocation':allocation,
            'complete_history_work':{'initial':initial_energy,'final':final_energy,
                                     'coarse_detector':coarse,'terminal_and_native_readout':readout,'residual':residual},
            'complete_preparation_and_work':{'input':'native source0, clock0, geometry ground and photon vacuum; idle reference flags',
                'input_energy_exact':'0','preparation':'declared finite unitary sends the pure zero-energy input to the stated source-ground/clock/packet/reference state; extends on its orthogonal complement',
                'mean_preparation_work':initial_energy,
                'mean_all_detector_and_readout_work':coarse+readout,
                'mean_total_work':initial_energy+coarse+readout,
                'final_mean_energy':final_energy,
                'total_work_minus_final_energy':initial_energy+coarse+readout-final_energy}}


def comparisons(a,b):
    eps=max(a['complete_history_TV_error_upper'],b['complete_history_TV_error_upper'])
    delay=b['conditional_arrival_mean_parameter']-a['conditional_arrival_mean_parameter']
    de=sum(12*eps/(x['arrival_success_probability']-eps)for x in (a,b))
    ta,tb=a['terminal_joint_records'][1],b['terminal_joint_records'][1]
    shift=tb['conditional_photon_position_mean'][1]-ta['conditional_photon_position_mean'][1]
    angle=tb['separate_complete_Fourier_instrument']['direction_angle_xy']-ta['separate_complete_Fourier_instrument']['direction_angle_xy']
    me=4*np.sqrt(2)*np.pi*eps
    ae=sum(me/(np.linalg.norm(t['separate_complete_Fourier_instrument']['mean_momentum'][:2])-me)for t in (ta,tb))
    phases=[]
    for j in range(2):
        local=[]
        for x in(a,b):
            t=x['terminal_joint_records'][0];theta=t['native_clock_phases'][j]
            modulus=t['native_clock_visibilities'][j]
            # Undo the1/4 setting weights conservatively before reconstructing
            # the two quadratures of the native Bloch coherence.
            delta=8*eps+1e-11;phase_error=delta/(modulus-delta)+1e-11
            assert theta-phase_error>0 and theta+phase_error<np.pi
            local.append((theta-phase_error,theta+phase_error))
        lo=local[1][0]/local[0][1]-1e-11;hi=local[1][1]/local[0][0]+1e-11
        phases.append([lo,hi])
    return {'conditional_arrival_delay':delay,'arrival_delay_error_upper':de,
            'transverse_position_difference':shift,'position_difference_error_upper':24*eps,
            'separate_Fourier_direction_deflection':angle,'direction_error_upper':ae,
            'same_joint_native_clock_phase_ratios':[b['terminal_joint_records'][0]['native_clock_phases'][j]/a['terminal_joint_records'][0]['native_clock_phases'][j]for j in range(2)],
            'same_joint_native_clock_phase_ratio_intervals':phases}


def run():
    trials=[execute(s)for s in (0,1)];comparison=comparisons(*trials)
    controls=[]
    for gamma in (-1,0):
        pair=[execute(s,gamma)for s in (0,1)];c=comparisons(*pair)
        if gamma==-1:
            assert abs(c['conditional_arrival_delay'])<1e-10 and abs(c['transverse_position_difference'])<1e-10
        controls.append({'gamma':gamma,'comparison':c,'complete_executions':pair})
    polarization=[execute(s,1,1)for s in (0,1)]
    d=field()
    assert comparison['conditional_arrival_delay']>comparison['arrival_delay_error_upper']
    assert comparison['transverse_position_difference']+comparison['position_difference_error_upper']<0
    assert comparison['separate_Fourier_direction_deflection']+comparison['direction_error_upper']<0
    assert all(hi<1 for lo,hi in comparison['same_joint_native_clock_phase_ratio_intervals'])
    # The same complete instruments are executed for controls. Their receipt
    # keeps every coarse/clock/work output and the native/geometry marginal,
    # avoiding duplicate bulk probability arrays.
    def compact(value):
        if isinstance(value,dict):
            out={}
            for k,v in value.items():
                if k=='all4_setting4_native2_geometry_vacuum_and_every_edge_probabilities':
                    out['complete_native_geometry_marginal']=np.array(v).sum(axis=-1).tolist()
                    out['executed_full_fine_table_shape']=[4,4,2,883]
                    out['complete_fine_table_marginalized_for_receipt']=True
                else:out[k]=compact(v)
            return out
        if isinstance(value,list):return[compact(v)for v in value]
        return value
    return {'schema':'oph.common_transverse_Maxwell_native_clock_geometry.v1',
            'verdict':'PASS_COMPLETE_TRANSVERSE_PHOTON_CLOCK_GEOMETRY_RECORDS',
            'source':{'complete_spatial_cube_cells':6,'vertices':343,'directed_edges':882,
                      'transverse_one_photon_dimension':540,'photon_vacuum_retained':True,
                      'native_source_and_clock_dimensions':[8,8],'geometry_dimension':2,
                      'complete_vacuum_one_photon_native_dimension':8*8*2*541,
                      'full_bosonic_source':'dGamma(sqrt(P S curl*curl S P)|ker grad*) with all photon-number sectors',
                      'S':'exp((1+gamma)Q F_edge), edge-average pinned profile',
                      'Maxwell_identification':'complete divergence-free electric-link modes, exact curl/gradient complex and both transverse polarizations; declared number-preserving positive-frequency constitutive source',
                      'canonical_metric_QED_Bogoliubov_transfer_claimed':False,
                      'all_boundary_Gauss_constraints_retained':True,'harmonic_one_forms':0,
                      'native_lapse':'exp(QF)','q_exact':'1/8','omega_exact':'1',
                      'spatial_response_gamma':1,'same_attractive_geometry_and_clock_source':True,
                      'edge_profile_exact':d['edge_profile_exact'],'packet_integer_coordinates':packet()[1],
                      'source_time_is_fundamental_global_clock':False},
            'full_source_executions':trials,'same_source_photon_comparison':comparison,
            'spatial_metric_controls':compact(controls),
            'second_physical_polarization_control':{'comparison':comparisons(*polarization),'complete_executions':compact(polarization)},
            'vacuum_sector_all_seven_native_clock_ratios':clock_ratios(),'vacuum_sector_weak_probe_universality':weak_clock_bound(),
            'operational_clock_position_records':[{'source_parameter':r['source_parameter'],
                'principal_native_clock_phase':r['native_clock_phases'][0],
                'native_clock_phase_coordinate_modulo_period':r['native_clock_phases'][0]/level_values()[1],
                'phase_unwrapped_by_external_clock':False,
                'conditional_photon_position_mean':r['conditional_photon_position_mean']}
                for r in trials[1]['terminal_joint_records']],
            'root_and_transverse_certificates':d['certificates'],
            'independent_energy_length_attachment':unit_comparison(),
            'scope':{'physical_Newton_strength_derived':False,'Einstein_continuum_proved':False,
                     'full_strength_clocks_exactly_universal':False,'physical_seconds_or_lensing_angle_identified':False,
                     'finite_cavity_transverse_EM_and_native_clock_records_share_one_Hamiltonian':True}}
