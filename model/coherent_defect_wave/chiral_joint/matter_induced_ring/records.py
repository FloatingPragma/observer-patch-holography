"""Complete source-clock loop-arrival records from the unchanged Hamiltonian."""
import sympy as s
from .model import physical_split,NORM_BOUND


def complete_records():
    block=physical_split();h,g,p=block['h'],block['g'],block['p']
    initial=s.zeros(72,1);initial[p[0]]=1
    assert (initial.T*g*initial)[0]==1
    targets={}
    for name,j in (('singlet_loop',2),('doublet_loop',4)):
        vector=s.zeros(72,1);vector[p[j]]=block['normalization'][j,j]
        targets[name]=vector
    effects={name:vector*vector.T*g for name,vector in targets.items()}
    effects['all_other_states']=s.eye(72)-sum(effects.values(),s.zeros(72))
    assert sum(effects.values(),s.zeros(72))==s.eye(72)
    for name,effect in effects.items():
        assert effect*effect==effect and effect.T*g==g*effect
        for other,f in effects.items():
            if name!=other:assert effect*f==s.zeros(72)
    moment_rows=[];power=initial;moments={name:[] for name in targets}
    for n in range(5):
        moment_rows.append({'power':n,'full_state_norm_squared':str((power.T*g*power)[0])})
        for name,target in targets.items():moments[name].append((target.T*g*power)[0])
        if n<4:power=h*power
    assert moments['singlet_loop']==[0,0,s.Rational(1,16),s.Rational(3,16),s.Rational(329,128)]
    assert moments['doublet_loop']==[0,0,0,0,s.Rational(3,8)]
    norm4=(power.T*g*power)[0]
    assert norm4==s.Rational(238254041,8192)
    bound4=s.Integer(s.ceiling(s.sqrt(norm4)));assert bound4**2>=norm4
    time=s.Rational(1,1<<14)
    error=NORM_BOUND*bound4*time**5/(s.factorial(5)*(1-NORM_BOUND*time))
    assert NORM_BOUND*time<1
    lower_ec=time**2/32-s.Rational(329,3072)*time**4-error
    upper_ec=time**2/32+time**3/32+s.Rational(329,3072)*time**4+error
    lower_L=time**4/64-error;upper_L=time**4/64+error
    assert 0<lower_ec<upper_ec and 0<lower_L<upper_L
    bounds={'singlet_loop':(lower_ec**2,upper_ec**2),'doublet_loop':(lower_L**2,upper_L**2)}
    bounds['all_other_states']=(1-bounds['singlet_loop'][1]-bounds['doublet_loop'][1],
                                1-bounds['singlet_loop'][0]-bounds['doublet_loop'][0])
    assert 0<bounds['all_other_states'][0]<bounds['all_other_states'][1]<1
    rows=[{'outcome':name,'effect_rank':int(s.trace(effect)),
           'probability_interval':[str(value) for value in bounds[name]],
           'postmeasurement_state':'E_'+name+' exp(-itH) |initial>; retained without renormalization'}
          for name,effect in effects.items()]
    assert sum(bounds[name][0] for name in bounds)<=1<=sum(bounds[name][1] for name in bounds)
    return {'initial_raw_key':[64,0,0,0,0],'initial_state_norm_squared':'1',
            'preparation':'declared pure physical source state with static compensator (0,-1)',
            'preparation_reference':'opposite hypercharge, lepton number and spin projection; matching fermion parity',
            'ground_preparation_used':False,'time_in_unchanged_source_clock':str(time),
            'complete_three_outcome_projective_instrument':True,'outcome_pointer_states':3,
            'source_block_dimension':72,'effects_outside_prepared_block':'zero success and identity other-states effect',
            'full_H_power_norms':moment_rows,
            'singlet_loop_full_H_moments':[str(v) for v in moments['singlet_loop']],
            'doublet_loop_full_H_moments':[str(v) for v in moments['doublet_loop']],
            'fourth_order_Taylor_state_error_upper_bound':str(error),
            'full_H_fourth_power_norm_upper_bound':str(bound4),
            'error_formula':'30*bound4*t^5/(120*(1-30*t))',
            'outcomes':rows,'both_loop_arrival_probability_lower_bounds_strictly_positive':True,
            'failure_and_all_other_states_retained':True,'success_conditioning_applied':False,
            'effective_resolvent_kernel_used_as_time_generator':False}
