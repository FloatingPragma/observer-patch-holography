"""Exact occupation records, isospin transport and certified radial moments."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import json
from pathlib import Path
from model.openwave_common_hyperfine_spectrum import source as native
from model.openwave_axial_state.certificate import exact_action
from model.openwave_calibrated_hadrons.source import integral
from model.openwave_baryon_weak_formfactors.native import currents
from model.openwave_string_positive_cone.source import sqrt_interval

HERE=Path(__file__).resolve().parent
INPUTS=json.loads((HERE/'inputs.json').read_text())


def sinc_interval(x,order=18):
    """Taylor's theorem, with rational endpoints, including x=0."""
    x=F(x)
    value=sum((-1)**j*x**(2*j)/F(factorial(2*j+1)) for j in range(order+1))
    error=abs(x)**(2*order+2)/F(factorial(2*order+3))
    return value-error,value+error


def radial_records(modes,basis,vector):
    norm=sum(x*x for x in vector)
    assert norm>0
    profile={}
    for mask,amplitude in zip(basis,vector,strict=True):
        if not amplitude:continue
        occupied=[m for i,m in enumerate(modes)if mask>>i&1]
        assert len(occupied)==3
        center=F(sum(m[1]for m in occupied),3)
        for _,x,spin,flavour,colour in occupied:
            radius=abs(F(x)-center)
            profile[radius]=profile.get(radius,F(0))+F(amplitude*amplitude,norm)*(F(2,3)if flavour==0 else -F(1,3))
    return profile


def charge_form_factor(profile,momentum_GeV,length_fm,sine_bound):
    if not momentum_GeV:
        q=sum(profile.values());return q,q
    k=F(momentum_GeV)*F(length_fm)/F(INPUTS['hbar_c_GeV_fm'])
    low=high=F(0)
    for radius,charge in profile.items():
        a,b=sinc_interval(k*radius)
        low+=charge*(a if charge>=0 else b)
        high+=charge*(b if charge>=0 else a)
    # The charge density has norm at most two on every N3 occupation.
    return low-4*sine_bound,high+4*sine_bound


@lru_cache(None)
def certificate():
    # Consume the certified parent vector, then independently apply this source
    # and its current. A fresh floating eigensolver proposal is unnecessary.
    parent=json.loads((HERE.parent/'openwave_common_hyperfine_spectrum/results.json').read_text())
    selected=parent['native_axial_ground'];data=native.components('baryon')
    neutron=list(map(int,selected['integer_native_trial']))
    norm=sum(x*x for x in neutron);sine=F(selected['ground_vector_sine_bound'])
    raising=integral(sum(currents()['full'][x][0][0]for x in range(2)))
    proton=exact_action(raising,neutron)
    assert sum(x*x for x in proton)==norm
    p=native.parameters()
    delta=p['delta'];unit=F(native.INPUTS['electron_unit_GeV'])
    commutators=[]
    for key,A in data['parts'].items():
        if key not in native.weights(p):continue
        residual=A@raising-raising@A
        if key=='D':residual=residual+2*raising
        assert residual.nnz==0,key
        commutators.append(key)
    def action(v):
        out=[F(0)for _ in v]
        for key,weight in native.weights(p).items():
            for i,x in enumerate(exact_action(data['parts'][key],v)):out[i]+=weight*x
        return out
    hn=action(neutron);hp=action(proton)
    # Raising is integral, but the Hamiltonian vector is rational.
    rn=[sum(int(raising.data[k])*hn[raising.indices[k]]for k in range(raising.indptr[i],raising.indptr[i+1]))for i in range(len(hn))]
    assert hp==[x-delta*y for x,y in zip(rn,proton,strict=True)]
    attachment=parent['common_length_attachment'];length=F(attachment['a_fm'])
    # Every occupation row belongs to one exact colour-epsilon column.
    R=data['Binteger'].tocsr();gram=list(map(int,data['Gram']))
    diagonal=[]
    for i in range(len(data['basis'])):
        terms=list(range(R.indptr[i],R.indptr[i+1]));assert len(terms)<=1
        diagonal.append(sum((F(int(R.data[k])**2,gram[R.indices[k]])for k in terms),F(0)))
    assert sum(diagonal)==120
    rows=[];radii={}
    for name,charge,vector,hv in [('neutron',0,neutron,hn),('proton',1,proton,hp)]:
        profile=radial_records(data['modes'],data['basis'],vector)
        assert sum(profile.values())==charge
        if name=='neutron':
            assert profile[F(1,3)]>4*sine and profile[F(2,3)]<-4*sine
        moment=sum(weight*r*r for r,weight in profile.items())
        error=F(16,9)*sine
        interval=((moment-error)*length**2,(moment+error)*length**2)
        radii[name]=interval
        mean=sum(x*y for x,y in zip(vector,hv,strict=True))/norm
        variance=sum((x-mean*y)**2 for x,y in zip(hv,vector,strict=True))/norm
        assert variance==F(selected['trial_variance'])
        probabilities=[F(x*x,norm)for x in vector]
        assert sum(probabilities)==1 and all(0<=x<=d for x,d in zip(probabilities,diagonal,strict=True))
        original=[{'native_mask':str(mask),'success':str(prob/120),'failure':str((d-prob)/120)}
                  for mask,prob,d in zip(data['basis'],probabilities,diagonal,strict=True)if d]
        assert sum(F(r['success'])+F(r['failure'])for r in original)==1
        ff=[]
        for momentum in map(F,('0','.05','.1','.2','.4','.6','1')):
            lo,hi=charge_form_factor(profile,momentum,length,sine)
            ff.append({'momentum_GeV':str(momentum),'exact_ground_form_factor_interval':list(map(str,(lo,hi))),'interval_diagnostic':[float(lo),float(hi)]})
        row={'assigned_channel':name,'charge_exact':str(charge),'trial_rest_energy_GeV':str(mean*unit),
             'energy_role':'same seven-input strong calibration; no weak-contact or continuum pole transfer',
             'ground_vector_sine_bound':str(sine),'native_signed_radial_charge':[
                 {'radius_graph':str(r),'radius_fm':str(r*length),'trial_charge_exact':str(q),
                  'ground_charge_interval':list(map(str,(q-4*sine,q+4*sine)))}for r,q in sorted(profile.items())],
             'mean_square_charge_radius_graph_exact':str(moment),
             'mean_square_charge_radius_fm2_interval':list(map(str,interval)),
             'mean_square_charge_radius_fm2_diagnostic':[float(x)for x in interval],
             'orientation_averaged_form_factors':ff,'complete_native_records':original,
             'unlisted_native_occupation_records':'Both outcomes have probability zero; all2024 native masks and the external Fock/Gauss complement are retained.'}
        if name=='proton':
            bounds=(sqrt_interval(interval[0])[0],sqrt_interval(interval[1])[1])
            row['rms_charge_radius_fm_interval']=list(map(str,bounds));row['rms_charge_radius_fm_diagnostic']=[float(x)for x in bounds]
        rows.append(row)
    assert F(rows[0]['trial_rest_energy_GeV'])-F(rows[1]['trial_rest_energy_GeV'])==delta*unit
    # A common length cancels exactly in this independent comparison.
    nl,nh=radii['neutron'];pl,ph=radii['proton'];ratio=(-nh/ph,-nl/pl)
    rp=F(INPUTS['proton_rms_radius_fm']);rp_error=F(INPUTS['proton_radius_uncertainty_fm'])
    rn=F(INPUTS['neutron_mean_square_radius_fm2']);rn_error=F(INPUTS['neutron_radius_uncertainty_fm2'])
    observed=(-(rn+rn_error)/(rp+rp_error)**2,-(rn-rn_error)/(rp-rp_error)**2)
    assert ratio[0]>observed[1]
    return {'schema':'oph.native_baryon_form_factors.v1','verdict':'PASS_NATIVE_SIGNED_PROFILES_AND_SHARED_LENGTH_RADIUS_TEST',
            'source':'selected seven-input hyperfine strong member; exact two-cell Gauss transfer',
            'physical_length_fm_exact':str(length),'length_ancestry':attachment['ancestry'],'unused_radius_comparisons':INPUTS,
            'native_states':rows,'mass_ordering_GeV_exact':str(delta*unit),
            'exact_full_N3_current_commutator_checks':commutators,
            'complete_physical_instrument':{'source_input':'maximally mixed on all120 N3 colour singlets, with charge-balancing reference correlated by sector',
                'native_probe_effects':'E_i=B*|i><i|B; square-root Kraus operators; all2024 labels retained',
                'filter_success_exact':'1/120','filter_failure_exact':'119/120',
                'full_unconditional_record_TV_to_ground_upper':str(2*sine/120),
                'arbitrary_reference_half_channel_bound':str(4*sine),
                'nonselective_preparation_plus_probe_source_work_exact':'0 on the declared maximally mixed physical input',
                'reason':'Every colour-epsilon occupation row occupies one orbital column; both filter and compressed probe are unital and preserve charge.',
                'all_other_sectors':'Full complementary identity outcome; zero on this declared input.'},
            'shared_length_comparison':{'minus_neutron_over_proton_radius_squared_interval':list(map(str,ratio)),
                'observed_one_sigma_box_interval':list(map(str,observed)),
                'no_single_length_matches_both_radii_in_comparison_box':True,
                'physical_radius_target_used_to_fit_source_or_length':False},
            'scope':'Complete native charge profiles, same-source calibrated GeV energies and orientation-averaged graph form factors. The two compared charge radii are not matched by any common length in this source. Rotational averaging is a declared observation law, not a relativistic Sachs-current, 3D spatial dynamics or physical nucleon pole theorem.'}
