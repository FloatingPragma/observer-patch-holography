"""Integer/Fraction admission of complete block spectra and neutral response."""
from fractions import Fraction as F
from math import comb
import numpy as np
from . import source,proposal
from ..centered_particle_spectrum import certificate as prior,packet as prior_packet


def eigensystem(s,e,k,kind,witness):
    prefix=proposal.key(s,e,k,kind)
    columns=np.asarray(witness[prefix+'_columns'])
    if columns.ndim!=1 or columns.dtype!=np.dtype(np.int64):raise ValueError('exact column indices required')
    b=source.block(s,e,k,kind,tuple(int(i) for i in columns));n=b['dimension']
    ei=prior._integer_array(witness[prefix+'_energies'],(n,),'energy')
    vi=prior._integer_array(witness[prefix+'_vectors'],(n,n),'vector')
    energies=tuple(F(int(v),1<<proposal.ENERGY_BITS) for v in ei)
    if any(a>bb for a,bb in zip(energies,energies[1:])):raise ValueError('ordered spectrum required')
    if not n:return {'block':b,'energies':(),'radius':F(0)}
    y=b['basis'].astype(object)@vi
    ty=np.zeros(y.shape,dtype=object)
    for i,row in enumerate(b['coefficient_source']):
        for j in np.flatnonzero(row):ty[i]+=int(row[j])*y[j]
    residual=ty*(1<<proposal.ENERGY_BITS)-16*y*ei[None,:]
    gram=y.T@(b['sizes'][:,None]*y);unit=1<<(2*proposal.BITS)
    norms=[F(int(gram[i,i]),unit) for i in range(n)]
    defect=gram-unit*np.eye(n,dtype=object)
    delta=F(max(sum(abs(int(v)) for v in row) for row in defect),unit)
    if delta>=F(1,10**6):raise ArithmeticError('basis Gram defect too large')
    squares=[F(sum(int(d)*int(v)**2 for d,v in zip(b['sizes'],col,strict=True)),1<<(2*(proposal.BITS+proposal.ENERGY_BITS+4))) for col in residual.T]
    radius=prior.sqrt_bounds(sum(squares,F(0))/(1-delta))[1]
    if radius>=F(1,10**8):raise ArithmeticError('complete block enclosure failed')
    errors=[];intervals=[]
    for i,en in enumerate(energies):
        r=prior.sqrt_bounds(squares[i]/norms[i])[1]
        sep=min((abs(en-v) for j,v in enumerate(energies) if j!=i),default=F(1))-radius
        isolated=sep>2*radius and r<sep
        if kind=='probe' and not isolated:raise ArithmeticError('probe lines not individually resolved')
        errors.append(2*r/sep if isolated else None)
        intervals.append((en-(r if isolated else radius),en+(r if isolated else radius)))
    return dict(block=b,energies=energies,radius=radius,delta=delta,vectors=vi,orbit_vectors=y,
                norms=norms,errors=errors,intervals=intervals)


def complete_neutral(witness):
    summaries=[];systems={};coverage=[]
    for s,e,k in source.labels():
        row=eigensystem(s,e,k,'global',witness);b=row['block'];n=b['dimension']
        coverage.append((s,e,k,n));systems[(s,e,k)]=row
        summaries.append({'spin':s,'eta':e,'lowest_eta_translation_representative':k,
                          'N8_translation_character':k^(7 if e%2 else 0),
                          'character_multiplicity':comb(3,k.bit_count()),'dimension':n,
                          'complete_enclosure_radius':str(row['radius']),
                          'eigenvalue_intervals':[prior.decimal_bounds(x-row['radius'],x+row['radius'],12) for x in row['energies']]})
    covered=source.coverage(coverage)
    g=systems[(0,0,0)];gl,gu=g['intervals'][0]
    candidates=[]
    for key,row in systems.items():
        for i,en in enumerate(row['energies']):
            if key==(0,0,0) and i==0:continue
            candidates.append((en-row['radius'],en+row['radius'],key,i))
    low,high,key,index=min(candidates,key=lambda x:x[1])
    if key!=(1,0,7) or index!=0:raise ArithmeticError('lowest neutral sector changed')
    if min(x[0] for x in candidates)<=gu:raise ArithmeticError('ground not isolated from entire neutral spectrum')
    if any(x[0]<=high for x in candidates if x[2:]!=(key,index)):
        raise ArithmeticError('first neutral multiplet not isolated from other blocks')
    gap=(low-gu,high-gl)
    if not F(0)<gap[0]<=gap[1]<F(1,8):raise ArithmeticError('neutral gap must lie below proven charged floor')
    return systems,{'coverage':covered,'blocks':summaries,'ground_energy_interval':prior.decimal_bounds(gl,gu,12),
        'neutral_gap_interval':prior.decimal_bounds(*gap,places=12),
        'full_Fock_gap_interval':prior.decimal_bounds(*gap,places=12),
        'imported_other_charge_gap_strict_floor':'1/8','lowest_excitation_spin':1,'lowest_excitation_eta':0,
        'lowest_excitation_translation_character':7,'lowest_excitation_multiplicity':3,
        'gap_certifies_full_65536_state_source':True,'continuum_gap_claimed':False}


def spatial_moments(system):
    b=system['block'];q=b['q'].tocoo();y=np.zeros((q.shape[0],len(system['energies'])),dtype=object)
    y[q.row]=q.data[:,None]*system['orbit_vectors'][q.col]
    moments=[]
    for shift in range(8):
        rows,sign=source.parent.fock_action(b['up'],b['down'],(0,1,2),shift)
        vals=np.sum(y*(sign[:,None]*y[rows]),axis=0)
        moments.append([F(int(v),1<<(2*proposal.BITS))/norm for v,norm in zip(vals,system['norms'],strict=True)])
    ranks=[]
    for i,error in enumerate(system['errors']):
        val=sum((m[i]**2 for m in moments),F(0));tol=32*error
        allowed=[d for d in (1,3) if val-tol<=F(8,d)<=val+tol]
        if len(allowed)!=1:raise ArithmeticError('neutral spatial rank not resolved')
        ranks.append(allowed[0])
    return ranks,moments


def neutral_response(spin,eta,witness,ground):
    row=eigensystem(spin,eta,0,'probe',witness);b=row['block'];n=b['dimension']
    name=proposal.key(spin,eta,0,'probe')
    z=prior._integer_array(witness[name+'_dark'],(n,),'dark vector')
    a=prior._integer_array(witness[name+'_pullback'],(ground['block']['dimension'],),'dark pullback')
    den=prior._integer_array(witness[name+'_pullback_denominator'],(1,),'dark denominator')[0]
    if den<=0 or not np.any(z):raise ArithmeticError('invalid exact dark witness')
    t=source.transition(spin,ground['block'],b).astype(object)
    if np.any((b['numerator'].astype(object)+32*b['gram'].astype(object))@z):
        raise ArithmeticError('claimed dark vector is not an energy-minus-two eigenvector')
    if np.any(ground['block']['gram'].astype(object)@a-den*(t.T@z)):
        raise ArithmeticError('exact dark pullback mismatch')
    if np.any((ground['block']['numerator'].astype(object)+32*ground['block']['gram'].astype(object))@a):
        raise ArithmeticError('dark pullback does not lie in the energy-minus-two subspace')
    dark=[i for i,(lo,hi) in enumerate(row['intervals']) if lo<=-2<=hi]
    if len(dark)!=1:raise ArithmeticError('exact dark eigenline not isolated')
    amplitudes=row['vectors'].T@(t@ground['vector_numerators'][:,0])
    ranks,moments=spatial_moments(row)
    factor=F(1,2) if spin else F(2) # Sz versus S+; density n-1 versus pair annihilation.
    records=[];weightbounds=[];momentbounds=[]
    for i,rawnum in enumerate(amplitudes):
        lo,hi=row['intervals'][i];elo=lo-ground['energy_intervals'][0][1];ehi=hi-ground['energy_intervals'][0][0]
        if elo<=0:raise ArithmeticError('positive neutral response energy not proved')
        rootlo,roothi=prior.sqrt_bounds(row['norms'][i]*ground['norms'][0])
        raw=F(abs(int(rawnum)),1<<(proposal.BITS+prior.VECTOR_BITS))
        err=row['errors'][i]+ground['vector_errors'][0]
        if i in dark:wlo=whi=F(0)
        else:
            alo,ahi=raw/roothi,raw/rootlo
            if alo<=err:raise ArithmeticError('neutral weight unresolved: '+str((spin,i,alo,err)))
            wlo,whi=factor*(alo-err)**2,factor*(ahi+err)**2
        weightbounds.append((wlo,whi));momentbounds.append((wlo*elo,whi*ehi))
        spatial=[]
        for x in range(8):
            entries=[]
            for y in range(8):
                m=moments[x^y][i];lower=max(F(-1),m-2*row['errors'][i]);upper=min(F(1),m+2*row['errors'][i])
                # Raising eta from pair m=-1 to density m=0 introduces staggered signs.
                if not spin and (x^y).bit_count()%2:lower,upper=-upper,-lower
                values=(wlo*lower,wlo*upper,whi*lower,whi*upper)
                entries.append(prior.decimal_bounds(min(values),max(values),32))
            spatial.append(entries)
        records.append({'index':i,'energy_interval':prior.decimal_bounds(elo,ehi,14),
            'weight_interval':prior.decimal_bounds(wlo,whi,32),'strictly_positive':i not in dark,
            'exactly_dark':i in dark,'candidate_spatial_rank':ranks[i],
            'visible_spatial_rank':ranks[i] if i not in dark else 0,
            'spatial_residue_intervals':spatial,
            'current_divergence_weight_interval':prior.decimal_bounds(wlo*elo**2,whi*ehi**2,32)})
    positive=[r for r in records if r['strictly_positive']]
    return {'probe':'local hidden Sz' if spin else 'local density n-1','spin':spin,'eta':eta,
        'candidate_lines':n,'positive_lines':len(positive),'exact_dark_lines':len(dark),
        'rows':records,'total_weight_interval':prior.decimal_bounds(sum(x[0] for x in weightbounds),sum(x[1] for x in weightbounds),14),
        'first_moment_interval':prior.decimal_bounds(sum(x[0] for x in momentbounds),sum(x[1] for x in momentbounds),14),
        'fixed_component_visible_dimension':sum(r['visible_spatial_rank'] for r in positive),
        'minimum_visible_energy':positive[0]['energy_interval'],'maximum_visible_energy':positive[-1]['energy_interval'],
        'unresolved_weight':'0','current_probe':'outgoing current divergence; not independent transverse edge currents'}


def verify(witness):
    expected={proposal.key(s,e,k,kind)+'_'+suffix
              for s,e,k,kind in [(*x,'global') for x in source.labels()]+[(1,0,0,'probe'),(0,1,0,'probe')]
              for suffix in ('columns','energies','vectors')}
    expected.update(proposal.key(s,e,0,'probe')+'_'+suffix for s,e in ((1,0),(0,1))
                    for suffix in ('dark','pullback','pullback_denominator'))
    if set(witness)!=expected:raise ValueError('incomplete or extra neutral witness arrays')
    _,global_result=complete_neutral(witness)
    old=prior_packet.load_witness()
    ground=prior.eigensystem('ground',old['ground_columns'],old['ground_energies'],old['ground_vectors'])
    responses=[neutral_response(s,e,witness,ground) for s,e in ((1,0),(0,1))]
    # Exact local identity Sz^2+(n-1)^2/4=I/4 on the complete four-state site.
    low=F(responses[0]['total_weight_interval'][0])+F(responses[1]['total_weight_interval'][0])/4
    high=F(responses[0]['total_weight_interval'][1])+F(responses[1]['total_weight_interval'][1])/4
    if not low<=F(1,4)<=high:raise ArithmeticError('spin-density local sum rule failed')
    spin_first=tuple(map(F,responses[0]['first_moment_interval']))
    density_first=tuple(map(F,responses[1]['first_moment_interval']))
    if not density_first[0]-4*spin_first[1]<=0<=density_first[1]-4*spin_first[0]:
        raise ArithmeticError('exact density/spin first-moment ratio failed')
    reachable=3*responses[0]['fixed_component_visible_dimension']+responses[1]['fixed_component_visible_dimension']
    return {'global':global_result,'responses':responses,
            'all_spin_components_and_density_reachable_N8_dimension':reachable,
            'N8_dimension_not_reached_by_these_probes':comb(16,8)-reachable,
            'excited_N8_dimension_not_reached_by_these_probes':comb(16,8)-reachable-1,'local_sum_rule':'<Sz^2>+<(n-1)^2>/4=1/4',
            'first_moment_sum_rule':'M1(density)=4 M1(Sz)',
            'finite_gap_is_relativistic_mass':False,'physical_species_identified':False,
            'native_detector_executed':False,'continuum_mass_gap_or_gravity_proved':False}
