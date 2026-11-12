from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
from math import gcd,comb
import numpy as np
import scipy.sparse as sp
from model.openwave_six_quark.source import colour48,exact_product,native_Fierz,op
from model.openwave_calibrated_hadrons.source import parameters,INPUTS


@lru_cache(None)
def local_census(number):
    native_Fierz();basis=tuple(sum(1<<i for i in ids)for ids in combinations(range(12),number))
    C=colour48(basis,(0,));assert np.all(C.data%16==0);C=C.copy();C.data//=16
    colours=[sum(1<<i for i in range(12)if i%3==c)for c in range(3)]
    keep=[i for i,b in enumerate(basis)if all((b&m).bit_count()==number//3 for m in colours)]
    C=C.tocsr()[keep][:,keep].tocsr();I=sp.eye(len(keep),dtype=np.int64,format='csr');P=I.copy();den=1
    for value in ((9,18)if number==3 else(9,18,24,36,54)):
        P=exact_product(value*I-C,P);den*=value;factor=gcd(den,*map(int,P.data));P.data//=factor;den//=factor
    assert (P-P.T).nnz==0 and(exact_product(P,P)-den*P).nnz==0 and exact_product(C,P).nnz==0
    weights=defaultdict(F)
    for k,index in enumerate(keep):
        b=basis[index];up=sum(bool(b&(1<<i))for i in range(6));nu=sum(bool(b&(1<<i))for i in (0,1,2,6,7,8))
        weights[2*up-number,2*nu-number]+=F(int(P[k,k]),den)
    rows=[]
    for s,i in product(range(number%2,number+1,2),repeat=2):
        count=weights[s,i]-weights[s+2,i]-weights[s,i+2]+weights[s+2,i+2]
        assert count.denominator==1 and count>=0
        if count:rows.append({'twice_spin':s,'twice_isospin':i,'multiplicity':int(count),'dimension':int(count)*(s+1)*(i+1)})
    expected={(1,1):1,(3,3):1}if number==3 else{(2,0):1,(0,2):1,(6,0):1,(0,6):1,(4,2):1,(2,4):1}
    assert {(r['twice_spin'],r['twice_isospin']):r['multiplicity']for r in rows}==expected
    dimension=sum(r['dimension']for r in rows);assert dimension==int(P.diagonal().sum())//den
    axial=None
    if number==3:
        src=tuple(basis[k]for k in keep);spins=[];weak=[];ax=[]
        for f,c in product(range(2),range(3)):spins.append((3*f+c,6+3*f+c,1))
        for s,c in product(range(2),range(3)):
            weak.append((6*s+c,6*s+3+c,1));ax.append((6*s+c,6*s+3+c,1 if s==0 else -1))
        Sp=op(spins,src,src);Ip=op(weak,src,src);Jp=op(ax,src,src)
        Zs=np.array([2*sum(bool(b&(1<<i))for i in range(6))-3 for b in src],dtype=np.int64)
        Zi=np.array([2*sum(bool(b&(1<<i))for i in(0,1,2,6,7,8))-3 for b in src],dtype=np.int64)
        S4=sp.diags(Zs*Zs+2*Zs,dtype=np.int64)+4*(Sp.T@Sp);I4=sp.diags(Zi*Zi+2*Zi,dtype=np.int64)+4*(Ip.T@Ip)
        PN=exact_product(15*I-S4,exact_product(15*I-I4,P));dn=den*144
        PN=exact_product(PN,sp.diags(((Zs==1)&(Zi==-1)).astype(np.int64),dtype=float,format='csr'))
        factor=gcd(dn,*map(int,PN.data));PN.data//=factor;dn//=factor
        assert(PN-PN.T).nnz==0 and(exact_product(PN,PN)-dn*PN).nnz==0 and PN.diagonal().sum()==dn
        axial=F(int(exact_product(PN,exact_product(Ip.T,Jp)).diagonal().sum()),dn);assert axial==F(5,3)
    return {'local_occupation_dimension':len(basis),'zero_Cartan_dimension':len(keep),'projector_denominator':den,
       'exact_colour_singlet_dimension':dimension,'complete_spin_isospin_multiplets':rows,
       'exact_native_local_nucleon_axial_ratio':str(axial)if axial is not None else None,
       'exact_weight_multiplicities':[{'twice_spin_z':s,'twice_isospin_z':i,'weight':str(v)}for(s,i),v in sorted(weights.items())if v]}


def casimir(j):return j*(j+1)
def minimum_square_occupancy(number,sites):
    q,r=divmod(number,sites);return r*(q+1)**2+(sites-r)*q*q


def graph(name,sites,edges):
    p=parameters();unit=F(INPUTS['electron_unit_GeV']);degree=[sum(x in e for e in edges)for x in range(sites)];z=max(degree)
    baryon=[]
    for twice in (1,3):
        j=F(twice,2);levels=[p['A']+(p['Js']+p['Jw'])*casimir(j)+(p['eta']*casimir(j)+3*p['t'])*d for d in degree]
        baryon.append({'twice_spin':twice,'twice_isospin':twice,'PHP_internal_levels':levels})
    bar_max=max(x for row in baryon for x in row['PHP_internal_levels']);bar_floor=F(8,3)*p['U']+p['A']*minimum_square_occupancy(3,sites)/9
    bar_gap=bar_floor-bar_max;assert bar_gap>0;bar_b=3*p['t']*z;bar_error=bar_b**2/bar_gap
    nuclear=[];local=local_census(6)['complete_spin_isospin_multiplets'];total=0
    for S,I in product(range(4),repeat=2):
        rows=[]
        if any(r['twice_spin']==2*S and r['twice_isospin']==2*I for r in local):
            for x in range(sites):
                energy=4*p['A']+p['Js']*casimir(S)+p['Jw']*casimir(I)+(p['eta']*casimir(I)+6*p['t'])*degree[x]
                rows.append({'occupation':'six at one cell','cells':[x],'PHP_internal_energy':energy})
        for x,y in combinations(range(sites),2):
            adjacent=(x,y)in edges or(y,x)in edges
            for a,b in product((F(1,2),F(3,2)),repeat=2):
                if not(abs(a-b)<=S<=a+b and abs(a-b)<=I<=a+b):continue
                energy=2*p['A']+(p['Js']+p['Jw'])*(casimir(a)+casimir(b))+3*p['t']*(degree[x]+degree[y])
                energy+=p['eta']*(casimir(a)*degree[x]+casimir(b)*degree[y]+int(adjacent)*(casimir(I)-casimir(a)-casimir(b)))
                rows.append({'occupation':'three at each of two cells','cells':[x,y],'local_twice_spin_isospin':[int(2*a),int(2*b)],'PHP_internal_energy':energy})
        rows.sort(key=lambda row:row['PHP_internal_energy']);total+=len(rows)*(2*S+1)*(2*I+1)
        nuclear.append({'spin':S,'isospin':I,'per_magnetic_component_low_dimension':len(rows),'PHP_states':rows})
    assert total==50*sites+400*comb(sites,2)
    nmax=max(x['PHP_internal_energy']for row in nuclear for x in row['PHP_states'])
    floor=F(8,3)*p['U']+p['A']*minimum_square_occupancy(6,sites)/9;gap=floor-nmax;assert gap>0
    coupling=6*p['t']*z;epsilon=coupling**2/gap
    for row in nuclear:
        row['complete_low_energy_intervals_relative_to_6mbar_minus_delta_Iz']=[[str(x['PHP_internal_energy']-epsilon),str(x['PHP_internal_energy'])]for x in row['PHP_states']]
        row['next_eigenvalue_internal_floor']=str(floor)
        success=F(row['per_magnetic_component_low_dimension'],comb(12*sites,6));tail=(coupling/gap)**2
        row['complete_original_preparation_and_spectral_records']={
          'preparation_success':str(success),'preparation_failure':str(1-success),
          'success_and_low_spectrum_interval':[str(success*(1-tail)),str(success)],
          'success_and_high_spectrum_interval':['0',str(success*tail)],
          'retained_failure_state':'Full complementary output of the exact local-colour/spin/isospin filter on the uniform full N6 occupation input; identity outcome on every other Fock sector.',
          'reference':'Initially correlated opposite native charge on each charge block; unchanged by the charge-preserving complete instruments.',
          'total_mean_work':'0 on the uniform full sector for preparation, and0 for the subsequent complete energy-band readout.'}
        for x in row['PHP_states']:x['PHP_internal_energy']=str(x['PHP_internal_energy'])
    for row in baryon:
        row['complete_low_energy_intervals_relative_to_3mbar_minus_delta_Iz']=[[str(x-bar_error),str(x)]for x in sorted(row['PHP_internal_levels'])]
        row['PHP_internal_levels']=list(map(str,row['PHP_internal_levels']));row['next_eigenvalue_internal_floor']=str(bar_floor)
    bmin=min(F(x)for x in baryon[0]['PHP_internal_levels']);bindings=[]
    for S,I in ((1,0),(0,0),(0,1)):
        row=next(x for x in nuclear if x['spin']==S and x['isospin']==I);emin=min(F(x['PHP_internal_energy'])for x in row['PHP_states'])
        lower=2*(bmin-bar_error)-emin;upper=2*bmin-emin+epsilon
        bindings.append({'spin':S,'isospin':I,'binding_interval':[str(lower),str(upper)],
          'binding_GeV_interval':[str(unit*lower),str(unit*upper)],'strictly_bound':lower>0,'excluded_below_two_baryons':upper<0})
    return {'graph':name,'vertices':sites,'edges':list(map(list,edges)),'degree':degree,
      'same_fitted_parameters_rational':{k:str(v)for k,v in p.items()},'baryon_complete_low_bands':baryon,'six_quark_complete_low_bands':nuclear,
      'nuclear_total_low_dimension_including_magnetic_components':total,
      'full_N6_occupation_dimension':comb(12*sites,6),'all_other_Fock_and_global_colour_sectors':'Retained; the stated low-spectrum theorem concerns the complete global-colour-singlet N6 sector, resolved by exact S,I and Iz.',
      'nuclear_complement_internal_floor':str(floor),'nuclear_off_block_norm_upper':str(coupling),'uniform_nuclear_energy_enclosure_error':str(epsilon),
      'baryon_uniform_energy_error':str(bar_error),'binding_channels':bindings,
      'complete_low_or_high_spectral_instrument':{'prepared_density':'normalized projector onto the complete locally-colour-singlet subspace in each stated spin/isospin/charge sector',
        'low_outcome_probability_interval':[str(1-(coupling/gap)**2),'1'],'high_or_complement_probability_interval':['0',str((coupling/gap)**2)],
        'all_times_and_external_quantum_references':True,'mean_spectral_readout_work':'0; exact spectral projectors commute with H',
        'reference':'idle opposite-charge reference; no discarded outcome or conditional renormalization',
        'execution_scope':'Exact native local CAR census and rational spectral/record bounds; no full high-dimensional wavefunction evolution is claimed'}}


@lru_cache(None)
def certificate():
    census={str(n):local_census(n)for n in(3,6)};rows=[]
    for L in(2,3,4):rows.append(graph(f'path{L}',L,tuple((x,x+1)for x in range(L-1))))
    for L in(3,4):rows.append(graph(f'cycle{L}',L,tuple((x,(x+1)%L)for x in range(L))))
    for row in rows:
        assert next(x for x in row['binding_channels']if(x['spin'],x['isospin'])==(1,0))['strictly_bound']
        assert next(x for x in row['binding_channels']if(x['spin'],x['isospin'])==(0,0))['strictly_bound']
    return {'schema':'oph.native_graph_nuclear_low_spectrum.v1','verdict':'PASS_COMPLETE_NATIVE_LOW_BANDS_GRAPH_REFINEMENT_WITH_ADDITIONAL_BOUND_SINGLET',
       'exact_native_local_CAR_census':census,'same_source_graphs':rows,
       'finite_theorem_shortcut':'P H P is diagonal on every locally-colour-singlet3+3 or6 occupation multiplet. Q H Q has a proved native-colour floor. With floor C>max PHP=a and ||PHQ||<=b, all dim(P) low energies lie in[lambda_j(PHP)-b²/(C-a),lambda_j(PHP)], and the next energy is at least C.',
       'additional_unassigned_bound_channel':'S0,I0 is strictly bound on every stated graph. Removing the S0,I1 exchange partner does not remove this different native state. Its physical orbital/spacetime dictionary is unestablished, so it is neither a proved new particle nor a proof of a forbidden physical nucleus.',
       'axial_scalar_deformation_obstruction':'The local native3Q nucleon irrep(S,I)=(1/2,1/2) has multiplicity1. Any colour/spin/weak scalar deformation is scalar on it by Schur, so it cannot change the fixed native axial matrix element5/3. Off-site, sea, orbital or current changes are distinct routes.',
       'scope':'All low nuclear spin/isospin bands below an explicit coloured-complement floor on paths and cycles of2–4 cells, with fixed calibration and complete low/high records. Graph variation is a controlled finite experiment, not a continuum or physical length calibration.'}
