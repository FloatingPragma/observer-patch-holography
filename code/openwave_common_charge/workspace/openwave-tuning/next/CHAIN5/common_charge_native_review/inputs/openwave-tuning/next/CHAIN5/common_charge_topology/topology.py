"""Local electric curvature, exact Gauss Chern spectrum and global controls."""
from bridge import *
from fractions import Fraction
import math

def sphere_complex(degree=6):
    faces=[];N=degree;S=degree+1
    for j in range(degree):faces.extend([(N,j,(j+1)%degree),(S,(j+1)%degree,j)])
    edges=sorted({tuple(sorted((x,y))) for f in faces for x,y in zip(f,f[1:]+f[:1])})
    edge_index={e:j for j,e in enumerate(edges)};boundary=s.zeros(len(edges),len(faces));incidence=s.zeros(degree+2,len(edges))
    for j,(a,b) in enumerate(edges):incidence[a,j]=-1;incidence[b,j]=1
    for j,f in enumerate(faces):
        for a,b in zip(f,f[1:]+f[:1]):boundary[edge_index[tuple(sorted((a,b)))],j]+=1 if a<b else -1
    assert incidence*boundary==s.zeros(degree+2,len(faces))
    assert boundary*s.ones(len(faces),1)==s.zeros(len(edges),1)
    assert len(boundary.nullspace())==1 and boundary.nullspace()[0]==s.ones(len(faces),1)
    # Every coefficient in the fundamental cycle is one. Integrality of
    # the source period uses this primitive sphere, without a degree cover.
    return {'vertices':degree+2,'edges':edges,'oriented_faces':faces,'Euler_characteristic':degree+2-len(edges)+len(faces),
            'boundary':[[int(x) for x in row] for row in boundary.tolist()],'primitive_fundamental_cycle':[1]*len(faces)}

def sampled_clutching(flux):
    flux=list(map(F,flux));N=max(16,math.ceil(16*max(map(abs,flux),default=F()))+8)
    phases=[];prefix=F()
    for value in flux:
        for j in range(N):
            t=F(j,N);cdf=10*t**3-15*t**4+6*t**5
            phases.append(2*math.pi*float(prefix+value*cdf))
        prefix+=value
    phases.append(2*math.pi*float(prefix));z=np.exp(1j*np.array(phases))
    increments=np.angle(z[1:]/z[:-1]);assert max(abs(increments))<math.pi/2
    winding=float(sum(increments)/(2*math.pi))
    return {'flux_sum_exact':str(prefix),'clutching_endpoint_real':float(z[-1].real),'clutching_endpoint_imag':float(z[-1].imag),
        'sampled_winding':winding,'subdivisions_per_native_edge_patch':N,'maximum_phase_increment':float(max(abs(increments)))}

def run():
    pins=verify();d=exact_geometry();B,A,H,Q=d['B'],d['A'],d['harmonics'],d['Q']
    embedding=read(HERE/'EMBEDDING.json');den=embedding['coordinate_denominator']
    starts=s.Matrix(embedding['oriented_edge_start_lifts']);ends=s.Matrix(embedding['oriented_edge_end_lifts'])
    spatial_period=(ends-starts).T/(2*den)
    port_part=s.zeros(3,96)
    for vertex,point in enumerate(embedding['vertex_integer_coordinates']):
        cell=d['graph'].cell_coordinates[vertex//12]
        for axis in range(3):port_part[axis,vertex]=-s.Rational(point[axis]-den*cell[axis],2*den)
    assert spatial_period-d['periods']==port_part*B
    assert spatial_period*A==s.zeros(3,208)
    mesh=sphere_complex();assert mesh['Euler_characteristic']==2
    examples=[];maxerror=0.
    # Rational samples lie in the full physical affine Gauss fibers. The
    # exact operator theorem covers every continuous canonical coordinate.
    transverse=A[:,0]/3-A[:,17]/5+A[:,77]/7+H[:,0]/11-H[:,1]/13+H[:,2]/17
    assert B*transverse==s.zeros(96,1)
    for cell in range(8):
        E=d['particulars'][:,cell]+transverse
        assert B*E==d['charges'][:,cell]
        for vertex in range(96):
            flux=[F(B[vertex,e]*E[e]) for e in range(288) if B[vertex,e]]
            assert len(flux)==6
            result=sampled_clutching(flux);expected=int(d['charges'][vertex,cell])
            assert F(result['flux_sum_exact'])==expected and abs(result['sampled_winding']-expected)<2e-11
            assert abs(result['clutching_endpoint_real']-1)<1e-12 and abs(result['clutching_endpoint_imag'])<2e-12
            maxerror=max(maxerror,abs(result['sampled_winding']-expected))
            if vertex in (0,12*cell):examples.append({'cell':cell,'surface_vertex':vertex,'patch_fluxes_exact':list(map(str,flux)),**result})
    # Exact physical U lies in ker B through KQ Lambda^-1. Q itself does not.
    BQ=B*Q;entry=next((i,j) for i in range(BQ.rows) for j in range(BQ.cols) if BQ[i,j])
    v,j=entry;amplitude=1/(2*BQ[v,j]);assert BQ[v,j]*amplitude==s.Rational(1,2)
    rounded={'vertex':v,'mode':j,'exact_BQ_entry':str(BQ[v,j]),'coefficient_giving_half_flux':str(amplitude),
             'fractional_flux_exact':'1/2','control':'Replacing U by rounded Q as the exact physical field violates Gauss and local Chern integrality. Q is only a separately bounded apparatus approximation.'}
    # Half charge obeys a real continuity equation but has nonperiodic
    # transition function. It is outside the integer occupation fiber.
    # Keep the declared unit background and total particle charge fixed.
    # This mutant has sharp half charges at two distinct daughter roots.
    half=(d['particulars'][:,1]+d['particulars'][:,2])/2
    half_boundary=B*half
    assert half_boundary[0]==1 and half_boundary[12]==half_boundary[24]==s.Rational(-1,2)
    half_flux=[F(B[12,e]*half[e]) for e in range(288) if B[12,e]]
    half_control=sampled_clutching(half_flux)
    half_control['declared_background_preserved']=True
    half_control['two_daughter_fluxes_exact']=['-1/2','-1/2']
    half_control['total_relative_charge_exact']='1'
    assert F(half_control['flux_sum_exact'])==F(-1,2) and half_control['clutching_endpoint_real']<-.999999
    # Global torus periods are not constrained by the same local Gauss law.
    harmonic_projection=H*(H.T*H).inv()*H.T
    global_half=harmonic_projection*d['winding'][:,0]/2
    assert B*global_half==s.zeros(96,1) and A.T*global_half==s.zeros(208,1)
    periods=d['periods']*global_half
    assert periods==s.Matrix([s.Rational(1,2),0,0])
    # The exact regional-flux observable is diagonal in the occupied cell.
    regions={'particle_singleton':[12],'reference_singleton':[0],
        'x_positive_roots':[12*c for c in range(4,8)],'daughter_A':[12,24,36],
        'daughter_B':[48,60,72,84],'all_ports':list(range(96))}
    operators={}
    for name,vertices in regions.items():
        indicator=s.zeros(96,1)
        for vertex in vertices:indicator[vertex]=1
        cut=B.T*indicator
        exact=indicator.T*d['charges']
        assert cut.T*d['particulars']==exact and cut.T*H==s.zeros(1,3) and cut.T*A==s.zeros(1,208)
        # cut^T U=0 follows from B U=0, independently of frame rounding.
        values=[int(x) for x in exact]
        assert all(x in (-1,0,1) for x in values)
        operators[name]={'vertices':vertices,'Gauss_reference_in_region':int(0 in vertices),'flux_eigenvalues_by_cell':values,
                         'physical_relative_charge_eigenvalues_by_cell':[-x for x in values],
                         'oriented_boundary_edges':[[e,int(cut[e])] for e in range(288) if cut[e]]}
    assert operators['x_positive_roots']['flux_eigenvalues_by_cell']==[0]*4+[-1]*4
    # A particle's charge cannot be present in both disjoint daughters.
    qa=[-x for x in operators['daughter_A']['flux_eigenvalues_by_cell']]
    qb=[-x for x in operators['daughter_B']['flux_eigenvalues_by_cell']]
    assert all(a*b==0 and a+b in (0,1) for a,b in zip(qa,qb))
    out={'schema':'oph.common_charge.local_electric_Chern.v1','checks':{'actual_degree6_dual_sphere':True,
        'primitive_oriented_fundamental_cycle':True,'all768_local_gauss_curvature_periods_integer':True,
        'exact_U_transverse_not_rounded_Q':True,'half_charge_clutching_rejected':True,'global_harmonic_bundle_claim_rejected':True,
        'regional_flux_equals_actual_occupation_operator':True,'disjoint_half_daughters_excluded':True},
        'sphere_complex':mesh,'examples':examples,'maximum_sampled_winding_error':maxerror,
        'curvature_formula':'On six labelled edge patches, F=2pi f(phi) w(theta) dtheta wedge dphi, integral_patch F=2pi B_ve E_e. f is the normalized 30t^2(1-t)^2 patch bump; w has unit integral and vanishes near the poles.',
        'connection_formula':'A_N=2pi f(phi) W(theta) dphi; A_S=2pi f(phi)(W(theta)-1)dphi; transition exp(i chi(phi)), chi=2pi integral_0^phi f. Its endpoint is exp(2pi i sum_e B_ve E_e).',
        'topological_charge':'c1_v=(B E)_v=G_v-n_v; physical relative source charge n_v-G_v=-c1_v.',
        'quantum_scope':'Joint commuting electric spectral measure on each exact Gauss fiber. The local Chern-class observable is the integer-spectrum operator G_R-n_R, not an assigned definite premeasurement field.',
        'rounded_Q_counterexample':rounded,'half_charge_counterexample':half_control,
        'embedding_period_join':'Actual lifted-edge displacement cocycle differs from source periods by the exact vertex coboundary port_part*B. Every closed-current period and every source-face null period agrees.',
        'global_control':{'harmonic_periods_exact':list(map(str,periods)),'local_chern_numbers':[0]*96,
            'global_line_bundle_exists_for_this_curvature':False,'local_bounding_surface_bundles_exist':True,
            'nonbounding_half_flux_is_not_half_charge':True},'regions':operators,
        'half_expectation_control':{'uniform_source_plane_mean_flux':'-1/2','sharp_flux_values':[0,-1],'weights':['1/2','1/2']},
        'source_g':1,'extra_history_register':False,'extra_field_dynamics':False,'official_rows_changed':[]}
    assert verify()==pins;save('TOPOLOGY.json',out)
if __name__=='__main__':run()
