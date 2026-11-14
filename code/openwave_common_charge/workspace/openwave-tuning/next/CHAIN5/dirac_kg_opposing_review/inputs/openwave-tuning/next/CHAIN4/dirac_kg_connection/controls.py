"""Reject wrong rail, residual Coulomb, scalar-source and flavor transfers."""
from connection import *

def run():
    pins=verify()
    # An always-positive spatial increment is not the signed Dirac rail.
    wrong=[]
    for j in reversed(range(2)):wrong+=multi_x({k:1 for k in range(j)},j)
    W,Wf,Wg,Wb=execute(events(gates_word(wrong)),3,np.eye(8))
    right=signed_shift_gates(2,(0,1),2);R,Rf,Rg,Rb=execute(events(gates_word(right)),3,np.eye(8))
    rail_error=float(np.linalg.norm(W-R,2));assert rail_error>1.9
    # Disabling only transverse current coupling leaves a position-dependent
    # Coulomb potential. It does not yield the claimed free Fourier factor.
    energies=np.array([float(x['coulomb']) for x in g.data()['sectors']]);uniform=np.ones(8)/math.sqrt(8)
    charged=np.exp(-1j*energies)*uniform;leak=float(np.linalg.norm(charged-uniform*np.vdot(uniform,charged)))
    assert leak>1e-3
    # The scalar acoustic velocities are rank two, not Clifford generators.
    old=[]
    for axis in range(3):
        V=np.zeros((4,4));V[0,axis+1]=V[axis+1,0]=1;old.append(V)
    scalar_defect=max(float(np.linalg.norm(V@V-np.eye(4),2)) for V in old);assert scalar_defect==1
    # Replacing this source's c=1 with 1/sqrt48 changes a native allowed mode.
    h=F(1,8);k=np.array([math.pi/4,0,0]);mass=MASSES[0]
    U=bloch(k,h,mass);bad=bloch(k/math.sqrt(48),h,mass)
    speed_error=float(np.linalg.norm(U-bad,2));assert speed_error>.6
    # The flavor padding has zero mass, but carries the defined kinetic
    # source action. Treating it as a projected-away output loses a block.
    padding=bloch(k,h,F());padding_defect=float(np.linalg.norm(padding-np.eye(4),2));assert padding_defect>.7
    # Strong, state-controlled g -> 0 connection on the same finite L2 source.
    # No operator-norm smallness is claimed for unbounded oscillator currents.
    d=g.data();eta=F(d['spectral']['actual_curvature_operator_error_upper'])
    omega_lower=g.interval_sqrt(g.A*g.B*(min(d['lambda'])-eta),100)[0]
    vco=g.A/(2*omega_lower);coulomb=max(row['coulomb'] for row in d['sectors'])
    coupled=[]
    for n in (1,2,4):
        qvar=vco+F(1,2)+F(n*n,8)
        kick=g.interval_sqrt(5*qvar,100)[1]
        for coupling in (F(1,100000),F(1,1000),F(1)):
            error=min(F(2),3*n*abs(coupling)*kick+n*coupling**2*coulomb)
            coupled.append({'macros':n,'coupling_exact':str(coupling),'state_error_upper_exact':str(error),'state_error_upper':float(error)})
    assert verify()==pins
    out={'all_checks_pass':True,'wrong_native_rail_operator_separation':rail_error,
        'retained_Coulomb_plane_leakage':leak,'scalar_velocity_Clifford_defect':scalar_defect,
        'wrong_scalar_clock_speed_operator_separation':speed_error,'discarded_flavor_padding_kinetic_defect':padding_defect,
        'actual_native_rail_complete_Gram_residuals':[Rg,Wg],
        'finite_size2_small_coupling_connection':coupled,
        'small_coupling_domain':'Normalized true coexact vacuum and normalized harmonic Gaussians; arbitrary matter/flavor/reference state. Exact-source telescoping evaluated on free prefixes, not an unbounded-generator operator-norm assertion.',
        'source_domain':'g=0 requires both g current and g² Coulomb scaling; no scalar-source or full90-Wilson equality is admitted.'}
    (HERE/'CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':run()
