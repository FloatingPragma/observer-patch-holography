"""Exact spin-current, Ward and reflection-positivity frame controls."""
import sympy as s


def comm(a,b):return a*b-b*a


def certificate():
    h=s.Rational(1,2); eye=s.eye(2)
    spin=(s.Matrix([[0,h],[h,0]]),s.Matrix([[0,-s.I*h],[s.I*h,0]]),s.diag(h,-h))
    left=[s.kronecker_product(x,eye) for x in spin]
    right=[s.kronecker_product(eye,x) for x in spin]
    bond=sum((a*b for a,b in zip(left,right)),s.zeros(4))-s.eye(4)/4
    f0,f1,g0,g1=s.symbols('f0 f1 g0 g1',real=True)
    rows=[]
    for a in range(3):
        b,c=(a+1)%3,(a+2)%3
        current=left[b]*right[c]-left[c]*right[b]
        assert s.simplify(s.I*comm(bond,left[a])+current)==s.zeros(4)
        assert s.simplify(s.I*comm(bond,right[a])-current)==s.zeros(4)
        assert comm(bond,left[a]+right[a])==s.zeros(4)
        probe=(f0+s.I*g0)*left[a]+(f1+s.I*g1)*right[a]
        difference=(f0-f1)**2+(g0-g1)**2
        transverse=left[b]*right[b]+left[c]*right[c]
        assert s.simplify(comm(probe.conjugate().T,comm(bond,probe))+difference*transverse)==s.zeros(4)
        assert max(abs(v) for v in transverse.eigenvals())==h
        assert s.simplify(comm(left[b],left[c])-s.I*left[a])==s.zeros(4)
        rows.append({'axis':a,'current_norm':str(max(abs(v) for v in current.eigenvals())),
                     'transverse_bond_norm':'1/2','conserved_total_spin':True})
    # Sublattice rotation about axis3 flips the transverse exchange signs.
    u=s.kronecker_product(eye,s.I*2*spin[2])
    rotated=u.conjugate().T*bond*u
    target=-left[0]*right[0]-left[1]*right[1]+left[2]*right[2]-s.eye(4)/4
    assert rotated==target
    # The distinct DLS reflection frame: real components 1,3 negative,
    # imaginary component 2 positive. It is not identified with the PF frame.
    reflection=s.kronecker_product(eye,s.I*2*spin[1])
    assert reflection.conjugate().T*bond*reflection == -left[0]*right[0]+left[1]*right[1]-left[2]*right[2]-s.eye(4)/4
    # At infinitesimal positive staggered field along axis1, it becomes a
    # uniform negative off-diagonal single-spin field after the same rotation.
    assert s.simplify(u.conjugate().T*(left[0]-right[0])*u-left[0]-right[0])==s.zeros(4)
    # Broken symmetry alone does not imply positive-energy gaplessness.
    free_spin_counterexample=-right[2]
    assert sorted(free_spin_counterexample.eigenvals().items())==[(-h,2),(h,2)]
    return {'exact_spin_current_and_moment_axes':rows,'sublattice_rotation_checked':True,'DLS_real_imaginary_reflection_frame_checked':True,
            'staggered_field_becomes_uniform_transverse_field':True,
            'free_spin_counterexample':{'ground_multiplicity':2,'gap_above_ground_space':'1'},
            'current_is_neutral_internal_spin_current':True}
