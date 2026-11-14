"""Actual signed-current transfer of local charge on the coupled source."""
from bridge import *

def run():
    pins=verify();d=exact_geometry();B,A,H=d['B'],d['A'],d['harmonics']
    floating=source.data();rows=[];maxerror=0.;drop=[];reverse=[]
    for cell in range(8):
        for axis in range(3):
            for sign in (-1,1):
                destination,path,current=matter.route(cell,axis,sign)
                compiled=g.current(cell,axis,sign)
                assert compiled['destination']==destination and tuple(compiled['path'])==path
                assert list(compiled['integer_current'])==list(current)
                J=s.Matrix(list(map(int,current)))
                difference=d['particulars'][:,destination]-d['particulars'][:,cell]
                delta_charge=d['charges'][:,destination]-d['charges'][:,cell]
                assert B*J==delta_charge and B*difference==delta_charge
                transverse=J-difference
                assert B*transverse==s.zeros(96,1)
                assert A.T*difference==s.zeros(208,1) and H.T*difference==s.zeros(3,1)
                assert difference.dot(transverse)==0
                # Source axis_step uses these exact signed routes and shifts
                # every retained canonical momentum by frame^T J.
                matrix=floating['frame'];reconstructed=np.array(difference,float).reshape(-1)+matrix@(matrix.T@np.array(current,float))
                error=float(np.linalg.norm(reconstructed-np.array(current,float)));maxerror=max(maxerror,error)
                assert error<2e-11
                bad_drop=float(s.sqrt(transverse.dot(transverse)))
                bad_reverse=float(s.sqrt((2*delta_charge).dot(2*delta_charge)))
                drop.append(bad_drop);reverse.append(bad_reverse)
                # Every source branch moves one primitive local charge.
                assert sorted(x for x in delta_charge if x)==[-1,1]
                rows.append({'cell':cell,'axis':axis,'signed_rail':sign,'destination':int(destination),
                    'native_path':list(map(int,path)), 'integer_current_nonzero':[[j,int(x)] for j,x in enumerate(current) if x],
                    'exact_charge_changes':[[j,int(x)] for j,x in enumerate(delta_charge) if x],
                    'source_full_field_shift_reconstruction_residual':error,
                    'drop_current_translation_operator_shift_error':bad_drop,
                    'reverse_current_Gauss_error':bad_reverse})
    assert min(drop)>1 and min(reverse)>2
    assert verify()==pins
    save('DYNAMICS.json',{'all_checks_pass':True,'signed_native_routes':len(rows),'g':1,
        'maximum_full_field_current_shift_residual':maxerror,'routes':rows,
        'exact_transfer':'E_after−E_before=delta_particular+P_kerB J=J; B J=delta(G−n). Free quadratic field, Coulomb, mass and spin terms preserve the Gauss identity.',
        'numerical_frame_boundary':'Full-shift numerical residual is a diagnostic. Exact longitudinal/transverse decomposition and polar-frame range prove the identity without using the rounded frame as physical E.',
        'no_retained_current_history':True,'dropped_current_control_minimum_error':min(drop),
        'reversed_current_control_minimum_Gauss_error':min(reverse)})
if __name__=='__main__':run()
