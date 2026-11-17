"""Remove the installed spin cost, retaining the actual positive hopping.

The point is to test a distinct operator choice without scanning a target.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import audit
HERE=Path(__file__).resolve().parent


def run():
    t=F(1,8);U=F(1);eta=F(1,8)
    P_floor=3*t
    Q_floor=8*U/3
    coupling_squared=6*t*t
    uniform_floor=min(P_floor-eta,Q_floor-coupling_squared/eta)
    assert uniform_floor==F(1,4)
    # Verify the exact P/Q coupling norm bound on every finite frozen path,
    # including all spin/weak basis vectors of both relevant projections.
    controls=[]
    for L in [2,3,4,5]:
        for w in [1,2]:
            d=audit.source.data(L,F(1),w,J=F())
            P=[i for i,x in enumerate(d['span']) if x==0]
            Q=[i for i,x in enumerate(d['span']) if x!=0]
            # Different colocated source columns have disjoint one-hop images:
            # in a (2,1) pattern the majority cell identifies its source.
            support=[]
            for i in P:
                columns={j:z for (j,k),z in d['entries'].items() if k==i and j in Q}
                norm={}
                for z in columns.values():norm=audit.source.add(norm,audit.source.mul(z,z))
                assert set(norm)<={1} and norm.get(1,F())<=coupling_squared
                support.append(set(columns))
            assert all(not a&b for i,a in enumerate(support) for b in support[i+1:])
            controls.append(dict(L=L,weakup=w,colocated_columns=len(P),disjoint_hop_images=True))
    d=audit.source.data(2,F(1),1,J=F())
    interval=audit.ground_interval(*audit.rational_matrix(d,audit.extension_values(d,'none')))
    massless_interval=[x-3 for x in interval]
    assert massless_interval[0]>uniform_floor
    inputs=json.loads((HERE/'RAW_COEFFICIENT_INPUTS.json').read_text())
    anchor=F(inputs['electron_anchor_GeV']);comparison=F('0.93827208943')
    rows=[]
    for name in ['432','3600']:
        elo,ehi=map(F,inputs[name]['light_interval'])
        lower=3*anchor+uniform_floor*anchor/ehi
        exact_interval=[3*anchor+massless_interval[0]*anchor/ehi,
                        3*anchor+massless_interval[1]*anchor/elo]
        assert lower>comparison
        rows.append(dict(analyzer=name,uniform_mass_lower_GeV=str(lower),
            uniform_mass_lower_GeV_decimal=float(lower),
            two_cell_mass_interval_GeV=list(map(str,exact_interval)),
            two_cell_mass_interval_GeV_decimal=list(map(float,exact_interval))))
    out=dict(schema='oph.obstruction_repair.hadrons.Jzero.v1',
        amendment='Remove the separately installed local spin-square term, J=0, retain t=1/8,U=1,kappa=1 and the full Gauss-chain source',
        new_calibration=False,all_path_lengths_L_at_least_2_bound='H_kinetic+colour+string >= 1/4 in the three-quark singlet domain',
        P_floor=str(P_floor),Q_floor=str(Q_floor),offdiagonal_squared_upper=str(coupling_squared),
        completing_square_eta=str(eta),uniform_floor=str(uniform_floor),
        finite_exact_controls=controls,two_cell_massless_energy_interval=list(map(str,massless_interval)),
        rows=rows,physical_mass_solution=False,
        conclusion='Removing the imposed spin coefficient is insufficient. The positive graph kinetic and colour localization cost already exceeds the physical proton at either actual raw-electron coefficient. A signed Dirac/current source with its actual dressed vacuum is a different construction and is not excluded by this bound.')
    (HERE/'KINETIC_FLOOR_RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['uniform_floor','two_cell_massless_energy_interval','rows','physical_mass_solution']},indent=2))


if __name__=='__main__':run()
