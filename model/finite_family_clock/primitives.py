"""Minimal FR family-register primitives, with no external path dependency.

Layout and construction are attributed in PROVENANCE.md to repair-law/work/FR.
The actual native seam event is independently compared in family.py.
"""
import sympy as s

Z=s.Symbol('z',nonzero=True)


def index(load,record,orientation):
    return 4*load+2*record+orientation


def repair():
    matrix=s.eye(12)
    for orientation in range(2):
        a,b=(0,1) if orientation==0 else (1,0)
        vector=s.zeros(12,1)
        for load,record,sign in ((a,0,1),(b,0,-1),(a,1,-1),(b,1,-1)):
            vector[index(load,record,orientation)]=sign/s.sqrt(2)
        matrix-=vector*vector.T
    return matrix


def address_reflection():
    matrix=s.zeros(12)
    for load in range(3):
        for record in range(2):
            for orientation in range(2):
                matrix[index(load,record,1-orientation),index(load,record,orientation)]=1
    return matrix


def oriented_record():
    return s.diag(*[-1 if (i//2)%2==1 and i%2==1 else 1 for i in range(12)])


def state(tail,head,idle,record=0,address='uniform'):
    vector=s.zeros(12,1)
    for load,amplitude in enumerate((tail,head,idle)):
        orientations=range(2) if address=='uniform' else (address,)
        for orientation in orientations:
            vector[index(load,record,orientation)]=amplitude/(s.sqrt(2) if address=='uniform' else 1)
    return vector


def shift():
    return s.Matrix(3,3,lambda i,j:int(i==(j+1)%3))


def circulant(z=Z):
    r=shift()
    return s.eye(3)/s.sqrt(6)+(z*r+r.T/z)/s.sqrt(12)


def flags():
    idle=state(0,0,1)
    agree=state(1/s.sqrt(2),1/s.sqrt(2),0)
    plus=state(1/s.sqrt(2),1/s.sqrt(2),0,1,0)
    minus=state(1/s.sqrt(2),1/s.sqrt(2),0,1,1)
    return {'C_idle':idle,'E_idle':idle,'E_agree':agree,'C_rep+':plus,'C_rep-':minus}


def analyzer():
    f=flags()
    return (s.kronecker_product(f['C_idle'],f['E_agree'])
            +s.kronecker_product(f['C_rep+'],f['E_idle'])
            +s.kronecker_product(f['C_rep-'],f['E_idle']))/s.sqrt(3)


def compression(phi,after,z=Z):
    matrix=s.zeros(3)
    for ic in range(12):
        record,orientation=(ic//2)%2,ic%2
        step=(1 if orientation==0 else -1) if record else 0
        for ie in range(12):
            coefficient=s.conjugate(phi[12*ic+ie])*after[ic]*after[ie]*z**step
            if coefficient:
                for j in range(3):
                    matrix[(j+step)%3,j]+=coefficient
    return matrix.applyfunc(s.simplify)
