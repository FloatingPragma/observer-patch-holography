import json,math
M4=155+(5739/256)*math.pi
C=M4/(3*math.pi)
rows=[]
for z in [.3,.6,1.]:
    for gap in [16,32,64,128,1000]:
        a=gap/2;b=gap/(2*math.exp(abs(z)))
        raw=C*(a**-3+b**-3)
        rows.append(dict(rapidity=z,gap_over_R=gap,raw_analytic_upper_bound=raw,operator_norm_upper_bound=min(1.,raw)))
with open('BOUND.json','w') as f:json.dump(dict(status='ANALYTIC_FULL_CONTINUUM_ONE_DIMENSIONAL_OPERATOR_BOUND',M4_upper_bound=M4,tail_coefficient=C,rows=rows,three_dimensional_bound='M6/(6 pi^2) (R/r)^3, with source-defined M6 finite but not yet numerically evaluated',no_modes_discarded=True),f,indent=2)
print(json.dumps(dict(M4=M4,C=C,selected=rows[8]),indent=2))
