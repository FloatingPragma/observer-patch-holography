import numpy as np
from .analysis import influence,autocorrelation
from model.openwave_scale_coupling.coupling import creutz,tree_coefficients

def test_influence_is_actual_nonlinear_gradient():
 x=np.array([[.59,.32,.19,.065,.023],[.61,.33,.20,.067,.025]])
 mean=x.mean(axis=0);eps=1e-6;expected=[]
 for row in x:
  direction=row-mean
  expected.append((creutz(mean+eps*direction)-creutz(mean-eps*direction))/(2*eps*tree_coefficients(10)))
 assert np.max(abs(influence(x,10)-expected))<1e-8

def test_constant_record_and_correlated_control():
 assert autocorrelation(np.ones(400))['effective_records']==400
 rng=np.random.default_rng(151);noise=rng.normal(size=4000);series=np.zeros(4000)
 for j in range(1,len(series)):series[j]=.9*series[j-1]+noise[j]
 assert autocorrelation(series)['tau_integrated']>3
 assert autocorrelation(series)['effective_records']<800

def test_reordering_correlated_records_changes_uncertainty_diagnostic():
 rng=np.random.default_rng(91);x=np.repeat(rng.normal(size=100),20)
 shuffled=x.copy();rng.shuffle(shuffled)
 assert autocorrelation(x)['tau_integrated']>5*autocorrelation(shuffled)['tau_integrated']
