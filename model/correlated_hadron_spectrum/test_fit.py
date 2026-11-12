import unittest

import numpy as np

from .fit import covariance_metric, fit_mass, fit_window, mean_covariance, stream_means, gev_row


class CorrelatedFitTests(unittest.TestCase):
    def test_known_correlated_cosh_and_exponential(self):
        t=np.arange(3,9);n=32;amplitude=1.7;mass=.43;extent=24
        corr=.72**abs(t[:,None]-t[None,:])
        cov=np.outer(np.exp(-mass*t),np.exp(-mass*t))*corr*1e-6
        for kind in ['exp','cosh']:
            shape=np.exp(-mass*t)
            if kind=='cosh':shape+=np.exp(-mass*(extent-t))
            y=amplitude*shape
            # Orthogonal zero-mean stream directions give the chosen mean
            # covariance exactly, without a lucky random noise realization.
            basis=np.sqrt(2/n)*np.cos(np.pi*(np.arange(n)[:,None]+.5)*np.arange(1,len(t)+1)[None,:]/n)
            samples=y+np.sqrt(n*(n-1))*basis@np.linalg.cholesky(cov).T
            jk=(samples.sum(0)-samples)/(n-1)
            got=fit_window(y,jk,np.arange(len(t)),extent-2*t[0],kind)
            # Shifted coordinates change the amplitude convention, not m.
            self.assertAlmostEqual(got['central']['mass'],mass,places=7)
            direct=fit_mass(y,cov,t,extent,kind)
            self.assertAlmostEqual(direct['mass'],mass,places=7)
            self.assertAlmostEqual(direct['amplitude'],amplitude,places=6)
            self.assertEqual(got['replica_count'],32)
            self.assertTrue(got['central_covariance_equals_original_jackknife'])

    def test_singular_and_ill_conditioned_covariance_rejected(self):
        with self.assertRaisesRegex(ValueError,'singular'):
            covariance_metric(np.ones((4,4)))
        with self.assertRaises(ValueError):
            covariance_metric(np.array([[1.,1-1e-10],[1-1e-10,1.]]))
        with self.assertRaises(ValueError):
            covariance_metric(np.diag([1.,0.,1.]))

    def test_stream_deletion_and_common_scale(self):
        x=np.arange(120,dtype=float).reshape(20,6)+np.sin(np.arange(120)).reshape(20,6)
        mean=x.mean(0);jk=(x.sum(0)-x)/19
        self.assertTrue(np.allclose(stream_means(mean,jk),x))
        cov=mean_covariance(x)
        j=jk-jk.mean(0)
        self.assertTrue(np.allclose(cov,19/20*j.T@j))
        row={'mrho_phys':.6,'mN_phys':.9,'sqrt_sigma':.4}
        a=gev_row(row,.178,.77526,.44)
        b=gev_row(row,.178,2*.77526,2*.44)
        self.assertTrue(all(abs(b[k]-2*a[k])<1e-14 for k in a))
        self.assertEqual(a['pion_rho_anchor'],.178*.77526)


if __name__=='__main__':unittest.main()
