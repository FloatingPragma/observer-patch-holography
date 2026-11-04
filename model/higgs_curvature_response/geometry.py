"""Direct coordinate curvature and covariant derivatives for exact controls."""
import sympy as s


def clean(value):
    return s.cancel(s.expand(value))


class Geometry:
    """Four-dimensional nondegenerate metric; no field equation is imposed."""

    def __init__(self, coordinates, metric):
        self.coordinates = tuple(coordinates)
        self.metric = s.ImmutableMatrix(metric)
        if len(self.coordinates) != 4 or self.metric.shape != (4, 4) or self.metric != self.metric.T:
            raise ValueError('four coordinates and a symmetric4x4 metric required')
        if self.metric.det() == 0 or any(x.has(s.Float) for x in self.metric):
            raise ValueError('nondegenerate exact metric required')
        self.inverse = s.ImmutableMatrix(self.metric.inv().applyfunc(clean))
        q, g, gi = self.coordinates, self.metric, self.inverse
        self.gamma = [[[clean(sum(gi[a, d]*(s.diff(g[d, c], q[b])+
                                     s.diff(g[d, b], q[c])-s.diff(g[b, c], q[d]))
                                  for d in range(4))/2)
                       for c in range(4)] for b in range(4)] for a in range(4)]
        G = self.gamma
        self.ricci = s.ImmutableMatrix(4, 4, lambda a, b: clean(
            sum(s.diff(G[c][a][b], q[c])-s.diff(G[c][a][c], q[b])+
                sum(G[c][c][d]*G[d][a][b]-G[c][b][d]*G[d][a][c]
                    for d in range(4)) for c in range(4))))
        self.scalar = clean(sum(gi[a, b]*self.ricci[a, b] for a in range(4) for b in range(4)))
        self.einstein = s.ImmutableMatrix((self.ricci-g*self.scalar/2).applyfunc(clean))

    def hessian(self, value):
        q, G = self.coordinates, self.gamma
        return s.ImmutableMatrix(4, 4, lambda a, b: clean(s.diff(value, q[a], q[b])-
                     sum(G[c][a][b]*s.diff(value, q[c]) for c in range(4))))

    def box(self, value):
        hessian = self.hessian(value)
        return clean(sum(self.inverse[a, b]*hessian[a, b] for a in range(4) for b in range(4)))

    def divergence(self, tensor):
        tensor = s.ImmutableMatrix(tensor)
        if tensor.shape != (4, 4):
            raise ValueError('covariant4x4 tensor required')
        q, G, gi = self.coordinates, self.gamma, self.inverse
        return s.ImmutableMatrix([clean(sum(gi[a, c]*(s.diff(tensor[a, b], q[c])-
                       sum(G[d][c][a]*tensor[d, b]+G[d][c][b]*tensor[a, d]
                           for d in range(4))) for a in range(4) for c in range(4)))
                       for b in range(4)])

    def curvature_response(self, F):
        """Half the metric Euler tensor from sqrt(-g)F R."""
        return s.ImmutableMatrix((F*self.einstein+self.metric*self.box(F)-self.hessian(F)).applyfunc(clean))
