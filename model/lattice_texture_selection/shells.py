"""Exact oriented cubic shells and rational regular-value degree counts."""
from collections import Counter
from fractions import Fraction as F
from itertools import product


def integer(value, *, positive=False):
    if type(value) is not int or (positive and value < 1):
        raise ValueError('an exact integer in the stated domain is required')
    return value


def vector(value):
    if type(value) not in (tuple, list) or len(value) != 3:
        raise ValueError('three exact rational components required')
    if any(type(x) not in (int, F) for x in value):
        raise ValueError('three exact rational components required')
    answer = tuple(F(x) for x in value)
    if not any(answer):
        raise ValueError('a nonzero direction required')
    return answer


def dot(a, b):
    return sum(x*y for x, y in zip(a, b, strict=True))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def det(a, b, c):
    return dot(a, cross(b, c))


def edge(a, b):
    return tuple(sorted((a, b)))


def shell(radius, center=(0, 0, 0)):
    radius = integer(radius, positive=True)
    if type(center) is not tuple or len(center) != 3 or any(type(x) is not int for x in center):
        raise ValueError('an exact integer lattice center required')
    triangles, sides = [], set()
    for normal in range(3):
        a, b = [axis for axis in range(3) if axis != normal]
        ea, eb = tuple(int(i == a) for i in range(3)), tuple(int(i == b) for i in range(3))
        orientation = cross(ea, eb)[normal]
        for sign in (-1, 1):
            for u, v in product(range(-radius, radius), repeat=2):
                corners = []
                for da, db in ((0, 0), (1, 0), (0, 1), (1, 1)):
                    point = list(center)
                    point[normal] += sign*radius
                    point[a] += u+da
                    point[b] += v+db
                    corners.append(tuple(point))
                p, q, r, s = corners
                for triangle in ((p, q, r), (s, r, q)):
                    if orientation != sign:
                        triangle = (triangle[0], triangle[2], triangle[1])
                    triangles.append(triangle)
                sides.update(edge(x, y) for x, y in ((p, q), (p, r), (s, q), (s, r)))
    counted = Counter(edge(t[0], v) for t in triangles for v in t[1:])
    if set(counted) != sides or set(counted.values()) != {2}:
        raise ArithmeticError('surface source-edge count is not exactly two')
    if any(sum(abs(x-y) for x, y in zip(a, b, strict=True)) != 1 for a, b in sides):
        raise ArithmeticError('an interpolation diagonal entered the source-energy edge set')
    if len(triangles) != 48*radius**2 or len(sides) != 48*radius**2:
        raise ArithmeticError('cubic-shell face or edge count changed')
    return tuple(triangles), frozenset(sides)


def admissible_directions(triangles, directions):
    oriented, links = Counter(), {}
    for triangle in triangles:
        if len(triangle) != 3 or len(set(triangle)) != 3:
            raise ValueError('three distinct vertices required in every triangle')
        for i in range(3):
            apex, a, b = triangle[i], triangle[(i+1) % 3], triangle[(i+2) % 3]
            oriented[(apex, a)] += 1
            links.setdefault(apex, {}).setdefault(a, set()).add(b)
            links[apex].setdefault(b, set()).add(a)
    if not oriented or any(count != 1 or oriented[(b, a)] != 1
                           for (a, b), count in oriented.items()):
        raise ValueError('a complete consistently oriented closed triangulation is required')
    vertices = set(links)
    if len(vertices)-len(oriented)//2+len(triangles) != 2:
        raise ValueError('sphere Euler characteristic required')
    for adjacency in links.values():
        if any(len(neighbors) != 2 for neighbors in adjacency.values()):
            raise ValueError('each vertex link must be a cycle')
        reached, pending = set(), [next(iter(adjacency))]
        while pending:
            point = pending.pop()
            if point not in reached:
                reached.add(point)
                pending.extend(adjacency[point]-reached)
        if reached != set(adjacency):
            raise ValueError('each vertex link must be one connected cycle')
    reached, pending = set(), [next(iter(vertices))]
    while pending:
        point = pending.pop()
        if point not in reached:
            reached.add(point)
            pending.extend(set(links[point])-reached)
    if reached != vertices:
        raise ValueError('one connected shell required')
    if set(directions) != vertices:
        raise ValueError('directions must cover every shell vertex exactly')
    values = {v: vector(directions[v]) for v in vertices}
    for triangle in triangles:
        u, v, w = (values[p] for p in triangle)
        if min(dot(u, v), dot(v, w), dot(w, u)) < 0:
            raise ValueError('pairwise nonnegative directions are required on every triangle')
    return values


def _regular(target, pairs):
    for a, b in pairs:
        normal = cross(a, b)
        if any(normal):
            if dot(normal, target) == 0:
                return False
        elif not any(cross(a, target)):
            return False
    return True


def degree(triangles, directions, target=None):
    """Signed regular-target preimage count on every normalized triangle."""
    values = admissible_directions(triangles, directions)
    pairs = [(values[t[i]], values[t[j]]) for t in triangles for i, j in ((0, 1), (1, 2), (2, 0))]
    if target is None:
        # A nonzero plane meets this moment curve in at most two points.
        target = next((tuple(map(F, (1, k, k*k+1))) for k in range(1, 2*len(pairs)+3)
                       if _regular(tuple(map(F, (1, k, k*k+1))), pairs)), None)
        if target is None:
            raise ArithmeticError('finite regular-target search exhausted its polynomial bound')
    else:
        target = vector(target)
        if not _regular(target, pairs):
            raise ValueError('the target meets an edge span and is not admitted as regular')
    total, preimages, degenerate = 0, 0, 0
    for triangle in triangles:
        a, b, c = (values[p] for p in triangle)
        determinant = det(a, b, c)
        if determinant == 0:
            degenerate += 1
            continue
        coefficients = (det(target, b, c)/determinant,
                        det(a, target, c)/determinant,
                        det(a, b, target)/determinant)
        if min(coefficients) > 0:
            total += 1 if determinant > 0 else -1
            preimages += 1
    return {'degree': total, 'regular_target': [str(x) for x in target],
            'preimage_count': preimages, 'degenerate_triangle_count': degenerate,
            'all_oriented_triangles_retained': len(triangles)}


def radial_directions(triangles, center=(0, 0, 0), orientation=1):
    if type(orientation) is not int or orientation not in (-1, 1):
        raise ValueError('orientation must be exactly plus or minus one')
    return {p: tuple(orientation*(x-y) for x, y in zip(p, center, strict=True))
            for t in triangles for p in t}


def certificate():
    rows, used = [], set()
    for radius in (1, 2, 3, 4):
        triangles, edges = shell(radius)
        if used.intersection(edges):
            raise ArithmeticError('concentric shell source edges overlap')
        used.update(edges)
        for orientation in (-1, 1):
            result = degree(triangles, radial_directions(triangles, orientation=orientation))
            if result['degree'] != orientation:
                raise ArithmeticError('oriented radial shell has the wrong exact degree')
            rows.append({'radius': radius, 'orientation': orientation,
                         'source_edge_count': len(edges), **result})
    triangles, _ = shell(2)
    uniform = degree(triangles, {v: (1, 0, 0) for t in triangles for v in t})
    if uniform['degree'] != 0:
        raise ArithmeticError('uniform-reference degree is nonzero')
    first = shell(2, (-3, 0, 0))[1]
    second = shell(2, (3, 0, 0))[1]
    if first.intersection(second):
        raise ArithmeticError('separated opposite-charge shell supports overlap')
    return {'exact_radial_shells': rows, 'uniform_reference': uniform,
            'concentric_edges_disjoint': True, 'separated_center_distance': 6,
            'opposite_charge_outer_radius': 2, 'separated_edges_disjoint': True,
            'floating_angles_or_eigenvectors_used': False,
            'degree_is_microscopic_superselection_operator': False}
