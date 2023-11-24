from math import sin, cos, sqrt, pi

import matplotlib.pyplot as plt
import numpy as np

from Euclidian import *


def init():
    figure, ax = plt.subplots()
    ax.set_aspect(1)
    square = plt.Rectangle((-2, -2), 4, 4, color='lightskyblue')
    ax.add_artist(square)
    unit_circle = plt.Circle((0, 0), 1, color='pink')
    ax.add_artist(unit_circle)
    unit_boundary = plt.Circle((0, 0), 1, fill=False, color='red')
    ax.add_artist(unit_boundary)
    ax.scatter(0, 0, color='red')
    ax.axis([-2, 2, -2, 2])
    return ax


class Mobius:
    def __init__(self, u: complex):
        self.u = u

    def __call__(self, obj):
        return obj.mobius(self)

    def f(self, z):
        return (self.u + z) / (1 - (self.u.conjugate() * z) / _a2)

    @classmethod
    def movePointToOrigin(cls, p: 'Point'):
        return cls(-p.toComplex())

    def transform(self, *objs):
        return list(map(lambda x: x.mobius(self), objs))


class Point:
    def __init__(self, r, phi, name=None, show=True):
        self.r = r
        self.phi = phi
        self.name = name
        self.x = r * cos(phi)
        self.y = r * sin(phi)
        if show: _plot_list.append(self)

    def valid(self):
        if _hyperbolic:
            return self.r <= 1
        else:
            return True

    def toComplex(self):
        return self.x + 1j * self.y

    @classmethod
    def fromComplex(cls, z: complex):
        return cls(abs(z), np.angle(z), show=False)

    @classmethod
    def fromXY(cls, xy: tuple):
        return cls.fromComplex(complex(xy[0], xy[1]))

    def plot(self, ax):
        ax.scatter(self.x, self.y, color='mediumslateblue')
        if self.name is not None:
            ax.text(self.x, self.y, self.name, fontsize=14)

    def setName(self, name):
        self.name = name
        return self

    def setPlot(self):
        _plot_list.append(self)
        return self

    def mobius(self, mob: Mobius):
        res = Point.fromComplex(mob.f(self.toComplex()))
        res.name = self.name
        return res


class GeodesicTwoPoints:
    def __init__(self, p1: Point, p2: Point, show=True):
        self.p1 = p1
        self.p2 = p2
        x1, y1 = p1.x, p1.y
        x2, y2 = p2.x, p2.y
        self.lamda = (_a2 + x1 * x2 + y1 * y2) / (x2 * y1 - x1 * y2)
        self.Xc = -(self.lamda * (y1 - y2) - (x1 + x2)) / 2
        self.Yc = -(self.lamda * (x2 - x1) - (y1 + y2)) / 2
        self.Radius = 1 / 2 * sqrt(((x1 - x2) ** 2 + (y1 - y2) ** 2) * (1 + self.lamda ** 2))
        if show: _plot_list.append(self)

    def plot(self, ax):
        draw_circle = plt.Circle((self.Xc, self.Yc), self.Radius, fill=False, color='blue')
        ax.add_artist(draw_circle)

    def mobius(self, mob: Mobius):
        p1m, p2m = self.p1.mobius(mob), self.p2.mobius(mob)
        try:
            res = GeodesicTwoPoints(p1m, p2m, show=False)
        except:
            res = None

        # Check if the circle is mapped to a line
        if res is not None and res.Radius < 1e6:
            return res
        else:
            if p1m.r > 1e-6:
                return GeodesicThroughOrigin.fromPoint(p1m, show=False)
            else:
                return GeodesicThroughOrigin.fromPoint(p2m, show=False)

    def toEuclidian(self):
        return self.Xc, self.Yc, self.Radius

    def intersect(self, other):
        if isinstance(other, GeodesicTwoPoints):
            points = find_intersections_two_circles(self.toEuclidian(), other.toEuclidian())
            lst = [Point.fromXY(p) for p in points]
            lst = list(filter(lambda x: x.valid(), lst))
            return lst
        if isinstance(other, GeodesicThroughOrigin):
            points = find_intersections_circle_line(self.toEuclidian(), other.toEuclidian())
            lst = [Point.fromXY(p) for p in points]
            lst = list(filter(lambda x: x.valid(), lst))
            return lst


class GeodesicThroughOrigin:
    def __init__(self, angle, show=True):
        self.angle = angle
        if show: _plot_list.append(self)

    @classmethod
    def fromPoint(cls, p: Point, show=True):
        return cls(p.phi, show)

    def points(self):
        return (0, 0), (cos(self.angle), sin(self.angle))

    def toEuclidian(self):
        return -sin(self.angle), cos(self.angle), 0

    def plot(self, ax):
        plt.axline(*self.points(), color='blue', linewidth=1)

    def intersect(self, other):
        if isinstance(other, GeodesicTwoPoints):
            points = find_intersections_circle_line(other.toEuclidian(), self.toEuclidian())
            lst = [Point.fromXY(p) for p in points]
            lst = list(filter(lambda x: x.valid(), lst))
            if len(lst) > 0:
                return lst[0]
            else:
                raise ValueError
        if isinstance(other, GeodesicThroughOrigin):
            return Point.fromXY((0, 0))


# global config
_hyperbolic = False
_a2 = -1 if _hyperbolic else 1
ax = init()
_plot_list = []


def show_all():
    for obj in _plot_list:
        obj.plot(ax)
    plt.show()


def clear_all():
    _plot_list.clear()


def plot(*objs):
    _plot_list.extend(list(objs))


"""
if __name__ == '__main__':
    a = Point(1, -pi / 6, 'A')
    b = Point(1, pi / 3, 'B')
    c = Point(1, pi, 'C')
    ab = GeodesicTwoPoints(a, b)
    bc = GeodesicTwoPoints(b, c)
    ca = GeodesicTwoPoints(c, a)
    jab = GeodesicThroughOrigin(pi / 4, show=False)
    jbc = GeodesicThroughOrigin(8 * pi / 5, show=False)
    jca = GeodesicThroughOrigin(-pi / 3, show=False)
    p = jab.intersect(ab).setName('P').setPlot()
    q = jbc.intersect(bc).setName('Q').setPlot()
    r = jca.intersect(ca).setName('R').setPlot()
    pq = GeodesicTwoPoints(p, q)
    qr = GeodesicTwoPoints(q, r)
    rp = GeodesicTwoPoints(r, p)

    clear_all()
    mob = Mobius.movePointToOrigin(p)
    lst = mob.transform(a, b, c, p, ab, bc, ca, p, q, r, pq, qr, rp)
    plot(*lst)
    '''
    clear_all()
    mob = Mobius.movePointToOrigin(q)
    lst = mob.transform(a, c, p, ca, p, q, r, pc, qa, rb)
    plot(*lst)
    '''
    show_all()
"""
if __name__ == '__main__':
    a = Point(0.5, -pi / 6, 'P')
    b = Point(0.6, pi / 3, 'Q')
    c = Point(0.6, pi, 'R')
    o = Point(0, 0, show=False)
    ab = GeodesicTwoPoints(a, b)
    bc = GeodesicTwoPoints(b, c)
    ca = GeodesicTwoPoints(c, a)
    # unit = GeodesicTwoPoints(Point(1, 1, show=False), Point(1, 2, show=False))
    p = ab.intersect(ca)[1].setName('A').setPlot()
    q = bc.intersect(ab)[1].setName('B').setPlot()
    r = ca.intersect(bc)[1].setName('C').setPlot()
    '''
    pq = GeodesicTwoPoints(p, q)
    qr = GeodesicTwoPoints(q, r)
    rp = GeodesicTwoPoints(r, p)
    '''
    clear_all()
    mob = Mobius.movePointToOrigin(a)
    lst = mob.transform(a, b, c, ab, bc, ca, q, r)
    plot(*lst)

    '''
    clear_all()
    mob = Mobius.movePointToOrigin(q)
    lst = mob.transform(a, c, p, ca, p, q, r, pc, qa, rb)
    plot(*lst)
    '''
    show_all()
