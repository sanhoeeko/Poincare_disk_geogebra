import random as rd
from math import sin, cos, sqrt, pi, tan, atan, tanh, atanh

import matplotlib.pyplot as plt
import numpy as np

from Euclidean import *


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


# global config
hyperbolic = True
a2 = -1 if hyperbolic else 1
ax = init()
plot_list = []


def setHyperbolic(tf: bool):
    global hyperbolic, a2
    hyperbolic = tf
    a2 = -1 if hyperbolic else 1


def safe_atanh(x):
    if x >= 1:
        return math.inf
    else:
        return atanh(x)


def tg(x):
    if hyperbolic:
        return tanh(x)
    else:
        return tan(x)


def atg(x):
    if hyperbolic:
        return safe_atanh(x)
    else:
        return atan(x)


class Mobius:
    def __init__(self, u: complex):
        self.u = u

    def __call__(self, obj):
        return obj.mobius(self)

    def f(self, z):
        return (self.u + z) / (1 - (self.u.conjugate() * z) / a2)

    @classmethod
    def movePointToOrigin(cls, p: 'Point'):
        return cls(-p.toComplex())

    def transform(self, *objs):
        return list(map(lambda x: x.mobius(self), objs))

    def inv(self):
        return Mobius(-self.u)


class Point:
    def __init__(self, r, phi, name=None, show=True):
        self.r = r
        self.phi = phi
        self.name = name
        self.x = r * cos(phi)
        self.y = r * sin(phi)
        if show: plot_list.append(self)

    def valid(self):
        if hyperbolic:
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
        plot_list.append(self)
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
        self.lamda = (a2 + x1 * x2 + y1 * y2) / (x2 * y1 - x1 * y2)
        self.Xc = -(self.lamda * (y1 - y2) - (x1 + x2)) / 2
        self.Yc = -(self.lamda * (x2 - x1) - (y1 + y2)) / 2
        self.Radius = 1 / 2 * sqrt(((x1 - x2) ** 2 + (y1 - y2) ** 2) * (1 + self.lamda ** 2))
        if show: plot_list.append(self)

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

    def intersect(self, other) -> list[Point]:
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

    def infinityPoints(self) -> (Point, Point):
        points = find_intersections_two_circles(self.toEuclidian(), (0, 0, 1))
        return tuple([Point.fromXY(p) for p in points])

    def proportionalDivisionPoint(self, prop):
        return ProportionalDivisionPoint(self.p1, self.p2, prop)


class GeodesicThroughOrigin:
    def __init__(self, angle, show=True):
        self.angle = angle
        if show: plot_list.append(self)

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
                raise ValueError('No intersection: hyper parallel.')
        if isinstance(other, GeodesicThroughOrigin):
            return Point.fromXY((0, 0))


class Polygon:
    p: list[Point]
    edges: list[GeodesicTwoPoints]

    def __init__(self, *points: tuple, names=None):
        if names is not None:
            names = names.split(' ')
        else:
            names = [None] * len(points)
        assert len(names) == len(points)
        self.p = [Point(points[i][0], points[i][1], names[i]) for i in range(len(points))]
        self.edges = []
        for i in range(len(self.p)):
            self.edges.append(GeodesicTwoPoints(self.p[i], self.p[(i + 1) % len(self.p)]))

    @classmethod
    def fromPoints(cls, *points: Point):
        lst = [(p.r, p.phi) for p in points]
        return cls(*lst)


class Circle:
    def __init__(self, center: Point, radius, show=True):
        self.center = center
        rho = tg(radius)
        self.rho = rho
        R2 = center.r ** 2
        self.Xc = a2 * (a2 + rho ** 2) / (1 - R2 * rho ** 2) * center.x
        self.Yc = a2 * (a2 + rho ** 2) / (1 - R2 * rho ** 2) * center.y
        self.Radius = (a2 + R2) * rho / abs((1 - R2 * rho ** 2))
        if show: plot_list.append(self)

    def plot(self, ax):
        draw_circle = plt.Circle((self.Xc, self.Yc), self.Radius, fill=False, color='green')
        ax.add_artist(draw_circle)

    def point(self, angle):
        phi = angle
        p = Point(self.rho, phi, show=False)
        mob = Mobius.movePointToOrigin(self.center).inv()
        p = p.mobius(mob)
        return p

    def arbitraryPoint(self, angle_range=None):
        if angle_range is None:
            angle_range = (0, 2 * pi)
        phi = rd.uniform(*angle_range)
        return self.point(phi)

    def arbitraryPoints(self, names: str, angle_range=None):
        if angle_range is None:
            angle_range = (0, 2 * pi)
        names = names.split(' ')
        ranges = np.linspace(*angle_range, len(names) + 1)
        lst = []
        for i in range(len(names)):
            rg = (ranges[i], ranges[i + 1])
            lst.append(self.arbitraryPoint(rg).setName(names[i]))
        return lst


def ProportionalDivisionPoint(a: Point, b: Point, prop):
    mob = Mobius.movePointToOrigin(a)
    b_prime = b.mobius(mob)
    r, phi = b_prime.r, b_prime.phi
    k = prop
    c = Point(tg(k * atg(r)), phi, show=False)
    return c.mobius(mob.inv())


def Midpoint(a: Point, b: Point):
    return ProportionalDivisionPoint(a, b, 0.5)


def Links(dic: dict) -> list[GeodesicTwoPoints]:
    lines = []
    for k in dic.keys():
        for end in dic[k]:
            lines.append(GeodesicTwoPoints(k, end))
    return lines


def show_all():
    for obj in plot_list:
        obj.plot(ax)
    plt.show()


def clear_all():
    plot_list.clear()


def plot(*objs):
    plot_list.extend(list(objs))
