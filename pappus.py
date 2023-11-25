from Poincare2 import *

setHyperbolic(False)
Q = Point(0.1, pi / 4, 'Q', show=False)
cir = Circle(Q, 0.6)
a, b, c = cir.arbitraryPoints('A B C', angle_range=(0, pi))
f, e, d = cir.arbitraryPoints('F E D', angle_range=(-pi, 0))
plot(a, b, c, d, e, f)
lines = Links({a: [e, f], b: [d, f], c: [d, e]})
plot(*lines)
p = lines[0].intersect(lines[2])[0].setName('P').setPlot()
q = lines[1].intersect(lines[4])[0].setName('Q').setPlot()
r = lines[3].intersect(lines[5])[0].setName('R').setPlot()
pq = GeodesicTwoPoints(p, q)
show_all()
