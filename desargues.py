from Poincare2 import *

setHyperbolic(True)
# ABC = Polygon((0.5, 2.1), (0.8, -2.3), (0.7, 0.2), names='A B C') # this data for the sphere
ABC = Polygon((0.3, 2.1), (0.2, -2.3), (0.4, 0.2), names='A B C')  # this data for the hyperbolic plane
A, B, C = ABC.p
AB, BC, CA = ABC.edges
p = Point(0.8, -0.8, 'P')
pa, pb, pc = Links({p: (A, B, C)})
# 1.5, 1.7, 2 for the sphere
a = pa.proportionalDivisionPoint(1.5).setName('a').setPlot()
b = pb.proportionalDivisionPoint(0.9).setName('b').setPlot()
c = pc.proportionalDivisionPoint(0.3).setName('c').setPlot()
abc = Polygon.fromPoints(a, b, c)
ab, bc, ca = abc.edges
R = AB.intersect(ab)[0].setName('R').setPlot()
S = BC.intersect(bc)[0].setName('S').setPlot()
T = CA.intersect(ca)[0].setName('T').setPlot()
'''
r = AB.intersect(ab)[1].setName('r').setPlot()
s = BC.intersect(bc)[1].setName('s').setPlot()
t = CA.intersect(ca)[1].setName('t').setPlot()
'''
l = GeodesicTwoPoints(S, T)
show_all()
