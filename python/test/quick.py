from parametric_curves import polynomial, curve_constraints, forcecurve, spline

from numpy import matrix

cc = curve_constraints()
cc.init_vel = matrix(range(3)).T

fc = forcecurve()
sp = spline()
pn = polynomial(matrix('1 2 3;4 5 6;7 8 9'))

assert pn.min() == 0.0
assert pn.max() == 1.0
