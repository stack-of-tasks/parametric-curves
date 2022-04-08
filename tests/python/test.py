import unittest

from numpy import matrix
from numpy.linalg import norm

from spline import (
    bezier,
    bezier6,
    curve_constraints,
    exact_cubic,
    polynom,
    spline_deriv_constraint,
)


class ParametricCurvesTests(unittest.TestCase):
    def test_all(self):
        waypoints = matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).transpose()
        waypoints6 = matrix(
            [[1.0, 2.0, 3.0, 7.0, 5.0, 5.0], [4.0, 5.0, 6.0, 4.0, 5.0, 6.0]]
        ).transpose()
        time_waypoints = matrix([0.0, 1.0])

        # testing bezier curve
        a = bezier6(waypoints6)
        a = bezier(waypoints, -1.0, 3.0)

        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)

        prim = a.compute_primitive(1)

        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue((a(t) == prim.derivate(t, 1)).all())
        self.assertTrue((prim(0) == matrix([0.0, 0.0, 0.0])).all())

        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.0
            self.assertTrue((a(t) == prim.derivate(t, 2)).all())
        self.assertTrue((prim(0) == matrix([0.0, 0.0, 0.0])).all())

        # testing bezier with constraints
        c = curve_constraints()
        c.init_vel = matrix([0.0, 1.0, 1.0])
        c.end_vel = matrix([0.0, 1.0, 1.0])
        c.init_acc = matrix([0.0, 1.0, -1.0])
        c.end_acc = matrix([0.0, 100.0, 1.0])

        a = bezier(waypoints, c)
        self.assertLess(norm(a.derivate(0, 1) - c.init_vel), 1e-10)
        self.assertLess(norm(a.derivate(1, 2) - c.end_acc), 1e-10)

        # testing polynom function
        a = polynom(waypoints)
        a = polynom(waypoints, -1.0, 3.0)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue(((a.derivate(0.4, 0) == a(0.4)).all()))
        a.derivate(0.4, 2)

        # testing exact_cubic function
        a = exact_cubic(waypoints, time_waypoints)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue(((a.derivate(0.4, 0) == a(0.4)).all()))
        a.derivate(0.4, 2)

        # testing spline_deriv_constraints
        c = curve_constraints()
        c.init_vel
        c.end_vel
        c.init_acc
        c.end_acc

        c.init_vel = matrix([0.0, 1.0, 1.0])
        c.end_vel = matrix([0.0, 1.0, 1.0])
        c.init_acc = matrix([0.0, 1.0, 1.0])
        c.end_acc = matrix([0.0, 1.0, 1.0])

        a = spline_deriv_constraint(waypoints, time_waypoints)
        a = spline_deriv_constraint(waypoints, time_waypoints, c)


if __name__ == "__main__":
    unittest.main()
