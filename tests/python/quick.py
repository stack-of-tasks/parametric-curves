import unittest

from numpy import array
from parametric_curves import curve_constraints, forcecurve, polynomial, spline


class ParametricCurvesQuickTests(unittest.TestCase):
    def test_quick(self):
        cc = curve_constraints()
        cc.init_vel = array([0, 1, 2])

        forcecurve()
        spline()
        pn = polynomial(array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))

        self.assertEqual(pn.min(), 0.0)
        self.assertEqual(pn.max(), 1.0)


if __name__ == "__main__":
    unittest.main()
