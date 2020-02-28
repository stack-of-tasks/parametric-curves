import unittest

import eigenpy
from numpy import matrix

from libparametric_curves_pywrap import curve_constraints, forcecurve, polynomial, spline

eigenpy.switchToNumpyMatrix()


class ParametricCurvesQuickTests(unittest.TestCase):
    def test_quick(self):
        cc = curve_constraints()
        cc.init_vel = matrix(range(3)).T

        forcecurve()
        spline()
        pn = polynomial(matrix('1 2 3;4 5 6;7 8 9'))

        self.assertEqual(pn.min(), 0.0)
        self.assertEqual(pn.max(), 1.0)


if __name__ == '__main__':
    unittest.main()
