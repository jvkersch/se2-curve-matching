import sys
sys.path.append('..')

import unittest

import math
import numpy as np
import numpy.random
from curves import FigureEight, LogSpiral
from cfunctionals import energy_functional, energy_functional_slow


class EnergyValueTest(unittest.TestCase):
    def setUp(self):
        self.N = 100

        self.c0 = FigureEight(N=self.N).points
        self.c1 = LogSpiral(N=self.N).points

        self.m = 2
        self.h = 1./(self.N-1)

    def test_value(self):
        thetas = numpy.random.rand(self.N)
        X = np.array([np.cos(thetas), np.sin(thetas)]).T
    
        E0 = energy_functional_slow(thetas, self.c0, self.c1, self.m, self.h)
        E1 = energy_functional(X, self.c0, self.c1, self.m, self.h)

        print "Exact energy:", E0
        print "Optimized energy:", E1
        
        rel_diff = (E0-E1)/E0
        print "Relative difference:", rel_diff

        self.assertTrue(abs(rel_diff) < 1e-3)

if __name__ == '__main__':
    unittest.main()
