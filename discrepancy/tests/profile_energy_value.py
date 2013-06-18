"""
Note to self: before running this, enable global option for profiling in 
cfunctionals.pyx

"""

import sys
sys.path.append('..')

import pstats, cProfile

import numpy as np
import numpy.random
from curves import FigureEight, LogSpiral
from cfunctionals import energy_functional

def setup(N=100):
    c0 = FigureEight(N=N).points
    c1 = LogSpiral(N=N).points

    m = 2
    h = 1./(N-1)

    thetas = numpy.random.rand(N)
    X = np.array([np.cos(thetas), np.sin(thetas)]).T

    return X, c0, c1, m, h


if __name__ == '__main__':

    X, c0, c1, m, h = setup()
    cProfile.runctx("energy_functional(X, c0, c1, m, h)", 
                    globals(), locals(), "Profile.prof")

    s = pstats.Stats("Profile.prof")
    s.strip_dirs().sort_stats("time").print_stats()
