import timeit

setup = """
import sys
sys.path.append('..')

import numpy as np
import numpy.random
from curves import FigureEight, LogSpiral
from cfunctionals import energy_functional, energy_functional_angles

N = 100

c0 = FigureEight(N=N).points
c1 = LogSpiral(N=N).points

m = 2
h = 1./(N-1)

thetas = numpy.random.rand(N)
X = np.array([np.cos(thetas), np.sin(thetas)]).T

"""

if __name__ == '__main__':
    n = 100000
    t = timeit.timeit('energy_functional(X, c0, c1, m, h)', 
                      number=n, setup=setup)

    print "Time to run %e iterations (direct): %f." % (n, t)

    t = timeit.timeit('energy_functional_angles(thetas, c0, c1, m, h)', 
                      number=n, setup=setup)

    print "Time to run %e iterations (using angles): %f." % (n, t)
