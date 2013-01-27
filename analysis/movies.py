import os
import shutil
import sys
sys.path.append('../src')

from subprocess import call

from matching import MatchingProblem
from curves import DiscreteCurve, FigureEight, LogSpiral
from interpolation import LiegroupInterpolationMethod
from plotting import LineImage

TMP_FOLDER = '/Users/Joris/tmp/'

def make_movie(c0, c1, theta0, movie_name, 
               m=2, dims=(-1.5, 1.5, -1.5, 1.5), numframes=100):
    """
    Create linear interpolation movie.
    """
    xmin, xmax, ymin, ymax = dims

    mp = MatchingProblem(c0, c1, m)
    g, res, n = mp.match(theta0)

    interp = LiegroupInterpolationMethod(c0, c1, g)

    # Generate frames
    for n in xrange(0, numframes):

        epsilon = float(n)/(numframes-1)
    
        c0_transformed = interp.transform(epsilon)
        im = LineImage(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax); 
        im.plot_curve(c0_transformed)
        im.save(TMP_FOLDER + 'still%04d.png' % n)
        im.close()

    # Compile movie
    cwd = os.getcwd()
    os.chdir(TMP_FOLDER)
    os.system('/Users/Joris/localsoft/bin/ffmpeg -i still%%04d.png %s' % movie_name)
    shutil.move(movie_name, cwd)
    os.chdir(cwd)

