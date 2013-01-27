"""
Computes the discrepancy between a set of handwritten digits, in files
'digit-name_raw.mat' (where digit-name ranges from 'zero' to 'nine').

"""

import sys
sys.path.append('../src')

import numpy as np

from matching import MatchingProblem
from curves import DiscreteCurve
from tools import separate_energies
from movies import make_movie

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

from options import DATA_FOLDER

########## Preprocessing ###########

from curves import _load_raw_data, sample_discrete_curve
from plotting import LineImage

def load_and_normalize(filename):
    """
    Load curve from file, center on mean, scale such that height is unity.

    """    
    x, y = _load_raw_data(DATA_FOLDER + filename, 'c0')
    c0 = sample_discrete_curve(x, y)
    c0.translate(-c0[0, :])
    w, h = c0.dimensions()
    c0.scale(1/h)
    c0.translate(-c0.mean())
    
    return c0

def make_figure(c0, figname):
    "Make standard figure of curve."
    
    fig = plt.figure(figsize=(3,3))
    ax = fig.add_subplot(111)

    ax.plot(c0.points[:, 0], c0.points[:, 1])
    ax.set_aspect('equal')

    ax.set_xlim((-.6, .6))
    ax.set_ylim((-.6, .6))

    plt.savefig(figname)
    plt.close()

def sample_from_curve(c0):
    "Uniformize number of sample points."
    c1 = sample_discrete_curve(c0.points[:, 0], c0.points[:, 1], N=100)
    return c1

def preprocess_digits():
    """
    Preprocess every handwritten digit and save result.

    """
    for digit in xrange(0, 10):
        c0 = load_and_normalize('%d_raw.mat' % digit)
        make_figure(c0, '%d_fig.pdf' % digit)
        sample_from_curve(c0).save(DATA_FOLDER + '/%d_processed.mat' % digit)


########## Matching ###########

def digits_pairwise_match():
    """
    Find discrepancy and minimizing angle for each pair of digits.

    """

    discrepancy = np.empty((10, 10))
    theta_min = np.empty((10, 10))

    for m in xrange(0, 10):
        c0 = DiscreteCurve()
        c0.load(DATA_FOLDER + '/%d_processed.mat' % m)
        c0.translate(-c0[0,:])
    
        print "*** SOURCE: %d ***" % m
    
        for n in xrange(0, 10):
            c1 = DiscreteCurve()
            c1.load(DATA_FOLDER + '/%d_processed.mat' % n)
            c1.translate(-c1[0,:])
        
            print "*** TARGET: %d ***" % n
        
            mp = MatchingProblem(c0, c1, 2)
            energies = mp.sweep_theta_range(40)
        
            E_conv, E_div = separate_energies(energies)
            if len(E_conv) == 0:
                print "*** ERROR: No convergence for %d to %d ***" % (m, n)
        
            argmin = E_conv[:, 1].argmin()
            theta, E_min = E_conv[argmin, 0:2]
        
            print "* theta_min = %f, E_min = %f *" % (theta, E_min)
        
            discrepancy[m, n] = E_min
            theta_min[m, n] = theta

    return theta_min, discrepancy


########## Generate movies ###########

def digits_one_movie(m, n, theta):
    """
    Create movie of interpolation between two digits.

    Parameters
    ----------

    m, n: integer
          Index of source and target digit.

    theta: real
           Initial angle for matching problem.

    """

    c0 = DiscreteCurve()
    c0.load(DATA_FOLDER + '/%d_processed.mat' % m)
    c0.translate(-c0[0,:])

    print "*** SOURCE: %d ***" % m

    c1 = DiscreteCurve()
    c1.load(DATA_FOLDER + '/%d_processed.mat' % n)
    c1.translate(-c1[0,:])

    print "*** TARGET: %d ***" % n

    movie_name = '%d_to_%d.avi' % (m, n)
    make_movie(c0, c1, theta, movie_name, 
               m=2, dims=(-1, 1.5, -1.5, 1), numframes=100)

    print "*** COMPILED MOVIE: %s ***" % movie_name

def create_digits_movies():
    """
    Create movies for interpolation of each pair of handwritten digits.

    """
    theta_min = np.load(DATA_FOLDER + 'digits_theta_min.npy')

    for m in xrange(0, 10):
        for n in xrange(0, 10):
            theta = theta_min[m, n]
            digits_one_movie(m, n, theta)


if __name__ == '__main__':
    """
    Preprocess digits, compute discrepancies, and generate movies.

    """
    print "Preprocessing digits..."
    preprocess_digits()

    print "Matching..."
    #theta_min, discrepancy = digits_pairwise_match()

    #np.save(DATA_FOLDER + 'digits_theta_min.npy', theta_min)
    #np.save(DATA_FOLDER + 'digits_discrepancy.npy', discrepancy)

    print "Creating movies..."
    create_digits_movies()
