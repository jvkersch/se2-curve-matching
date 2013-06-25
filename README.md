Overview
--------

This repository provides a new algorithm for matching/registration of planar, parametrized curves using translations and rotations (i.e. elements of the special Euclidean group SE(2)).

How to use this package
-----------------------

It's recommended to create a virtual environment for this package, e.g.

    virtualenv curve_matching
	source curve_matching/bin/activate

Then install the package as follows:

    pip install -e git+git@github.com:jvkersch/se2-curve-matching.git@master#egg=Discrepancy

If you've downloaded the package as a zip file, install it instead as follows:

	pip install <path-to-zip-file>/se2-curve-matching-master.zip

This will add the discrepancy code to your site-packages. The directory `curve_matching/src/discrepancy/worksheets` contains a bunch of ipython worksheets which illustrate working with the package. To fire up ipython, do 

    cd curve_matching/src/discrepancy/worksheets
	ipython notebook --pylab inline	

To re-do the simulations in the paper, just run the four worksheets. This will take some time, and will put the figures of the paper in the `Figures` subfolder and some (temporary) data in `data`.
	
Further reading
---------------

Lyle Noakes, Darryl D. Holm, Joris Vankerschaver: _Relative Geodesics in the Special Euclidean Group_, [arXiv:1305.3572](http://arxiv.org/abs/1305.3572)
