Overview
--------

This repository provides a new algorithm for matching/registration of planar, parametrized curves using translations and rotations (i.e. elements of the special Euclidean group SE(2)).

How to use this package
-----------------------

It's recommended to install this package in a virtual environment:

    virtualenv curve_matching
	source curve_matching/bin/activate
	pip install -e git+git@github.com:jvkersch/se2-curve-matching.git@master#egg=Discrepancy
	
This will add the discrepancy code to your site-packages. 
	
Further reading
---------------

Lyle Noakes, Darryl D. Holm, Joris Vankerschaver: _Relative Geodesics in the Special Euclidean Group_, [arXiv:1305.3572](http://arxiv.org/abs/1305.3572)
