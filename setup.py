from distutils.core import setup

setup(
    name='Discrepancy',
    version='0.1.0',
    author='Joris Vankerschaver',
    author_email='joris.vankerschaver@gmail.com',
    packages=['discrepancy'],
    scripts=['tools/draw_curves.py'],
    url='https://github.com/jvkersch/se2-curve-matching',
    license='LICENSE.txt',
    description='Discrepancy of planar curves',
    long_description=open('README.md').read(),
)
