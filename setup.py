# from setuptools import setup
from numpy.distutils.core import Extension, setup
# from distutils.core import setup

setup(name='pygs',
        version = '0.14b',
        description = 'Package for geospatial analysis and plotting',
        url = 'https://ishwarmbhat.github.io/pygs.html',
        author = 'Ishwar M',
        author_email = 'ishwarm.bhat@gmail.com',
        license = 'MIT',
        packages = ['pygs'],
        package_data = {'pygs': ['data/landsea.nc']},
        ext_modules = [Extension('__vibeta__', ['src/vibeta_dp.f'])],
        zip_safe = False)
