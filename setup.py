from setuptools import setup

setup(name='pygs',
        version = '0.13',
        description = 'Package for geospatial analysis and plotting',
        url = 'https://irondrummer.github.io/pygs.html',
        author = 'Ishwar M',
        author_email = 'ishwarm.bhat@gmail.com',
        license = 'MIT',
        packages = ['pygs'],
        package_data={'pygs': ['data/landsea.nc']},
        zip_safe = False)
