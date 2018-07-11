from setuptools import setup
    
setup(name='kimmy',
      version='0.1',
      description='Chemical evolution in galaxies',
      author='Jo Bovy',
      author_email='bovy@astro.utoronto.ca',
      license='MIT',
      url='http://github.com/jobovy/kimmy',
      package_dir = {'kimmy/': ''},
      packages=['kimmy'],
      package_data={"": ["README.md","LICENSE"]},
      include_package_data=True,
      install_requires=['numpy','astropy']
      )
