from setuptools import setup
    
long_description= ''
previous_line= ''
with open('README.md') as dfile:
    for line in dfile:
        if '[!' in line: continue
        long_description+= line

setup(name='kimmy',
      version='0.1',
      description='Chemical evolution in galaxies',
      long_description=long_description,
      long_description_content_type='text/markdown',     
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
