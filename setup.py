from setuptools import setup, find_packages

setup(name='visflux',
      version='0.3',
      description='A d3.js-based metabolic visualization tool for cobra models. A variant version of original d3flux.',
      url='https://github.com/ikspike/visflux',
      download_url='https://github.com/ikspike/visflux',
      author='Peter St. John, Haoyu Zhang',
      author_email='peter.stjohn@nrel.gov, haoyu_z@outlook.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['pandas', 'cobra', 'jinja2', 'ipython', 'csscompressor'],
      package_data={'visflux': ['templates/*', 'main/*', 'include/*', 'vendor/*']},
      )
