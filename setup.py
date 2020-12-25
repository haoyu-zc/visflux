from setuptools import setup, find_packages

setup(name='d3flux_vis',
      version='0.3',
      description='A d3.js-based metabolic visualization tool for cobra models. A variant version of original d3flux.',
      url='https://github.com/ikspike/d3flux_vis',
      download_url='https://github.com/ikspike/d3flux_vis',
      author='Peter St. John, Haoyu Zhang',
      author_email='peter.stjohn@nrel.gov, haoyu_z@outlook.com',
      license='MIT',
      packages=find_packages(),
      install_requires=['pandas', 'cobra', 'jinja2', 'ipython', 'csscompressor'],
      package_data={'d3flux_vis': ['templates/*', 'main/*', 'include/*', 'vendor/*']},
      )
