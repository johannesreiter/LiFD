import setuptools

with open('README.md', 'r', encoding='utf8') as f:
    long_description = f.read()

__version__ = 'unknown'
exec(open('lifd/version.py').read())

setuptools.setup(
      name='lifd',                                  # package name
      # packages=setuptools.find_packages(),
      packages=['lifd', 'lifd.databases', 'lifd.predictors', 'lifd.test'],
      version=__version__,
      description='LiFD is a two-phase algorithm that predicts likely functional driver (LiFD) mutations. ',
      long_description=long_description,
      long_description_content_type='text/markdown',
      install_requires=['numpy', 'scipy', 'pandas', 'statsmodels'],
      url='https://github.com/johannesreiter/LiFD_dev',
      author='Johannes Reiter',
      author_email='johannes.reiter@stanford.edu',
      license='GNUv3',
      classifiers=[
        'Programming Language :: Python :: 3.6',
      ],
      test_suite='lifd.test'
)
