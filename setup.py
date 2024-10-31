from setuptools import setup


setup(name='run_enterprise',
      version=1.0,
      description='Run "enterprise" on a single pulsar',
      author='Michael Keith',
      author_email='mkeith@pulsarastronomy.net',
      packages=['run_enterprise','run_enterprise.model_components','run_enterprise.samplers'],
      install_requires=['enterprise-pulsar','libstempo','astropy','dynesty','corner','emcee','matplotlib','numpy'],
      scripts=['bin/run_enterprise','utilities/make_pulsar_plots.py','chop_tim.py','tidy_glitches.py'],
      )
