#!/usr/bin/env python3
from setuptools import setup
from setuptools.command.install import install
import os
import sys

__version__ = '0.1.34'

def readme():
    with open('README.md') as f:
        return f.read()

def check_dir_write_permission(directory):
    if os.path.isdir(directory) and not os.access(directory, os.W_OK):
        sys.exit('Error: no write permission for ' + directory + '  ' +
                 'Perhaps you need to use sudo?')

class AmpliconTyperInstall(install):

    def run(self):
        check_dir_write_permission(self.install_lib)
        install.run(self)



setup(name='AmpliconTyper',
      version=__version__,
      description='AmpliconTyper',
      long_description=readme(),
      python_requires='>=3.10',
      classifiers=['Development Status :: Beta',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research'],
      keywords='PCR amplicon primers',
      url='https://github.com/AntonS-bio/AmpliconTyper.git',
      author='Anton Spadar',
      author_email='',
      packages=['scripts'],
      include_package_data=True,
      entry_points={'console_scripts': ['classify = classify:main', 'train = train:main', 'genotyper_utilities = genotyper_utilities:main']},
      scripts=[
          'scripts/check_for_update.py',
          'scripts/classifier_report.py',
          'scripts/classify.py',
          'scripts/data_classes.py',
          'scripts/genotyper_utilities.py',
          'scripts/hierarchy_utils.py',
          'scripts/html_head.txt',
          'scripts/input_processing.py',
          'scripts/inputs_validation.py',
          'scripts/map.py',
          'scripts/model_manager.py',
          'scripts/read_classifier.py',
          'scripts/read_quality_check.py',       
          'scripts/reporting_classes.py',
          'scripts/train.py'
      ],
      cmdclass={'install': AmpliconTyperInstall}
)
