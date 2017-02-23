from setuptools import setup

import dingo

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='dingo',
      version=dingo.__version__,
      description=dingo.__description__,
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GPLv3',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics',
      url=dingo.__url__,
      author=dingo.__author__,
      author_email=dingo.__author_email__,
      license=dingo.__license__,
      packages=['dingo'],
      install_requires=[
        'click',
        'pandas',
        'progressbar2'
      ],
      test_suite='nose.collector',
      tests_require=[],
      entry_points={
          'console_scripts': ['dingo=dingo.dingo:main'],
      },
      include_package_data=True,
      zip_safe=False)
