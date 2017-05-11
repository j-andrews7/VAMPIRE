"""
VENUSAR is a set of tools to identify and prioritize variants that alter TF binding, affect gene
expression, and change regulatory element activity.
"""
from setuptools import find_packages, setup

dependencies = ['click>=6', 'pyfaidx', 'scipy', 'biopython']

setup(
    name='venusar',
    version='0.1.0',
    url='https://github.com/j-andrews7/VENUSAR',
    license='MIT',
    author='Payton Lab',
    author_email='jared.andrews07@gmail.com',
    description=('VENUSAR is a set of tools to identify and prioritize variants that alter TF binding,'
                 ' affect gene expression, and change regulatory element activity'),
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    keywords='variant genomics biology bioinformatics epigenetic chromatin transcription motif research',
    entry_points={
        'console_scripts': [
            'venusar = venusar.venusar:cli',
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
