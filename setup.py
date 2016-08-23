"""
Driven is a full-fledged bioinformatics suite geared towards easy 
integration of multiple types of big data to make reasonable and interesting biological conclusions.
"""
from setuptools import find_packages, setup

dependencies = ['click', 'pyfaidx', 'numpy', 'scipy', 'biopython']

setup(
    name='driven',
    version='0.1.0',
    url='https://github.com/j-andrews7/DRIVEN',
    license='MIT',
    author='Jared Andrews',
    author_email='jared.andrews07@gmail.com',
    description='Driven is a full-fledged bioinformatics suite geared towards easy integration of multiple types of big data to make reasonable and interesting biological conclusions.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'driven = driven.cli:main',
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
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
        'Operating System :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
