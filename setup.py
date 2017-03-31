import setuptools


LONG_DESCRIPTION = """
**aospy**: automated gridded climate data analysis and management

A framework that enables automated calculations using gridded climate data.
Following some basic description of where your data lives and defining any
functions of variables stored in that data you want to compute, aospy enables
you to fire off an arbitrary number of calculations using that data.

Important links
---------------
- HTML documentation: http://aospy.readthedocs.io/en/latest
- Mailing list: https://groups.google.com/d/forum/aospy
- Issue tracker: https://github.com/spencerahill/aospy/issues
- Source code: https://github.com/spencerahill/aospy
"""


setuptools.setup(
    name="aospy",
    version="0.1.2",
    packages=setuptools.find_packages(),
    author="aospy Developers",
    author_email="aospy@googlegroups.com",
    description="Automated gridded climate data analysis and management",
    long_description=LONG_DESCRIPTION,
    install_requires=['numpy >= 1.7',
                      'scipy >= 0.16',
                      'pandas >= 0.15.0',
                      'netCDF4 >= 1.2',
                      'toolz >= 0.7.2',
                      'dask >= 0.12',
                      'xarray >= 0.9.1',
                      'cloudpickle >= 0.2.1'],
    tests_require=['pytest >= 2.7.1',
                   'pytest-catchlog >= 1.0'],
    package_data={'aospy': ['test/data/netcdf/*.nc']},
    scripts=['aospy/examples/aospy_main.py',
             'aospy/examples/example_obj_lib.py',
             'aospy/examples/tutorial.ipynb'],
    license="Apache",
    keywords="climate science netcdf xarray",
    url="https://github.com/spencerahill/aospy",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Atmospheric Science'
    ]
)
