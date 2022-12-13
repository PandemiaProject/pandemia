import pathlib
from setuptools import setup, find_packages

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')
version = {}
with open(here / "src/pandemia/version.py") as fp:
    exec(fp.read(), version)
setup(
    name='pandemia',
    version=version['VERSION'],
    description='An individual-based pandemic simulator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/?',
    author='J. Thompson, & S. Wattam',
    author_email='j.thompson22@imperial.ac.uk',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3 :: Only',
    ],
    keywords='epidemiology, model',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.8, <3.11',
    install_requires=[
                      'tqdm~=4.64',
                      'numpy~=1.23',
                      'scipy~=1.8',
                      'matplotlib~=3.5',
                      'pyproj~=3.3',
                      'joblib~=1.1',
                      'pyshp~=2.3',
                      'pygame~=2.1',
                      'Shapely~=1.8',
                      'pygad~=2.17',
                      'pyyaml'
                     ],

    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    # Similar to `install_requires` above, these must be valid existing
    # projects.
    extras_require={  # Optional
        # 'dev': ['check-manifest'],
        'test': ['pandas',
                 'pdoc',
                 'pytest',
                 'pytest-cov',
                 "icecream"
                 ],
    },
    # package_data={'pandemia': ["C/build/*"]},
    include_dirs=["C/build/"],
    entry_points={
        'console_scripts': [
            'pandemia=pandemia.__main__:main'
        ],
    }
)
