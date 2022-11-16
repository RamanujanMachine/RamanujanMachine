from setuptools import setup

setup(
    name="ramanujan",
    version="0.0.1",
    # due to Amazon EC2 issues, we probably can't go beyond Python3.9
    # (Python3.10.8 seems to have issues with pip). As such, I don't care if
    # "this is an outdated version of python" until one of 2 things happens:
    #     1. We no longer need EC2
    #     2. Someone manages a working EC2 instance with a later python version
    # Even then, as of now numba (including llvmlite) and ortools do not support python3.11, see respectively:
    # https://github.com/numba/numba/issues/8304 ; https://github.com/google/or-tools/issues/3515
    # Because of this, upgrading to Python3.11 will also require these packages to update.
    python_requires=">=3.8.10, <3.9",
    description="Ramanujan Machine",
    packages=['ramanujan'],
    install_requires=[
        'cycler>=0.10.0',
        'kiwisolver>=1.1.0',
        'matplotlib>=3.3.3',
        'mpmath>=1.2.1',
        'numpy>=1.19.5',
        'ordered-set>=3.1.1',
        'ortools>=7.4.7247',
        'pandas>=1.0.1',
        'protobuf>=3.11.3',
        'psutil>=5.5.1',
        'psycopg2>=2.8.6',
        'pybloom-live',
        'PyLaTeX>=1.3.1',
        'pyparsing>=2.4.6',
        'pytest>=6.2.4',
        'python-dateutil>=2.8.1',
        'pytz>=2019.3',
        'scipy>=1.6.0',
        'simplejson>=3.16.0',
        'six>=1.14.0',
        'sqlacodegen>=2.3.0',
        'sympy>=1.5.1',
        'xlrd>=2.0.1',
    ]
)
