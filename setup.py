import setuptools

setuptools.setup(
    name="ramanujan", # Replace with your own username
    version="0.0.1",
    description="Ramanujan Machine",
    packages=['ramanujan'],
    install_requires=[
        'cycler>=0.10.0',
        'kiwisolver>=1.1.0',
        'matplotlib>=3.2.0',
        'mpmath>=1.1.0',
        'numpy>=1.18.1',
        'ordered-set>=3.1.1',
        'pandas>=1.0.1',
        'protobuf>=3.11.3',
        'PyLaTeX>=1.3.1',
        'pyparsing>=2.4.6',
        'python-dateutil>=2.8.1',
        'pytz>=2019.3',
        'six>=1.14.0',
        'sympy>=1.5.1',
        'ortools>=7.4.7247',
        'pybloom-live'
    ]
)