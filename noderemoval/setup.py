from setuptools import setup, find_packages

setup(
    name='NodeImpact',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'networkx',
        'numpy',
        'giotto-tda',
        'scipy'
    ],
    author='Sina Roshandel',
    description='A package for processing sequences and analyzing persistence diagrams'
)