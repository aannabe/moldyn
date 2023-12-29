#!/usr/bin/env python

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='moldyn',
    version='0.1.0',
    description='Molecular Dynamics using Lennard-Jones Potential',
    long_description=readme,
    author='Abdulgani Annaberdiyev',
    author_email='annaberdiyev@gmail.com',
    url='https://github.com/aannabe/moldyn',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

