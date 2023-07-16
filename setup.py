import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'rga',
    version = '0.0.1',
    author='Roman Sermyagin',
    author_email='roman.ags@gmail.com',
    description=('Read CG-6 or CG-5 data file and proc by relative gravity adjustment'),
    keywords='scintrex cg-6 gravimeter',
    packages=['src', 'longman'],
    long_description=read('README.md'),
    install_requires=[],
    scripts=['rga.py']
)