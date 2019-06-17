#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    # TODO: put package requirements here
]

setup_requirements = [
    'pytest-runner',
    # TODO(konrad): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'pytest',
    # TODO: put package test requirements here
]

setup(
    name='gffpandas',
    version='1.2.0',
    description="GFF annotations in panda data frames",
    long_description=readme + '\n\n' + history,
    author="Vivian A. Monzon, Konrad U. Förstner",
    author_email='vivian.monzon@stud-mail.uni-wuerzburg.de, konrad@foerstner.org',
    url='https://github.com/foerstner-lab/gffpandas',
    packages=find_packages(include=['gffpandas']),
    include_package_data=True,
    install_requires=requirements,
    license="ISC license",
    zip_safe=False,
    keywords='gffpandas',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
