#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Dr. Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2021 Dr. Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CGRdb.
#
#  CGRdb is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from os import listdir
from pathlib import Path
from setuptools import setup, find_packages


version = '0.0.1'

setup(
    name='CGRdb',
    version=version,
    packages=find_packages(),
    url='https://github.com/stsouko/CGRdb',
    license='LGPLv3',
    author=['Dr. Timur Gimadiev', 'Dr. Ramil Nugmanov'],
    author_email=['timur.gimadiev@gmail.com', 'nougmanoff@protonmail.com'],
    python_requires='>=3.9.0',
    install_requires=['CGRtools>4.1.20,<4.2',  'StructureFingerprint>=1.24',
                      'pony>=0.7.14,<0.8', 'psycopg2-binary>=2.8.6', 'tqdm>=4.61.0',
                      'datasketch>=1.5.3'],
    extras_require={'autocomplete': ['argcomplete'],
                    'index': ['pyroaring>=0.2.9', 'aiohttp>=3.7', 'datasketch>=1.5.3', 'tqdm>=4.61.0']},
    #long_description=(Path(__file__).parent / 'README.md').open().read(),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Software Development :: Libraries :: Python Modules',
                 'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.9',
                 ]
)
