# -*- coding: utf-8 -*-
#
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from pony.orm import PrimaryKey, Required, Set, Optional
from .config import Entity


class ModificationStamp(Entity):
    id = PrimaryKey(int, auto=True)
    fileName = Optional(str)
    time = Required(str, unique=True)


class DataSource(Entity):
    id = PrimaryKey(int, auto=True)
    name = Required(str, unique=True)
    molecules = Set('MoleculeSourceID')
    reactions = Set('ReactionSourceID')


class MoleculeSourceID(Entity):
    source = Required('DataSource')
    molecule = Required('Molecule')
    id = PrimaryKey(source, molecule)


class ReactionSourceID(Entity):
    source = Required('DataSource')
    reaction = Required('Reaction')
    id = PrimaryKey(source, reaction)




__all__ = ['ModificationStamp', 'DataSource', 'MoleculeSourceID', 'ReactionSourceID']
