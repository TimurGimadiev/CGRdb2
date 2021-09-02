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
from datasketch import MinHashLSH
organic_set = {'B', 'C', 'N', 'O', 'Si', 'P', 'S', 'Se', 'F', 'Cl', 'Br', 'I'}


def validate_molecule(mol):
    if int(mol):
        return False
    for n, a in mol.atoms():
        if a.atomic_symbol not in organic_set:
            return False
    return True


class MinHashLSH(MinHashLSH):

    def _byteswap(self, hs):
        return int.from_bytes(hs.byteswap().data, 'big')

    def _hashed_byteswap(self, hs):
        return self.hashfunc(int.from_bytes(hs.byteswap().data, 'big'))


__all__ = ['validate_molecule', 'MinHashLSH']
